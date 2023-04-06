#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Recent Update: 10/23/2021
@author: Miguel Guardado Miguel.Guardado@ucsf.edu
"""
import networkx as nx
import pandas as pd
import numpy as np
import subprocess
def add_contig(vcf_file, contig_length=None):
    """
    This function will add a Contig Line to the VCF file since Slim appears to not have this included in their
    output vcf function. If Slim has the ability to do this and you are reading this, please contact me!!

    :param vcf_file: vcf file, specifically after slim simulations are done
    :param contig_length: Length of the genome this file was simulated under.
    :return:
    """

    awk_cmd = """'/^#CHROM/ { printf("##contig=<ID=1,length=195471971>\\n");} {print;}'"""
    cmd = f"awk {awk_cmd} {vcf_file} > tmp.vcf"
    subprocess.run([cmd], shell=True)

    mv_cmd = f'mv tmp.vcf {vcf_file}'
    subprocess.run([mv_cmd], shell=True)

def update_vcf_header(vcf_filepath, fam_graph):
    """
    This function will take the simulated genomes(vcf file) from our software, and utilizes bcftools to
    reindex the file to the list of IDs used in the family pedigree graph. My software outputs the individuals genomes
    as the indexes of thier ID's from lowest to highest. Utilizing bcftools reheader function we rewrite the
    node id as a list of individual
    :param vcf_filepath: output vcf file path from ped_slim.
    :param sorted_fam_list: sorted list of individual ID based on family pedigree graph.
    :return: vcf file

    comment 02/24/2022:  Currently, bcftools does not rewrite a vcf file in place, instead we have to create a temporary
    file, remove the original file, and then rename the new file the original files names. Could be a big cost on
    larger datasets. Meaning
    """
    #  First we check if the length of the list inputted is the correct length found in the vcf file

    find_length_cmd = f"bcftools query -l {vcf_filepath} | wc -l"
    process = subprocess.Popen([find_length_cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    vcf_indiv_length, err = process.communicate()
    vcf_indiv_length = int(vcf_indiv_length.decode('ascii'))

    di_family = nx.read_edgelist(fam_graph, create_using=nx.DiGraph())

    # Create a numpy array of the list of nodes (individuals) in our family pedigree.
    fam_list = np.array(di_family.nodes)
    fam_list = np.sort(fam_list.astype(int)).astype(str)

    if len(fam_list) == vcf_indiv_length:
        np.savetxt('sorted_sampleid.txt', fam_list, fmt='%s')

        reheader_cmd = f'bcftools reheader -s sorted_sampleid.txt -o tmp.vcf {vcf_filepath}'
        subprocess.run([reheader_cmd], shell=True)

        rm_cmd = f'rm sorted_sampleid.txt {vcf_filepath}'
        subprocess.run([rm_cmd], shell=True)

        mv_cmd = f'mv tmp.vcf {vcf_filepath}'
        subprocess.run([mv_cmd],shell=True)


def convert_networkx_to_ped_wprofiles(networkx_file, output_prefix, profiles_file):
    """
    This function will be used to convert a networkx based family pedigree into a .ped file that is often used for other
    tools such as plink. This function will take in an additional file for individuals sex and birth year  profiles.

    :param networkx_file: File path where the networkx graph is found.
    :param networkx_file: File path where the profiles file is found
    :param output_prefix: Directory  where you want to output the file in

    :return: will return a ped file of the networkx graph that was inputted in.
    """
    sex_map = {'male':1, 'female':2}


    ped_dir_edgelist = nx.read_edgelist(networkx_file, create_using=nx.DiGraph())

    ped_undir_edgelist = nx.read_edgelist(networkx_file, create_using=nx.Graph())
    sub_fams = list(nx.connected_components(ped_undir_edgelist))

    # we will check of the profiles:
    profiles = pd.read_csv(profiles_file, sep='\t')

    if len(np.intersect1d(profiles.columns, ['Sex', 'Birth_Year'])) != 2:
        print("profiles entered do not contain 'Birth_Year' and 'Sex', please check profile file")
        exit(0)

    fam_id = 1

    # Load in this as a pandas dataframe so you can add a header to the text file!
    for sub_family in sub_fams:
        ped_file = []
        # Subset the current family we are looking at
        ped_dir_edgelist.sub_fam_graph = ped_dir_edgelist.subgraph(sub_family)

        for indiv in ped_dir_edgelist.sub_fam_graph:

            print(indiv)
            sex = profiles[profiles['ID'] == int(indiv)]['Sex'].to_list()[0]
            sex = sex_map[sex]

            birthyear = profiles[profiles['ID'] == int(indiv)]['Birth_Year'].to_list()[0]

            cur_node_pred = list(ped_dir_edgelist.sub_fam_graph.predecessors(indiv))

            if len(cur_node_pred) == 2:

                p1_sex = profiles[profiles['ID'] == int(cur_node_pred[0])]['Sex'].to_list()[0]
                sp1_sex = sex_map[p1_sex]

                if sp1_sex == 1:
                    print(p1_sex, int(cur_node_pred[0]), int(cur_node_pred[1]))
                    line = [fam_id, int(indiv), int(cur_node_pred[0]), int(cur_node_pred[1]), sex, -9, birthyear]
                    print(line)
                else:
                    print(p1_sex, int(cur_node_pred[0]), int(cur_node_pred[1]))
                    line = [fam_id, int(indiv), int(cur_node_pred[1]), int(cur_node_pred[0]), sex, -9, birthyear]
                    print(line)

                ped_file.append(line)

            if len(cur_node_pred) == 1:

                p1_sex = profiles[profiles['ID'] == int(cur_node_pred[0])]['Sex'].to_list()[0]
                sp1_sex = sex_map[p1_sex]

                if sp1_sex == 1:
                    line = [fam_id, int(indiv), int(cur_node_pred[0]), 0, sex, -9, birthyear]
                else:
                    line = [fam_id, int(indiv), 0, int(cur_node_pred[0]), sex, -9, birthyear]

                ped_file.append(line)

            if len(cur_node_pred) == 0:
                line = [fam_id, int(indiv), 0, 0, sex,-9, birthyear]
                ped_file.append(line)

        fam_id += 1

    ped_filepath = output_prefix + ".ped"

    column_names=['#Family_ID', "Indiv_ID", "Paternal_ID", "Maternal_ID", "Sex", "Phenotype", "Birth_Year"]

    ped_file = pd.DataFrame(ped_file,columns=column_names)
    ped_file.sort_values(by=['Indiv_ID'], inplace=True)

    ped_file.to_csv(ped_filepath, sep=" ", index=False, header=False)


def convert_networkx_to_ped(networkx_file, output_prefix):
    """
    This function will be used to convert a networkx based family pedigree into a .ped file that is often used for other
    tools such as plink.

    :param networkx_file: File path where the networkx graph is found.
    :param output_prefix: Directory  where you want to output the file in

    :return: will return a ped file of the networkx graph that was inputted in.
    """
    ped_dir_edgelist = nx.read_edgelist(networkx_file, create_using=nx.DiGraph())

    ped_undir_edgelist = nx.read_edgelist(networkx_file, create_using=nx.Graph())
    sub_fams = list(nx.connected_components(ped_undir_edgelist))


    fam_id = 1
    ped_file = []
    # Load in this as a pandas dataframe so you can add a header to the text file!
    for sub_family in sub_fams:
        # Subset the current family we are looking at
        ped_dir_edgelist.sub_fam_graph = ped_dir_edgelist.subgraph(sub_family)

        for indiv in ped_dir_edgelist.sub_fam_graph:


            cur_node_pred=list(ped_dir_edgelist.sub_fam_graph.predecessors(indiv))

            if len(cur_node_pred) == 2:
                line = [fam_id, int(indiv), int(cur_node_pred[0]), int(cur_node_pred[1]), 0,-9]
                ped_file.append(line)

            if len(cur_node_pred) == 1:
                line = [fam_id, int(indiv), int(cur_node_pred[0]), 0, 0, -9]
                ped_file.append(line)

            if len(cur_node_pred) == 0:
                line = [fam_id, int(indiv), 0, 0, 0, -9]
                ped_file.append(line)

        fam_id += 1

    ped_filepath = output_prefix + ".ped"

    column_names=['#Family_ID', "Indiv_ID", "Paternal_ID", "Maternal_ID", "Sex", "Phenotype"]

    ped_file = pd.DataFrame(ped_file,columns=column_names)
    ped_file.sort_values(by=['Indiv_ID'],inplace=True)

    ped_file.to_csv(ped_filepath, sep=" ", index=False, header=False)


def convert_ped_to_networkx(ped_file, output_prefix):
    """


    :param ped_file:
    :param output_dir:
    :return:

    TO DO, add a check to see if the pedigree file is readable, for networkx conversion.
    """
    #

    networkx_pedigree = nx.DiGraph()

    ped_file = open(f"{ped_file}", "r")

    #We ignore the first line since that be the header
    ped_file_lines = ped_file.readlines()[1:]

    for line in ped_file_lines:
        split_line = line.split()
        split_line = list(map(int, split_line))

        if split_line[2] != 0:
            networkx_pedigree.add_edge(split_line[2],split_line[1])
        if split_line[3] != 0:
            networkx_pedigree.add_edge(split_line[3],split_line[1])

    output_filepath = f'{output_prefix}.nx'
    nx.write_edgelist(networkx_pedigree, f"{output_filepath}")

def filter_vcf_for_slim(vcf_filepath):
    '''
        This method will be called to update the user's inputted vcf file, this file will only be called under certain
        conditions of the vcf file.
        1. The vcf file will be filter for any multi-allelic site, only bi-allelic sites are allowed for current
        simulations.
        2. Will remove any empty sites found in the vcf file.
        3. We will update the AA col found in the info column, to match with SLiM's standard upper case A/C/T/G input.
        SliM simulations that require an input vcf MUST have the AA snp specified in the INFO column.

        TO DO
        I am actually right, you only need to update the AA column for nucleotide specific simulations, should I modify
        to add an extra flag if the preprocessing needs to be done?

        03/28/23
        Inputting the AA allele is only needed for nucleotide specific simulations, might just put it on the user to
        make sure the info column has info/AA information.


        filter_vcf_for_slim will only preform steps 1. and 2. to make sure vcf files only contain bi-allelic snps.
        :return:
        self.founder_genomes - assigned to the output founder genome vcf file that is able to be read in by SLiM.

    :param vcf_filepath:
    :return:
    '''

    vcf_prefix = vcf_filepath.split('.')[0]

    shell_cmd = f"bcftools view -m2 -M2 -v snps {vcf_filepath} -O v -o tmp_snps.vcf"
    subprocess.run([shell_cmd], shell=True)

    # This will filter any sites that empty, Minor Allel Count == 0
    shell_cmd = "bcftools filter -e 'MAC == 0' tmp_snps.vcf -O v -o tmp_only_snps.vcf"
    subprocess.run([shell_cmd], shell=True)

    mv_cmd = f'mv tmp_only_snps.vcf {vcf_prefix}_slim_fil.vcf'
    subprocess.run([mv_cmd], shell=True)

    shell_cmd = f'rm tmp_snps.vcf'
    subprocess.run([shell_cmd], shell=True)


def find_founders(networkx_file, shell_output=False):
    """
    This function is used inside ped_sim to check the founder information of a networkx represented family pedigree.
    This will return the number of founder, both implicit and explicit. We additionally will return information on the
    number of desendants (children of founders) found from the simulation, as well as the number of indivs found from
    the simulation. Note this will only print the results out, and should be used to aid the pedigree simulations


    ***Update function to return the number for one output and print the output for "find_founder" -t feature ***
    :param networkx_file: Networkx representation of the family pedigree.
    :return:
    """
    ped_dir_edgelist = nx.read_edgelist(networkx_file, create_using=nx.DiGraph())
    # nodes_num_pred = []

    explicit_founders = []
    implicit_founders =[]
    descendants = []
    for node in ped_dir_edgelist.nodes():
        num_pred = len(list(ped_dir_edgelist.predecessors(node)))
        # nodes_num_pred.append(num_pred)
        if(num_pred==0):
            explicit_founders.append(node)
        elif(num_pred==1):
            implicit_founders.append(node)
        else:
            descendants.append(node)

    if shell_output:
        print(f'Number of explicit founders (known founders): {len(explicit_founders)}')
        #print(explicit_founders)
        print(f'Number of impicit founders (unknown founders): {len(implicit_founders)}')
        #print(implicit_founders)
        print(f'Number of desendants: {len(descendants)}')
        #print(descendants)
        print(f'Total individuals present: {len(ped_dir_edgelist.nodes)}')
        return

    return len(explicit_founders) + len(implicit_founders)