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
            sex = profiles[profiles['ID'] == int(indiv)]['Sex'].to_list()[0]
            sex = sex_map[sex]

            birthyear = profiles[profiles['ID'] == int(indiv)]['Birth_Year'].to_list()[0]

            cur_node_pred=list(ped_dir_edgelist.sub_fam_graph.predecessors(indiv))
            if len(cur_node_pred)==2:
                line = [fam_id, int(indiv),int(cur_node_pred[0]), int(cur_node_pred[1]), sex,-9, birthyear]
                ped_file.append(line)

            if len(cur_node_pred)==1:
                line = [fam_id, int(indiv), int(cur_node_pred[0]), 0, sex,-9, birthyear]
                ped_file.append(line)

            if len(cur_node_pred)==0:
                line = [fam_id, int(indiv), 0, 0, sex,-9, birthyear]
                ped_file.append(line)

        fam_id+=1

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
                line = [fam_id,int(indiv),int(cur_node_pred[0]),int(cur_node_pred[1]), 0,-9]
                ped_file.append(line)

            if len(cur_node_pred) == 1:
                line = [fam_id,int(indiv),int(cur_node_pred[0]),0,0,-9]
                ped_file.append(line)

            if len(cur_node_pred) == 0:
                line = [fam_id,int(indiv), 0, 0, 0,-9]
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
        print(explicit_founders)
        print(f'Number of impicit founders (unknown founders): {len(implicit_founders)}')
        print(implicit_founders)
        print(f'Number of desendants: {len(descendants)}')
        print(descendants)
        print(f'Total individuals present: {len(ped_dir_edgelist.nodes)}')
        return

    return(len(explicit_founders)+len(implicit_founders))