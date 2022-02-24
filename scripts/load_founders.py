#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Recent Update: 10/23/2021
@author: Miguel Guardado Miguel.Guardado@ucsf.edu
"""

import numpy as np
import subprocess
from scripts import convert_pedigree
from scripts import util


class load_founders:

    """
    This function is used to run the pedigree simulations in where the user provides a vcf file to load the founders
    genetic information. This will allow us to not simulate the founders genomes and go directly to the main simulation.
    This simulation will initialize the founders randomly from the inputted vcf file.

    Additionally, if the positions of the variants in your simulation matter, then you can provide a fasta file to have
    as the reference genomes prior to when the simulations start.

    This function will....
    1) Convert and sort the family pedigree in a way that SLiM will be able to read in.
    2) Run the pedigree simulations where we initialize the genetic information based off the user inputted vcf file.

    :return:
    This will results in a vcf file that will hold all the genetic variants(snps) found for the individuals in the pedigree
    {output_prefix}_founder_file.txt #Slim converted pedigree for simulation
    {output_prefix}_genomes.vcf #Genotype file for the simulation

    """
    def __init__(self, networkx_file, cur_dir, out_pref, vcf_file,
                 mutation_rate, recomb_rate, seed_number, fasta_file=None):
        self.networkx_file = networkx_file
        self.num_founder = int(util.find_founders(self.networkx_file))
        self.cur_dir = cur_dir
        self.output_prefix = out_pref
        self.vcf_file = vcf_file
        self.founder_genomes = ''
        self.implicit_founders_vcf_filepath = ''
        self.explicit_founders_vcf_filepath = ''
        self.genome_length = 0
        self.mutation_rate = mutation_rate
        self.recomb_rate = recomb_rate
        self.seed_number = seed_number
        self.fasta_file = fasta_file
        self.is_nuc_seq = (self.fasta_file is None)
        self.run_simulation()



    def check_vcf(self):

        find_length_cmd = f"bcftools query -l {self.vcf_file} | wc -l"
        process = subprocess.Popen([find_length_cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        vcf_indiv_length, err = process.communicate()
        vcf_indiv_length = int(vcf_indiv_length.decode('ascii'))

        if vcf_indiv_length < self.num_founder:
            print("Vcf file does not have enough individuals inside to run the simulations, please check your vcf file")
            print(f"Number of individuals inside vcf file: {vcf_indiv_length}")
            print(f"Number of founders inside family pedigree: {self.num_founder}")
            exit(0)



    def extract_founders(self):
        """
        Internal function, this will extract the list of sample id to be used for the founders genome,
        load_founder will assign the founders randomly across the vcf file.
        :return:
        self.founder_genomes - initally assigned to the sample subset of the vcf file that is inputted by the user.
        """
        sample_query_cmd = f"bcftools query -l {self.vcf_file} > user_vcf_sample_id.txt"
        subprocess.run([sample_query_cmd], shell=True)

        #  Load in generated list into python
        present_sampleid = np.loadtxt('user_vcf_sample_id.txt', dtype='str')
        founder_samples = np.random.choice(present_sampleid, self.num_founder, replace=False)
        np.savetxt('founder_sampleid.txt', founder_samples, fmt='%s')
        #

        self.founder_genomes = f'{self.output_prefix}_founder_genomes.vcf'
        subset_vcf_cmd = f'bcftools view -S founder_sampleid.txt {self.vcf_file} -O v -o {self.founder_genomes}'
        subprocess.run([subset_vcf_cmd], shell=True)

        rm_cmd = "rm user_vcf_sample_id.txt founder_sampleid.txt"
        subprocess.run([rm_cmd], shell=True)

    def filter_vcf_for_slim(self):
        """
        This method will be called to update the user's inputted vcf file, this file will only be called under certain
        conditions of the vcf file.
        1. The vcf file will be filter for any multi-allelic site, only bi-allelic sites are allowed for current
        simulations.
        2. Will remove any empty sites found in the vcf file.
        3. We will update the AA col found in the info column, to match with SLiM's standard upper case A/C/T/G input
        :return:
        self.founder_genomes - assigned to the output founder genome vcf file that is able to be read in by SLiM.
        """
        #This will filter any sites that are not bi-allelic
        shell_cmd = f"bcftools view -m2 -M2 -v snps {self.founder_genomes} -O v -o tmp_snps.vcf"
        subprocess.run([shell_cmd], shell=True)

        #This will filter any sites that empty, Minor Allel Count == 0
        shell_cmd = "bcftools filter -e 'MAC == 0' tmp_snps.vcf -O v -o tmp_only_snps.vcf"
        subprocess.run([shell_cmd], shell=True)

        #Extract a list of each snps infomration for ancestral allele info correction
        shell_cmd = "bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%REF\n' tmp_only_snps.vcf | bgzip -c > annot.txt.gz"
        subprocess.run([shell_cmd], shell=True)

        shell_cmd = "tabix -s1 -b2 -e2 annot.txt.gz"
        subprocess.run([shell_cmd], shell=True)

        self.founder_genomes = f"{self.output_prefix}_slim_corrected.vcf"

        shell_cmd = "bcftools annotate -a annot.txt.gz -c CHROM,POS,REF,ALT,INFO/AA tmp_only_snps.vcf -O v -o " \
                    f"{self.output_prefix}_slim_corrected.vcf"
        subprocess.run([shell_cmd], shell=True)

        shell_cmd = "rm tmp_snps.vcf tmp_only_snps.vcf annot.txt.gz annot.txt.gz.tbi"
        subprocess.run([shell_cmd], shell=True)

    def split_founders(self, founder_vcf_filepath, num_explicit, num_implicit):
        """
        This function split the founder vcf file into implicit vs explict vcf files. This function will only be called
        if implicit founders are inputted.

        founder_vcf_filepath(str) - filepath to where the vcf file is found, containing all the explicit and implicit founders
        num_explicit(int) - number of explicit founders
        num_implicit(int) - number of implicit founders
        :return:

        self.explicit_founders_vcf_filepath
        self.implicit_founders_vcf_filepath

        Note
        {founder_vcf_filepath}.gz and .csi files are created from bcftools index/compression input requirement, should
        delete outside the function. For some reason if you delete them at any point before the end of the simulation
        my code will break, even though all forms of founder_vcf_filepath is not used at any point downsteam in the
        implicit route of the code. So I am not sure why this happens but should figure out one day..... but not today.

        """
        #Compress and index the vcf file so we can use bcftools downsteam to seperate implicit/explicit founder
        subprocess.run([f'bcftools view {founder_vcf_filepath} -O z -o {founder_vcf_filepath}.gz'], shell=True)
        subprocess.run([f'bcftools index {founder_vcf_filepath}.gz'], shell=True)
        founder_vcf_filepath = f'{founder_vcf_filepath}.gz'


        #This will return the id list of every founder simulated inside the genetic model, based of compressed vcf file
        proc = subprocess.Popen(f'bcftools query -l {founder_vcf_filepath} ', stdout=subprocess.PIPE, shell=True)
        sample_id = proc.communicate()[0].decode("utf-8").split("\n")
        sample_id.pop()


        #Randomly sample explicit and founder id, based the count inputted by num_explicit and num_implicit
        explicit_founders_id = np.random.choice(sample_id, num_explicit, replace=False)

        sample_id = list(set(sample_id) - set(explicit_founders_id))

        implicit_founders_id = np.random.choice(sample_id, num_implicit, replace=False)

        explicit_founders_filepath = self.output_prefix + '_explicit_founders.txt'
        implicit_founders_filepath = self.output_prefix + '_implicit_founders.txt'

        np.savetxt(explicit_founders_filepath, explicit_founders_id, fmt='%s')
        np.savetxt(implicit_founders_filepath, implicit_founders_id, fmt='%s')

        self.explicit_founders_vcf_filepath = f"'{self.output_prefix}_explicit_founders.vcf'"
        self.implicit_founders_vcf_filepath = f"'{self.output_prefix}_implicit_founders.vcf'"


        subprocess.Popen(
            f'bcftools view -S {explicit_founders_filepath} {founder_vcf_filepath} > '
            f'{self.explicit_founders_vcf_filepath}',
            stdout=subprocess.PIPE, shell=True).wait()

        subprocess.Popen(
            f'bcftools view -S {implicit_founders_filepath} {founder_vcf_filepath} > '
            f'{self.implicit_founders_vcf_filepath}',stdout=subprocess.PIPE, shell=True).wait()

        rm_cmd = f"rm {implicit_founders_id} {explicit_founders_filepath}"
        subprocess.run([rm_cmd], shell=True)



    def find_genome_length(self):
        """
        Internal function inside the class to initalize self.genome_length based the last position found inside the
        vcf file.

        :return: self.genome_length
        """

        pos_query_cmd = f"bcftools view {self.founder_genomes} | tail -n 1"
        process = subprocess.Popen([pos_query_cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

        vcf_indiv_line, err = process.communicate()
        vcf_indiv_line = vcf_indiv_line.decode('ascii').split("\t")
        self.genome_length = int(vcf_indiv_line[1])





    def run_simulation(self):
        """
        Main function for the class. Where all the magic happens.
        """
        #  This will check the user inputted vcf file.
        self.check_vcf()

        self.extract_founders()

        self.filter_vcf_for_slim()

        self.find_genome_length()
        #call on convert_pedigree object to convert networkx pedigree into a SLiM readable readable
        ped_converter = convert_pedigree.convert_pedigree(ped_filepath=self.networkx_file,
                                                          output_prefix=self.output_prefix)

#       double quote syntax is required for SLiM to have command line input parameters
        self.output_vcf = f"'{self.output_prefix}_genomes.vcf'"
        ped_converter.slim_filepath = f"'{ped_converter.slim_filepath}'"
        ped_converter.founder_filepath = f"'{ped_converter.founder_filepath}'"


#       If the user provides a fasta file input, then we will run the slim script that is nucleotide specific.
#       for all other uses, where the nucleotide positions are not desired, we will use a simulations that is not
#       nucleotide sepecific.

        if (self.is_nuc_seq):
            if (ped_converter.num_implicit > 0):

                self.split_founders(self.founder_genomes, ped_converter.num_explicit, ped_converter.num_implicit)
                subprocess.run([f'slim -d pedigree_filepath="{ped_converter.slim_filepath}"'
                                f' -d founder_filepath="{ped_converter.founder_filepath}"'
                                f' -d exp_founder_vcf_filepath="{self.explicit_founders_vcf_filepath}"'
                                f' -d imp_founder_vcf_filepath="{self.implicit_founders_vcf_filepath}"'
                                f' -d genome_length="{self.genome_length}"'
                                f' -d mu_rate="{self.mutation_rate}"'
                                f' -d recomb_rate="{self.recomb_rate}"'
                                f' -s {self.seed_number}'
                                f' -d output_filename="{self.output_vcf}" scripts/simulate_pedigree.slim'], shell=True)
            else:

                self.founder_genomes = f"'{self.founder_genomes}'"
                subprocess.run([f'slim -d pedigree_filepath="{ped_converter.slim_filepath}"'
                                f' -d founder_filepath="{ped_converter.founder_filepath}"'
                                f' -d exp_founder_vcf_filepath="{self.founder_genomes}"'
                                f' -d genome_length="{self.genome_length}"'
                                f' -d mu_rate="{self.mutation_rate}"'
                                f' -d recomb_rate="{self.recomb_rate}"'
                                f' -s {self.seed_number}'
                                f' -d output_filename="{self.output_vcf}" scripts/simulate_pedigree.slim'],
                               shell=True)
        else:
            self.fasta_file = f"'{self.fasta_file}'"
            if (ped_converter.num_implicit > 0):
                self.split_founders(self.founder_genomes, ped_converter.num_explicit, ped_converter.num_implicit)
                subprocess.run([f'slim -d pedigree_filepath="{ped_converter.slim_filepath}"'
                                f' -d founder_filepath="{ped_converter.founder_filepath}"'
                                f' -d exp_founder_vcf_filepath="{self.explicit_founders_vcf_filepath}"'
                                f' -d imp_founder_vcf_filepath="{self.implicit_founders_vcf_filepath}"'
                                f' -d genome_length="{self.genome_length}"'
                                f' -d mu_rate="{self.mutation_rate}"'
                                f' -d recomb_rate="{self.recomb_rate}"'
                                f' -s {self.seed_number}'
                                f' -d fasta_file="{self.fasta_file}"'
                                f' -d output_filename="{self.output_vcf}" scripts/simulate_pedigree_wnuc.slim'], shell=True)
            else:
                self.founder_genomes = f"'{self.founder_genomes}'"
                subprocess.run([f'slim -d pedigree_filepath="{ped_converter.slim_filepath}"'
                                f' -d founder_filepath="{ped_converter.founder_filepath}"'
                                f' -d exp_founder_vcf_filepath="{self.founder_genomes}"'
                                f' -d genome_length="{self.genome_length}"'
                                f' -d mu_rate="{self.mutation_rate}"'
                                f' -d recomb_rate="{self.recomb_rate}"'
                                f' -s {self.seed_number}'
                                f' -d fasta_file="{self.fasta_file}"'
                                f' -d output_filename="{self.output_vcf}" scripts/simulate_pedigree_wnuc.slim'],
                               shell=True)
        # Now that the simulation is complete, this will reindex the vcf files sample ID to match the
        # id found in the family pedigree graph
        util.update_vcf_header(self.output_vcf, self.networkx_file)


        #  Now that the simulation is done, we will delete all files not desired by user, feel free to undelete these
        # if you want these outputs.

        rm_cmd = f"rm {self.founder_genomes}* {ped_converter.founder_filepath} {self.output_prefix}_founder_genomes.vcf" \
                 f" {self.explicit_founders_vcf_filepath} {self.implicit_founders_vcf_filepath} &> /dev/null"
        subprocess.run([rm_cmd], shell=True)
