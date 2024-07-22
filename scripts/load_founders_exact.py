#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Recent Update: 10/23/2021
@author: Miguel Guardado Miguel.Guardado@ucsf.edu
"""

import numpy as np
import pandas as pd
import subprocess
from scripts import convert_pedigree
from scripts import util

class load_founders_exact:

    """
    This function is used to run the pedigree simulations in where the user provides a vcf file to load the founders
    genetic information. This will allow us to not simulate the founders genomes and go directly to the main simulation.
    This simulation will initialize the founders using an exact id text file inputted by the user.

    Additionally, if the positions of the variants in your simulation matter, then you can provide a fasta file to have
    as the reference genomes prior to when the simulations start.

    This function will....
    1) Convert and sort the family pedigree in a way that SLiM will be able to read in.
    2) Run the pedigree simulations where we initialize the genetic information based off the user inputted vcf file.

    :return:
    This will results in a vcf file that will hold all the genetic variants(snps) found for the individuals in the pedigree
    {output_prefix}_founder_file.txt #Slim converted pedigree for simulation
    {output_prefix}_genomes.vcf #Genotype file for the simulation


    TO DO:
    Should check on this function to make sure the exact_founder_id list has correct values for each column


    """
    def __init__(self, networkx_file, cur_dir, out_pref, vcf_file,
                 mutation_rate, recomb_rate, seed_number, exact_founder_id, fasta_file=None):
        self.networkx_file = networkx_file
        self.cur_dir = cur_dir
        self.output_prefix = out_pref
        self.vcf_file = vcf_file
        self.output_vcf = 0
        self.founder_genomes = ''
        self.exact_founder_id = exact_founder_id
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
        """

        :return:
        """

        founder_length = int(util.find_founders(self.networkx_file))
        find_length_cmd = f"bcftools query -l {self.vcf_file} | wc -l"
        process = subprocess.Popen([find_length_cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        vcf_indiv_length, err = process.communicate()
        vcf_indiv_length = int(vcf_indiv_length.decode('ascii'))

        if vcf_indiv_length < founder_length:
            print("Vcf file does not have enough individuals inside to run the simulations, please check your vcf file")
            print(f"Number of individuals inside vcf file: {vcf_indiv_length}")
            print(f"Number of founders inside family pedigree: {founder_length}")
            exit(0)

    def extract_founders(self):
        """
        This function is used to assign individuals in the user inputted vcf file to the founders in the inputted
        pedigree. Since founders are based in ascending order of the pedigree id's.

        1) We will sort the exact founder table to have an ordered list of ID's.
        2) Subset those individuals via bcftools, the resulting vcf file will not be ordered based on the user inputted table.

        :return:
        """

        exact_founder_filepath = self.exact_founder_id
        exact_founder = pd.read_csv(exact_founder_filepath, header=None, index_col=None, sep=' ')
        exact_founder.columns = ['vcf_sample_id','network_sample_id']
        exact_founder.sort_values(by=['network_sample_id'], inplace=True)

        np.savetxt(f'{self.output_prefix}_sampleid.txt', exact_founder['vcf_sample_id'], fmt="%s")

        subset_vcf_cmd = f'bcftools view -S {self.output_prefix}_sampleid.txt {self.vcf_file} -O v -o ' \
                         f'{self.output_prefix}_founder_genomes.vcf'
        subprocess.run([subset_vcf_cmd], shell=True)

        #vcf file will now point to the updated vcf subset
        self.founder_genomes = f"{self.output_prefix}_founder_genomes.vcf"

        rm_cmd = f"rm {self.output_prefix}_sampleid.txt"
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
        main function for the class. Where all the magic happens.
        """
        # Check vcf file if SLiM readable
        self.check_vcf()

        # This function will idenitfy the founders in the vcf file and align them to user input.
        self.extract_founders()

        # Next we find the genome length.
        self.find_genome_length()

        # Call on convert_pedigree object to convert networkx pedigree into a something SLiM can read.
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
                                f' -d output_filename="{self.output_vcf}" scripts/simulate_pedigree.slim &> /dev/null'], shell=True)
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
        util.update_vcf_header(f'{self.output_prefix}_genomes.vcf', self.networkx_file)
        #Now that the simulation is done, we will delete all files not desired by user, feel free to undelete these
        # if you want these outputs.

        rm_cmd = f"rm {self.founder_genomes}* {ped_converter.founder_filepath} {self.output_prefix}_founder_genomes.vcf" \
                 f" {self.explicit_founders_vcf_filepath} {self.implicit_founders_vcf_filepath} &> /dev/null"
        subprocess.run([rm_cmd], shell=True)
