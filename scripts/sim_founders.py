#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Recent Update: 10/23/2021
@author: Miguel Guardado Miguel.Guardado@ucsf.edu
"""

import subprocess
import numpy as np
from scripts import convert_pedigree
from scripts import util

class sim_founders:
    """
    This function is used to run the pedigree simulations in where we simulate the founders genomes. This function will...
    1. convert and sort the family pedigree in a way that SLiM will be able to read in.
    2. Run the populations burn in to add genetic variance onto the founder.
    3. Seperate the vcf file into the number of implicit and explicit founders to input for the pedigree simulation.
    4. Run the pedigree simulations where we initialize the founders genomes from the previous simulation(step 2)


    :return:
    This will results in a vcf file that will hold all the genetic variants(snps) found for the individuals in the pedigree
    {output_prefix}_founder_file.txt #Slim converted pedigree for simulation
    {output_prefix}_genomes.vcf #Genotype file for the simulation
    """

    def __init__(self, networkx_file, cur_dir, out_prefix, genome_length, mutation_rate, recomb_rate,
                 seed_number, num_of_indivs, num_of_gens):
        self.networkx_file = networkx_file
        self.cur_dir = cur_dir
        self.output_prefix = out_prefix
        self.num_of_indivs = num_of_indivs
        self.num_of_gens = num_of_gens
        self.genome_length = genome_length
        self.mutation_rate = mutation_rate
        self.recomb_rate = recomb_rate
        self.seed_number = seed_number
        self.sim_founders_genomes_filepath = ""
        self.output_vcf = ""
        self.run_simulation()

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
            stdout=subprocess.PIPE, shell=True)

        subprocess.Popen(
            f'bcftools view -S {implicit_founders_filepath} {founder_vcf_filepath} > '
            f'{self.implicit_founders_vcf_filepath}',stdout=subprocess.PIPE, shell=True)


    def run_simulation(self):
        """
        main function for the class. Where all the magic happens.
        """
        #Using a networkx inputted family pedigree, we will input this genome into
        ped_converter = convert_pedigree.convert_pedigree(ped_filepath=self.networkx_file,
                                                          output_prefix=self.output_prefix)

        #For input into slim via command line, we must "'double quote'" are variable names
        founder_filepath = f"'{ped_converter.founder_filepath}'"
        self.sim_founders_genomes_filepath = f"'{self.output_prefix}_founders_genomes.vcf'"
        ped_converter.slim_filepath = f"'{ped_converter.slim_filepath}'"
        ped_converter.founder_filepath = f"'{ped_converter.founder_filepath}'"


        # This will call the slim script to simulate the founders genomes for the main pedigree simulation.
        subprocess.run([f'slim -d file_path="{founder_filepath}"'
                        f' -d n_indiv="{self.num_of_indivs}"'
                        f' -d genome_length="{self.genome_length}"'
                        f' -d gen_stop="{self.num_of_gens}"'
                        f' -d mu_rate="{self.mutation_rate}"'
                        f' -d recomb_rate="{self.recomb_rate}"'
                        f' -s {int(self.seed_number)}'
                        f' -d output_filename="{self.sim_founders_genomes_filepath}" scripts/simulate_founders.slim &> /dev/null'],
                       shell=True)
        self.output_vcf = f"'{self.output_prefix}_genomes.vcf'"

        if (ped_converter.num_implicit > 0):


            self.split_founders(f"{self.output_prefix}_founders_genomes.vcf",
                                ped_converter.num_explicit, ped_converter.num_implicit)


            #       this will call upon the slim script to take in the vcf files of founders + number of individuals at home.
            subprocess.run([f'slim -d pedigree_filepath="{ped_converter.slim_filepath}"'
                            f' -d founder_filepath="{ped_converter.founder_filepath}"'
                            f' -d exp_founder_vcf_filepath="{self.explicit_founders_vcf_filepath}"'
                            f' -d imp_founder_vcf_filepath="{self.implicit_founders_vcf_filepath}"'
                            f' -d genome_length="{self.genome_length}"'
                            f' -d mu_rate="{self.mutation_rate}"'
                            f' -d recomb_rate="{self.recomb_rate}"'
                            f' -s {self.seed_number}'
                            f' -d output_filename="{self.output_vcf}" scripts/simulate_pedigree.slim &> /dev/null'], shell=True)

            #
            rm_cmd = f'rm {ped_converter.slim_filepath} {self.explicit_founders_vcf_filepath} ' \
                     f'{self.implicit_founders_vcf_filepath} {self.output_prefix}_explicit_founders.txt ' \
                     f'{self.output_prefix}_implicit_founders.txt {self.output_prefix}_founders_genomes.vcf*'
            subprocess.run([rm_cmd], shell=True)

        else:

            subprocess.run([f'slim -d pedigree_filepath="{ped_converter.slim_filepath}"'
                            f' -d founder_filepath="{ped_converter.founder_filepath}"'
                            f' -d exp_founder_vcf_filepath="{self.sim_founders_genomes_filepath}"'
                            f' -d genome_length="{self.genome_length}"'
                            f' -d mu_rate="{self.mutation_rate}"'
                            f' -d recomb_rate="{self.recomb_rate}"'
                            f' -s {self.seed_number}'
                            f' -d output_filename="{self.output_vcf}" scripts/simulate_pedigree.slim &> /dev/null'], shell=True)

            rm_cmd = f'rm {ped_converter.slim_filepath} {self.sim_founders_genomes_filepath} '
            subprocess.run([rm_cmd], shell=True)

        util.update_vcf_header(self.output_vcf, self.networkx_file)
'''
list of files that get deleted, and the variable/process they are attached to.
{output_prefix}_explicit_founders.txt - random sample for explicit founders
{output_prefix}_explicit_founders.vcf - vcf file input from explicit founders
{output_prefix}_implicit_founders.txt - random sample for implicit founders
{output_prefix}_implicit_founders.vcf - vcf file input from implicit founders


{output_prefix}_founder_genomes.vcf - simulated founder genomes from simulate_founders.slim lines 85-93
{output_prefix}_explicit_founders.vcf.gz -compressed founder genomes for bcftools manipulation
{output_prefix}_explicit_founders.txt.gz.csi -compressed index founder genomes for bcftools manipulation
 

{output_prefix}_genomes.vcf -output of the simulation, not deleted.
{output_prefix}_founder_file.txt -output of the simulation, not deleted.
'''
