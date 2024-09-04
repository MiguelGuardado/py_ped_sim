#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Recent Update: 10/23/2021
@author: Miguel Guardado Miguel.Guardado@ucsf.edu
"""
import pandas as pd
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
                 mutation_rate, recomb_rate, seed_number, fasta_file=None, recomb_map=None, exact_founder_id=None):
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
        self.exact_founder_id = exact_founder_id
        self.recomb_map = recomb_map
        self.is_nuc_seq = (self.fasta_file is not None) ## logical to determine if a nucleotide specific simualtion will be run.
        self.is_recomb_map = (self.recomb_map is not None) ## logical to determine if a recombination map is provided.
        self.is_exact_sim = (self.exact_founder_id is not None)  ## logical to determine if a recombination map is provided.
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
        sample_query_cmd = f"bcftools query -l {self.vcf_file} > {self.output_prefix}_sample_id.txt"
        subprocess.run([sample_query_cmd], shell=True)

        #  Load in generated list into python
        present_sampleid = np.loadtxt(f'{self.output_prefix}_sample_id.txt', dtype='str')
        founder_samples = np.random.choice(present_sampleid, self.num_founder, replace=False)
        np.savetxt(f'{self.output_prefix}_founder_sampleid.txt', founder_samples, fmt='%s')


        self.founder_genomes = f'{self.output_prefix}_founder_genomes.vcf'
        subset_vcf_cmd = f'bcftools view -S {self.output_prefix}_founder_sampleid.txt {self.vcf_file} ' \
                         f'-O v -o {self.founder_genomes}'
        subprocess.run([subset_vcf_cmd], shell=True)

        rm_cmd = f"rm {self.output_prefix}_founder_sampleid.txt {self.output_prefix}_sample_id.txt"
        subprocess.run([rm_cmd], shell=True)

    def extract_founders_exact(self):
        """
               This function is used to assign individuals in the user inputted vcf file to the founders in the inputted
               pedigree. Since founders are based in ascending order of the pedigree id's.

               1) We will sort the exact founder table to have an ordered list of ID's.
               2) Subset those individuals via bcftools, the resulting vcf file will not be ordered based on the user inputted table.

               :return:
               """

        exact_founder_filepath = self.exact_founder_id
        exact_founder = pd.read_csv(exact_founder_filepath, header=None, index_col=None, sep=' ')
        exact_founder.columns = ['vcf_sample_id', 'network_sample_id']
        exact_founder.sort_values(by=['network_sample_id'], inplace=True)

        np.savetxt(f'{self.output_prefix}_sampleid.txt', exact_founder['vcf_sample_id'], fmt="%s")

        subset_vcf_cmd = f'bcftools view -S {self.output_prefix}_sampleid.txt {self.vcf_file} -O v -o ' \
                         f'{self.output_prefix}_founder_genomes.vcf'
        subprocess.run([subset_vcf_cmd], shell=True)

        # vcf file will now point to the updated vcf subset
        self.founder_genomes = f"{self.output_prefix}_founder_genomes.vcf"

        rm_cmd = f"rm {self.output_prefix}_sampleid.txt"
        subprocess.run([rm_cmd], shell=True)

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
            f'{self.implicit_founders_vcf_filepath}',
            stdout=subprocess.PIPE, shell=True).wait()

        rm_cmd = f"rm {implicit_founders_filepath} {explicit_founders_filepath}"
        subprocess.run([rm_cmd], shell=True)

    def find_genome_length(self):
        """
        Internal function inside the class to initalize self.genome_length based the last position found inside the
        vcf file.

        at some point I could round this up to the closest MB, one day......

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

        # This function will identify the individuals in the vcf file to use as founders genomes
        if self.is_exact_sim:
            self.extract_founders_exact()
        else:
            self.extract_founders()

        # This will initialize the genome length from the vcf file.
        self.find_genome_length()

        # This willn call the convert_pedigree object to convert the pedigree into a format SLiM can read to preform genome simualtions.
        ped_converter = convert_pedigree.convert_pedigree(ped_filepath=self.networkx_file,
                                                          output_prefix=self.output_prefix)

        # Double quote syntax is required for SLiM to have command line input parameters. ANNOYING!!!
        self.output_vcf = f"'{self.output_prefix}_genomes.vcf'"
        ped_converter.slim_filepath = f"'{ped_converter.slim_filepath}'"
        ped_converter.founder_filepath = f"'{ped_converter.founder_filepath}'"

        # We have many potenital genomic simulations of pedigrees to run, in the event users want nucleotide specific
        # simulation, implicit vs explicit founder initialization simulations, or if a user inputs a recombination map

        ################################################################################################################
        # Nucleotide specific simulation, with implicit founders, and a recombination map.
        ################################################################################################################
        if self.is_nuc_seq and ped_converter.num_implicit > 0 and self.is_recomb_map:
            print('Nucleotide specific simulation, with implicit founder, and recombination map')

            # Initalize implicit and explicit founders
            self.split_founders(self.founder_genomes, ped_converter.num_explicit, ped_converter.num_implicit)

            # Correct fasta and add notation for SLiM
            self.fasta_file = util.check_fasta(self.fasta_file)
            self.fasta_file = f"'{self.fasta_file}'"
            self.recomb_map = f"'{self.recomb_map}'"

            # Run SLiM Command
            subprocess.run([f'slim -d pedigree_filepath="{ped_converter.slim_filepath}"'
                            f' -d founder_filepath="{ped_converter.founder_filepath}"'
                            f' -d exp_founder_vcf_filepath="{self.explicit_founders_vcf_filepath}"'
                            f' -d imp_founder_vcf_filepath="{self.implicit_founders_vcf_filepath}"'
                            f' -d genome_length="{self.genome_length}"'
                            f' -d mu_rate="{self.mutation_rate}"'
                            f' -d recomb_map="{self.recomb_map}"'
                            f' -s {self.seed_number}'
                            f' -d fasta_file="{self.fasta_file}"'
                            f' -d output_filename="{self.output_vcf}" scripts/simulate_pedigree_wnuc_wrecomb.slim'], shell=True)

        ################################################################################################################
        # Nucleotide specific simulation, no implicit founders, and a recombination map
        ################################################################################################################
        if self.is_nuc_seq and ped_converter.num_implicit == 0 and self.is_recomb_map:
            print('Nucleotide specific simulation, no implicit founder, and recombination map')


            # Correct fasta and add notation for SLiM
            self.fasta_file = util.check_fasta(self.fasta_file)
            self.fasta_file = f"'{self.fasta_file}'"
            self.founder_genomes = f"'{self.founder_genomes}'"
            self.recomb_map = f"'{self.recomb_map}'"

            print(f'slim -d pedigree_filepath="{ped_converter.slim_filepath}"'
                            f' -d founder_filepath="{ped_converter.founder_filepath}"'
                            f' -d exp_founder_vcf_filepath="{self.founder_genomes}"'
                            f' -d genome_length="{self.genome_length}"'
                            f' -d mu_rate="{self.mutation_rate}"'
                            f' -d recomb_map="{self.recomb_map}"'
                            f' -s {self.seed_number}'
                            f' -d fasta_file="{self.fasta_file}"'
                            f' -d output_filename="{self.output_vcf}" scripts/simulate_pedigree_wnuc_wrecomb.slim')

            # Run SLiM Command
            subprocess.run([f'slim -d pedigree_filepath="{ped_converter.slim_filepath}"'
                            f' -d founder_filepath="{ped_converter.founder_filepath}"'
                            f' -d exp_founder_vcf_filepath="{self.founder_genomes}"'
                            f' -d genome_length="{self.genome_length}"'
                            f' -d mu_rate="{self.mutation_rate}"'
                            f' -d recomb_map="{self.recomb_map}"'
                            f' -s {self.seed_number}'
                            f' -d fasta_file="{self.fasta_file}"'
                            f' -d output_filename="{self.output_vcf}" scripts/simulate_pedigree_wnuc_wrecomb.slim'], shell=True)

        ################################################################################################################
        # Non-Nucleotide specific simulation, with implicit founders, and a recombination map
        ################################################################################################################
        if not self.is_nuc_seq and ped_converter.num_implicit > 0 and self.is_recomb_map:
            print('Non-Nucleotide specific simulation, with implicit founder, and recombination map')

            # Idenitfy implicit and explicit founders
            self.split_founders(self.founder_genomes, ped_converter.num_explicit, ped_converter.num_implicit)

            # Get variables ready to input into SLiM
            self.founder_genomes = f"'{self.founder_genomes}'"
            self.recomb_map = f"'{self.recomb_map}'"

            # Run SLiM Command
            subprocess.run([f'slim -d pedigree_filepath="{ped_converter.slim_filepath}"'
                            f' -d founder_filepath="{ped_converter.founder_filepath}"'
                            f' -d exp_founder_vcf_filepath="{self.explicit_founders_vcf_filepath}"'
                            f' -d imp_founder_vcf_filepath="{self.implicit_founders_vcf_filepath}"'
                            f' -d genome_length="{self.genome_length}"'
                            f' -d mu_rate="{self.mutation_rate}"'
                            f' -d recomb_map="{self.recomb_map}"'
                            f' -s {self.seed_number}'
                            f' -d output_filename="{self.output_vcf}" scripts/simulate_pedigree_wrecomb.slim'], shell=True)


        ################################################################################################################
        # Non-Nucleotide specific simulation, no implicit founders, and a recombination map
        ################################################################################################################
        if not self.is_nuc_seq and ped_converter.num_implicit == 0 and self.is_recomb_map:
            print('Non-Nucleotide specific simulation, no implicit founder, and recombination map')


            # Get variables ready to input into SLiM
            self.founder_genomes = f"'{self.founder_genomes}'"
            self.recomb_map = f"'{self.recomb_map}'"

            subprocess.run([f'slim -d pedigree_filepath="{ped_converter.slim_filepath}"'
                            f' -d founder_filepath="{ped_converter.founder_filepath}"'
                            f' -d exp_founder_vcf_filepath="{self.founder_genomes}"'
                            f' -d genome_length="{self.genome_length}"'
                            f' -d mu_rate="{self.mutation_rate}"'
                            f' -d recomb_map="{self.recomb_map}"'
                            f' -s {self.seed_number}'
                            f' -d output_filename="{self.output_vcf}" scripts/simulate_pedigree_wrecomb.slim'],
                           shell=True)

        ################################################################################################################
        # Nucleotide specific simulation, with implicit founders, and no recombination map.
        ################################################################################################################
        if self.is_nuc_seq and ped_converter.num_implicit > 0 and not self.is_recomb_map:
            print('Nucleotide specific simulation, with implicit founder, and no recombination map')

            # Initalize implicit and explicit founders
            self.split_founders(self.founder_genomes, ped_converter.num_explicit, ped_converter.num_implicit)

            self.founder_genomes = f"'{self.founder_genomes}'"

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
        ################################################################################################################
        # Nucleotide specific simulation, no implicit founders, and no recombination map
        ################################################################################################################
        if self.is_nuc_seq and ped_converter.num_implicit == 0 and not self.is_recomb_map:
            print('Nucleotide specific simulation, no implicit founder, and no recombination map')

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
        ################################################################################################################
        # Non-Nucleotide specific simulation, with implicit founders, and no recombination map
        ################################################################################################################
        if not self.is_nuc_seq and ped_converter.num_implicit > 0 and not self.is_recomb_map:
            print('Non-Nucleotide specific simulation, no implicit founder, and no recombination map')
            # Initalize implicit and explicit founders
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

        ################################################################################################################
        # Non-Nucleotide specific simulation, no implicit founders, and no recombination map
        ################################################################################################################
        if not self.is_nuc_seq and ped_converter.num_implicit == 0 and not self.is_recomb_map:
            print('Non-Nucleotide specific simulation, no implicit founder, and recombination map')

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


        # Now that the simulation is complete, this will reindex the vcf files sample ID to match the
        # id found in the family pedigree graph
        util.update_vcf_header(f'{self.output_prefix}_genomes.vcf', self.networkx_file)

        # Add CONTIG header to the VCF file
        util.add_contig(f'{self.output_prefix}_genomes.vcf')

        # Finally we will correct the chr name based on the chr number used in the vcf output
        util.correct_chr_in_vcf(f'{self.output_prefix}_genomes.vcf', self.vcf_file)

        # #  Now that the simulation is done, we will delete all files not desired by user, feel free to undelete these
        # #  if you want these outputs.
        # rm_cmd = f"rm {self.founder_genomes}* {ped_converter.founder_filepath}" \
        #          f" {self.explicit_founders_vcf_filepath} {self.implicit_founders_vcf_filepath}"
        # subprocess.run([rm_cmd], shell=True)

        # if (self.is_nuc_seq):
        #     if (ped_converter.num_implicit > 0):
        #
        #         self.split_founders(self.founder_genomes, ped_converter.num_explicit, ped_converter.num_implicit)
        #
        #         subprocess.run([f'slim -d pedigree_filepath="{ped_converter.slim_filepath}"'
        #                         f' -d founder_filepath="{ped_converter.founder_filepath}"'
        #                         f' -d exp_founder_vcf_filepath="{self.explicit_founders_vcf_filepath}"'
        #                         f' -d imp_founder_vcf_filepath="{self.implicit_founders_vcf_filepath}"'
        #                         f' -d genome_length="{self.genome_length}"'
        #                         f' -d mu_rate="{self.mutation_rate}"'
        #                         f' -d recomb_rate="{self.recomb_rate}"'
        #                         f' -s {self.seed_number}'
        #                         f' -d output_filename="{self.output_vcf}" scripts/simulate_pedigree.slim'], shell=True)
        #     else:
        #
        #         self.founder_genomes = f"'{self.founder_genomes}'"
        #
        #         subprocess.run([f'slim -d pedigree_filepath="{ped_converter.slim_filepath}"'
        #                         f' -d founder_filepath="{ped_converter.founder_filepath}"'
        #                         f' -d exp_founder_vcf_filepath="{self.founder_genomes}"'
        #                         f' -d genome_length="{self.genome_length}"'
        #                         f' -d mu_rate="{self.mutation_rate}"'
        #                         f' -d recomb_rate="{self.recomb_rate}"'
        #                         f' -s {self.seed_number}'
        #                         f' -d output_filename="{self.output_vcf}" scripts/simulate_pedigree.slim'],
        #                        shell=True)
        # else:
        #     self.fasta_file = util.check_fasta(self.fasta_file)
        #     self.fasta_file = f"'{self.fasta_file}'"
        #
        #     if (ped_converter.num_implicit > 0):
        #
        #         self.split_founders(self.founder_genomes, ped_converter.num_explicit, ped_converter.num_implicit)
        #
        #         subprocess.run([f'slim -d pedigree_filepath="{ped_converter.slim_filepath}"'
        #                         f' -d founder_filepath="{ped_converter.founder_filepath}"'
        #                         f' -d exp_founder_vcf_filepath="{self.explicit_founders_vcf_filepath}"'
        #                         f' -d imp_founder_vcf_filepath="{self.implicit_founders_vcf_filepath}"'
        #                         f' -d genome_length="{self.genome_length}"'
        #                         f' -d mu_rate="{self.mutation_rate}"'
        #                         f' -d recomb_rate="{self.recomb_rate}"'
        #                         f' -s {self.seed_number}'
        #                         f' -d fasta_file="{self.fasta_file}"'
        #                         f' -d output_filename="{self.output_vcf}" scripts/simulate_pedigree_wnuc.slim'],
        #                        shell=True)
        #     else:
        #         self.founder_genomes = f"'{self.founder_genomes}'"
        #         subprocess.run([f'slim -d pedigree_filepath="{ped_converter.slim_filepath}"'
        #                         f' -d founder_filepath="{ped_converter.founder_filepath}"'
        #                         f' -d exp_founder_vcf_filepath="{self.founder_genomes}"'
        #                         f' -d genome_length="{self.genome_length}"'
        #                         f' -d mu_rate="{self.mutation_rate}"'
        #                         f' -d recomb_rate="{self.recomb_rate}"'
        #                         f' -s {self.seed_number}'
        #                         f' -d fasta_file="{self.fasta_file}"'
        #                         f' -d output_filename="{self.output_vcf}" scripts/simulate_pedigree_wnuc.slim'],
        #                        shell=True)

