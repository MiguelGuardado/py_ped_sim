#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Recent Update: 02/15/2022
@author: Miguel Guardado Miguel.Guardado@ucsf.edu
"""

import argparse
import os
import subprocess
import numpy as np

from scripts import convert_pedigree
from scripts import util
from scripts import load_founders
from scripts import load_founders_exact

parser = argparse.ArgumentParser()

def load_args():
    """
    This function is used to read in and initialize the user parameters
    entered from the command lines. Some of these parameters are predefined
    for the founders burn in.
    This wont check the parameters but just reads them in. Each of the individual python modules

    USER MUST INPUT THESE PARAMETERS to make any of the simulations run:
        -o output prefix
        -p pedigree file
        -t type of simulation to run.
        -v #if t= load_founder


    :return: returns the users command line inputs for each of the parameters defined
    """
    parser.add_argument('-t', '--type_of_sim', default=0, type=str, required=True)
    parser.add_argument('-v', '--vcf_file', type=str)
    parser.add_argument('-p', '--ped_file', type=str)
    parser.add_argument('-pr', '--profiles_file', type=str)
    parser.add_argument('-p1', '--prob_1', type=float, default=0.01)
    parser.add_argument('-p2', '--prob_2', type=float, default=0.50)
    parser.add_argument('-y', '--years_to_sample', nargs='+', type=int)
    parser.add_argument('-n', '--networkx_file', type=str)
    parser.add_argument('-c', '--census_filepath', type=str)
    parser.add_argument('-f', '--fasta_file', type=str)
    parser.add_argument('-e', '--exact_founder_id', type=str)
    parser.add_argument('-o', '--output_prefix', type=str)
    parser.add_argument('-l', '--genome_length', default=99999, type=int)
    parser.add_argument('-mu', '--mutation_rate', default='1e-7', type=str)
    parser.add_argument('-r', '--recomb_rate', default='1e-8', type=str)
    parser.add_argument('-rm', '--recomb_map', type=str)
    parser.add_argument('-n_gen,', '--number_of_gens', default=12000, type=int)
    parser.add_argument('-n_indiv', '--number_of_indivs', default=1000, type=int)
    parser.add_argument('-s', '--seed_number', default=np.random.randint(100000, size=1)[0], type=int)
    parser.add_argument('-sf', '--sample_file', type=str)
    parser.add_argument('-Ne', '--population_size', type=int, default=10000)
    parser.add_argument('-Nf', '--num_founder', type=int, default=5000)
    parser.add_argument('-n1', '--family1', type=str)
    parser.add_argument('-n2', '--family2', type=str)
    parser.add_argument('-pr1', '--profiles1', type=str)
    parser.add_argument('-pr2', '--profiles2', type=str)
    parser.add_argument('-cf', '--chosen_founder', type = str, default = 'empty')
    parser.add_argument('-cs', '--chosen_sub', type = str, default = 'empty')
    parser.add_argument('-mo', '--main_family_output_prefix', type=str, default = 'main_family')
    parser.add_argument('-mf', '--main_family')
    
    return parser.parse_args()


def raise_filepath_error(filepath):
    raise Exception(f"Filepath not found: '{filepath}'")

def check_exp_input(exp_expression):
    """
    helper function for check_params() to make sure the user entered an exponential number correctly. This will be
    necessary for inputting the parameter into SLiM correctly.

    :return: boolean value to represent if the number inputted was inputted correctly as an exponential.
    """

    try:
        exp_split = exp_expression.split('e')
#       Should crash if there is no e included for base/exp seperation

#       Should crash if the assignment are integer values, but can still be inputted as str for input into SLIM
        base = int(exp_split[0])
        exp = int(exp_split[1])

#       If the split produced more values than the base and exp, then it will return false.
#       and throw error in raise_error()
        if(len(exp_split) == 2):
            return (True)
    except:
        return(False)
    return(False)

def check_and_abs_path(file_path, raise_error=True):
    """Helper function to check file existence and convert to absolute path."""
    if file_path is not None:
        if raise_error and not os.path.isfile(file_path):
            raise_filepath_error(file_path)
        return os.path.abspath(file_path)
    return file_path

def check_output_prefix(required=True, family_broadening=False):
    """Helper function to check if output_prefix is required and convert to absolute path."""
    if required:
        if not args.output_prefix:
            raise ValueError("Output prefix is required but not provided.")
        args.output_prefix = os.path.abspath(args.output_prefix)
    elif args.output_prefix:
        # If optional and provided, convert to absolute path
        args.output_prefix = os.path.abspath(args.output_prefix)

    # This code will evaluate the main family output prefix found in run_family_broadening_full
    if family_broadening:
        if not args.main_family_output_prefix:
            raise ValueError("Output prefix is required but not provided.")
        args.main_family_output_prefix = os.path.abspath(args.main_family_output_prefix)



def check_params():
    """
    This function checks the user parameters and ensures they are correctly formatted.
    It also converts file paths to absolute paths to avoid local directory reference issues.
    """
    if args.type_of_sim == 'sim_founders':
        check_output_prefix()

    elif args.type_of_sim == 'sim_genomes':
        check_output_prefix()
        args.vcf_file = check_and_abs_path(args.vcf_file)
        args.networkx_file = check_and_abs_path(args.networkx_file)
        args.fasta_file = check_and_abs_path(args.fasta_file, raise_error=False)
        args.recomb_map = check_and_abs_path(args.recomb_map, raise_error=False)

    elif args.type_of_sim == 'sim_genomes_exact':
        check_output_prefix()
        args.vcf_file = check_and_abs_path(args.vcf_file)
        args.networkx_file = check_and_abs_path(args.networkx_file)
        args.exact_founder_id = check_and_abs_path(args.exact_founder_id)
        args.fasta_file = check_and_abs_path(args.fasta_file, raise_error=False)
        args.recomb_map = check_and_abs_path(args.recomb_map, raise_error=False)

    elif args.type_of_sim == 'enur_fam':
        # Output prefix is optional in this case
        check_output_prefix(required=False)
        args.networkx_file = check_and_abs_path(args.networkx_file)
        args.sample_file = check_and_abs_path(args.sample_file, raise_error=False)

    elif args.type_of_sim == 'sim_ped':
        check_output_prefix()
        args.census_filepath = check_and_abs_path(args.census_filepath, raise_error=False) or 'scripts/ipumps_sibship_dist.txt'
        args.years_to_sample = ' '.join(str(x) for x in args.years_to_sample)

    elif args.type_of_sim == 'sim_map':
        check_output_prefix()
        args.networkx_file = check_and_abs_path(args.networkx_file)
        args.profiles_file = check_and_abs_path(args.profiles_file)

    elif args.type_of_sim == 'ped_to_networkx':
        check_output_prefix()
        args.ped_file = check_and_abs_path(args.ped_file)

    elif args.type_of_sim == 'networkx_to_ped':
        check_output_prefix()
        args.networkx_file = check_and_abs_path(args.networkx_file)
        args.profiles_file = check_and_abs_path(args.profiles_file, raise_error=False)

    elif args.type_of_sim == 'filter_vcf':
        check_output_prefix(required=False)
        args.vcf_file = check_and_abs_path(args.vcf_file)

    elif args.type_of_sim == 'check_founders':
        check_output_prefix()
        args.networkx_file = check_and_abs_path(args.networkx_file)

    elif args.type_of_sim == 'fill_ped':
        check_output_prefix()
        args.networkx_file = check_and_abs_path(args.networkx_file)

    elif args.type_of_sim == 'convert_pedigree':
        check_output_prefix()
        args.networkx_file = check_and_abs_path(args.networkx_file)

    elif args.type_of_sim == 'run_single_family_broadening':
        check_output_prefix()
        args.networkx_file = check_and_abs_path(args.networkx_file)
        args.family1 = check_and_abs_path(args.family1)
        args.family2 = check_and_abs_path(args.family2)
        args.profile1 = check_and_abs_path(args.profile1)
        args.profile2 = check_and_abs_path(args.profile2)

    elif args.type_of_sim == 'run_full_family_broadening':
        check_output_prefix(required=False)
        args.years_to_sample = ' '.join(str(x) for x in args.years_to_sample)
        args.census_filepath = check_and_abs_path(args.census_filepath, raise_error=False) or 'scripts/ipumps_sibship_dist.txt'
        args.networkx_file = check_and_abs_path(args.networkx_file)
        args.main_family = check_and_abs_path(args.main_family)
        check_output_prefix(args.main_family_output_prefix, family_broadening=True)


#MAIN CLASS: this is where the ped_sim code starts.
if __name__ == '__main__':
    cur_user_dir = os.getcwd()

#   Load in the users arguments
    args = load_args()

#   Check in the user's arguments are valid.
    check_params()

#   We will also set the internal wd to be inside the ped_sim directory
    ped_sim_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(ped_sim_dir)

#   Run feature provided by the user in the -t parameter, will exit if user inputs non supported feature.'
    if args.type_of_sim == 'sim_founders': # Simulate founders genomes
        sim_founder_cmd = f'python scripts/sim_fam_msprime.py -Nf {args.num_founder} -Ne {args.population_size} ' \
              f'-l {args.genome_length} -s {args.seed_number} -mu {args.mutation_rate} -r {args.recomb_rate} ' \
              f'-o {args.output_prefix}'
        subprocess.run([sim_founder_cmd], shell=True)


    elif args.type_of_sim == 'sim_genomes' or args.type_of_sim == 'sim_genomes_exact':
        load_founders.load_founders(
            networkx_file=args.networkx_file,
            cur_dir=cur_user_dir,
            out_pref=args.output_prefix,
            vcf_file=args.vcf_file,
            mutation_rate=args.mutation_rate,
            recomb_rate=args.recomb_rate,
            seed_number=args.seed_number,
            fasta_file=args.fasta_file,
            recomb_map=args.recomb_map,
            exact_founder_id=args.exact_founder_id
        )

    elif args.type_of_sim == 'enur_fam':
        enur_fam_cmd = f'python scripts/enur_fam.py -n {args.networkx_file} -sf {args.sample_file} -o {args.output_prefix}'
        subprocess.run([enur_fam_cmd], shell=True)

    elif args.type_of_sim == 'sim_ped':
        sim_ped_cmd = f'python scripts/sim_pedigree.py -y {args.years_to_sample} -c {args.census_filepath} ' \
                      f'-o {args.output_prefix} -s {args.seed_number}'
        subprocess.run([sim_ped_cmd], shell=True)

    elif args.type_of_sim == 'sim_map':
        sim_ped_cmd = f'python scripts/sim_map.py -n {args.networkx_file} -pr {args.profiles_file} ' \
                      f'-p1 {args.prob_1} -p2 {args.prob_2} -o {args.output_prefix}'
        subprocess.run([sim_ped_cmd], shell=True)

    elif args.type_of_sim == 'ped_to_networkx':
        util.convert_ped_to_networkx(ped_file=args.ped_file, output_prefix=args.output_prefix)

    elif args.type_of_sim == 'networkx_to_ped':
        if args.profiles_file is not None:
            util.convert_networkx_to_ped_wprofiles(networkx_file=args.networkx_file,
                                                   output_prefix=args.output_prefix,
                                                   profiles_file=args.profiles_file)
        else:
            util.convert_networkx_to_ped(networkx_file=args.networkx_file, output_prefix=args.output_prefix)

    elif args.type_of_sim == 'check_founders':
        util.find_founders(networkx_file=args.networkx_file, shell_output=True)

    elif args.type_of_sim == 'filter_vcf':
        util.filter_vcf_for_slim(vcf_file=args.vcf_file, output_prefix=args.output_prefix)

    elif args.type_of_sim == 'fill_ped':
        util.fill_ped(networkx_file=args.networkx_file, output_prefix=args.output_prefix)

    elif args.type_of_sim == 'convert_pedigree':
        convert_pedigree.convert_pedigree(ped_filepath=args.networkx_file, output_prefix=args.output_prefix)

    elif args.type_of_sim == 'run_single_family_broadening':
        fb_cmd = f'python scripts/run_single_family_broadening.py -n1 {args.family1} -n2 {args.family2} -pr1 {args.profiles1} -pr2 {args.profiles2}' \
                      f'-o {args.output_prefix}'
        subprocess.run([fb_cmd], shell = True)

    elif args.type_of_sim == 'run_full_family_broadening':
        fb_cmd = f'python scripts/run_full_family_broadening.py -c {args.census_filepath} -y {args.years_to_sample} ' \
                      f'-mf {args.main_family} -mo {args.main_family_output_prefix} -o {args.output_prefix}'
        subprocess.run([fb_cmd], shell = True)

    else:
        exit('No input entered, please use -t to specify the action you want preformed')