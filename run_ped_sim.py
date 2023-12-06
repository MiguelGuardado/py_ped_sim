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
    parser.add_argument('-c', '--census_filepath', type=str, default='scripts/ipumps_sibship_dist.txt')
    parser.add_argument('-f', '--fasta_file', type=str)
    parser.add_argument('-e', '--exact_founder_id', type=str)
    parser.add_argument('-o', '--output_prefix', type=str)
    parser.add_argument('-l', '--genome_length', default=99999, type=int)
    parser.add_argument('-mu', '--mutation_rate', default='1e-7', type=str)
    parser.add_argument('-r', '--recomb_rate', default='1e-8', type=str)
    parser.add_argument('-n_gen,', '--number_of_gens', default=12000, type=int)
    parser.add_argument('-n_indiv', '--number_of_indivs', default=1000, type=int)
    parser.add_argument('-s', '--seed_number', default=np.random.randint(100000, size=1), type=int)

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

def check_params():
    """
    This function will be used to check in all the user parameters inputted by the user. This function will return
    nothing if the parameters are inputted correctly, and will go onto running the simulations. If the user input
    parameter are off, then this function will throw and error and tell you which parameters needs to be fixed.

    We will additionally find the absolute path of each file prior to running the simulations so we dont have to
    reference the users local directory .


    It might be annoying but we check the input for each simulation type reqested.
    """
    if args.output_prefix is not None:
        args.output_prefix = os.path.abspath(f"{args.output_prefix}")

#   We need to simulate the founders genomes, so this will check if SLIM input parameters are in the correct format
    if args.type_of_sim == 'sim_founders':

        if not os.path.isfile(args.networkx_file):
            raise_filepath_error(args.networkx_file)

        if not check_exp_input(args.recomb_rate):
            raise_filepath_error(args.recomb_rate)

        if not check_exp_input(args.mutation_rate):
            raise_filepath_error(args.mutation_rate)
        #Convert rel path to full path
        args.networkx_file = os.path.abspath(args.networkx_file)

    elif args.type_of_sim == 'load_founders':

        if not os.path.isfile(args.vcf_file):
            raise_filepath_error(args.vcf_file)

        if not os.path.isfile(args.networkx_file):
            raise_filepath_error(args.networkx_file)

        args.vcf_file = os.path.abspath(args.vcf_file)
        args.networkx_file = os.path.abspath(args.networkx_file)

        if args.fasta_file is not None:
            if not os.path.isfile(args.fasta_file):
                raise_filepath_error(args.fasta_file)
            args.fasta_file = os.path.abspath(args.fasta_file)


    elif args.type_of_sim == 'load_founders_exact':
        if not os.path.isfile(args.vcf_file):
            raise_filepath_error(args.vcf_file)
        if not os.path.isfile(args.networkx_file):
            raise_filepath_error(args.networkx_file)
        if not os.path.isfile(args.exact_founder_id):
            raise_filepath_error(args.exact_founder_id)

        # Correct rel path to full path
        args.vcf_file = os.path.abspath(f"{args.vcf_file}")
        args.networkx_file = os.path.abspath(f"{args.networkx_file}")
        args.output_prefix = os.path.abspath(f"{args.output_prefix}")
        args.exact_founder_id = os.path.abspath(f"{args.exact_founder_id}")
        if args.fasta_file is not None:
            if not os.path.isfile(args.fasta_file):
                raise_filepath_error("-fasta_file")
            args.fasta_file = os.path.abspath(args.fasta_file)

    elif args.type_of_sim == 'ped_to_networkx':
        if not os.path.isfile(args.ped_file):
            raise_filepath_error(args.ped_file)
        args.ped_file = os.path.abspath(f"{args.ped_file}")

    elif args.type_of_sim == 'networkx_to_ped':
        if not os.path.isfile(args.networkx_file):
            raise_filepath_error(args.networkx_file)

        if args.profiles_file is not None and os.path.isfile(args.profiles_file):
            args.profiles_file = os.path.abspath(f"{args.profiles_file}")
        args.networkx_file = os.path.abspath(f"{args.networkx_file}")

    elif args.type_of_sim == 'filter_vcf':
        if not os.path.isfile(args.vcf_file):
            raise_filepath_error(args.vcf_file)
        args.vcf_file = os.path.abspath(args.vcf_file)

    elif args.type_of_sim == 'check_founders':
        if not os.path.isfile(args.networkx_file):
            raise_filepath_error(args.networkx_file)
        args.networkx_file = os.path.abspath(args.networkx_file)

    elif args.type_of_sim == 'fill_ped':
        if not os.path.isfile(args.networkx_file):
            raise_filepath_error(args.networkx_file)
        args.networkx_file = os.path.abspath(args.networkx_file)

    elif args.type_of_sim == 'convert_pedigree':
        if not os.path.isfile(args.networkx_file):
            raise_filepath_error(args.networkx_file)
        args.networkx_file = os.path.abspath(args.networkx_file)

    elif args.type_of_sim == 'enur_fam':
        if not os.path.isfile(args.networkx_file):
            raise_filepath_error(args.networkx_file)
        args.networkx_file = os.path.abspath(args.networkx_file)

    elif args.type_of_sim == 'sim_ped':
        if not os.path.isfile(args.census_filepath):
            raise_filepath_error(args.census_filepath)
        args.census_filepath = os.path.abspath(args.census_filepath)
        args.years_to_sample = ' '.join(str(x) for x in args.years_to_sample)

    elif args.type_of_sim == 'sim_map':
        if not os.path.isfile(args.networkx_file):
            raise_filepath_error(args.networkx_file)

        if not os.path.isfile(args.profiles_file):
            raise_filepath_error(args.profiles_file)

        args.networkx_file = os.path.abspath(args.networkx_file)
        args.profiles_file = os.path.abspath(args.profiles_file)

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
        exit('sim founder under development, dont look here!!')


    elif args.type_of_sim == 'load_founders':

        load_founders.load_founders(
            networkx_file=args.networkx_file,
            cur_dir=cur_user_dir,
            out_pref=args.output_prefix,
            vcf_file=args.vcf_file,
            mutation_rate=args.mutation_rate,
            recomb_rate=args.recomb_rate,
            seed_number=args.seed_number,
            fasta_file=args.fasta_file
        )

    elif args.type_of_sim == 'load_founders_exact':

        load_founders_exact.load_founders_exact(
            networkx_file=args.networkx_file,
            exact_founder_id=args.exact_founder_id,
            cur_dir=cur_user_dir,
            out_pref=args.output_prefix,
            vcf_file=args.vcf_file,
            mutation_rate=args.mutation_rate,
            recomb_rate=args.recomb_rate,
            seed_number=args.seed_number,
            fasta_file=args.fasta_file
        )

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
        util.filter_vcf_for_slim(vcf_file=args.vcf_file)

    elif args.type_of_sim == 'fill_ped':
        util.fill_ped(networkx_file=args.networkx_file, output_prefix=args.output_prefix)

    elif args.type_of_sim == 'convert_pedigree':
        convert_pedigree.convert_pedigree(ped_filepath=args.networkx_file, output_prefix=args.output_prefix)

    elif args.type_of_sim == 'enur_fam':
        if args.output_prefix is None:
            enur_fam_cmd = f'python scripts/enur_fam.py -n {args.networkx_file}'
        else:
            enur_fam_cmd = f'python scripts/enur_fam.py -n {args.networkx_file} -o {args.output_prefix}'
        subprocess.run([enur_fam_cmd], shell=True)

    elif args.type_of_sim == 'sim_ped':
        sim_ped_cmd = f'python scripts/sim_pedigree.py -y {args.years_to_sample} -c {args.census_filepath} ' \
                      f'-o {args.output_prefix} -s {args.seed_number}'
        subprocess.run([sim_ped_cmd], shell=True)

    elif args.type_of_sim == 'sim_map':
        sim_ped_cmd = f'python scripts/sim_map.py -n {args.networkx_file} -pr {args.profiles_file} ' \
                      f'-p1 {args.prob_1} -p2 {args.prob_2} -o {args.output_prefix}'
        subprocess.run([sim_ped_cmd], shell=True)

    else:
        exit('No input entered, please use -t to specify the action you want preformed')