#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Recent Update: 10/23/2021
@author: Miguel Guardado Miguel.Guardado@ucsf.edu
"""

import argparse
import subprocess
import os
import errno
import numpy as np
from scripts import convert_pedigree
from scripts import sim_founders
from scripts import util
from scripts import load_founders
from scripts import load_founders_exact

parser = argparse.ArgumentParser()

def load_args():
    """
    This function is used to read in and initalize the user parameters
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
    parser.add_argument('-t', '--type_of_sim', default=0)
    parser.add_argument('-v', '--vcf_file')
    parser.add_argument('-p', '--ped_file')
    parser.add_argument('-n', '--networkx_file')
    parser.add_argument('-o', '--output_prefix')
    parser.add_argument('-l', '--genome_length', default=99999)
    parser.add_argument('-mu', '--mutation_rate', default='1e-7')
    parser.add_argument('-r', '--recomb_rate', default='1e-8')
    parser.add_argument('-n_gen,', '--number_of_gens', default=12000)
    parser.add_argument('-n_indiv', '--number_of_indivs', default=1000)
    parser.add_argument('-s', '--seed_number', default=1)
    parser.add_argument('-f', '--fasta_file')
    parser.add_argument('-e', '--exact_founder_id')
    parser.add_argument('-full_output', default=False)

    return (parser.parse_args())


def raise_error(flag):
    raise Exception(f"Improper parameter '{flag}' specified, please review docs and fix parameter")

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
    """
    args.output_prefix = os.path.abspath(f"{args.output_prefix}")
#   We need to simulate the founders genomes, so this will check if SLIM input parameters are in the correct format
    if(str(args.type_of_sim) == 'sim_founders'):

        if(not os.path.isfile(str(args.networkx_file))):
            raise_error("-n")
        if (not isinstance(int(args.genome_length), int)):
            raise_error("-l")
        if (not check_exp_input(args.recomb_rate)):
            raise_error("-r")
        if (not check_exp_input(args.mutation_rate)):
            raise_error("-mu")
        if (not isinstance(int(args.number_of_gens), int)):
            raise_error("-N_gen")
        if (not isinstance(int(args.number_of_indivs), int)):
            raise_error("-N_indiv")

        args.networkx_file = os.path.abspath(f"{args.networkx_file}")

    elif(str(args.type_of_sim) == 'load_founders'):
        if (not os.path.isfile(f'{args.vcf_file}')):
            raise_error("-vcf_file")
        if (not os.path.isfile(f'{args.networkx_file}')):
            raise_error("-networkx_file")
        if (not os.path.isfile(f'{args.fasta_file}') and args.fasta_file != None):
            raise_error("-fasta_file")

        if(args.fasta_file == None):
            args.vcf_file = os.path.abspath(f"{args.vcf_file}")
            args.networkx_file = os.path.abspath(f"{args.networkx_file}")
        else:
        #Correct rel path to full path
            args.vcf_file = os.path.abspath(f"{args.vcf_file}")
            args.networkx_file = os.path.abspath(f"{args.networkx_file}")
            args.fasta_file = os.path.abspath(f"{args.fasta_file}")

    elif(str(args.type_of_sim) == 'load_founders_exact'):
        if (not os.path.isfile(f'{args.vcf_file}')):
            raise_error("-vcf_file")
        if (not os.path.isfile(f'{args.networkx_file}')):
            raise_error("-networkx_file")
        if (not os.path.isfile(f'{args.exact_founder_id}')):
            raise_error("-exact_founder_id")
        if (not os.path.isfile(f'{args.fasta_file}') and args.fasta_file != None):
            raise_error("-fasta_file")

        if (args.fasta_file == None):
            args.vcf_file = os.path.abspath(f"{args.vcf_file}")
            args.networkx_file = os.path.abspath(f"{args.networkx_file}")
            args.output_prefix = os.path.abspath(f"{args.output_prefix}")
            args.exact_founder_id = os.path.abspath(f"{args.exact_founder_id}")
        else:
            # Correct rel path to full path
            args.vcf_file = os.path.abspath(f"{args.vcf_file}")
            args.networkx_file = os.path.abspath(f"{args.networkx_file}")
            args.fasta_file = os.path.abspath(f"{args.fasta_file}")
            args.output_prefix = os.path.abspath(f"{args.output_prefix}")
            args.exact_founder_id = os.path.abspath(f"{args.exact_founder_id}")

    elif(str(args.type_of_sim) == 'ped_to_networkx'):
        if (not os.path.isfile(str(args.ped_file))):
            raise_error("-p")
        args.ped_file = os.path.abspath(f"{args.ped_file}")

    elif(str(args.type_of_sim) == 'networkx_to_ped'):
        if (not os.path.isfile(str(args.networkx_file))):
            raise_error("-n")
        args.networkx_file = os.path.abspath(f"{args.networkx_file}")

    elif(str(args.type_of_sim) == 'check_founders'):
        if (not os.path.isfile(str(args.networkx_file))):
            raise_error("-n")
        args.networkx_file = os.path.abspath(f"{args.networkx_file}")

    elif(str(args.type_of_sim) == 'convert_pedigree'):
        if (not os.path.isfile(str(args.networkx_file))):
            raise_error("-n")
        args.networkx_file = os.path.abspath(f"{args.networkx_file}")

    else:
        print(args.type_of_sim)
        raise Exception("Improper type of simulation '-t' specified")

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

#   Determine what function the user wants preformed, currently there are 5 options supported
    if(str(args.type_of_sim) == 'sim_founders'): # Simulate founders genomes

        sim_founders.sim_founders(
            networkx_file=args.networkx_file,
            cur_dir=cur_user_dir,
            out_prefix=args.output_prefix,
            genome_length=args.genome_length,
            mutation_rate=args.mutation_rate,
            recomb_rate=args.recomb_rate,
            seed_number=args.seed_number,
            num_of_indivs=args.number_of_indivs,
            num_of_gens=args.number_of_gens
        )


    elif(str(args.type_of_sim) == 'load_founders'):  # Load vcf file as founders genomes

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

    elif(str(args.type_of_sim) == 'load_founders_exact'):

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

    elif(str(args.type_of_sim) == 'ped_to_networkx'):
        util.convert_ped_to_networkx(ped_file=args.ped_file, output_prefix=args.output_prefix)

    elif(str(args.type_of_sim) == 'networkx_to_ped'):
        util.convert_networkx_to_ped(networkx_file=args.networkx_file, output_prefix=args.output_prefix)

    elif(str(args.type_of_sim) == 'check_founders'):
        util.find_founders(networkx_file=args.networkx_file, shell_output=True)

    elif (str(args.type_of_sim) == 'convert_pedigree'):
        convert_pedigree.convert_pedigree(ped_filepath=args.networkx_file, output_prefix=args.output_prefix)


    else:
        print("Did not give a correct action to take, please check -t parameters")
