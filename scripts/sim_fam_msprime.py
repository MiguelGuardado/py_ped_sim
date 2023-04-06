import msprime
import argparse
import numpy.random as npr


'''


'''

def load_args():
    """
    This function is used to read in and initalize the user parameters
    entered from the command lines. Some of these parameters are predefined for the founders burn in.
    This wont check the parameters but just reads them in. Each of the individual python modules

    :return: returns the users command line inputs for each of the parameters defined


    """
    parser = argparse.ArgumentParser()

    parser.add_argument('-n', '--n_founder', required=True, type=int)
    parser.add_argument('-Ne', '--population_size', type=int, default=10000)
    parser.add_argument('-l', '--genome_length', required=True, type=int)
    parser.add_argument('-s', '--seed_num', type=int, default=npr.randint(1000000000, size=1)[0])
    parser.add_argument('-o', '--output_vcf', required=True, type=str)
    parser.add_argument('-mu', '--mu_rate', type=str, default='1e-7')
    parser.add_argument('-r', '--recomb_rate', type=str, default='1e-8')

    return parser.parse_args()



if __name__ == '__main__':

        user_args = load_args()

        founder_size = user_args.n_founder
        seed_num = user_args.seed_num
        seq_length = user_args.genome_length
        Ne = user_args.population_size
        output_fp = user_args.output_vcf

        # Simulate the ancestral populations of a set of founders.
        ts = msprime.sim_ancestry(samples=founder_size,
                recombination_rate=user_args.recomb_rate,
                sequence_length=seq_length,
                population_size=Ne,
                random_seed=seed_num)

        # Add mutations onto tree sequence.
        mutated_ts = msprime.sim_mutations(ts, rate=user_args.mu_rate, random_seed=seed_num)

        # Write mutated tree sequence as filepath
        with open(output_fp, "w") as vcf_file:
                mutated_ts.write_vcf(output=vcf_file)