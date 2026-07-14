import networkx as nx
import subprocess
import argparse
import os
import sys
import glob
import numpy as np
import pandas as pd
import time
from run_single_family_broadening import *
from networkx.relabel import convert_node_labels_to_integers
from networkx.algorithms.operators.binary import disjoint_union
from networkx.algorithms.traversal.depth_first_search import dfs_predecessors

"""
run_full_family_broadening.py takes in a main family or simulates a main family.
It then iterates through each non-root founder and simulates a family with a
similar number of generations, attaching it onto that founder.

If a main family is supplied (-mf), the root-founder generation must be omitted
from the -y flag.

Outputs:
  -o / --output_prefix            : the FINAL, fully broadened pedigree
                                    ({prefix}.nx + {prefix}_profiles.txt)
  -mo / --main_family_output_prefix : the simulated main family, written ONLY when
                                    no main family was supplied via -mf

ex: python run_ped_sim.py -t run_full_family_broadening -c ipumps_sibship_dist.txt -y 1850 1860 1870
    python run_full_family_broadening.py -c ipumps_sibship_dist.txt -y 1850 1860 1870
"""

# retry ceiling shared by every "simulate a family until it joins" loop
MAX_JOIN_TRIES = 10

# record start time
start = time.time()


def load_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-y', '--years_to_sample', nargs='+', type=int)
    parser.add_argument('-c', '--census_filepath', type=str)
    parser.add_argument('-o', '--output_prefix', type=str, default='joint_family_output')
    parser.add_argument('-mo', '--main_family_output_prefix', type=str, default='main_family')
    parser.add_argument('-mf', '--main_family')
    return parser.parse_args()


def fam_check_gen(profiles_path):
    '''
    Read a family profile file and return its number of generations
    (unique Gen values minus the root generation).
    '''
    df = pd.read_csv(profiles_path, sep='\t')
    return len(df['Gen'].unique()) - 1


def run(cmd, **kwargs):
    '''Thin wrapper so every shell call is invoked identically.'''
    return subprocess.run(cmd, shell=True, **kwargs)


def validate_years_in_census(census_path, years):
    '''
    sim_pedigree.py silently exits (writing no .nx) if any requested year is not
    in the census file. Catch that here so the run fails with a clear message
    instead of looping on satellite families that never get created.
    '''
    census_years = set(pd.read_csv(census_path, sep='\t')['Year'].astype(int))
    missing = [y for y in years if y not in census_years]
    if missing:
        sys.exit(
            f'ERROR: requested year(s) {missing} are not in the census file '
            f'{census_path}.\nAvailable years: {sorted(census_years)}'
        )


def validate_generations(profiles_path, years):
    '''
    run_single_family_broadening matches a founder to a satellite member by exact
    Gen value, and Gen IS the sample year (e.g. 1850), not a 0/1/2 index. Validate
    against EVERY generation present in the main family (all Gen values in its
    profiles, excluding the NA root generation) -- not just the generations that
    founders happen to occupy -- so a family generation like 1910 with no founders
    isn't silently dropped. Any family generation missing from -y is a fatal
    mismatch, so fail fast and say which.
    '''
    def _norm(v):
        # the NA root generation makes pandas read Gen as float (1850.0); normalize
        # every value to a plain integer-year string so comparisons are exact
        s = str(v).strip()
        if s in ('NA', 'nan', ''):
            return None
        try:
            return str(int(float(s)))
        except ValueError:
            return s

    prof = pd.read_csv(profiles_path, sep='\t')
    fam_gens = {g for g in (_norm(v) for v in prof['Gen']) if g is not None}
    year_set = {str(int(y)) for y in years}
    missing = sorted(fam_gens - year_set)
    if missing:
        sys.exit(
            f'ERROR: the main family has generation(s) {missing} not present in -y '
            f'{sorted(year_set)}.\nSet -y to every generation in the family '
            '(the Gen column of its _profiles.txt, excluding the NA root generation).'
        )


def simulate_main_family(years_str, census, out_prefix, n_years):
    '''
    Simulate a main family until it has at least as many generations as the user
    requested, then return its networkx graph. Writes {out_prefix}.nx and
    {out_prefix}_profiles.txt (this is the -mo output).
    '''
    num_gen = 0
    while num_gen < n_years:
        run(f'python scripts/sim_pedigree_v2.py {years_str} -c {census} -o {out_prefix}')
        graph = nx.read_edgelist(f'{out_prefix}.nx', create_using=nx.DiGraph())
        num_gen = fam_check_gen(f'{out_prefix}_profiles.txt')
        print('main family created')
    return graph


def simulate_and_join(years_str, census, scratch_prefix, fam_index,
                      n1_prefix, out_prefix, founder):
    '''
    Simulate one satellite family and attempt to broaden it onto `founder`.
    Retries with a fresh simulated family until the join succeeds or the retry
    ceiling is hit. On success the satellite's scratch files are removed.

      n1_prefix  : prefix of the pedigree we are attaching TO
                   (the main family for the first founder, then the growing
                    joint family for every founder after that)
      out_prefix : the joint-family output prefix (-o), always written here
    '''
    scratch_fam = f'{scratch_prefix}_fam{fam_index}'
    for attempt in range(1, MAX_JOIN_TRIES + 1):
        run(f'python scripts/sim_pedigree.py {years_str} -c {census} -o {scratch_fam}',
            capture_output=True)
        join_cmd = (
            f'python scripts/run_single_family_broadening.py '
            f'-n1 {n1_prefix}.nx -n2 {scratch_fam}.nx '
            f'-pr1 {n1_prefix}_profiles.txt -pr2 {scratch_fam}_profiles.txt '
            f'-o {out_prefix} -cf {founder}'
        )
        result = run(join_cmd, capture_output=True)
        if result.returncode != 1:
            print(f'connection successful for fam{fam_index} at founder {founder}')
            _safe_remove(f'{scratch_fam}.nx', f'{scratch_fam}_profiles.txt')
            return
        print(f'creating new fam (attempt {attempt}) for founder {founder}')
    # surface the underlying error so exhaustion is not silent
    last_err = (result.stderr or b'').decode(errors='replace').strip()
    sys.exit(
        f'ERROR: could not attach a family to founder {founder} after '
        f'{MAX_JOIN_TRIES} tries.\nLast error from run_single_family_broadening:\n{last_err}'
    )


def _safe_remove(*paths):
    for p in paths:
        if os.path.exists(p):
            os.remove(p)


def cleanup(scratch_prefix):
    '''Remove every scratch artifact; the -o and -mo outputs are left untouched.'''
    leftovers = glob.glob(f'{scratch_prefix}_fam*') + ['relabled_fam.nx']
    _safe_remove(*leftovers)


if __name__ == '__main__':
    user_args = load_args()

    if user_args.output_prefix in (None, 'None'):
        user_args.output_prefix = 'joint_family_output'

    u_years = user_args.years_to_sample
    u_census = user_args.census_filepath
    u_main_family = user_args.main_family
    u_output = user_args.output_prefix
    u_mo = user_args.main_family_output_prefix

    # intermediates live under their own prefix so -o only ever names the final result
    scratch_prefix = f'{u_output}_scratch'

    # build the -y terminal argument once
    years_str = '-y ' + ' '.join(str(y) for y in u_years)

    # every satellite family is simulated from these years; bail early if any are invalid
    validate_years_in_census(u_census, u_years)

    supplied_main_family = u_main_family not in (None, 'None')

    if supplied_main_family:
        # normalize once: strip .nx so {prefix}.nx / {prefix}_profiles.txt are always correct
        main_family_filepath = u_main_family.replace('.nx', '')
        main_family = nx.read_edgelist(f'{main_family_filepath}.nx', create_using=nx.DiGraph())

        # generation-based check (matches the check used on the simulated path)
        if fam_check_gen(f'{main_family_filepath}_profiles.txt') < len(u_years):
            sys.exit('Number of generations in submitted family (-mf) does not match '
                     'number of generations requested (-y)')

        if user_args.main_family_output_prefix != 'main_family':
            print('note: -mo is ignored when a main family is supplied via -mf')
    else:
        main_family_filepath = u_mo
        main_family = simulate_main_family(years_str, u_census, u_mo, len(u_years))

    # non-root founders to broaden upon (assumes the first two founders are the roots)
    founders = find_founders(main_family)
    roots, founders = founders[:2], founders[2:]
    assert len(roots) == 2, 'expected exactly two root founders at the head of find_founders()'
    print('Founders to extend breadth upon: ', founders)

    # every generation in the family (all years, not just founders') must be requested
    validate_generations(f'{main_family_filepath}_profiles.txt', u_years)

    try:
        # first founder: attach onto the main family, producing the joint family (-o)
        simulate_and_join(years_str, u_census, scratch_prefix, 0,
                          n1_prefix=main_family_filepath, out_prefix=u_output,
                          founder=founders[0])

        # remaining founders: attach onto the growing joint family so the main family is preserved
        for fam_index, founder in enumerate(founders[1:], start=1):
            simulate_and_join(years_str, u_census, scratch_prefix, fam_index,
                              n1_prefix=u_output, out_prefix=u_output,
                              founder=founder)
    finally:
        cleanup(scratch_prefix)

    print(f'\nDone in {time.time() - start:.1f}s')
    print(f'Final pedigree : {u_output}.nx, {u_output}_profiles.txt')
    if not supplied_main_family:
        print(f'Main family    : {u_mo}.nx, {u_mo}_profiles.txt')