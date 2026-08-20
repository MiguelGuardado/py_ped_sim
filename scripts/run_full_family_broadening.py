import networkx as nx
import subprocess
import argparse
import os
import sys
import glob
import pandas as pd
import time
from concurrent.futures import ThreadPoolExecutor
from run_single_family_broadening import find_founders

"""
run_full_family_broadening.py takes in a main family or simulates a main family.
It then iterates through each non-root founder and simulates a family with a
similar number of generations, attaching it onto that founder.

If a main family is supplied (-mf), the root-founder generation must be omitted
from the -y flag.

Parallelism: the satellite families are independent, so they are simulated
concurrently; the joins that follow stay sequential because each one reads the
joint family produced by the previous join. Concurrency is chosen automatically
(no flag): it uses the CPUs actually available to this process (affinity-aware,
so it respects an HPC scheduler's allocation). Override with the PED_SIM_JOBS or
SLURM_CPUS_PER_TASK environment variable if you want to cap or force it (set to 1
for fully serial behavior).

Outputs:
  -o / --output_prefix            : the FINAL, fully broadened pedigree
                                    ({prefix}.nx + {prefix}_profiles.txt)
  -mo / --main_family_output_prefix : the simulated main family, written ONLY when
                                    no main family was supplied via -mf
  -k / --keep_families True|False : keep the satellite families that were attached,
                                    as {prefix}_satellite_founder<ID>.nx (+ _profiles.txt).
                                    Defaults to False, i.e. satellites are deleted.

ex: python run_ped_sim.py -t run_full_family_broadening -c ipumps_sibship_dist.txt -y 1850 1860 1870
    python run_full_family_broadening.py -c ipumps_sibship_dist.txt -y 1850 1860 1870
"""

# retry ceiling shared by every "simulate a family until it joins" loop
MAX_JOIN_TRIES = 10

# Pass 1 gives every main-family founder two parents, but only ONE of those two
# (the satellite descendant) has parents of its own -- the other is the
# satellite's married-in founder, who has none. That shortfall propagates down:
# main-family individuals in the deepest generation end up with 6
# great-grandparents instead of 8, and the founders themselves get 2
# grandparents instead of 4. Pass 2 extends exactly those parentless parents --
# and nothing beyond them -- which closes both gaps at bounded cost.
# Set False to skip pass 2 and keep only the single broadening pass.
EXTEND_FOUNDER_PARENTS = True

SCRIPT_DIR = 'scripts'
SIM_PED = os.path.join(SCRIPT_DIR, 'sim_pedigree.py')
SIM_PED_V2 = os.path.join(SCRIPT_DIR, 'sim_pedigree_v2.py')
JOIN = os.path.join(SCRIPT_DIR, 'run_single_family_broadening.py')

# record start time
start = time.time()


def str2bool(value):
    '''
    Parse a True/False command-line value.

    argparse's `type=bool` cannot be used for this: it just calls bool() on the
    string, and every non-empty string is truthy -- so `-k False` would come back
    as True. This converter compares the text instead. 'None' maps to False so the
    flag survives being forwarded as an f-string from run_ped_sim.py.
    '''
    if isinstance(value, bool):
        return value
    v = str(value).strip().lower()
    if v in ('true', 't', 'yes', 'y', '1'):
        return True
    if v in ('false', 'f', 'no', 'n', '0', 'none', ''):
        return False
    raise argparse.ArgumentTypeError(f"expected True or False, got '{value}'")


def load_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-y', '--years_to_sample', nargs='+', type=int)
    parser.add_argument('-c', '--census_filepath', type=str)
    parser.add_argument('-o', '--output_prefix', type=str, default='joint_family_output')
    parser.add_argument('-mo', '--main_family_output_prefix', type=str, default='main_family')
    parser.add_argument('-mf', '--main_family')
    parser.add_argument('-k', '--keep_families', type=str2bool, nargs='?',
                        const=True, default=False,
                        help='True/False: keep the simulated satellite families that were '
                             'attached to the main family (default False, i.e. delete them '
                             'after joining). Bare -k also means True.')
    return parser.parse_args()


def run(args, **kwargs):
    '''
    Run a subprocess from an argument LIST with the current interpreter and no
    shell. Argument lists avoid shell quoting/globbing pitfalls and sys.executable
    (built into the command builders) ensures the child uses the same Python --
    both required for the code to behave the same on Linux, HPC, and Windows.
    '''
    return subprocess.run(args, **kwargs)


def sim_ped_cmd(years_args, census, out_prefix):
    return [sys.executable, SIM_PED, *years_args, '-c', census, '-o', out_prefix]


def sim_ped_v2_cmd(years_args, census, out_prefix):
    return [sys.executable, SIM_PED_V2, *years_args, '-c', census, '-o', out_prefix]


def join_cmd(n1_prefix, scratch_fam, out_prefix, founder):
    return [sys.executable, JOIN,
            '-n1', f'{n1_prefix}.nx', '-n2', f'{scratch_fam}.nx',
            '-pr1', f'{n1_prefix}_profiles.txt', '-pr2', f'{scratch_fam}_profiles.txt',
            '-o', out_prefix, '-cf', str(founder)]


def pick_n_jobs(n_tasks):
    '''
    Choose how many satellite simulations to run at once, with no user flag:

      1. an explicit PED_SIM_JOBS / SLURM_CPUS_PER_TASK env var, if set;
      2. otherwise the CPUs actually available to this process
         (os.sched_getaffinity respects an HPC scheduler's cgroup/cpuset
         allocation, so we don't oversubscribe a shared node);
      3. otherwise os.cpu_count() where affinity isn't available (macOS/Windows).

    Never exceeds the number of tasks, and is always at least 1.
    '''
    n = 0
    for var in ('PED_SIM_JOBS', 'SLURM_CPUS_PER_TASK'):
        val = os.environ.get(var)
        if val:
            try:
                n = int(val)
                break
            except ValueError:
                pass
    if n <= 0:
        try:
            n = len(os.sched_getaffinity(0))   # Linux / HPC: cores allocated to us
        except AttributeError:
            n = os.cpu_count() or 1            # macOS / Windows
    return max(1, min(n, n_tasks))


def fam_check_gen(profiles_path):
    '''
    Read a family profile file and return its number of generations
    (unique Gen values minus the root generation).
    '''
    df = pd.read_csv(profiles_path, sep='\t')
    return len(df['Gen'].unique()) - 1


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


def simulate_main_family(years_args, census, out_prefix, n_years):
    '''
    Simulate a main family until it has at least as many generations as the user
    requested, then return its networkx graph. Writes {out_prefix}.nx and
    {out_prefix}_profiles.txt (this is the -mo output).
    '''
    num_gen = 0
    while num_gen < n_years:
        run(sim_ped_v2_cmd(years_args, census, out_prefix))
        graph = nx.read_edgelist(f'{out_prefix}.nx', create_using=nx.DiGraph())
        num_gen = fam_check_gen(f'{out_prefix}_profiles.txt')
        print('main family created')
    return graph


def simulate_satellite(years_args, census, scratch_fam):
    '''Simulate one satellite family into {scratch_fam}.nx / _profiles.txt.'''
    run(sim_ped_cmd(years_args, census, scratch_fam), capture_output=True)


def frontier_founders(joint_prefix, original_founders):
    '''
    Return the parents of the original main-family founders that are themselves
    still parentless -- i.e. the married-in founder each satellite brought along
    in pass 1. These are exactly the individuals whose missing parents cause the
    deepest generation to have 6 great-grandparents instead of 8.

    Extending only this set keeps the cost bounded (one satellite per founder,
    ~7 on test_fam). Recursively extending every newly created founder instead
    does not converge -- measured 44 -> 696 -> 7,606 individuals over two rounds
    with 1,512 founders still outstanding.

    Individuals whose generation is NA (a satellite's root founders) are skipped:
    they sit at the top of the pedigree, with no earlier generation to sample from.
    '''
    G = nx.read_edgelist(f'{joint_prefix}.nx', create_using=nx.DiGraph())
    prof = pd.read_csv(f'{joint_prefix}_profiles.txt', sep='\t')
    gen_by_id = dict(zip(prof['ID'].astype(str), prof['Gen'].astype(str)))

    def has_real_gen(node):
        return gen_by_id.get(str(node), 'NA').strip() not in ('NA', 'nan', '')

    frontier = []
    for f in original_founders:
        if f not in G:
            continue
        for parent in G.predecessors(f):
            if (G.in_degree(parent) == 0 and has_real_gen(parent)
                    and parent not in frontier):
                frontier.append(parent)
    return frontier


def simulate_and_join(years_args, census, scratch_prefix, fam_index,
                      n1_prefix, out_prefix, founder, presimulated=False,
                      keep_families=False):
    '''
    Attach a satellite family onto `founder`, retrying with a fresh simulation
    until the join succeeds or the retry ceiling is hit. When `presimulated` is
    True the first attempt reuses the satellite already generated in the parallel
    pre-simulation phase; every retry simulates a new one. Scratch files are
    removed on success.

      n1_prefix  : prefix of the pedigree we are attaching TO
                   (the main family for the first founder, then the growing
                    joint family for every founder after that)
      out_prefix : the joint-family output prefix (-o), always written here
    '''
    scratch_fam = f'{scratch_prefix}_fam{fam_index}'
    result = None
    for attempt in range(1, MAX_JOIN_TRIES + 1):
        if not (attempt == 1 and presimulated):
            simulate_satellite(years_args, census, scratch_fam)
        result = run(join_cmd(n1_prefix, scratch_fam, out_prefix, founder),
                     capture_output=True)
        if result.returncode != 1:
            print(f'connection successful for fam{fam_index} at founder {founder}')
            if keep_families:
                for path in retain_satellite(scratch_fam, out_prefix, founder):
                    print(f'  kept satellite family: {path}')
            else:
                _safe_remove(f'{scratch_fam}.nx', f'{scratch_fam}_profiles.txt')
            return
        print(f'creating new fam (attempt {attempt}) for founder {founder}')
    # surface the underlying error so exhaustion is not silent
    last_err = (result.stderr or b'').decode(errors='replace').strip() if result else ''
    sys.exit(
        f'ERROR: could not attach a family to founder {founder} after '
        f'{MAX_JOIN_TRIES} tries.\nLast error from run_single_family_broadening:\n{last_err}'
    )


def _safe_remove(*paths):
    for p in paths:
        if os.path.exists(p):
            os.remove(p)


def retain_satellite(scratch_fam, out_prefix, founder):
    '''
    Promote a satellite family that was successfully joined out of the scratch
    namespace and into a stable, user-facing name (the -k / --keep_families path).

    Naming is by FOUNDER rather than by loop index: each founder is extended by
    exactly one satellite, in either pass, so the founder ID both is unique and
    records where the family was attached. Renaming also moves the files out of
    the {prefix}_scratch_fam* pattern, so the end-of-run cleanup leaves them alone.

    Only satellites that actually joined are kept -- families discarded during a
    retry were overwritten in place and never attached to anything.
    '''
    kept = []
    for suffix in ('.nx', '_profiles.txt'):
        src = f'{scratch_fam}{suffix}'
        if os.path.exists(src):
            dst = f'{out_prefix}_satellite_founder{founder}{suffix}'
            os.replace(src, dst)   # atomic, and works the same on Windows
            kept.append(dst)
    return kept


def cleanup(scratch_prefix, out_prefix=None, keep_families=False):
    '''
    Remove every scratch artifact; the -o and -mo outputs are left untouched.

    When keep_families is False we ALSO clear any {out_prefix}_satellite_founder*
    files. Those are only produced by a previous -k True run: retain_satellite()
    renames them out of the scratch namespace, so the scratch glob below can never
    match them and they would otherwise survive every later run with the same -o
    prefix. Without this, `-k False` appears not to delete the satellite families,
    when in fact it deleted its own and the user was looking at stale ones.
    '''
    leftovers = glob.glob(f'{scratch_prefix}_fam*') + ['relabled_fam.nx']
    _safe_remove(*leftovers)

    if not keep_families and out_prefix:
        stale = (glob.glob(f'{out_prefix}_satellite_founder*.nx')
                 + glob.glob(f'{out_prefix}_satellite_founder*_profiles.txt'))
        if stale:
            print(f'Removed {len(stale)} satellite file(s) left over from an '
                  f'earlier --keep_families run')
            _safe_remove(*stale)


if __name__ == '__main__':
    user_args = load_args()

    if user_args.output_prefix in (None, 'None'):
        user_args.output_prefix = 'joint_family_output'

    u_years = user_args.years_to_sample
    u_census = user_args.census_filepath
    u_main_family = user_args.main_family
    u_output = user_args.output_prefix
    u_mo = user_args.main_family_output_prefix
    u_keep = user_args.keep_families

    # intermediates live under their own prefix so -o only ever names the final result
    scratch_prefix = f'{u_output}_scratch'

    # -y as a list of tokens, spliced into each subprocess argument list
    years_args = ['-y'] + [str(y) for y in u_years]

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
        main_family = simulate_main_family(years_args, u_census, u_mo, len(u_years))

    # non-root founders to broaden upon (assumes the first two founders are the roots)
    founders = find_founders(main_family)
    roots, founders = founders[:2], founders[2:]
    assert len(roots) == 2, 'expected exactly two root founders at the head of find_founders()'
    print('Founders to extend breadth upon: ', founders)

    # every generation in the family (all years, not just founders') must be requested
    validate_generations(f'{main_family_filepath}_profiles.txt', u_years)

    # ---- parallel phase: satellite families are independent, so simulate them at once
    n_founders = len(founders)
    n_jobs = pick_n_jobs(n_founders)
    print(f'Pre-simulating {n_founders} satellite families with {n_jobs} parallel worker(s)')
    with ThreadPoolExecutor(max_workers=n_jobs) as pool:
        # list() forces evaluation so any worker exception propagates here
        list(pool.map(
            lambda i: simulate_satellite(years_args, u_census, f'{scratch_prefix}_fam{i}'),
            range(n_founders)))

    # ---- sequential phase: joins must be ordered (each reads the previous joint family)
    try:
        # first founder: attach onto the main family, producing the joint family (-o)
        simulate_and_join(years_args, u_census, scratch_prefix, 0,
                          n1_prefix=main_family_filepath, out_prefix=u_output,
                          founder=founders[0], presimulated=True, keep_families=u_keep)

        # remaining founders: attach onto the growing joint family so the main family is preserved
        for fam_index, founder in enumerate(founders[1:], start=1):
            simulate_and_join(years_args, u_census, scratch_prefix, fam_index,
                              n1_prefix=u_output, out_prefix=u_output,
                              founder=founder, presimulated=True, keep_families=u_keep)

        # ---- pass 2: extend the parentless parents pass 1 introduced, so the
        # deepest generation gets 8 great-grandparents (not 6) and the
        # main-family founders themselves get 4 grandparents (not 2)
        if EXTEND_FOUNDER_PARENTS:
            frontier = frontier_founders(u_output, founders)
            if frontier:
                print(f'Second pass: extending {len(frontier)} parentless founder '
                      f'parent(s) to complete the great-grandparent generation')
                offset = n_founders
                with ThreadPoolExecutor(max_workers=pick_n_jobs(len(frontier))) as pool:
                    list(pool.map(
                        lambda k: simulate_satellite(
                            years_args, u_census, f'{scratch_prefix}_fam{offset + k}'),
                        range(len(frontier))))
                for k, parent in enumerate(frontier):
                    simulate_and_join(years_args, u_census, scratch_prefix, offset + k,
                                      n1_prefix=u_output, out_prefix=u_output,
                                      founder=parent, presimulated=True,
                                      keep_families=u_keep)
    finally:
        cleanup(scratch_prefix, out_prefix=u_output, keep_families=u_keep)

    print(f'\nDone in {time.time() - start:.1f}s')
    print(f'Final pedigree : {u_output}.nx, {u_output}_profiles.txt')
    if not supplied_main_family:
        print(f'Main family    : {u_mo}.nx, {u_mo}_profiles.txt')
    if u_keep:
        kept = sorted(glob.glob(f'{u_output}_satellite_founder*.nx'))
        print(f'Satellites     : {len(kept)} kept as '
              f'{u_output}_satellite_founder<ID>.nx (+ _profiles.txt)')