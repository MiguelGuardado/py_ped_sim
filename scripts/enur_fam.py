import networkx as nx
import numpy as np
import pandas as pd
import argparse

relationship_code = {
    '1_1_direct': 'Parent Child',
    '2_2_direct': 'Grandparent',
    '3_3_direct': 'Great Grandparent',
    '4_4_direct': 'Great Great Grandparent',
    '5_5_direct': 'Great Great Great Grandparent',
    '0_2_full': 'Full sibling',
    '0_4_full': 'Full First Cousin',
    '0_6_full': 'Full Second Cousin',
    '0_8_full': 'Full Third Cousin',
    '0_10_full': 'Full Fourth Cousin',
    '0_12_full': 'Full Fifth Cousin',
    '1_3_full': 'Avuncular',
    '1_5_full': 'Avuncular Once Removed',
    '1_7_full': 'Avuncular Twice Removed',
    '1_9_full': 'Avuncular Thrice Removed',
    '1_11_full': 'Avuncular Forth Removed',
    '1_13_full': 'Avuncular Fifth Removed',
    '0_2_half': 'Half Sibling',
    '0_4_half': 'Half First Cousin',
    '0_6_half': 'Half Second Cousin',
    '0_8_half': 'Half Third Cousin',
    '0_10_half': 'Half Forth Cousin',
    '0_12_half': 'Half Fifth Cousin',
    '0_2_unknown': 'Unknown Sibling',
    '0_4_unknown': 'Unknown First Cousin',
    '0_6_unknown': 'Unknown Second Cousin'}


def find_relationship(gen, mc, rt):
    '''Map the (generation_depth, meiosis_count, relationship_type) triple to a
    named relationship, or 'N/A' if the combination is not catalogued.'''
    return relationship_code.get(f'{gen}_{mc}_{rt}', 'N/A')


def build_ancestor_distances(G):
    '''
    Precompute, once per node, the meiosis distance UP to every ancestor
    (including the node itself at distance 0). One BFS per node on the reversed
    graph replaces the per-pair ancestor/shortest-path traversals the original
    performed for every one of the O(n^2) pairs.

    Returns:
        up  : {node: {ancestor: distance}}   distance = # edges ancestor -> node
    '''
    G_rev = G.reverse(copy=False)
    return {n: nx.single_source_shortest_path_length(G_rev, n) for n in G}


def relationship_type(n_mrcas):
    '''Two shared lowest common ancestors -> full, one -> half. Three or more
    indicates a consanguinity loop the catalogue does not model.'''
    if n_mrcas == 2:
        return 'full'
    if n_mrcas == 1:
        return 'half'
    return 'N/A'


def classify_pair(node1, node2, d1, d2):
    '''
    Determine the relationship record for one ordered pair using only the
    precomputed ancestor->distance maps d1 (for node1) and d2 (for node2).
    Returns a row [Fam_ID, ID1, ID2, MC, GD, RT, RC] or None if unrelated.
    '''
    # Direct line: one individual is an ancestor of the other. Checked first so
    # ancestor/descendant pairs are labelled 'direct' (matches the original's
    # has_path precedence over the common-ancestor branch).
    if node1 in d2:            # node1 is an ancestor of node2
        dist = d2[node1]
        return ['1', node1, node2, dist, dist, 'direct',
                find_relationship(dist, dist, 'direct')]
    if node2 in d1:            # node2 is an ancestor of node1
        dist = d1[node2]
        return ['1', node1, node2, dist, dist, 'direct',
                find_relationship(dist, dist, 'direct')]

    # Otherwise look for common ancestors. Sorting the intersection reproduces
    # the original's np.intersect1d ordering, so ties pick the same MRCA.
    common = sorted(d1.keys() & d2.keys())
    if not common:
        return None

    sums = {c: d1[c] + d2[c] for c in common}
    best = min(sums.values())
    mrcas = [c for c in common if sums[c] == best]
    mrca = mrcas[0]

    meioses_count = d1[mrca] + d2[mrca]
    generation_depth = abs(d1[mrca] - d2[mrca])
    rt = relationship_type(len(mrcas))
    rc = find_relationship(generation_depth, meioses_count, rt)
    return ['1', node1, node2, meioses_count, generation_depth, rt, rc]


def find_pairwise_relationships(G, output_prefix, sample_file=None):
    assert nx.is_directed_acyclic_graph(G), 'Graph has to be acyclic and directed.'

    up = build_ancestor_distances(G)

    # weakly-connected component id per node: relationships only exist within a
    # component, so cross-component pairs are skipped without any ancestor work.
    comp_id = {}
    for cid, comp in enumerate(nx.weakly_connected_components(G)):
        for node in comp:
            comp_id[node] = cid

    family_list = list(G)
    family_list.sort(key=int)

    if sample_file is not None:
        fam_list_1 = list(np.atleast_1d(sample_file))
    else:
        fam_list_1 = family_list

    relationships = []
    for node_index, node1 in enumerate(fam_list_1):
        # sample mode compares each sampled node against everyone; default mode
        # walks the upper triangle to avoid duplicate (a,b)/(b,a) rows.
        fam_list_2 = family_list if sample_file is not None else family_list[node_index + 1:]
        d1 = up[node1]

        for node2 in fam_list_2:
            if node1 == node2 or comp_id[node1] != comp_id[node2]:
                continue
            row = classify_pair(node1, node2, d1, up[node2])
            if row is not None:
                relationships.append(row)

    output_df = pd.DataFrame(relationships,
                             columns=['Fam_ID', 'ID1', 'ID2', 'MC', 'GD', 'RT', 'RC'])
    output_df.to_csv(f'{output_prefix}_rel.csv', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--networkx_file', type=str, required=True)
    parser.add_argument('-sf', '--sample_file', type=str)
    parser.add_argument('-o', '--output_prefix', type=str, required=False, default=None)
    parser.add_argument('-t', '--true_half', action='store_true',
                        help='reserved: distinguish true half-sibs from single-known-ancestor pairs')
    user_args = parser.parse_args()

    if user_args.output_prefix in (None, 'None'):
        user_args.output_prefix = user_args.networkx_file.split('.nx')[0]

    fam_graph = nx.read_edgelist(user_args.networkx_file, create_using=nx.DiGraph())

    if user_args.sample_file in (None, 'None'):
        find_pairwise_relationships(fam_graph, output_prefix=user_args.output_prefix)
    else:
        sample_file = np.loadtxt(user_args.sample_file, dtype=str)
        find_pairwise_relationships(fam_graph, user_args.output_prefix, sample_file)