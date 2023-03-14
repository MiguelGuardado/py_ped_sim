import networkx as nx
import numpy as np
import pandas as pd
import argparse
import datetime

relationship_code = {
    '1_1_direct':'Parent Child',
    '2_2_direct':'Grandparent',
    '3_3_direct':'Great Grandparent',
    '4_4_direct':'Great Great Grandparent',
    '5_5_direct':'Great Great Great Grandparent',
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


def find_rt(node1, node2):
    global G_undir
    global true_half
    '''
    This function will determine the relationtype for a pair of individuals who share a common ancestor but are not
    considered direct connections.

    :param G:
    :param node1:
    :param node2:
    :return:
    '''
    # Returns all the shortest paths between a pair of two individuals
    sp = list(nx.all_shortest_paths(G_undir, node1, node2))

    if len(sp) == 2:
        return 'full'

    elif len(sp) == 1:
        if true_half:
            pass

        return 'half'

    else:
        return 'N/A'

def find_com_anc(G, node1, node2):
    '''
    Code block given to me by ChatGPT, very kind of them.
    :param G: Graph to use to find the most recent common ancestor
    :param node1: Individual 1
    :param node2: Individual 2
    :return: most recent common ancestor between a set of two nodes.
    '''
    try:
        # Find the set of ancestors for both nodes
        ancestors1 = set(nx.ancestors(G, node1))
        ancestors2 = set(nx.ancestors(G, node2))
        # Find the intersection of the two sets to find the common ancestors
        common_ancestors = list(ancestors1.intersection(ancestors2))
        if common_ancestors:

            # mrca is taking a long time so we can just instead look at the first intersecting ancestor.
            ca = common_ancestors[0]
            return ca
        else:
            return None
    except nx.NetworkXNoPath:
        return None

def find_mc_with_ca(G, node1, node2, com_anc):
    node1_path = nx.shortest_path(G, com_anc, node1)
    node2_path = nx.shortest_path(G, com_anc, node2)

    if len(node1_path) == len(node2_path):
        mc = len(node1_path) + len(node2_path) - 2
    else:
        mc = len(np.setdiff1d(node1_path, node2_path)) + len(np.setdiff1d(node2_path, node1_path))
    return mc


def find_relationship(gen, mc, rt):
    rc = f'{gen}_{mc}_{rt}'
    if rc in relationship_code.keys():
        return relationship_code[rc]
    else:
        return 'N/A'


def find_pairwise_relationships(G, networkx_prefix):
    relationships = []

    family_list = list(G)
    family_list.sort(key=int)

    node_index = 0
    for node1 in family_list:
        # print(node_index/len(family_list))
        for node2 in family_list[node_index + 1:]:
            if node1 == node2:
                continue
            # Function to return a single most recent common ancestor(mrca) between two individuals.
            ca = find_com_anc(G, node1, node2)

            if nx.has_path(G, node1, node2):
                # In this case, if a path exsist between two nodes in a graph then are a direct connection between these members.
                path = nx.shortest_path(G, node1, node2)
                meioses_count = len(path) - 1
                generation_depth = nx.shortest_path_length(G, node1, node2)

                rc = find_relationship(generation_depth, meioses_count, 'direct')
                relationships.append(['1', node1, node2, meioses_count, generation_depth, 'direct', rc])

            elif nx.has_path(G, node2, node1):
                # In this case, if a path exsist between two nodes in a graph then are a direct connection between these members.
                path = nx.shortest_path(G, node2, node1)
                meioses_count = len(path) - 1
                generation_depth = nx.shortest_path_length(G, node2, node1)

                rc = find_relationship(generation_depth, meioses_count, 'direct')
                relationships.append(['1', node1, node2, meioses_count, generation_depth, 'direct', rc])

            elif ca is not None:
                # If a common ancestor exist, then there is at least one common ancestor between a pair of two nodes.
                # meioses_count = nx.shortest_path_length(G, mrca, node1) + nx.shortest_path_length(G, mrca, node2)
                meioses_count = find_mc_with_ca(G, node1, node2, ca)
                generation_depth = np.abs(nx.shortest_path_length(G, ca, node1) - nx.shortest_path_length(G, ca, node2))
                rt = find_rt(node1, node2)

                rc = find_relationship(generation_depth, meioses_count, rt)
                relationships.append(['1', node1, node2, meioses_count, generation_depth, rt, rc])

            else:
                relationships.append(['1', node1, node2, 0, 0, 'NA', 'NA'])

        # Index to keep track of upper triangle of pairwise individuals in family_list
        node_index+=1

    # Convert list of pairwise relationships statistics into pandas dataframe to output.
    output_df = pd.DataFrame(relationships, columns=['Fam_ID', 'ID1', 'ID2', 'MC', 'GD', 'RT', "RC"])
    output_fn = f"{networkx_prefix}_rel.csv"
    output_df.to_csv(output_fn, index=False)

if __name__ == '__main__':

    # Load in command line user inputs
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--networkx_file', type=str, required=True)
    parser.add_argument('-o', '--output_prefix', type=str, required=False, default=None)
    parser.add_argument('-t', '--true_half', type=bool, default=False)
    user_args = parser.parse_args()

    if user_args.output_prefix is None:
        user_args.output_prefix = user_args.networkx_file.split('.nx')[0]

    # We convert true_half into a variable so we can make it a global variable for find_rt() to read in.
    true_half = user_args.true_half

    # Read in family pedigree represented as a directed acyclic graph
    G = nx.read_edgelist(user_args.networkx_file, create_using=nx.DiGraph())
    G_undir = G.to_undirected()

    # Feed into function to determine relationships statistics for the family pedigrees
    find_pairwise_relationships(G, networkx_prefix=user_args.output_prefix)



