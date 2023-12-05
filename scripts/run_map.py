import networkx as nx
import numpy as np
import argparse
import pandas as pd


# Example: python run_non_paternity_v2.py -f test_fam.nx -p test_fam_profile.txt -c .25

# Created flags for user input
def load_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output_prefix', type=str, default='NewFile')
    parser.add_argument('-n', '--pedigree_filepath', type=str)
    parser.add_argument('-p1', '--isnp_prob', type=float, default=0.01)  ## NP Event
    parser.add_argument('-p2', '--infam_prob', type=float, default=0.50)  ## NP Event happening inside or outside family
    parser.add_argument('-pr', '--profile_filepath', type=str)
    return parser.parse_args()


# Finds nodes that are male and within a 20-50 year age range of an individuals conception and puts them in a list.
def pot_parents(graph, u_indiv):
    sample_parents = []
    sex = nx.get_node_attributes(graph, name='sex')
    Gen = nx.get_node_attributes(graph, name='Gen')
    parents_indiv= list(graph.predecessors(u_indiv))[0]
    parents_gen = Gen[parents_indiv]
    for i in graph.nodes():
        if sex[i] != 'female':
            if Gen[i] == parents_gen:
                sample_parents.append(i)
    if parents_indiv in sample_parents:
        sample_parents.remove(parents_indiv)
    return sample_parents
# ERROR! Gen[i] should not be the same Gen[u_undiv]

def is_parents_connected(graph, indiv):
    '''
  This function will idenitfy an individuals father is the only connection to the main family. This will resolve an
  edge case in the event a np events happens and replaces the father from outside the family, this in turn can lose
  the connectiveness assumption of a family pedigrees and screw up other analysis downstream. If the result of this
  function returns false, then a non-paternity event must come from inside the family.


  We will check for two instances of this occuring in the family pedigree
  1) If the mom is connected to the main family, if so return true.

  :return:
  '''
    # Calculate individual parents
    indiv_pred = list(graph.predecessors(indiv))

    if len(indiv_pred) == 0:
        return np.nan
    # Ignore this command
    if len(indiv_pred) == 1:
        return False

    else:
        sex_dict = nx.get_node_attributes(graph, name='sex')

        # extract sex information for an infant parents
        sex_p1 = sex_dict[indiv_pred[0]]
        sex_p2 = sex_dict[indiv_pred[1]]

        # idenitfy the id of the female parent
        if sex_p1 == 'female':
            node_to_check = indiv_pred[0]
        else:
            node_to_check = indiv_pred[1]

        # Determine if the mom is connected to the main family
        female_pred = list(graph.predecessors(node_to_check))
        if len(female_pred) > 0:
            # if there is a connection then it will be to the main family.
            return True
        else:
            return False


# Iterates through each individual. each iteration has a chance (-c, default 1%) of experiencing a non_paternity.
# If a non-paternity happens, a new father is chosen from a pool of potential parents
def non_paternity(graph, prob1, prob2, output_name):
    '''
  non_paternity takes a .nx file to simulate the probility of non-paternity events.
  probability can be specified
  output file name can be specified
  '''
    count = 0  # count for non-paternity events

    rem_relations = []  # list of old parent child relations to remove
    add_relations = []  # list of new parent child relations to add
    new_pat_bydict = {}
    new_pat_sexdict = {}

    new_pat = max([int(i) for i in graph.nodes]) + 1

    for indiv in graph.nodes:  # For loop around each node(individual) in the graph pedigree
        cur_parents = list(graph.predecessors(indiv))  # initializes current parents of individual
        if len(cur_parents) == 0:
            continue
        potential_parents = pot_parents(graph, indiv)  # lists all nodes
        paternal_event = np.random.choice([0, 1], size=1, p=[1 - prob1, prob1])  # 0:true paternity 1:non paternity
        if len(potential_parents) == 0:
            continue
        # is a np event happens, we begin the process!
        if paternal_event == 1:
            # determines if the indivs np event will break the connectedness of the family
            is_conn = is_parents_connected(graph, indiv)

            # Calculates the prob if the np events happens inside the family. # 0 - draw in fam , 1 - new indiv
            prob_infam = np.random.choice([0, 1], size=1, p=[1 - prob2, prob2])[0]

            # If an indiv is not connected or has no parents (founder)
            if len(potential_parents) == 0 and np.isnan(is_conn):
                continue
            # If reassigning the father to a new would lose the connectedness of the family, draw the np event inside the fam
            elif is_conn == False:
                rep_pat = np.random.choice(potential_parents)
            # If potential parents is not zero and the father is not the only connection, then we let prob decided
            # if the np event happens inside or outside the family.
            elif prob_infam == 1:
                rep_pat = new_pat
            elif prob_infam == 0:
                rep_pat = np.random.choice(potential_parents)

            for i in cur_parents:

                if indiv in potential_parents:
                    potential_parents.remove(indiv)

                if i in potential_parents:
                    potential_parents.remove(i)  # removes parents from potential parents list

            #  First we will check if the individual is a founder (no parents found in nx.predecessor)
            if len(list(graph.predecessors(indiv))) == 0:
                continue

            elif len(list(graph.predecessors(indiv))) == 1:

                parent_1 = list(graph.predecessors(indiv))[0]
                if nx.get_node_attributes(graph, name='sex')[parent_1] == 'male':
                    add_relations.append((rep_pat, indiv))

                    if rep_pat == new_pat:
                        new_pat_bydict[new_pat] = graph.nodes[female_parent]["Gen"]
                        new_pat_sexdict[new_pat] = 'male'
                        new_pat += 1  # increments the newPat and keeps it as a string
                    count += 1
                    # potential_parents.remove(f'{rep_pat}')
            else:
                parent_1 = list(graph.predecessors(indiv))[0]
                parent_2 = list(graph.predecessors(indiv))[1]
                if nx.get_node_attributes(graph, name='sex')[parent_1] == 'male':
                    male_parent = parent_1
                    female_parent = parent_2
                else:
                    male_parent = parent_2
                    female_parent = parent_1
                if (paternal_event == 1):  # If false then we test one parent connection,
                    rem_relations.append((male_parent, indiv))
                    add_relations.append((rep_pat, indiv))
                    if rep_pat == new_pat:
                        new_pat_bydict[new_pat] = graph.nodes[female_parent]["Gen"]
                        new_pat_sexdict[new_pat] = 'male'
                        new_pat += 1  # increments the newPat and keeps it as a string
                    count += 1

                else:  # True paternity events means we fill the child and parent connection.
                    pass
            if rep_pat == new_pat:
                new_pat_bydict[new_pat] = graph.nodes[female_parent]["Gen"]
                new_pat_sexdict[new_pat] = 'male'
                new_pat += 1
        else:  ## no non-paternal event happening, move onto next individual
            pass

    print(f'number of non-paternity events: {count}')
    print(f'Probability of non-paternity event: {prob1 * 100}%')
    graph.remove_edges_from(rem_relations)
    graph.add_edges_from(add_relations)

    nx.set_node_attributes(graph, new_pat_sexdict, name="sex")
    nx.set_node_attributes(graph, new_pat_bydict, name="Gen")
    nx.write_edgelist(graph, f'{output_name}.nx')

    return count


if __name__ == '__main__':
    user_args = load_args()

    u_output = user_args.output_prefix
    u_prob = user_args.isnp_prob
    u_probinfam = user_args.infam_prob

    u_graph = nx.read_edgelist(F'{user_args.pedigree_filepath}', create_using=nx.DiGraph())
    profiles = pd.read_csv(f'{user_args.profile_filepath}', sep='\t')
    sex_dict = dict(zip(profiles['ID'].astype(str).to_numpy(), profiles['Sex'].to_numpy()))
    gen_dict = dict(zip(profiles['ID'].astype(str).to_numpy(), profiles['Gen'].to_numpy()))
    nx.set_node_attributes(u_graph, values=sex_dict, name="sex")
    nx.set_node_attributes(u_graph, values=gen_dict, name="Gen")

    non_paternity(u_graph, u_prob, u_probinfam, u_output)