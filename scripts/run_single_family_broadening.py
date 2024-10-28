import networkx as nx
import numpy as np
import argparse
import pandas as pd
import math
from networkx.relabel import convert_node_labels_to_integers
from networkx.algorithms.operators.binary import disjoint_union
from networkx.algorithms.traversal.depth_first_search import dfs_predecessors

'''
run_single_familly_broadening connects two families together at a single non-root founder
The founder can be specified. otherwise the founder is chosen at random

'''

# python run_single_family_broadening.py -n1 sym_family.nx -n2 sym_family.nx -p1 sym_family_profiles.txt -p2 sym_family_profiles.txt -o sym_joint
#example python run_single_family_broadening.py -n1 testfam_half.nx -n2 testfam_half.nx -p1 test_fam_profilev2.txt -p2 test_fam_profilev2.txt -s 3
# python run_single_family_broadening.py -n1 main_sim.nx -n2 fam0.nx -p1 main_sim_profiles.txt -p2 fam0_profiles.txt -o sym_joint -cf 9

def load_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n1', '--family1', type=str)
    parser.add_argument('-n2', '--family2', type=str)
    parser.add_argument('-pr1', '--profiles1', type=str)
    parser.add_argument('-pr2', '--profiles2', type=str)
    parser.add_argument('-o', '--output_prefix', type=str)
    parser.add_argument('-cf', '--chosen_founder', type = str, default = 'empty')
    parser.add_argument('-cs', '--chosen_sub', type = str, default = 'empty')
    return parser.parse_args()
    
    
def replace_nodes(family1, family2,  node1, node2):
  '''
  This fucntion takes 2 pedigrees on the same graph (family1 and family2) and 2 nodes (node1 and node2) 
  to make a connection 
  node1 replces nodes2 
  keep predecessors of node 1
  keep successors of node 2 but removes successors 
  '''
  #print(family2.nodes())
  n2_pred = list(family2.predecessors(node2))
  #print('predecessors of node2: ', n2_pred)
  replacement = []
  #print(replacement)
  remove = [node2]
  #add children to remove
  remove += list(family2.successors(node2))
  # adds all spouces to the list 'remove' 
  for i in list(family2.successors(node2)):
      remove += list(family2.predecessors(i))
  #set new list without duplicates  
  new_remove = [*set(remove)]    
  #print('removed node: ', new_remove)
  for i in n2_pred:
    replacement.append((i, node1))
  
  family1.add_edges_from(replacement)
  new_remove += list(nx.isolates(family1))
  print('Connections: ', replacement)
  family1.remove_nodes_from(new_remove)
  print('removed nodes: ', new_remove)

def find_founders(fam_pedigree):
    '''
    This function takes fam_pedigree as a networkx (.nx) file and returns a list of all individuals who
    are founders in the pedigree.
    '''
    founders = []
    for indiv in fam_pedigree.nodes():
        indiv_parents = list(fam_pedigree.predecessors(indiv))

        if len(indiv_parents) == 0:
            founders.append(indiv)

    return founders

def find_descendants(fam_pedigree):
    '''
    this fucntion takes fam_pedigree as a networkx (.nx) file and returns a list of all indivuals who 
    has at least 1 parent / not a founder.
    '''
    descendants = []
    for indiv in fam_pedigree.nodes():
        indiv_parents = list(fam_pedigree.predecessors(indiv))

        if len(indiv_parents) > 0:
            descendants.append(indiv)

    return descendants

def relabel_family(fam1, fam2):
    '''
    This fuction takes fam1 as a networkx (.nx) file, finda the indvidual with the highest id number
    and relables fam2 (.nx) ids by adding the highest id from fam1 to each individual in fam2
    '''
    mapping = {}
    last_indiv = np.max(list(map(int, fam1.nodes())))
    for node in fam2.nodes():
        mapping[node]= str(int(node) + (last_indiv))
    fam2_new = nx.relabel_nodes(fam2, mapping)
    
    return(fam2_new)

def set_attributes(family, attributes):
    '''
    This fucntion takes family as a networkx (.nx) file and attributes as a text file and assigns
    the attributes to the nodes bassed off of the attributes given to each individual
    '''
    sex_dict = dict(zip(attributes['ID'].astype(str).to_numpy(), attributes['Sex'].to_numpy()))
    gen_dict = dict(zip(attributes['ID'].astype(str).to_numpy(), attributes['Gen'].to_numpy()))
    nx.set_node_attributes(family, values=sex_dict, name="Sex")
    nx.set_node_attributes(family, values=gen_dict, name="Gen")
  
if __name__== '__main__':
    user_args = load_args()
      
    #Load in pedigrees and apply the atributes 
    main_family = nx.read_edgelist(f'{user_args.family1}', create_using = nx.DiGraph())
    sub_family = nx.read_edgelist(f'{user_args.family2}', create_using = nx.DiGraph())
    profile1 = pd.read_csv(f'{user_args.profiles1}', sep='\t')
    profile2 = pd.read_csv(f'{user_args.profiles2}', sep='\t')
    u_output = user_args.output_prefix
    set_attributes(main_family, profile1)
    set_attributes(sub_family, profile2)
    
    #loads in chosen founders and substitutes
    chosen_founder = user_args.chosen_founder
    chosen_sub = user_args.chosen_sub
    
    #finds founders and chooses one randomly if one is not specified
    main_founders = find_founders(main_family)
    
    #filters out root founders or founders with no generation
    for indiv in main_founders:
        remove_these = []

        if math.isnan(nx.get_node_attributes(main_family, name = 'Gen')[indiv]):
            remove_these.append(indiv)
        main_founders = [i for i in main_founders if i not in remove_these]
    
    #prints out chosen founders 
    founder_gen = int(nx.get_node_attributes(main_family, name = 'Gen')[chosen_founder])
    print('Chosen founder: ',chosen_founder, founder_gen)
    
    #takes sub family, relabels them, makes a list of successors and removes initial node
    new_sub_family = relabel_family(main_family, sub_family)
    
    #list of decendants 
    succ_sub_family = find_descendants(new_sub_family)
    print('All dencendants: ', succ_sub_family)
    
    # a list of all the decendents that fit with our desired pereameters 
    filtered_subs = []
    for i in succ_sub_family:
        sub_gen = int(nx.get_node_attributes(new_sub_family, name = 'Gen')[i])
        if sub_gen == founder_gen:
            filtered_subs.append(i)
    print('filtered subs: ', filtered_subs)
    
    #cases if chosen sub is not defined a sub will be chosen at random
    if chosen_sub == 'empty':
        chosen_sub = np.random.choice(filtered_subs)
    if chosen_sub not in succ_sub_family:
        print(f'individual {chosen_sub} is not a suitable substitute.')
        print('Please choose someone else')
        print(f'possible options: {succ_sub_family}')
        exit()    
    print('Substituted indiv: ', chosen_sub)
    
    #connect the two families onto one graph
    new_fam = nx.compose(main_family , new_sub_family)
  
    nx.write_edgelist(new_fam, 'relabled_fam.nx')
    
    #make connection
    replace_nodes(new_fam, new_sub_family, chosen_founder, chosen_sub)
    
    #creates output file name
    nx.write_edgelist(new_fam, f'{u_output}.nx')
    sex_dict = nx.get_node_attributes(new_fam, name='Sex')
    gen_dict = nx.get_node_attributes(new_fam, name='Gen')
    prof_df = pd.DataFrame(zip(sex_dict.keys(), sex_dict.values(), gen_dict.values()), columns=['ID', 'Sex', 'Gen'])
    prof_df.to_csv(f"{u_output}_profiles.txt", sep='\t', index=False)
    
    #filtering out the largest family
    out_fam = nx.read_edgelist(f'{u_output}.nx')
    largest_cc = max(nx.connected_components(out_fam), key=len)
    print('cc', largest_cc)