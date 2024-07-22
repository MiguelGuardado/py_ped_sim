import numpy as np
import networkx as nx
import argparse
import matplotlib.pyplot as plt
import pandas as pd

"""
Quick start: 
python 
"""

sex_map = {0:'male', 1:'female'}
# zero_event_list = []

def load_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-y', '--years_to_sample', nargs='+', type=int)
    parser.add_argument('-s', '--seed', nargs='+', type=int, default=np.random.randint(100000, size=1))
    parser.add_argument('-g', '--gen_threshold', type=int, default=0)
    parser.add_argument('-c', '--census_filepath', type=str, default='ipumps_sibship_dist.txt')
    parser.add_argument('-o', '--output_prefix', type=str)

    return parser.parse_args()

def sim_fam_recursive(graph, parent1, curGen, finalGenDepth):
    global curPer
    global years_to_sim
    global sibship_dist_df

    if(curGen >= finalGenDepth):
        return()

    # Subset pandas dataframe to give us the row for years' mean + sd to draw come
    curgen_row = sibship_dist_df[sibship_dist_df['Year'] == years_to_sim[curGen - 1]]

    # This will sample from a normal random distribution of the number of kids that will be simulated.
    # Note this will simulate a continuous number that will need to be rounded to a integer.
    # the mean and sd are defined based on the sibship file that the user enters.
    kid_num = np.random.normal(loc=curgen_row['Mean'], scale=curgen_row['SD'], size=1)[0]
    # If the random number sampled is below 0.5 then no kid will be simulated
    if kid_num <= 0.5:
        return()
    else:
        ## Convert and round to the nearest int
        kid_num = int(np.round(np.abs(kid_num)))

    # Assign the ID of parent 2 based on curPer
    parent2 = curPer + 1
    curPer = curPer + 1

    for i in range (0,kid_num):
        child = curPer + 1
        curPer = curPer + 1

        # Connect parent nodes to child. This simulates parent child relationships
        graph.add_edge(parent1, child)
        graph.add_edge(parent2, child)

        # We will check on the sex of the first parent to assign the sex of the second. In addition to assigning the
        # birth year to the child based on the by of the maternal parent.
        if graph.nodes.data()[parent1]['sex'] == 'male':
            nx.set_node_attributes(graph, values={parent2: 'female'}, name='sex')
            # by_child = randrange(by_parent2 + 15, by_parent2 + 35)

        else:
            nx.set_node_attributes(graph, values={parent2: 'male'}, name='sex')
            # by_child = randrange(by_parent1 + 15, by_parent1 + 35)

        # Assign the gen used to simulate the parent and
        nx.set_node_attributes(graph, values={parent2: years_to_sim[curGen - 2]}, name='gen')
        nx.set_node_attributes(graph, values={child: years_to_sim[curGen - 1]}, name='gen')

        # Randomly sample the sex of the child that was simulated
        sex_label = np.random.choice([0, 1], size=1, p=[0.5, 0.5])
        sex = sex_map[sex_label[0]]
        nx.set_node_attributes(graph, values={child: sex}, name='sex')

        # Get right back up there !!
        sim_fam_recursive(graph, child, curGen + 1, finalGenDepth)
    return

if __name__ == '__main__':
    # Load user arguments
    user_args = load_args()

    # Load in years to sample and
    years_to_sim = user_args.years_to_sample
    sibship_dist_df = pd.read_csv(user_args.census_filepath, sep='\t')

    # Set the random seed
    np.random.seed(user_args.seed)

    parent1 = 1
    curPer = 1
    fam_graph = nx.DiGraph()

    # Check to determine if the years imputed are able to be found in the census file
    years_isec = np.intersect1d(sibship_dist_df['Year'].astype(int), years_to_sim)
    if not set(years_to_sim) <= set(sibship_dist_df['Year'].astype(int)):
        print('EXIT: Years inputted does not match the years provided in the census file')
        exit(0)

    # Simulate birth year variation for initial individual in family pedigree simulation.
    # p1_by = randrange(years_to_sim[0] - 10, years_to_sim[0] + 10)

    # Add node and sex attribute for parent 1
    sex_label = np.random.choice([0, 1], size=1, p=[0.5, 0.5])
    sex = sex_map[sex_label[0]]
    fam_graph.add_node(1, sex=sex, gen=years_to_sim[0])

    # Begin simulating a family using recursive functions
    sim_fam_recursive(fam_graph, parent1, 1, len(years_to_sim) + 1)
    print(f'Number individuals simulated: {len(fam_graph.nodes)}, num generations:{len(years_to_sim)+1}')

    # Change the Gen of indivs 1 and 2 to NA since no generation is used to simulate the root founders
    nx.set_node_attributes(fam_graph, values={1:"NA", 2:"NA"}, name='gen')

    # save edge list
    nx.write_edgelist(fam_graph, f"{user_args.output_prefix}.nx")  # save edge list

    # getting the nodes to be saved in a pandas data frame for sex
    sex_dict = nx.get_node_attributes(fam_graph, name='sex')
    # birth_dict = nx.get_node_attributes(fam_graph, name='birth_year')
    gen_dict = nx.get_node_attributes(fam_graph, name='gen')
    prof_df = pd.DataFrame(zip(sex_dict.keys(), sex_dict.values(), gen_dict.values()), columns=['ID', 'Sex', 'Gen'])
    prof_df.to_csv(f"{user_args.output_prefix}_profiles.txt", sep='\t', index=False)