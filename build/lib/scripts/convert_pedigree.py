#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Recent Update: 10/23/2021
@author: Miguel Guardado Miguel.Guardado@ucsf.edu
"""

import numpy as np
import networkx as nx



class convert_pedigree:
    def __init__(self, ped_filepath='', output_prefix=''):
        self.ped_filepath = ped_filepath
        self.output_prefix = output_prefix
        self.ped_dir_edgelist = []
        self.ped_undir_edgelist = []
        self.sub_fams = []
        self.slim_ped = []
        self.founders = []
        self.num_implicit = 0
        self.run_conversion()

    def load_pedigree(self):
        '''
            Internal helper function used to load in the network based family pedigree.
            We will load in a directed and undirected edgelist. We load in an undirected edgelist to use
            connected_componets() function which will seperate families inside the inputted pedigree. In the case
            multiple families are inputted.

        :return:
        '''
        self.ped_dir_edgelist = nx.read_edgelist(self.ped_filepath, create_using=nx.DiGraph())
        self.ped_undir_edgelist = nx.read_edgelist(self.ped_filepath, create_using=nx.Graph())
        self.sub_fams = list(nx.connected_components(self.ped_undir_edgelist))

    def find_founders(self):
        """
        Internal function to identify founders from a subfamily. This will initialize self.sub_fams_founders and the
        self.sub_fam_num_gens variables to record.
        self.sub_fams_founders - dictionary of all founders, and thier deepest connection.
        seld.sub_fams_num_gens - integer value of the amount of generations this family contains.

        :return:
        """

        sub_fam_founders = {}  # Creates an empty dictionary to add founder nodes and their parent too.
        for node in self.sub_fam_graph.nodes:  # Iterates though all the nodes inside the subfamily
            if (len(list(self.sub_fam_graph.predecessors(node))) == 0):  # Checks if the node's predecessors is empty, meaning they are a founder
                node_depth = (nx.shortest_path_length(G=self.sub_fam_graph, source=node))  # Calculates the depth of the family
                keys = list(node_depth.keys())  # gets an array of all the key values
                sub_fam_founders[node] = node_depth[keys[-1]]

        self.sub_fams_founders = sub_fam_founders
        self.sub_fam_num_gens = max(sub_fam_founders.values())


    def build_subfam_depth(self):
        '''
        Internal Function used to find the generation of all the individuals inside the family pedigree.
        First we will start off identifying the eldest founders inside the family, and attach the generation to all the
        descendant. After we initialize all descendants from the eldest founders, we then will search though the
        remaining founders and find the generation of creation for remaining descendants. Since the family networkx is
        connected, each founder will have a descendants attached to the eldest family branch, where we can initialize the
        generation of creation for each descendants. (see publication for more details.)

        :return:
        '''

        #This will identify all the founders who are from the eldest generation
        eldest_founder = []
        for eldest_founder_id in self.sub_fams_founders:
            if (self.sub_fams_founders[eldest_founder_id] == self.sub_fam_num_gens):
                eldest_founder.append(eldest_founder_id)
        self.sub_fam_num_gens = (max(self.sub_fams_founders.values()))


        #  Now that we idenitfied the eldest founders, we will find the generations of their descendants creations by
        #  counting the number of directed edges found between them.
        self.family_generation = dict()
        for eldest_founder_id in eldest_founder:
            eldest_founder_desc = nx.shortest_path_length(self.sub_fam_graph, eldest_founder_id)
            for i in eldest_founder_desc:
                if (i not in self.family_generation.keys()):
                    self.family_generation[i] = eldest_founder_desc[i]
                else:
                    if self.family_generation[i] < eldest_founder_desc[i]:
                        self.family_generation[i] = eldest_founder_desc[i]
            #Remove each of the eldest founders root node.
            self.sub_fams_founders.pop(eldest_founder_id)


        #Now that we have found the generation time of all direct desendants of the eldest family memebers we will
        # look across the rest of the founders, and identify any connection to the root founders.
        for founder_id in self.sub_fams_founders:
            founder_desc = nx.shortest_path_length(self.sub_fam_graph, founder_id)

            for node in founder_desc:
                if int(node) == int(founder_id):
                    continue

                if node not in self.family_generation.keys():
                    # This means we have identified a individual who needs an assigned generation, but it not
                    # connected to the main family(eldest founder), meaning we need to find the closest known
                    # connection to the founder and identify how many generations they are apart.
                    known_connections = np.intersect1d(list(self.family_generation.keys()), list(founder_desc))[0]

                    known_connect_length = nx.shortest_path_length(self.sub_fam_graph, source=founder_id,
                                                                   target=known_connections)
                    founder_start_gen = known_connect_length - self.family_generation[known_connections]

                    founder_node_len = nx.shortest_path_length(self.sub_fam_graph, source=founder_id, target=node)

                    self.family_generation[node] = founder_start_gen + founder_node_len


    def run_conversion(self):
        """
        Main function for the class. Where all the magic happens.
        """
        self.load_pedigree()

        #self.ped_slim  and self.founders will hold the information for the slim pedigree / founder id files

        for sub_family in self.sub_fams:
            #Subset the currect family we are looking at
            self.sub_fam_graph = self.ped_dir_edgelist.subgraph(sub_family)

            #for this subfamily we will idenitfy all founders (no predessesors), and calulate the number of
            #generations present in the family.
            self.find_founders()

            #This will be used to create the depth list for all individuals based off the sub fams founders
            self.build_subfam_depth()

            for indiv in self.sub_fam_graph.nodes:
                # This will create the slim pedigree files with two known parents for the child
                if(len(list(self.sub_fam_graph.predecessors(indiv)))==2):

                    #This will create the slim pedigree files with two known parents for the child
                    indivs_parents = list(self.sub_fam_graph.predecessors(indiv))
                    cur_mating_line = [self.family_generation[indiv] + 1, indivs_parents[0], indivs_parents[1], indiv]
                    cur_mating_line = list(map(int, cur_mating_line))
                    self.slim_ped.append(cur_mating_line)

                elif(len(list(self.sub_fam_graph.predecessors(indiv)))== 1):

                    #This will create the slim pedigree with the known parents for the child.
                    self.num_implicit+=1
                    indivs_parents = list(self.sub_fam_graph.predecessors(indiv))
                    cur_mating_line = [self.family_generation[indiv] + 1, indivs_parents[0], 0, indiv]
                    cur_mating_line = list(map(int, cur_mating_line))
                    self.slim_ped.append(cur_mating_line)

                else:

                    #if there are no parents, then the person must be a founder for the family.
                    self.founders.append(int(indiv))

        self.num_explicit = len(self.founders)

        self.slim_ped.sort()
        self.founders.sort()

        self.founders.append(self.num_implicit)

        self.slim_filepath= self.output_prefix + "_slim_pedigree.txt"
        self.founder_filepath= self.output_prefix + "_founder_file.txt"


        np.savetxt(self.slim_filepath, self.slim_ped, fmt='%s', delimiter=" ")
        np.savetxt(self.founder_filepath, self.founders, fmt='%s', delimiter=",")


