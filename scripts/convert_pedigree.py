#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Recent Update: 10/23/2021
@author: Miguel Guardado Miguel.Guardado@ucsf.edu
"""

import numpy as np
import networkx as nx
import pandas as pd



"""
NOTE THIS SCRIPT IS USED TO CONVERT NETWORKX REPRESENTED PEDIGREES INTO SLIM READABLE PEDIGREES.
"""
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
            Internal function used to load in the network based family pedigree.
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
            if (len(list(self.sub_fam_graph.predecessors(
                    node))) == 0):  # Checks if the node's predecessors is empty, meaning they are a founder
                node_depth = (
                    nx.shortest_path_length(G=self.sub_fam_graph, source=node))  # Calculates the depth of the family
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

        # This will identify all the founders who are from the eldest generation
        eldest_founder = max(self.sub_fams_founders, key=self.sub_fams_founders.get)
        self.sub_fam_num_gens = self.sub_fams_founders[eldest_founder]

        #  Now that we idenitfied the eldest founders, we will find the generations of their descendants creations by
        #  counting the number of directed edges found between them.
        self.family_generation = nx.shortest_path_length(self.sub_fam_graph, eldest_founder)


        # Get a list of all individuals who are connected to the deepest founder + all the founders identified.
        all_indivs_in_fg = np.union1d(list(self.family_generation.keys()), list(self.sub_fams_founders.keys()))

        # Find all individuals who are connected to the main family
        indivs_left = np.setdiff1d(np.array(self.sub_fam_graph.nodes()), all_indivs_in_fg)

        # This begins the search for desc generation time not connected to the deepest connected founder.
        while len(indivs_left) > 0:

            for cur_founder in self.sub_fams_founders:
                founder_desc = nx.shortest_path_length(self.sub_fam_graph, cur_founder)

                for node in founder_desc:
                    if int(node) == int(cur_founder):
                        continue

                    if node not in self.family_generation.keys():
                        # print('node: ', node)

                        # This means we have identified a individual who needs an assigned generation, but it not
                        # connected to the main family (the deepest founder), meaning we need to find the closest known
                        # connection to the founder and identify how many generations they are apart.
                        known_connections = np.intersect1d(list(self.family_generation.keys()), list(founder_desc))

                        # print('known connections: ', known_connections)
                        if len(known_connections) > 0:
                            # Extracts first individual in known connection
                            kc_dict = dict((k, self.family_generation[k]) for k in known_connections if
                                           k in self.family_generation)
                            # Extract the shortest connection in known connection.
                            known_connections = min(kc_dict, key=kc_dict.get)

                            # Find the depth between the current founder and the known connection.
                            cur_founder_depth = nx.shortest_path_length(self.sub_fam_graph, source=cur_founder,
                                                                        target=known_connections)

                            # Difference between the depth (deepest->kc) and (cur_founder->kc)
                            founder_start_gen = self.family_generation[known_connections] - cur_founder_depth

                            #  Difference between the depth of current founder and the node.
                            founder_node_len = nx.shortest_path_length(self.sub_fam_graph, source=cur_founder,
                                                                       target=node)

                            if founder_start_gen < 0:
                                # If the difference is negative, then cur_founder is the node
                                founder_start_gen = np.abs(founder_start_gen)

                                self.family_generation = {k: v + founder_start_gen for k, v in
                                                          self.family_generation.items()}
                                self.family_generation[node] = founder_node_len

                            else:
                                # If the difference is positive, then add the depth between the cur_founder->node + difference
                                self.family_generation[node] = founder_node_len + founder_start_gen
                            # Remove the individual once a generation time is assigned to them
                            indivs_left = np.delete(indivs_left, np.where(indivs_left == node))

            # Recalculated the amount of individuals left after we sort the founders
            all_indivs_in_fg = np.union1d(list(self.family_generation.keys()), list(self.sub_fams_founders.keys()))
            indivs_left = np.setdiff1d(np.array(self.sub_fam_graph.nodes()), all_indivs_in_fg)


    def correct_family_gen(self):

        '''
        This Internal function is used after all descendants have their generation creation time calculated. In the case
        of realistic families, there are often cases in which parent and children are born in the same generation
        (0-30 years). While this biologically makes sense SLiM cannot make sense of it. We will check each individual
        ID and check if they have the same generation creation time as their parents. If they do, we will simply add
        a generation to the child so they are birthed after the parent. We will loop though each individual in
        self.family_generation until no relationships are corrected for.

        :return:
        '''

        num_gen_shift = 1
        # We will run this alg until no shifts are made. Costly operation as sample size increases
        while num_gen_shift > 0:

            num_gen_shift = 0
            for indiv in self.family_generation.keys():
                indiv_pred = list(self.sub_fam_graph.predecessors(indiv))
                indiv = str(indiv)  # Convert this to str to input inside self.family_generation

                # individual has two parents
                if len(indiv_pred) == 2:
                    # Check if first parent is a desc, not a founder
                    if indiv_pred[0] in self.family_generation:
                        if self.family_generation[indiv] <= self.family_generation[indiv_pred[0]]:
                            gen_diff = self.family_generation[indiv_pred[0]] - self.family_generation[indiv] + 1
                            self.family_generation[indiv] = self.family_generation[indiv] + gen_diff
                            num_gen_shift += 1

                    if indiv_pred[1] in self.family_generation:
                        # Check if second parent is a desc, not a founder.
                        if self.family_generation[indiv] <= self.family_generation[indiv_pred[1]]:
                            gen_diff = self.family_generation[indiv_pred[1]] - self.family_generation[indiv] + 1
                            self.family_generation[indiv] = self.family_generation[indiv] + gen_diff
                            num_gen_shift += 1

                # Individual has one parent
                elif len(indiv_pred) == 1:
                    # Check if the single parent is a desc.
                    if indiv_pred[0] in self.family_generation:
                        if self.family_generation[indiv] <= self.family_generation[indiv_pred[0]]:
                            gen_diff = self.family_generation[indiv_pred[0]] - self.family_generation[indiv] + 1
                            self.family_generation[indiv] = self.family_generation[indiv] + gen_diff
                            num_gen_shift += 1
                else:
                    continue


    def run_conversion(self):
        """
        Main function for the class. Where all the magic happens.
        """
        self.load_pedigree()

        # self.ped_slim  and self.founders will hold the information for the slim pedigree / founder id files

        for sub_family in self.sub_fams:
            # Subset the currect family we are looking at
            self.sub_fam_graph = self.ped_dir_edgelist.subgraph(sub_family)

            # For this subfamily we will idenitfy all founders (no predessesors), and calculated the number of
            # generations present in the family.
            self.find_founders()

            # This will be used to create the depth list for all individuals based off the sub fams founders
            self.build_subfam_depth()

            # Now that all generation times are added look into family_generation to correct any generation errors
            self.correct_family_gen()

            for indiv in self.sub_fam_graph.nodes:
                # This will create the slim pedigree files with two known parents for the child
                if len(list(self.sub_fam_graph.predecessors(indiv))) == 2:

                    # This will create the slim pedigree files with two known parents for the child
                    indivs_parents = list(self.sub_fam_graph.predecessors(indiv))
                    cur_mating_line = [self.family_generation[indiv] + 1, indivs_parents[0], indivs_parents[1], indiv]
                    cur_mating_line = list(map(int, cur_mating_line))
                    self.slim_ped.append(cur_mating_line)

                elif len(list(self.sub_fam_graph.predecessors(indiv))) == 1:
                    # This will create the slim pedigree with the known parents for the child.
                    self.num_implicit += 1
                    indivs_parents = list(self.sub_fam_graph.predecessors(indiv))
                    cur_mating_line = [self.family_generation[indiv] + 1, indivs_parents[0], 0, indiv]
                    cur_mating_line = list(map(int, cur_mating_line))

                    self.slim_ped.append(cur_mating_line)

                else:

                    # if there are no parents, then the person must be a founder for the family.
                    self.founders.append(int(indiv))

        self.num_explicit = len(self.founders)

        # Sort slim ped by generation time (1st col) and founders by thier ID
        self.slim_ped.sort()
        self.founders.sort()

        # Add the number of implicit founders to our founder pedigree
        self.founders.append(self.num_implicit)

        # Write and output slim and founders filepath.
        self.slim_filepath = self.output_prefix + "_slim_pedigree.txt"
        self.founder_filepath = self.output_prefix + "_founder_file.txt"

        print(self.slim_filepath)
        np.savetxt(self.slim_filepath, self.slim_ped, fmt='%s', delimiter=" ")
        np.savetxt(self.founder_filepath, self.founders, fmt='%s', delimiter=",")