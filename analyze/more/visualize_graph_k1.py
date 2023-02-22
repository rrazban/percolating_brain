"""
Visualize graph (nodes and edges) for targeted attack 
at an average degree of 1. 

Creates right graph of Figure S3 in the 
Supplement.

"""

import os, sys
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import numpy as np

sys.path.append('../')
from Pcurve_Fig2 import preprocess, get_k_and_P 


layout = nx.spring_layout


def plot_graph(structure, thresholds, which):

    print('Number of nodes: {0}'.format(len(structure)))
    for lim in thresholds:
        structure[structure<=lim] = 0	
        G = nx.from_numpy_array(structure)

        avg_degree, P_one = get_k_and_P(G)
        if avg_degree <= 1:
            print("Average degree: {0}".format(avg_degree))

            G = nx.convert_matrix.from_numpy_array(structure)
            pos = layout(G, iterations=200) #increasing iterations helps it find a better visual sometimes
            plt.title('Increasing Tract ' + r"$\bf{" + which.capitalize() + "}$" + " Targeted Attack", fontsize=16)  #default size is 12

            R = G.copy()
            R.remove_nodes_from(list(nx.isolates(R)))	#remove orphans from visualization
            nx.draw(R, pos, with_labels=False, node_size=10)

            # identify largest connected component
            Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
            G0 = G.subgraph(Gcc[0])
            nx.draw_networkx_edges(G0, pos, edge_color="r", width=6.0)

            # show other connected components
            for Gi in Gcc[1:]:
                if len(Gi) > 1:
                    nx.draw_networkx_edges(
                    G.subgraph(Gi), pos, edge_color="r", alpha=0.3, width=5.0)

            plt.tight_layout()
            plt.show()
            break


if __name__ == '__main__':
    which = 'density'    #tract length or tract density
    dataset = 'ukb'    #abcd or ukb

    filenames = ['../../sample_outputs/atlas/HarOx_6025360_20250_2_0_{0}.txt'.format(which)]

    for f,filename in enumerate(filenames):
        thresholds = np.logspace(-1, 4.2, num=200)
        adjacency_matrix = pd.read_csv(filename, delimiter=" ", header=None).values 
        adjacency_matrix = preprocess(adjacency_matrix)
        plot_graph(adjacency_matrix, thresholds, which)
