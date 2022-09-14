"""
Plot distribution of cluster sizes for 
for increating targeted attack at 
various average degrees. 

"""

import os, sys
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import numpy as np

sys.path.append('../')
from Pcurve import preprocess, get_k_and_P 



def plot_graph(structure, thresholds, which):

   # collection = [16, 4, 2, 1, 0.5]
    collection = [16, 4, 1, 0.25]

    print('Number of nodes: {0}'.format(len(structure)))
    n = len(structure) 
    indi = 0
    for lim in thresholds:
        structure[structure<=lim] = 0	
        G = nx.from_numpy_matrix(structure)

        avg_degree, P_one = get_k_and_P(G)

        if avg_degree <= collection[indi]:
            Gcc = nx.connected_components(G)
            cluster_sizes = [len(x)/n for x in Gcc]   #Gcc is destroyed after this operation

           # plt.hist(np.log(cluster_sizes), label=collection[indi], bins=n+1)
 #           weights = np.zeros(len(cluster_sizes))
#            weights.fill(1/n)
            plt.hist(cluster_sizes, label=collection[indi], bins=np.arange(0, n+1, 1)/n, alpha=0.5)#, weights=weights)


            indi+=1
            if indi==len(collection):
                break

    plt.title('Increasing Tract ' + r"$\bf{" + which.capitalize() + "}$" + " Targeted Attack", fontsize=16)  #default size is 12
    leg = plt.legend(prop={'size':12})
    leg.set_title(title='average degree',prop={'size':12}) 
    plt.xlabel('normalized cluster size', fontsize=14)
    plt.ylabel('frequency', fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    which = 'density'    #tract length or tract density
    dataset = 'ukb'    #abcd or ukb

    filenames = ['../../sample_outputs/atlas/HarOx_6025360_20250_2_0_{0}.txt'.format(which)]

    for f,filename in enumerate(filenames):
        thresholds = np.logspace(-1, 4.2, num=200)
        adjacency_matrix = pd.read_csv(filename, delimiter=" ", header=None).values 
        adjacency_matrix = preprocess(adjacency_matrix)
        plot_graph(adjacency_matrix, thresholds, which)
