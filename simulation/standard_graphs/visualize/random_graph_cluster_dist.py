"""
Plot distribution of cluster sizes for 
random graphs constructed at various 
average degrees. 

Creates left graph of Figure S4 in the 
Supplement.

"""

import os, sys
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

sys.path.append('../../../analyze')
from Pcurve_Fig2 import get_k_and_P 



def plot_graph():

    n=64
    collection = [16, 4, 1, 0.25]

    for k in collection:
        G = nx.binomial_graph(n, k/n)

        avg_degree, P_one = get_k_and_P(G)

        Gcc = nx.connected_components(G)
        cluster_sizes = [len(x)/n for x in Gcc]   #Gcc is destroyed after this operation

        plt.hist(cluster_sizes, label=k, bins=np.arange(0, n+1, 1)/n, alpha=0.5)


    plt.title(f"Random Graph Simulation", fontsize=16)
    leg = plt.legend(prop={'size':12})
    leg.set_title(title='average degree',prop={'size':12}) 
    plt.xlabel('normalized cluster size', fontsize=14)
    plt.ylabel('frequency', fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    plot_graph()
