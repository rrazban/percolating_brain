"""
Calculate small-world properties, such as 
average clustering coefficient and average
path length for our theory and compre to 
random graph results with the same number
of nodes and average degree.

"""


import sys
import matplotlib.pyplot as plt

import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns

from theory import make_graph

#sys.path.append('../analyze')
sys.path.append('/shared/home/rostam/percolating_brain/analyze')
from Pcurve import get_k_and_P, sample_equidistant


def brain_like_network(repeat, max_k, alpha, n):
    clus_coeffs = []
    path_lengths = []
    for r in range(repeat):
        print(r)	
        G, ks, result = make_graph(alpha, max_k, n)

        try:
            clus_coeffs.append(nx.average_clustering(G))
            path_lengths.append(nx.average_shortest_path_length(G))
        except: #dont crash if not fully connected
            print('here')
            pass#if P_one not equal to 1 then error raised for calculating path lenghts
    return clus_coeffs, path_lengths

def random_graph(repeat, max_k, n):
    p = max_k/n
    
    clus_coeffs = []
    path_lengths = []
    for r in range(repeat):
        print(r)	
        G = nx.binomial_graph(n, p) #random graph
        clus_coeffs.append(nx.average_clustering(G))
        path_lengths.append(nx.average_shortest_path_length(G))
    return clus_coeffs, path_lengths


def streamline(output, label, prop, data):
    for x in data:
        output.append([label,prop,x]) 
    return output

def setup_dataframe2(cc, pl, rg_cc, rg_pl, alpha):
    output = []
    output = streamline(output, 'theory, $\\alpha = ${0}'.format(alpha), 'average clustering coefficient', cc)
    output = streamline(output, 'random graph', 'average clustering coefficient', rg_cc)
    return output

def setup_dataframe(cc, pl, rg_cc, rg_pl, alpha):
    output = []
    output = streamline(output, 'theory, $\\alpha = ${0}'.format(alpha), 'average path length', pl)
    output = streamline(output, 'random graph', 'average path length', rg_pl)
    return output

def plotout(df_output):
    sns.set(font_scale=1.25)
    sns.boxplot(x="property", y="value", hue="graph", data=df_output)  #very little variability
    plt.legend()
    plt.xlabel('')
    plt.ylabel('')
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    repeat = 10#00

    n = 100
    max_k = 50  #maximum average degree of graph   
    alpha = 11 

    clus_coeffs, path_lengths = brain_like_network(repeat, max_k, alpha, n)
    rg_clus_coeffs, rg_path_lengths = random_graph(repeat, max_k, n)

    output = setup_dataframe(clus_coeffs, path_lengths, rg_clus_coeffs, rg_path_lengths, alpha)
    plotout(pd.DataFrame(output, columns=['graph', 'property', 'value']))

    output = setup_dataframe2(clus_coeffs, path_lengths, rg_clus_coeffs, rg_path_lengths, alpha)
    plotout(pd.DataFrame(output, columns=['graph', 'property', 'value']))

    #dont have real brain data here cuz the average degrees will be different
    ##clustering coefficient and path length very sensitive to k
