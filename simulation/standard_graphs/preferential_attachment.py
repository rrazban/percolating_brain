"""
Generate the percolation probability curve for
a scale-free graph (probability of edge formation
depends linearly on the multiplication of the 
two nodes' degrees)

"""

import sys
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import random

from random_graph import plotout, d_n_max_k

sys.path.append('../../analyze')
from Pcurve import get_k_and_P



def make_graph(n, k):
    Ps = []
    ks = []

    p = k/n

    G=nx.Graph()
    G.add_nodes_from(range(n))

    choose_from = list(range(n))
    avgdegree=0

    while avgdegree<k: 
        ind1, ind2 = np.random.choice(choose_from, size=2, replace=False)    #have it just choose two at once, no replacement
        potential_edge = ((ind1, ind2))

        if not G.has_edge(*potential_edge) and ind1!=ind2:  #still can get self-edge cuz choose_from reflects node degree
            if random.random() < p:
                G.add_edge(*potential_edge)
	
                avgdegree, P_one = get_k_and_P(G)
                ks.append(avgdegree)
                Ps.append(P_one)

                choose_from.append(ind1)  
                choose_from.append(ind2)

    return np.array(ks), np.array(Ps)




if __name__ == '__main__':
    repeat = 10

    atlas = 'Harvard-Oxford'
#    atlas = 'Talairach'
    n, max_k = d_n_max_k[atlas]

    collection = np.arange(0, int(max_k+0.5), 0.5)  #match range of atlas

    output = [[] for _ in collection]
    for r in range(repeat):
        ks, result = make_graph(n, max_k)
        for i, a0 in enumerate(collection):
            indi = (np.abs(ks-a0).argmin())
            output[i].append(result[indi])

    plotout(collection, output, 'preferential attachment')
