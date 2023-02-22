"""
Generate the percolation probability curve for
a general scale-free graph where the probability 
of edge formation depends on the multiplication 
of the two nodes' degrees to some power x

"""

import sys
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import random

from random_graph import plotout, d_n_max_k

sys.path.append('../../analyze')
from Pcurve_Fig2 import get_k_and_P



def make_graph(n, k, exp):
        Ps = []
        ks = []

        p = k/n

        G=nx.Graph()
        G.add_nodes_from(range(n))

        choose_from = list(range(n))
        avgdegree=0

        degrees = np.zeros(n)
        while avgdegree<k: 
            weights = np.array([degrees[ind]**(exp) for ind in choose_from]) + 1	#add one to make sure all have at least one 
            normalized_weights = weights/np.sum(weights)	#random choice does not do the normalization
            ind1, ind2 = np.random.choice(choose_from, size=2, p=normalized_weights, replace=False)    #have it just choose two at once, no replacement

            potential_edge = ((ind1, ind2)) 
            if not G.has_edge(*potential_edge):# and ind1!=ind2:
                if random.random() < p:
                    G.add_edge(*potential_edge)
                    degrees[ind1]+=1
                    degrees[ind2]+=1

                    avgdegree, P_one = get_k_and_P(G)

                    ks.append(avgdegree)
                    Ps.append(P_one)

        return np.array(ks), np.array(Ps)


if __name__ == '__main__':
    repeat = 10#00
    exp = 1	#set to one to have edge addition scale linearly in k
                #exp=1 reproduces preferential_attachment.py
            #N=64 can do 2, but N=727 exp=2 is too large;choosing the same edge over and over!

    atlas = 'Harvard-Oxford'
#    atlas = 'Talairach'
    n, max_k = d_n_max_k[atlas]

    collection = np.arange(0, int(max_k+0.5), 0.5)  #match range of atlas

    output = [[] for _ in collection]
    for r in range(repeat):
        ks, result = make_graph(n, max_k, exp)
        for i, a0 in enumerate(collection):
            indi = (np.abs(ks-a0).argmin())
            output[i].append(result[indi])

    plotout(collection, output, 'preferential attachment')
