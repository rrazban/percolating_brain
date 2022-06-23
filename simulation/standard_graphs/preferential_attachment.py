"""
Generate the percolation probability curve for
a scale-free graph (probability of edge formation
depends on node degree)

"""

import sys
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import random
import itertools

from random_graph import plotout

sys.path.append('../../analyze')
from Pcurve import get_k_and_P


n = 727 


def setup_edges(n):
    edges=itertools.combinations(range(n),2)
    edges = list(edges)

    edge_per_node = [[] for _ in range(n)]
    for e in edges:
        ind1=e[0]
        ind2=e[1]
        edge_per_node[ind1].append(e)
        edge_per_node[ind2].append((ind2, ind1))
    return edge_per_node


def make_graph():
    Ps = []
    ks = []

    k = 14
    p = k/n

    G=nx.Graph()
    G.add_nodes_from(range(n))
    edge_per_node = setup_edges(n)

    choose_from = list(range(n))

    while any(edge_per_node):
        ind1 = random.choice(choose_from)
        potential_edges = edge_per_node[ind1]
        potential_ind2s = [edge[1] for edge in potential_edges]

        weights = G.degree()

        weighted_potential_ind2s = []
        for ind2 in potential_ind2s:
            for _ in range(weights[ind2]):
                weighted_potential_ind2s.append(ind2)	
            weighted_potential_ind2s.append(ind2)	#make sure have at least one	

        ind2 = random.choice(weighted_potential_ind2s)
        potential_edge = ((ind1, ind2))

        edge_per_node[ind1].remove(potential_edge)
        edge_per_node[ind2].remove((ind2, ind1))

        if random.random() < p:
            G.add_edge(*potential_edge)
	
            avgdegree, P_one = get_k_and_P(G)
            ks.append(avgdegree)
            Ps.append(P_one)

            choose_from.append(ind1)
            choose_from.append(ind2)

        if len(edge_per_node[ind1])==0:
            choose_from = [x for x in choose_from if x!=ind1]	#cant just do a simple remove cuz multiple instances
        if len(edge_per_node[ind2])==0:
            choose_from = [x for x in choose_from if x!=ind2]

    return np.array(ks), np.array(Ps)




if __name__ == '__main__':
    repeat = 10#00

    collection = np.arange(0, 14, 0.5)
    output = [[] for _ in collection]

    for r in range(repeat):
        ks, result = make_graph()
        for i, a0 in enumerate(collection):
            indi = (np.abs(ks-a0).argmin())
            output[i].append(result[indi])

    plotout(collection, output, 'preferential attachment')
