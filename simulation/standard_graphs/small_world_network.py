"""
Generate the percolation probability curve for
a small-world network (implementation from 
Python's networkx package).

"""


import sys
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

from random_graph import plotout

sys.path.append('../../analyze')
from Pcurve import get_k_and_P


n = 727 #adjust to match number of nodes in parcellation 
pstrog = 0.1 


def make_graph(collection):
    Ps = []
    ks = []

    for k in collection:
        G = nx.watts_strogatz_graph(n, k, pstrog)

        try:
            avgdegree, P_one = get_k_and_P(G)
            ks.append(avgdegree)
            Ps.append(P_one)
        except: #k=0 and k=1 go here
            ks.append(k)
            Ps.append(0)
            
    return np.array(ks), np.array(Ps)


if __name__ == '__main__':
    repeat = 1000

    collection = np.arange(0, 8, 1)
    output = [[] for _ in collection]

    for r in range(repeat):
        ks, result = make_graph(collection)
        for i, a0 in enumerate(collection): #no need since input 'collection', but keep the same format as other graphs
            indi = (np.abs(ks-a0).argmin())
            output[i].append(result[indi])

    plotout(collection, output, 'small-world network, p={0}'.format(pstrog))


