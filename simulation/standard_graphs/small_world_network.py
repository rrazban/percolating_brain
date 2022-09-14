"""
Generate the percolation probability curve for
a small-world network (implementation from 
Python's networkx package).

"""
#note for some reason watts_strogatz_graph does not work with current decorator
#if i use previous, pip install decorator==5.0.9, it works
#https://stackoverflow.com/questions/66920533/networkx-shows-random-state-index-is-incorrect

import sys
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

from random_graph import plotout, d_n_max_k

sys.path.append('../../analyze')
from Pcurve import break_apart 

pstrog = 0.5 


def make_graph(collection, n):
    Ps = []
    ks = []

    G = nx.watts_strogatz_graph(n, int(collection[-1]), pstrog)
    coords = nx.circular_layout(G) #get coordinates to get distances
    edges = list(G.edges())

    distances = np.zeros((n, n))
    for i, j in edges:
        dist = np.linalg.norm(np.array(coords[i]) - np.array(coords[j]))	#probly faster to pick out terms rather than recalculate distance over and over	
        distances[i][j] = dist


    min_ = np.min(np.log10(distances[distances!=0].flatten()))
    thresholds = np.logspace(min_, np.log10(np.max(distances)), 500) #checked
    avg_degrees, P_ones = break_apart(distances, thresholds)

    return np.array(avg_degrees), np.array(P_ones)


if __name__ == '__main__':
    repeat = 10#00

    atlas = 'Harvard-Oxford'
#    atlas = 'Talairach'
    n, max_k = d_n_max_k[atlas]

    collection = np.arange(0, int(max_k+0.5) + 0.5, 0.5)  #final k must be an integer
    output = [[] for _ in collection]

    for r in range(repeat):
        ks, result = make_graph(collection, n)
        for i, a0 in enumerate(collection): #no need since input 'collection', but keep the same format as other graphs
            indi = (np.abs(ks-a0).argmin())
            output[i].append(result[indi])

    plotout(collection, output, 'small-world network, p={0}'.format(pstrog))
