"""
Generate the percolation probability curve for
a random graph (any edge as the same probability
of forming)

"""


import sys
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import networkx as nx
import numpy as np
import random
import itertools

sys.path.append('../../analyze')
from Pcurve import get_k_and_P



def add_label(label):    #need to manually add legend labels cuz violinplots not compatible
    labels = []
    labels.append((mpatches.Patch(color='#1f77b4'), label))
    return labels

def plotout(collection, output, label):
    plt.violinplot(output, positions = collection, showmeans=True, showmedians=True)#, showextrema=False)

    labels = add_label(label)    #set legend labels
    plt.legend(*zip(*labels), loc='lower right', prop={'size': 12})
    plt.xlim([-1, max(collection)+1])
    plt.ylim([-0.05, 1.05])
    plt.ylabel('probability in the giant cluster $P$', fontsize = 14 )
    plt.xlabel('average degree $\langle k \\rangle$', fontsize = 14 )
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    plt.tight_layout()
    plt.show()


def make_graph(n, max_k):
    Ps = []
    ks = []

    p = max_k/n

    #copied from gnp_random_graph source code
    #initialize graph
    G=nx.Graph()
    G.add_nodes_from(range(n))
    
    edges=itertools.combinations(range(n),2)
    edges = list(edges)
    random.shuffle(edges)

    for e in edges: #could equivalently make independent graph at certain k
        if random.random() < p:
            G.add_edge(*e)

            avgdegree, P_one = get_k_and_P(G)
            ks.append(avgdegree)
            Ps.append(P_one)

    return np.array(ks), np.array(Ps)

d_n_max_k = {'Harvard-Oxford': (64, 25.53), 'Talairach': (727, 29.93)}

if __name__ == '__main__':
    repeat = 10#00

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

    print(list(collection))
    print(output)   #copy and paste into saved_outputs.py, 1000 repeats
    plotout(collection, output, 'random graph')
