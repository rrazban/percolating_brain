"""
Generate the percolation probability curve based 
on theory that assumes no secondary cluster formation
and rescaled probability (alpha) of edge fromation 
between giant cluster and non-giant cluster nodes 

"""


import sys
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

import networkx as nx
import numpy as np
import random
from scipy.special import lambertw
from scipy.optimize import curve_fit

sys.path.append('../analyze')
from Pcurve import get_k_and_P, sample_equidistant

n = 100 #number of nodes	


def theory(ks, alpha):
	p_soln = (1 + lambertw(-np.exp(-(ks/alpha + 1)))).real
	p_soln[0]=0
	return p_soln

def make_graph(rescale_p_ngc, max_k):
    Ps, ks = [],[]

    G=nx.Graph()
    G.add_nodes_from(range(n))

    ind1,ind2 = np.random.choice(range(n), size=2, replace=False)
    potential_edge = ((ind1, ind2))
    G.add_edge(*potential_edge)

	
    p = max_k/n
    avgdegree = 0
    while avgdegree < max_k:
        Gcc = list(max(nx.connected_components(G), key=len))

        ind1 = random.choice(Gcc)	#one node must always be a part of the giant cluster
        ind2 = random.choice(list(range(n)))	#choose from any node
        potential_edge = ((ind1, ind2))

        if ind2 in Gcc:
            if random.random() < p:
                G.add_edge(*potential_edge)
                avgdegree, P_one = get_k_and_P(G)

                ks.append(avgdegree)
                Ps.append(P_one)
        else:
            if random.random() < p/rescale_p_ngc:
                G.add_edge(*potential_edge)
                avgdegree, P_one = get_k_and_P(G)

                ks.append(avgdegree)
                Ps.append(P_one)

    return np.array(ks), np.array(Ps)


def add_label(alpha, which):    #need to manually add legend labels cuz violinplots not compatible
    labels = []

    if alpha:
        labels.append((Line2D([0], [0], color='#1f77b4') , 'theory, $\\alpha$={:.1f}'.format(alpha)))	
        labels.append((Line2D([0], [0], color='#d62728') , 'simulation, $\\alpha$={:.1f}'.format(alpha)))	
    elif which:
        labels.append((Line2D([0], [0], color='#1f77b4') , 'theory'))	
        labels.append((mpatches.Patch(color='#d62728'), 'model simulation'))
    else:
        labels.append((Line2D([0], [0], color='#1f77b4') , 'theory'))	
        labels.append((Line2D([0], [0], color='#d62728'), 'simulation'))
            
    return labels

def plotout(collection, output, pred_output, alpha, which):
    plt.plot(collection, pred_output, label = 'theory') #label in legend not set here
    plt.plot(0,0)	#cycle through colors to get to red
    plt.plot(0,0)   #hard to set color of violinplot
    plt.violinplot(output, positions = collection, showmeans=True, showmedians=True)#, showextrema=False)

    labels = add_label(alpha, which)    #set legend labels
    plt.legend(*zip(*labels), loc='lower right', prop={'size': 12})
    plt.xlim([-1, max(collection)+1])
    plt.ylim([-0.05, 1.05])
    plt.ylabel('probability in the giant cluster $P$', fontsize = 14 )
    plt.xlabel('average degree $\langle k \\rangle$', fontsize = 14 )
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    if which:
        plt.title('Increasing Tract ' + r"$\bf{" + which.capitalize() + "}$" + " Targeted Attack", fontsize=16)

    plt.tight_layout()
    plt.show()



if __name__ == '__main__':
    repeat = 10#00

    max_k = 20  #maximum average degree of graph
    rescale_p_ngc_gc = 5
    alpha = (rescale_p_ngc_gc)*2 + 1

    collection = np.arange(0, max_k, 0.5)
    output = [[] for _ in collection]


    check_alpha_dist = False
    alphas = [] #to check theory alpha match with fitted alpha
    for r in range(repeat):
        print(r)	
        ks, result = make_graph(rescale_p_ngc_gc, max_k)
        for i, a0 in enumerate(collection): #save Ps closeset to target <k> for accurate error bars
            indi = (np.abs(ks-a0).argmin())
            output[i].append(result[indi])

        if check_alpha_dist:
            avgdegree, Pones = sample_equidistant(ks, result)
            popt, pcov = curve_fit(theory, avgdegree, Pones)
            alphas.append(popt[0])

    if check_alpha_dist:
        plt.hist(alphas)    #check whether theory alpha matches fitted alpha 
        plt.show()

    pred_output = theory(collection, alpha)
    plotout(collection, output, pred_output, alpha, "")
