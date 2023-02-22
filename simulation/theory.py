"""
Generate the percolation probability curve based 
on Giant Cluster Self Preference theory. Assumes 
no secondary cluster formation and rescaled 
probability (alpha) of edge fromation between 
giant cluster and non-giant cluster nodes.

Creates Figure S19 in the Supplement.

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
from Pcurve_Fig2 import get_k_and_P, sample_equidistant
from stochastic_numerics import numerics_theory, transform_p_to_matrix, get_mean_std 


#set the average degree such that all simulations have P=1 at the end
#based on the number of nodes
d_k_alpha11 = {100: 50, 727:100}


#Equation 2 of the main text. Corresponds to our Giant Cluster Self Preference theory
def theory(ks, a):
    p_soln = (1 + a/(a-2) * lambertw(-(a-2)/a * np.exp(-ks/a - 1 + 2/a))).real
    p_soln[0]=0
    return p_soln


def make_graph(rescale_p_ngc, max_k, n):
    Ps, ks = [],[]
    all_nodes = list(range(n)) 

    G=nx.Graph()
    G.add_nodes_from(all_nodes)

    ind1,ind2 = np.random.choice(all_nodes, size=2, replace=False)
    potential_edge = ((ind1, ind2))
    G.add_edge(*potential_edge)

    p = max_k/n
    avgdegree= 2/n
    P_one = 2* 1/n
    Gcc = [ind1, ind2]
    while avgdegree < max_k:

        ind1 = random.choice(Gcc)	#one node must always be a part of the giant cluster
        ind2 = random.choice(all_nodes)

        potential_edge = ((ind1, ind2))

        if ind2 in Gcc:
            if not G.has_edge(*potential_edge) and ind1!=ind2:  #no repeating edges, no self-loops
                if random.random() < p/2:
                    G.add_edge(*potential_edge)  
                    avgdegree+=2/n
        else:
            if random.random() < p/rescale_p_ngc:
                G.add_edge(*potential_edge)

                avgdegree+=2/n
                P_one+=1/n
                Gcc.append(ind2)
                ks.append(avgdegree)
                Ps.append(P_one)

#    print(avgdegree, P_one)
    return G, np.array(ks), np.array(Ps)


def add_label(alpha, which):    #need to manually add legend labels cuz violinplots not compatible
    labels = []

    if which:
        labels.append((Line2D([0], [0], color='#1f77b4') , 'theory, $\\alpha = ${0}'.format(alpha)))	
        labels.append((mpatches.Patch(color='#d62728'), 'model simulation, $\\alpha = ${0}'.format(alpha)))
    else:
        labels.append((Line2D([0], [0], color='#1f77b4') , 'numerical solution'))	
        labels.append((Line2D([0], [0], color='#ff7f0e') , 'analytical equation'.format(alpha)))	
        labels.append((mpatches.Patch(color='#d62728'), 'simulation'.format(alpha)))

    return labels

def plotout(collection, output, pred_output, alpha, which):
    
    plt.plot(collection, pred_output, label = 'theory', zorder=20) #label in legend not set here
    plt.plot(0,0)	#cycle through colors to get to red #hard to set color of violinplot
    if which:
        plt.plot(0,0)   
    plt.violinplot(output, positions = collection, showmeans=True, showmedians=True)#, showextrema=False)

    labels = add_label(alpha, which)    #set legend labels
    plt.legend(*zip(*labels), loc='lower right', prop={'size': 12})
    plt.xlim([-1, max(collection)+1])
    #plt.ylim([-0.05, 1.05])
    plt.ylabel('probability in the giant cluster $P$', fontsize = 14 )
    plt.xlabel('average degree $\langle k \\rangle$', fontsize = 14 )
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.title('$\\alpha = ${0}'.format(alpha), fontsize=14)
    if which:
        plt.title('Increasing Tract ' + r"$\bf{" + which.capitalize() + "}$" + " Targeted Attack", fontsize=16)

    plt.tight_layout()
    plt.show()


def no_stochasticity_numerics_theory(N_edges, a, n):
    Ps = []
    ks = []

    n_e =0
    P=1/n

    ks.append(n_e/n)
    Ps.append(P)
    while n_e<N_edges:
        pre_term = (P* n*(1-P))/( (P* n*(1-P)) + a*P*(n*P-1)/2 - a*n_e/(2*n) )
        P+= pre_term/n 
        Ps.append(P)
        n_e+=2 
        ks.append(n_e/n)

    return ks,Ps



if __name__ == '__main__':
    repeat = 5#00
    
    n = 727 
    alpha = 11
    max_k = d_k_alpha11[n] 

    if n==727:
        collection = np.arange(0, 30, 0.5)
    else:
        collection = np.arange(0, 20, 0.5)


    output = [[] for _ in collection]

    check_alpha_dist = False#True
    alphas = [] #to check theory alpha match with fitted alpha
    for r in range(repeat):
        print(r)	
        _, ks, result = make_graph(alpha, max_k, n)
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
#    ks, pred_output3 = no_stochasticity_numerics_theory(n*max(collection), alpha, n)

    E_tot = int(n*max(collection)/2)#divide by 2 to avoid overcounting
    d_p = numerics_theory(E_tot, alpha, n)
    M_mean, M_var = transform_p_to_matrix(d_p, n, E_tot)
    mean, std = get_mean_std(M_mean, M_var)

    freq = int(E_tot/(len(collection)-1))
    plt.errorbar(collection, mean[::freq], yerr=std[::freq], zorder=10) 

    plotout(collection, output, pred_output, alpha, "")
