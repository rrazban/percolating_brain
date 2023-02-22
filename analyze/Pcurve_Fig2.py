"""
Generate the percolation probability curve from
targeted attack of edges based on sequential order
of length or density.


Creates Figure 2 of the main text. Also responsibe for 
several supplementary figures: Figures S1, S2, S8, S14,
S15, S17.
"""

import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
from scipy.special import lambertw
from scipy.optimize import curve_fit



def beautify_figure(ax, avg_degrees, P_ones, alpha, exp_label):
    plt.title('Increasing Tract ' + r"$\bf{" + exp_label.capitalize() + "}$" + " Targeted Attack", fontsize=16)
    plt.xlabel('average degree $\langle k \\rangle$', fontsize = 14)
    plt.ylabel('probability in the giant cluster $P$', fontsize=14)
    plt.legend(loc = (0.63, 0.58), prop={'size': 12})

    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlim([-1, max(avg_degrees)+1])
    plt.ylim([-0.05, 1.05])

    #have inset to focus on early <k>
    axins = ax.inset_axes([0.5, 0.03, 0.47, 0.47])

    axins.scatter(avg_degrees, P_ones, color='r')
    axins.plot(avg_degrees, theory(avg_degrees, alpha))
    axins.plot(avg_degrees, random_graph(np.array(avg_degrees)))

    x1, x2, y1, y2 = -0.03, 2, -0.03, 0.51
    axins.set_xlim(x1, x2)
    axins.set_ylim(y1, y2)
    axins.set_xticklabels([])
    axins.set_yticklabels([])
    ax.indicate_inset_zoom(axins, edgecolor="black")

    plt.tight_layout()
    plt.show()

#find at which index average degrees (ks) becomes 0
def find_zero_index(ks):
    zero_i = -1
    if 0 in ks:
        zero_i = ks.index(0)
    return zero_i


#generate P curves with subsequent points equally distant from each other
#***make sure sampling from dense points so not too much interpolation
#code from https://stackoverflow.com/questions/19117660/how-to-generate-equispaced-interpolating-values
def sample_equidistant(x, y):    
    xd = np.diff(x)
    yd = np.diff(y)
    dist = np.sqrt(xd**2+yd**2)
    u = np.cumsum(dist)
    u = np.hstack([[0],u])
    t = np.linspace(0,u.max(),50)   #alpha values are the same for 50 or 500
    new_x = np.interp(t, u, x)
    new_y = np.interp(t, u, y)
    return new_x, new_y

#calculate the two graph properties: average degree and 
def get_k_and_P(G):
    degree = [val for (node, val) in G.degree()]

    Gcc = max(nx.connected_components(G), key=len)	
    P_one = len(Gcc)/len(G)
    return np.mean(degree), P_one
    

#perform targeted attack procedure
def break_apart(structure, thresholds):
    print_extra_info = False 

    avg_degrees = []
    P_ones = []

    for lim in thresholds:#reversed(thresholds):       
        structure[structure<lim] = 0	#switch < to > for reversed

        G = nx.from_numpy_array(structure)
        avg_degree, P_one = get_k_and_P(G)

        avg_degrees.append(avg_degree)
        P_ones.append(P_one)

    return np.array(avg_degrees), np.array(P_ones)

#remove from adjacency matrix the rows and columns of regions corresponding to background, white matter regions, or regions with no edges
def preprocess(r):
    print_extra_info = False 

    np.fill_diagonal(r, 0)	#if dont remove diagonals, shifts avg degree to the right!

    r2 = r[~np.all(r == 0, axis=1)]
    adjacency_matrix = r2[:, ~np.all(r2 == 0, axis=0)]

    if print_extra_info:
        print(list(np.where(~r.any(axis=1))[0]))    #indices of excluded regions
        print(r.shape)
        print(r2.shape)
        print(adjacency_matrix.shape)

    return adjacency_matrix


def get_experimental_data(filename):
    thresholds = np.logspace(-1, 4.2, num=200)  #thresholds for the targeted attack procedure
#    thresholds = np.logspace(-.7, 1.5, num=50)  #for mice densities    #    thresholds = np.logspace(3.8, 4.2, num=50)  #for mice distances (units are um)

    adjacency_matrix = pd.read_csv(filename, delimiter=" ", header=None).values 
    adjacency_matrix = preprocess(adjacency_matrix) #remove certain isolated regions (no edges) from adjacency matrix
    avg_degrees, P_ones = break_apart(adjacency_matrix, thresholds)

    return avg_degrees, P_ones


#Equation 1 of main text. Corresponds to well-known result from random graphs.
def random_graph(ks):
    constant = 1    #change critical point location
    p_solns = 1 + lambertw(-( ks/constant) * np.exp(-ks/constant))/(ks/constant)
    p_solns[-1] = 0  #lambertw is undefined at <k>=0
    return p_solns

#Equation 2 of the main text. Corresponds to our Giant Cluster Self Preference theory
def theory(ks, a):
    p_soln = (1 + a/(a-2) * lambertw(-(a-2)/a * np.exp(-ks/a - 1 + 2/a))).real
    p_soln[-1]=0    #lambertw is undefined at <k>=0
    return p_soln


if __name__ == '__main__':
    which = 'length'    #tract length or tract density

    filename = '../sample_outputs/standard/6025360_20250_2_0_{0}.txt'.format(which) #this subject is from the ukb dataset #all subjects in sample_output/ are from ukb

    ori_avg_degrees, ori_P_ones = get_experimental_data(filename)
    avg_degrees, P_ones = sample_equidistant(ori_avg_degrees, ori_P_ones)

    zero_i = find_zero_index(list(avg_degrees))   #solver has trouble with nan values obtained at k=0, remove these values
    popt, pcov = curve_fit(theory, avg_degrees[:zero_i], P_ones[:zero_i])

    fig, ax = plt.subplots()
    plt.plot(avg_degrees, theory(avg_degrees, popt[0]), label = 'theory, $\\alpha$={:.1f}'.format(popt[0]))
    plt.plot(avg_degrees, random_graph(np.array(avg_degrees)), label = 'random graph')
    plt.scatter(avg_degrees, P_ones, color='r', label='human subject')
    beautify_figure(ax, avg_degrees, P_ones, popt[0], which)
