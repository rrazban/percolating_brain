"""
Generate the percolation probability curve from
targeted attack of edges based on sequential order
of length or density.

"""


import sys, os
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import glob
import pandas as pd
from scipy.special import lambertw
from scipy.optimize import curve_fit



def beautify_figure(avg_degrees, P_ones, alpha, exp_label):
    plt.title('Increasing Tract ' + r"$\bf{" + exp_label.capitalize() + "}$" + " Targeted Attack", fontsize=16)  #default size is 12
    plt.xlabel('average degree $\langle k \\rangle$', fontsize = 14)   #default size is 10
    plt.ylabel('probability in the giant cluster $P$', fontsize=14)
    plt.legend(loc = (0.63, 0.58), prop={'size': 12})    #default is 10

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


def find_zero_index(ks):
    zero_i = -1
    if 0 in ks:
        zero_i = ks.index(0)
    return zero_i

def sample_equidistant(x, y):    #***make sure sampling from dense points so not too much interpolation
#code from https://stackoverflow.com/questions/19117660/how-to-generate-equispaced-interpolating-values
    xd = np.diff(x)
    yd = np.diff(y)
    dist = np.sqrt(xd**2+yd**2)
    u = np.cumsum(dist)
    u = np.hstack([[0],u])
    t = np.linspace(0,u.max(),50)   #alpha values are the same for 50 or 500
    new_x = np.interp(t, u, x)
    new_y = np.interp(t, u, y)
    return new_x, new_y


def get_k_and_P(G):
    degree = [val for (node, val) in G.degree()]

    Gcc = max(nx.connected_components(G), key=len)	
    P_one = len(Gcc)/len(G)

    return np.mean(degree), P_one
    
def break_apart(structure, thresholds):
    print_extra_info = False 

    avg_degrees = []
    P_ones = []

    for lim in thresholds:#reversed(thresholds):       
        structure[structure<lim] = 0	#switch < to > for reversed

        G = nx.from_numpy_matrix(structure)
        avg_degree, P_one = get_k_and_P(G)

        avg_degrees.append(avg_degree)
        P_ones.append(P_one)

        if lim==thresholds[0] and print_extra_info:
            print('average degree: {0:.2f}'.format(np.mean([deg for (node, deg) in G.degree()])))
            print('average clustering coefficient: {0:.2f}'.format(nx.average_clustering(G)))
            print('average shortest path length: {0}'.format(nx.average_shortest_path_length(G)))

        #   plt.hist(structure[np.nonzero(structure)])  #check out degree distribution
         #  plt.show()
     
    return np.array(avg_degrees), np.array(P_ones)

def preprocess(r):
    print_extra_info = False

    np.fill_diagonal(r, 0)	#if dont remove diagonals, shifts avg degree to the right!

    #remove background and white matter regions
    r2 = r[~np.all(r == 0, axis=1)]
    adjacency_matrix = r2[:, ~np.all(r2 == 0, axis=0)]

    if print_extra_info:
        print(list(np.where(~r.any(axis=1))[0]))    #indices of excluded regions
        print(r.shape)
        print(r2.shape)
        print(adjacency_matrix.shape)

    return adjacency_matrix

def get_experimental_data(filename):
    thresholds = np.logspace(-1, 4.2, num=200)
 #   thresholds = np.logspace(-.7, 1.5, num=50)  #for mice densities    #np.logspace(3.8, 4.2, num=50)  #for mice distances (units are um)

    adjacency_matrix = pd.read_csv(filename, delimiter=" ", header=None).values  #faster than load.txt
    adjacency_matrix = preprocess(adjacency_matrix)
    avg_degrees, P_ones = break_apart(adjacency_matrix, thresholds)

    return avg_degrees, P_ones


def random_graph(ks):
    constant = 1    #change critical point location
    p_solns = 1 + lambertw(-( ks/constant) * np.exp(-ks/constant))/(ks/constant) 
    p_solns[-1] = 0  #lambertw is undefined at <k>=0
    return p_solns

def theory(ks, alpha):
    p_solns = (1 + lambertw(-( 1) * np.exp(-(ks/alpha + 1)))).real
   # p_solns = 1 + lambertw(-(ks**2/2 + ks + 1) * np.exp(-(ks + 1))).real    #correct for coarser parcellation (smaller N)
    p_solns[-1] = 0  #lambertw is undefinedi at <k>=0
    return p_solns 



if __name__ == '__main__':
    which = 'distance'    #tract length or tract density
    dataset = 'ukb'    #abcd or ukb

    filenames = ['../sample_outputs/standard/6025360_20250_2_0_{0}.txt'.format(which)]
#   filenames = glob.glob('/shared/datasets/public/{0}/derivatives/*_conn.txt'.format(dataset))[:1] #fix conn name convention

    for f,filename in enumerate(filenames):
        print(filename)
        ori_avg_degrees, ori_P_ones = get_experimental_data(filename)
        avg_degrees, P_ones = sample_equidistant(ori_avg_degrees, ori_P_ones)

        zero_i = find_zero_index(list(avg_degrees))   #solver has trouble with nan values obtained at k=0
        popt, pcov = curve_fit(theory, avg_degrees[:zero_i], P_ones[:zero_i])

        fig, ax = plt.subplots()
        plt.plot(avg_degrees, theory(avg_degrees, popt[0]), label = 'theory, $\\alpha$={:.1f}'.format(popt[0]))
        plt.plot(avg_degrees, random_graph(np.array(avg_degrees)), label = 'random graph')
        plt.scatter(avg_degrees, P_ones, color='r', label='human subject')
        beautify_figure(avg_degrees, P_ones, popt[0], which)
