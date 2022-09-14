"""
Compare percolation probability curve with numerical
results from simulations of standard graphs: random 
graph and scale-free graph.

"""


import sys, os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

#sys.path.append('../')
sys.path.append('/shared/home/rostam/percolating_brain/analyze')
from Pcurve import preprocess, break_apart 

#sys.path.append('../../simulation/standard_graphs/')
sys.path.append('/shared/home/rostam/percolating_brain/simulation/standard_graphs')
import saved_outputs
from random_graph import d_n_max_k



def plotout(x, y, xlabel, ylabel, exp_label):
    fig, ax = plt.subplots()
    plt.scatter(x,y, color='r', label='human subject')

    plt.plot(0,0)   #blue is reserved for theory
    k_random, P_random = d_readin['{0} {1}'.format('random graph', N)]()
    plt.violinplot(P_random, positions = k_random, showmeans=True, showmedians=True, showextrema=False)

    k_pref, P_pref = d_readin['{0} {1}'.format('preferential attachment', N)]()
    plt.violinplot(P_pref, positions = k_pref, showmeans=True, showmedians=True, showextrema=False)

    labels = add_label()
    plt.legend(*zip(*labels), loc= (0.61, 0.52), prop={'size':12})

    plt.xlabel(xlabel, fontsize = 14)   #default size is 10
    plt.ylabel(ylabel, fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlim([-1, max(x)+1])
    plt.ylim([-0.05, 1.05])

    plt.title('Increasing Tract ' + r"$\bf{" + exp_label.capitalize() + "}$" + " Targeted Attack", fontsize=16)#.format(which))  #default size is 12
#    plt.title('Increasing Tract ' + r"$\bf{" + exp_label.capitalize() + "}$" + " Targeted Attack (H-O atlas)", fontsize=14)#.format(which))  #default size is 12


    #have inset to focus on early <k>
    axins = ax.inset_axes([0.5, 0.03, 0.47, 0.47])
    axins.scatter(x,y, color='r', label='human subject')

    axins.plot(0,0) 
    axins.violinplot(P_random, positions = k_random, showmeans=True, showmedians=True, showextrema=False)
    axins.violinplot(P_pref, positions = k_pref, showmeans=True, showmedians=True, showextrema=False)

    x1, x2, y1, y2 = -0.03, 2, -0.03, 0.51
    axins.set_xlim(x1, x2)
    axins.set_ylim(y1, y2)
    axins.set_xticklabels([])
    axins.set_yticklabels([])
    ax.indicate_inset_zoom(axins, edgecolor="black")

    plt.tight_layout()
    plt.show()


def get_experimental_data(filename):
    thresholds = np.logspace(-1, 4.2, num=200)

    adjacency_matrix = pd.read_csv(filename, delimiter=" ", header=None).values
    adjacency_matrix = preprocess(adjacency_matrix)

    avg_degrees, P_ones = break_apart(adjacency_matrix, thresholds)

    return avg_degrees, P_ones, len(adjacency_matrix)


def sample(ks, y, collection):

    new_y = []
    for i, a0 in enumerate(collection):
        indi = (np.abs(ks-a0).argmin())
        new_y.append(y[indi])

    return collection, np.array(new_y)


def add_label():
    labels = []

    labels.append((mpatches.Patch(color='#ff7f0e'), 'random graph')) #color exactly set to match violinplot light blue
    labels.append((mpatches.Patch(color='#2ca02c'), 'pref attachment'))
    labels.append((Line2D([0], [0], marker='o',color='white',markerfacecolor='r', markersize=8), 'human subject'))
    return labels



d_readin = {'roi random 727': saved_outputs.roi_random_727,'random graph 64': saved_outputs.random_64, 'random graph 727': saved_outputs.random_727, 'preferential attachment 64': saved_outputs.prefattach_64, 'preferential attachment 727': saved_outputs.prefattach_727, 'roi random 64': saved_outputs.roi_random_64}


if __name__ == '__main__':
    which = 'density'    #tract length or tract density

#    atlas = 'Harvard-Oxford'
    atlas = 'Talairach'
    n, max_k = d_n_max_k[atlas]
    collection = np.arange(0, int(max_k+0.5), 0.5)

    if atlas=='Harvard-Oxford':
        filenames = ['../../sample_outputs/atlas/HarOx_6025360_20250_2_0_{0}.txt'.format(which)]
    elif atlas=='Talairach':
        filenames = ['../../sample_outputs/standard/6025360_20250_2_0_{0}.txt'.format(which)]


    for f,filename in enumerate(filenames):
        print(filename)
        ori_avg_degrees, ori_P_ones, N = get_experimental_data(filename)
        avg_degrees, P_ones = sample(ori_avg_degrees, ori_P_ones, collection) 

        plotout(avg_degrees, P_ones, 'average degree $\langle k \\rangle$', 'probability in the giant cluster $P$', which)
