"""
Compare percolation probability curve with numerical
results from simulations of general preferential 
attachment models

"""


import sys, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

from Pcurve_simulations import get_experimental_data, sample

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
    k_random, P_random = d_readin['{0} {1}'.format('pref attach k0.5', N)]()
    plt.violinplot(P_random, positions = k_random, showmeans=True, showmedians=True, showextrema=False)

    k_pref, P_pref = d_readin['{0} {1}'.format('pref attach k1', N)]()
    plt.violinplot(P_pref, positions = k_pref, showmeans=True, showmedians=True, showextrema=False)
    plt.plot(0,0)   #dark red color looks to similiar to the red color for experiment
    k_random, P_random = d_readin['{0} {1}'.format('pref attach k1.5', N)]()
    plt.violinplot(P_random, positions = k_random, showmeans=True, showmedians=True, showextrema=False)



    labels = add_label()
#    plt.legend(*zip(*labels), loc= (0.61, 0.52), prop={'size':12})
    plt.legend(*zip(*labels), prop={'size':12})

    plt.xlabel(xlabel, fontsize = 14)   #default size is 10
    plt.ylabel(ylabel, fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlim([-1, max(x)+1])
    plt.ylim([-0.05, 1.05])

    plt.title('Increasing Tract ' + r"$\bf{" + exp_label.capitalize() + "}$" + " Targeted Attack", fontsize=16)#.format(which))  #default size is 12
#    plt.title('Increasing Tract ' + r"$\bf{" + exp_label.capitalize() + "}$" + " Targeted Attack (H-O atlas)", fontsize=14)#.format(which))  #default size is 12
    plt.tight_layout()
    plt.show()


def add_label():
    labels = []

    labels.append((mpatches.Patch(color='#ff7f0e'), 'pref attach, $p_{ij} \propto (k_i k_j)^{1/2}$')) #color exactly set to match violinplot light blue
    labels.append((mpatches.Patch(color='#2ca02c'), 'pref attach, $p_{ij} \propto (k_i k_j)^{1}$'))
    labels.append((mpatches.Patch(color='#9467bd'), 'pref attach, $p_{ij} \propto (k_i k_j)^{3/2}$'))
    labels.append((Line2D([0], [0], marker='o',color='white',markerfacecolor='r', markersize=8), 'human subject'))
    return labels



d_readin = {'pref attach k1 727': saved_outputs.prefattach1_727, 'pref attach k0.5 727': saved_outputs.prefattach05_727,'pref attach k1.5 727': saved_outputs.prefattach15_727,'pref attach k1 64': saved_outputs.prefattach1_64, 'pref attach k0.5 64': saved_outputs.prefattach05_64, 'pref attach k2 64': saved_outputs.prefattach2_64}


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
