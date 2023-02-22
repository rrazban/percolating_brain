"""
Compare percolation probability curve with numerical
results from simulations of small-world graphs with 
different reorganization probabilities and random 
graph constrained to atlas coordinates. Both graph 
types are limited to distance targeted attack.

Creates Figures S9 and S10 in the Supplement.

"""


import sys, os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

from Pcurve_simulations import sample, get_experimental_data

sys.path.append('../')
from Pcurve_Fig2 import preprocess, break_apart 

sys.path.append('../../simulation/standard_graphs/')
import saved_outputs
from random_graph import d_n_max_k


def plotout(x, y, xlabel, ylabel, exp_label):
    fig, ax = plt.subplots()
    plt.scatter(x,y, color='r', label='human subject')

    plt.plot(0,0)   #blue is reserved for theory
    k_roi, P_roi = d_readin['{0} {1}'.format('roi random', N)]()
    plt.violinplot(P_roi, positions = k_roi, showmeans=True, showmedians=True, showextrema=False)
 
    k_pdot1, P_pdot1 = d_readin['{0} {1}'.format('pdot1', N)]()
    plt.violinplot(P_pdot1, positions = k_pdot1, showmeans=True, showmedians=True, showextrema=False)

    plt.plot(0,0)   #dark red color looks to similiar to the red color for experiment
    k_pdot5, P_pdot5 = d_readin['{0} {1}'.format('pdot5', N)]()
    plt.violinplot(P_pdot5, positions = k_pdot5, showmeans=True, showmedians=True, showextrema=False)

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
    axins.violinplot(P_roi, positions = k_roi, showmeans=True, showmedians=True, showextrema=False)
    axins.violinplot(P_pdot1, positions = k_pdot1, showmeans=True, showmedians=True, showextrema=False)
    axins.plot(0,0)
    axins.violinplot(P_pdot5, positions = k_pdot5, showmeans=True, showmedians=True, showextrema=False)

    x1, x2, y1, y2 = -0.03, 2, -0.03, 0.51
    axins.set_xlim(x1, x2)
    axins.set_ylim(y1, y2)
    axins.set_xticklabels([])
    axins.set_yticklabels([])
    ax.indicate_inset_zoom(axins, edgecolor="black")

    plt.tight_layout()
    plt.show()


def add_label():
    labels = []

    labels.append((mpatches.Patch(color='#ff7f0e'), 'random on atlas'))
    labels.append((mpatches.Patch(color='#2ca02c'), 'SWN, p=0.1'))
    labels.append((mpatches.Patch(color='#9467bd'), 'SWN, p=0.5'))
    labels.append((Line2D([0], [0], marker='o',color='white',markerfacecolor='r', markersize=8), 'human subject'))
    return labels


d_readin = {'pdot1 64': saved_outputs.watts_strog_pdot1_64, 'pdot5 64': saved_outputs.watts_strog_pdot5_64, 'pdot1 727': saved_outputs.watts_strog_pdot1_727, 'pdot5 727': saved_outputs.watts_strog_pdot5_727, 'roi random 64': saved_outputs.random_64, 'roi random 727': saved_outputs.random_727} 


if __name__ == '__main__':
    which = 'distance'    #only works for distance because random ROI and SWN limited to distance targeted attack

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
