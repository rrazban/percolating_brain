"""
Compare tract density vs tract distance. 

"""


import sys, os
import numpy as np
import matplotlib.pyplot as plt
import glob
import pandas as pd
from scipy.stats import spearmanr, pearsonr

sys.path.append('../')
from Pcurve import preprocess


def plotout(densities, dists):
    rho, pval = spearmanr(densities, dists)
    plt.scatter(dists, densities, label = '$\\rho=$ {0:.2f} ({1:.2E})\n$N=${2}'.format(rho, pval, len(densities)))

    plt.xlabel('tract distance (mm)', fontsize=14)
    plt.ylabel('tract density', fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.title('Human Subject', fontsize=16)
    plt.legend(prop={'size':12}, loc='upper left')

#   print(pearsonr(density, dist))
    print(np.poly1d(np.polyfit(densities, dists, 1)))

    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    dataset = 'ukb'    #abcd or ukb

    density_fnames = ['../../sample_outputs/standard/6025360_20250_2_0_{0}.txt'.format('density')]

    for density_fname in density_fnames:
        distance_fname = "{0}_distance.txt".format(density_fname[:-12])
        M_density = pd.read_csv(density_fname, delimiter=" ", header=None).values
        M_dist = pd.read_csv(distance_fname, delimiter=" ", header=None).values

        M_density = preprocess(M_density)
        M_dist = preprocess(M_dist)

        densities = M_density[M_density!=0]
        distances = M_dist[M_dist!=0]

        plotout(densities, distances)
