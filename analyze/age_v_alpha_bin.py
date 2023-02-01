"""
Plot alphas as a function of ages. Data is binned
for ease in visualization, however, correlations
are calculated by considering data points

"""


import sys, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import spearmanr, binned_statistic


def plotout(xs, ys, which, dataset):
    rho, pval = spearmanr(xs, ys)

    means, bins, binnums = binned_statistic(xs, ys, statistic='mean', bins=5)
    std, bins, binnums = binned_statistic(xs, ys, statistic='std', bins=5)
    bins = np.diff(bins)/2.+bins[:-1]
    plt.errorbar(bins, means, yerr = std, fmt='o-',markersize=8, capsize=6, color='r', label='$\\rho=$ {0:.2f} ({1:.2E})\n$N=$ {2}'.format(rho, pval, len(xs)))

    if dataset=='ukb':
        dataset_title = 'UK Biobank'
    elif dataset=='abcd':
        dataset_title = 'ABCD Study'
    elif dataset=='dhcp':
        dataset_title = 'dHCP'


    plt.title('Increasing Tract ' + r"$\bf{" + which.capitalize() + "}$" + " Targeted Attack ({0})".format(dataset_title), fontsize=14)  #default size is 12
#    plt.title('Increasing Tract ' + r"$\bf{" + which.capitalize() + "}$" + " Targeted Attack (UKB, dti)".format(dataset_title), fontsize=14)  #default size is 12
    plt.ylabel('growth parameter $\\alpha$', fontsize=14)

    if dataset=='dhcp':
        plt.xlabel('gestational age in weeks', fontsize=14)
    else:
        plt.xlabel('age in years', fontsize=14)

    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(prop={'size':12})

    plt.tight_layout()
    plt.show()



if __name__ == '__main__':
    dataset = 'ukb' #ukb, abcd, dhcp
    which = 'density'    #length or density

    alpha_file = 'database_output/{0}_dti.csv'.format(dataset)
    alpha_output = pd.read_csv(alpha_file)
    if which=='density':
        d_alpha = dict(zip(alpha_output.id, alpha_output.alpha_density))
    elif which=='length':
        d_alpha = dict(zip(alpha_output.id, alpha_output.alpha_length)) 

    phenotypes = pd.read_csv('./phenotypes/{0}/phenotypes.csv'.format(dataset))

    ages = []
    alphas = []
    for eid, status in zip(phenotypes['id'], phenotypes['age']):
        if eid not in d_alpha:
            continue

        age = status
        if dataset=='abcd':
            age/=12     #ABCD reports age in months

        ages.append(age)
        alphas.append(d_alpha[eid])
        d_alpha.pop(eid, None)  #seems like we got reps in phenotypes.csv for ABCD

    plotout(ages, alphas, which, dataset)
