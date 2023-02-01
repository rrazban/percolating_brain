"""
Plot alphas as a function of ages. Data is binned
for ease in visualization, however, correlations
are calculated by considering all data points.
All three databases are shown in one figure.

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import spearmanr, binned_statistic


def plotout(xs, ys, which, dataset, ax):
    rho, pval = spearmanr(xs, ys)

    means, bins, binnums = binned_statistic(xs, ys, statistic='mean', bins=5)
    std, bins, binnums = binned_statistic(xs, ys, statistic='std', bins=5)
    bins = np.diff(bins)/2.+bins[:-1]
    ax.errorbar(bins, means, yerr = std, fmt='o-',markersize=8, capsize=6, color='r', label='$\\rho=$ {0:.2f} ({1:.2E})\n$N=$ {2}'.format(rho, pval, len(xs)))


    ax.tick_params(labelsize=14)
    if dataset=='ukb':
        dataset_title = 'UK Biobank'
        ax.tick_params(labelsize=0, axis='y')
    elif dataset=='abcd':
        dataset_title = 'ABCD Study'
        ax.tick_params(labelsize=0, axis='y')
    elif dataset=='dhcp':
        dataset_title = 'dHCP'
        ax.set_ylabel('growth parameter $\\alpha$', fontsize=14)

    ax.set_title(dataset_title, fontsize=14)
 

    if dataset=='dhcp':
        ax.set_xlabel('gestational age in weeks', fontsize=14)
    else:
        ax.set_xlabel('age in years', fontsize=14)

    if which=='length':
        ax.set_ylim((10.6, 17))
    elif which=='density':
        ax.set_ylim((9.6, 13.75))

    ax.legend(prop={'size':12})



if __name__ == '__main__':
    which = 'length'    #length or density

    fig, axes = plt.subplots(1, 3, figsize=(10,4.4))

    for shift, dataset in enumerate(['dhcp', 'abcd', 'ukb']):
        alpha_file = 'database_output/{0}.csv'.format(dataset)
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

        plotout(ages, alphas, which, dataset,axes[shift])
    plt.suptitle('Increasing Tract ' + r"$\bf{" + which.capitalize() + "}$" + " Targeted Attack", fontsize=16)  #default size is 12
    plt.tight_layout()
    plt.show()
