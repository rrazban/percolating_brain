"""
Check that the alpha vs age relationship is not
biased by the final average degree (kf) of the 
network.

Creates the bottom graphs of Figure S18 in the 
Supplement.

"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import binned_statistic


def plotout(xs, ys, which, dataset, thresh):

    set_bins = [45.0, 52.0, 59.0, 66.0, 73.0, 80.0] #taken from those determined by age_v_alpha_bin.py 
    means, bins, binnums = binned_statistic(xs, ys, statistic='mean', bins=set_bins)
    std, bins, binnums = binned_statistic(xs, ys, statistic='std', bins=set_bins)
    bins = np.diff(bins)/2.+bins[:-1]

    plt.errorbar(bins, means, yerr = std, fmt='o-',markersize=8, capsize=6, label='{0} ($N=$ {1})'.format(thresh, len(xs)))

    plt.ylabel('growth parameter $\\alpha$', fontsize=14)
    plt.xlabel('age in years', fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    dataset_title = 'UK Biobank'
    plt.title('Increasing Tract ' + r"$\bf{" + which.capitalize() + "}$" + " Targeted Attack".format(dataset_title), fontsize=16)  #default size is 12

if __name__ == '__main__':
    which = 'length'    #length or density

    for thresh in [(25, 27), (27,29), (29,31)]: 
        dataset='ukb'   #only implemented for UK Biobank, so far

        alpha_file = 'database_output/{0}.csv'.format(dataset)
        alpha_output = pd.read_csv(alpha_file)
        if which=='density':
            d_alpha = dict(zip(alpha_output.id, alpha_output.alpha_density))
        elif which=='length':
            d_alpha = dict(zip(alpha_output.id, alpha_output.alpha_length)) 

        basic_file = 'database_output/{0}_basic.csv'.format(dataset)
        basic_output = pd.read_csv(basic_file)

        match_first_criteria = basic_output[basic_output.kf<=thresh[1]]
        eids = list(match_first_criteria[match_first_criteria.kf>thresh[0]].id) #second criteria

        phenotypes = pd.read_csv('./phenotypes/{0}/phenotypes.csv'.format(dataset))
        ages = []
        alphas = []
        for eid, status in zip(phenotypes['id'], phenotypes['age']):
            if eid not in d_alpha or eid not in eids:
                continue

            ages.append(status)
            alphas.append(d_alpha[eid])
            d_alpha.pop(eid, None)  #seems like we got reps in phenotypes.csv for ABCD

        plotout(ages, alphas, which, dataset, thresh)

    legend = plt.legend(prop={'size':12}, title='average degree')
    legend.get_title().set_fontsize('14')

    plt.tight_layout()
    plt.show()
