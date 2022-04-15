"""
Compare alpha distribution of individuals with 
some disease compared to healthy individuals.

"""


import sys, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import ks_2samp 



def group_by_disease(phenotypes, disease):
#diabetes code: https://biobank.ndph.ox.ac.uk/ukb/coding.cgi?id=100291
#bipolar/dep code: https://biobank.ndph.ox.ac.uk/ukb/coding.cgi?id=100695

    ingroup = []
    outgroup = []
    for eid, status in zip(phenotypes['id'], phenotypes[disease]):
        if eid not in d_alpha:  #technically do not need cuz do this when make phenotype file
            continue

        if status>0:    #for diabetes: -1 = do not know, -3 = prefer not to answer  #no one in dataset has been diagnosed with diabetes at birth (status=0)
                        #for bipolar/depression: 0 = none
            ingroup.append(d_alpha[eid])
        else:   #assume empty, as well as -1 and -3 correspond to no disorder
            outgroup.append(d_alpha[eid])

    return ingroup, outgroup


def plotout(ingroup, outgroup, disease, which):
    plt.hist(ingroup, density = True, alpha = 0.3, label='{0} ($N$={1})'.format(disease, len(ingroup)))
    plt.hist(outgroup, density = True, alpha = 0.3, label='healthy ($N$={0})'.format(len(outgroup)))

    plt.title('Increasing Tract ' + r"$\bf{" + which.capitalize() + "}$" + " Targeted Attack (UK Biobank)", fontsize=14)#.format(which))  #default size is 12
    plt.xlabel('modularity $\\alpha$', fontsize=14)
    plt.ylabel('frequency', fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    ks, pval = (ks_2samp(ingroup, outgroup))
    leg = plt.legend(title='K-S = {0:.2f} ({1:.2E})'.format(ks, pval), prop={'size':12})
    leg.set_title(title = 'K-S = {0:.2f} ({1:.2E})'.format(ks, pval), prop={'size':12})
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    which = 'density'
    disease = 'diabetes' #bipolar/depression or diabetes

    alpha_file = 'ukb.csv' #this analysis is only for ukb
    alphas = pd.read_csv(alpha_file)

    if which=='density':
        d_alpha = dict(zip(alphas.id, alphas.alpha_density))
    elif which=='length':
        d_alpha = dict(zip(alphas.id, alphas.alpha_length))


    phenotypes = pd.read_csv('./phenotypes/ukb/phenotypes.csv')
    ingroup, outgroup = group_by_disease(phenotypes, disease)
    plotout(ingroup, outgroup, disease, which)
