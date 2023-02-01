"""
Compare alpha distribution of individuals with 
some categorical value compared to another 
categorical value for UK Biobank. Categories include:
    female vs male
    diabetes vs healthy
    bipoloar/depression vs healthy
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import ks_2samp 



def group_by_category(phenotypes, cat):
#diabetes code: https://biobank.ndph.ox.ac.uk/ukb/coding.cgi?id=100291
#bipolar/dep code: https://biobank.ndph.ox.ac.uk/ukb/coding.cgi?id=100695

    group1 = []
    group2 = []
    for eid, status in zip(phenotypes['id'], phenotypes[cat]):
        if eid not in d_alpha:  #technically do not need cuz do this when make phenotype file
            continue

        if status>0:    #for diabetes: -1 = do not know, -3 = prefer not to answer  #no one in dataset has been diagnosed with diabetes at birth (status=0)
                        #for bipolar/depression: 0 = none
                        #for sex: 0=female, 1=male
            group1.append(d_alpha[eid])
        else:   #assume empty, as well as -1 and -3 correspond to no disorder
            group2.append(d_alpha[eid])

    return group1, group2


def plotout(group1, group2, cat, which):
    if cat=='bipolar/depression':
        plt.hist(group1, density = True, alpha = 0.3, label='bipolar or depression ($N$={0})'.format(len(group1)))
        plt.hist(group2, density = True, alpha = 0.3, label='healthy ($N$={0})'.format(len(group2)))
    elif cat=='sex':
        plt.hist(group1, density = True, alpha = 0.3, label='male ($N$={0})'.format(len(group1)))
        plt.hist(group2, density = True, alpha = 0.3, label='female ($N$={0})'.format(len(group2)))
    else:
        plt.hist(group1, density = True, alpha = 0.3, label='{0} ($N$={1})'.format(cat, len(group1)))
        plt.hist(group2, density = True, alpha = 0.3, label='healthy ($N$={0})'.format(len(group2)))

    plt.title('Increasing Tract ' + r"$\bf{" + which.capitalize() + "}$" + " Targeted Attack (UK Biobank)", fontsize=14)#.format(which))  #default size is 12
    plt.xlabel('growth parameter $\\alpha$', fontsize=14)
    plt.ylabel('frequency', fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    ks, pval = (ks_2samp(group1, group2))
    leg = plt.legend(title='K-S = {0:.2f} ({1:.2E})'.format(ks, pval), prop={'size':12})
    leg.set_title(title = 'K-S = {0:.2f} ({1:.2E})'.format(ks, pval), prop={'size':12})
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    which = 'length'
    category = 'sex' #bipolar/depression, diabetes or sex 

    alpha_file = 'database_output/ukb.csv' #this analysis is only implemented for ukb
    alphas = pd.read_csv(alpha_file)

    if which=='density':
        d_alpha = dict(zip(alphas.id, alphas.alpha_density))
    elif which=='length':
        d_alpha = dict(zip(alphas.id, alphas.alpha_length))


    phenotypes = pd.read_csv('./phenotypes/ukb/phenotypes.csv')
    group1, group2 = group_by_category(phenotypes, category)
    plotout(group1, group2, category, which)
