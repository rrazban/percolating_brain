"""
Fit alpha values to experimental P curves for
each individual. Write out alpha vales to a 
text file named by the respective database,
UK Biobank or ABCD Study.

"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
from multiprocessing import Pool
from scipy.optimize import curve_fit

from Pcurve import get_experimental_data, find_zero_index, theory, sample_equidistant
from scipy.stats import spearmanr
from datetime import datetime 


def get_ages_abcd(filenames):
    d_age = parse_phenotype('abcd')

    ages = []
    eids = []
    not_present = []

    for filename in filenames:
        pre_eid = filename.split('/')[-1]
        first = pre_eid[pre_eid.index('-')+1:pre_eid.index('_')]
        new_eid = pre_eid[pre_eid.index('_')+1:]
        second = new_eid[new_eid.index('-')+1:new_eid.index('_')]
        fname='{0}_{1}'.format(first, second)
        if fname in d_age:
            age = d_age[fname]
            ages.append(age)
            eids.append(fname)
        else:
            print("{0} not present in fmriresults01.txt".format(fname))
            not_present.append(filename)
    return eids, ages, not_present


def parse_phenotype(dataset):
    pheno_file = './phenotypes/{0}/phenotypes.csv'.format(dataset) 
    phenotypes = pd.read_csv(pheno_file)
	
    d_age = dict(zip(phenotypes.id, phenotypes.age))
    return d_age

def get_ages_ukb(filenames):
    d_age = parse_phenotype('ukb')

    ages = []
    not_present = []
    eids = []

    for filename in filenames:
        pre_eid = filename.split('/')[-1]
        eid = int(pre_eid[:pre_eid.index('_')])
        if eid in d_age:
            age = d_age[eid]
            ages.append(age)
            eids.append(eid)
        else:
            print("{0} not present in phenotypes.csv".format(eid))
            not_present.append(filename)
    return eids, ages, not_present


def run(filename):
   
    avg_degrees, P_ones = get_experimental_data(filename)
    avg_degrees, P_ones = sample_equidistant(avg_degrees, P_ones) 

    zero_i = find_zero_index(list(avg_degrees))

    popt, pcov = curve_fit(theory, avg_degrees[:zero_i], P_ones[:zero_i])
    return popt[0]

    #check goodness of fit
 #   pred_P_ones = theory(avg_degrees, popt[0])
#    r, pval = (spearmanr(pred_P_ones[:zero_i], P_ones[:zero_i]))
#    return -np.log10(pval) 




if __name__ == '__main__':
    start_time = datetime.now()    
    print("Start time: {0}".format(start_time))

    dataset = 'abcd'    #abcd or ukb

    collect = []
    for which in ['density', 'length']: 
        filenames = sorted(glob.glob('/shared/datasets/public/{0}/derivatives/*_{1}.txt'.format(dataset, which)))

        print("working on {0} dataset's {1}.txt files, N={2}".format(dataset, which, len(filenames)))

        if dataset == 'ukb':
            eids, ages, not_present = get_ages_ukb(filenames)
        elif dataset == 'abcd':
            eids, ages, not_present = get_ages_abcd(filenames)

        for fname in not_present:
            filenames.remove(fname)

        with Pool(28) as p: #28 is the number of processors 
            alphas = p.map(run, filenames)
        
        collect.append(alphas)

    df = pd.DataFrame(list(zip(*[eids, collect[0], collect[1]])))
    df.columns = ['id', 'alpha_density', 'alpha_length']
    df.to_csv('{0}.csv'.format(dataset), index=False)

    end_time = datetime.now()
    print("End time: {0}".format(end_time))
    print("Duration of analyses: {0}".format(end_time-start_time))

