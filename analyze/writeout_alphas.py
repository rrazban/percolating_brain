"""
Fit alpha values to experimental P curves for
each individual. Write out alpha vales to a 
text file named by the respective database,
UK Biobank, ABCD Study or dHCP.

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


def get_ids_abcd(filenames):
    eids = []

    for filename in filenames:
        pre_eid = filename.split('/')[-1]
        first = pre_eid[pre_eid.index('-')+1:pre_eid.index('_')]
        new_eid = pre_eid[pre_eid.index('_')+1:]
        second = new_eid[new_eid.index('-')+1:new_eid.index('_')]
        fname='{0}_{1}'.format(first, second)
        eids.append(fname)
    return eids


def get_ids_ukb(filenames):
    eids = []

    for filename in filenames:
        pre_eid = filename.split('/')[-1]
        eid = int(pre_eid[:pre_eid.index('_')])
        eids.append(eid)
    return eids

def get_ids_dhcp(filenames):
    eids = []

    for filename in filenames:
        pre_eid = filename.split('/')[-1]
        eid = pre_eid[pre_eid.index('-')+1:pre_eid.index('_')]

        full_id='{0}_{1}'.format(first, second)
        eids.append(eid)
    return eids


def run(filename):
   
    avg_degrees, P_ones = get_experimental_data(filename)
    avg_degrees, P_ones = sample_equidistant(avg_degrees, P_ones) 

    zero_i = find_zero_index(list(avg_degrees))

    popt, pcov = curve_fit(theory, avg_degrees[:zero_i], P_ones[:zero_i])
    return popt[0]

    #check goodness of fit
#    pred_P_ones = theory(avg_degrees, popt[0])
 #   r, pval = (spearmanr(pred_P_ones[:zero_i], P_ones[:zero_i]))
#    return -np.log10(pval) 
#    return pcov[0][0]**0.5  #transform variance to std #better than r and pval cuz essentially 1 and 0



if __name__ == '__main__':
    start_time = datetime.now()    
    print("Start time: {0}".format(start_time))

    dataset = 'dhcp'    #abcd, ukb, dhcp

    collect = []
    for which in ['density', 'length']: 
        filenames = sorted(glob.glob('/shared/datasets/public/{0}/derivatives/*_{1}.txt'.format(dataset, which)))#[:2]

        print("working on {0} dataset's {1}.txt files, N={2}".format(dataset, which, len(filenames)))

        if dataset == 'ukb':
            eids = get_ids_ukb(filenames)
        elif dataset == 'abcd':
            eids = get_ids_abcd(filenames)
        elif dataset == 'dhcp':
            eids = get_ids_abcd(filenames)

        with Pool(28) as p: #28 is the number of processors 
            alphas = p.map(run, filenames)
        
        collect.append(alphas)

    df = pd.DataFrame(list(zip(*[eids, collect[0], collect[1]])))
    df.columns = ['id', 'alpha_density', 'alpha_length']
    df.to_csv('{0}.csv'.format(dataset), index=False)

    end_time = datetime.now()
    print("End time: {0}".format(end_time))
    print("Duration of analyses: {0}".format(end_time-start_time))

