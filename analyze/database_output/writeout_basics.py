"""
Calculate basic network properties, number 
of nodes, average degree and length-density
correlation. Write out vales to a text file 
named by the respective database, 
UK Biobank, ABCD Study or dHCP.

"""


import numpy as np
import pandas as pd
import glob
from multiprocessing import Pool, set_start_method
import networkx as nx

from writeout_alphas import get_ids_abcd, get_ids_ukb, get_ids_dhcp
from Pcurve import get_experimental_data, preprocess
from scipy.stats import spearmanr
from datetime import datetime 


set_start_method("spawn", force=True)   #not sure if this helps
#https://pythonspeed.com/articles/python-multiprocessing/


def run(filename):
   
    M_density = pd.read_csv(filename, delimiter=" ", header=None).values
    M_density = preprocess(M_density)

    len_fname = "{0}_length.txt".format(filename[:-12])
    M_length = pd.read_csv(len_fname, delimiter=" ", header=None).values
    M_length = preprocess(M_length)

    densities = M_density[M_density!=0]
    lens = M_length[M_length!=0]
    rho, pval = spearmanr(densities, lens)

    N = len(M_density)

    G = nx.from_numpy_array(M_density)
    degree = [val for (node, val) in G.degree()]
    avg_degree = np.mean(degree)

    return N, avg_degree, rho 


if __name__ == '__main__':
    start_time = datetime.now()    
    print("Start time: {0}".format(start_time))

    dataset = 'abcd'    #abcd, ukb, dhcp

    filenames = sorted(glob.glob('/shared/datasets/public/{0}/derivatives/*_{1}.txt'.format(dataset, 'density')))#[:500]#[:2]

    print("working on {0} dataset, N={1}".format(dataset,len(filenames)))

    with Pool(30) as p:
        alphas = p.map(run, filenames)
        p.close()   #not sure if helps, but standard practice to have close and join
        p.join()

    if dataset == 'ukb':
        eids = get_ids_ukb(filenames)
    elif dataset == 'abcd':
        eids = get_ids_abcd(filenames)
    elif dataset == 'dhcp':
        eids = get_ids_abcd(filenames)


    collect = []
    for eid,alpha in zip(eids, alphas):
        new_tuple=(eid,)+ alpha
        collect.append(new_tuple)

    df = pd.DataFrame(collect)
    df.columns = ['id', 'N','kf', 'rho_length_density']
    df.to_csv('{0}_basic.csv'.format(dataset), index=False)

    end_time = datetime.now()
    print("End time: {0}".format(end_time))
    print("Duration of analyses: {0}".format(end_time-start_time))

