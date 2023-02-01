"""
Plot histograms of basic newtork properties, 
number of nodes, average degree and length-
density correlation. All three databases are
shown in one figure.

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



def remove_replicates(d_alpha, dataset):    #only ABCD needs cuz has ~1000 people who have 2 scans
    phenotypes = pd.read_csv('./phenotypes/{0}/phenotypes.csv'.format(dataset))

    outs = []
    for eid, status in zip(phenotypes['id'], phenotypes['age']):
        if eid not in d_alpha:
            continue
        outs.append(d_alpha[eid])
        d_alpha.pop(eid, None)  #seems like we got reps in phenotypes.csv for ABCD
    return outs

if __name__ == '__main__':
    which = 'length'    #length or density

    fig, axes = plt.subplots(3, 3, figsize=(10,2*4.4))

    for i,dataset in enumerate(['ukb', 'abcd','dhcp']):
        alpha_file = 'database_output/{0}_basic.csv'.format(dataset)
        alpha_output = pd.read_csv(alpha_file)

        d_N = dict(zip(alpha_output.id, alpha_output.N))
        d_kf = dict(zip(alpha_output.id, alpha_output.kf))
        d_rho = dict(zip(alpha_output.id, alpha_output.rho_length_density))
        for j, d_alpha in enumerate([d_N, d_kf, d_rho]):
            outs = remove_replicates(d_alpha, dataset)
#            print(np.mean(outs), np.std(outs), min(outs))
            axes[i][j].hist(outs)#, density=True) 

            if i==2 and j==0:
                axes[i][j].set_ylabel('dHCP', fontsize=14)
            if i==1 and j==0:
                axes[i][j].set_ylabel('ABCD Study', fontsize=14)
            if i==0 and j==0:
                axes[i][j].set_ylabel('UK Biobank', fontsize=14)

            if i==2:
                if j==0:
                    axes[2][j].set_xlabel('number of nodes', fontsize=14)
                elif j==1:
                    axes[2][j].set_xlabel('average degree', fontsize=14)
                elif j==2:
                    axes[2][j].set_xlabel('$\\rho$(length, density)', fontsize=14)
    plt.show()
