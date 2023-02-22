"""
Pcurve_Fig2.py across ages from neonatal (dHCP), 
adolescent (ABCD Study) to adult (UK Biobank). 

Creates the top portion of Figure 3 in the main 
text.

Note that this Figure cannot be reproduced from 
output files in this GitHub directory. dHCP and 
ABCD output files are not provided and must be 
generated from process/dmri2adjacency_matrix.py 
after acquiring dMRI scans from respective dataset.

"""


import numpy as np
import matplotlib.pyplot as plt

from Pcurve_Fig2 import get_experimental_data


def beautify_figure(ax, avg_degrees, P_ones, alpha, exp_label):
    plt.title('Increasing Tract ' + r"$\bf{" + exp_label.capitalize() + "}$" + " Targeted Attack", fontsize=16)
    plt.xlabel('average degree $\langle k \\rangle$', fontsize = 14)
    plt.ylabel('probability in the giant cluster $P$', fontsize=14)
    legend=plt.legend(title='subject age', prop={'size': 14})

    legend.get_title().set_fontsize('14')

    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlim([-1, max(avg_degrees)+1])
    plt.ylim([-0.05, 1.05])

    plt.tight_layout()
    plt.show()



def sample(ks, y):  #same as Pcurve_simulations.py
    collection = np.arange(0, 30, 0.5) 

    new_y = []
    for i, a0 in enumerate(collection):
        indi = (np.abs(ks-a0).argmin())
        new_y.append(y[indi])

    return collection, np.array(new_y)





if __name__ == '__main__':
    which = 'length'    #tract length or tract density

    ages = ['36 weeks (dHCP)', '10 years (ABCD)', '51 years (UKB)', '80 years (UKB)']
    filenames = ['/shared/datasets/public/dhcp/derivatives/sub-CC00063AN06_ses-15102_desc-preproc_dwi_{0}.txt'.format(which), '/shared/datasets/public/abcd/derivatives/sub-NDARINVNVF8N71U_ses-baselineYear1Arm1_run-01_dwi_{0}.txt'.format(which),  '/shared/datasets/public/ukb/dMRI_derivatives/qball/6025360_20250_2_0_{0}.txt'.format(which), '/shared/datasets/public/ukb/dMRI_derivatives/qball/4482035_20250_2_0_{0}.txt'.format(which)]   #must generate these output files for yourself, see comment at the top of this script 

    fig, ax = plt.subplots()
    for age,filename in zip(ages, filenames):
        print(filename)
        ori_avg_degrees, ori_P_ones = get_experimental_data(filename)
        avg_degrees, P_ones = sample(ori_avg_degrees, ori_P_ones)

        plt.plot(avg_degrees, P_ones, label=age)
    beautify_figure(ax, avg_degrees, P_ones, 1, which)
