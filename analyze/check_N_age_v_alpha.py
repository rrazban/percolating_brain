"""
Check that the alpha vs age relationship is not
biased by the number of the nodes (N) of the 
network.


"""


import pandas as pd
import matplotlib.pyplot as plt

from check_kf_age_v_alpha import plotout


if __name__ == '__main__':
    which = 'density'    #length or density

    for thresh in [724, 727, 730]:
        dataset='ukb'

        alpha_file = 'database_output/{0}.csv'.format(dataset)
        alpha_output = pd.read_csv(alpha_file)
        if which=='density':
            d_alpha = dict(zip(alpha_output.id, alpha_output.alpha_density))
        elif which=='length':
            d_alpha = dict(zip(alpha_output.id, alpha_output.alpha_length)) 


        basic_file = 'database_output/{0}_basic.csv'.format(dataset)
        basic_output = pd.read_csv(basic_file)

        eids = list(basic_output[basic_output.N==thresh].id)

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

    legend = plt.legend(prop={'size':12}, title='number of nodes')
    legend.get_title().set_fontsize('14')

    plt.tight_layout()
    plt.show()
