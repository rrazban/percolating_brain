"""
Parse raw dHCP phenotype file and write out 
into a more easily readable file for scripts.

"""


import sys, os
import numpy as np
import pandas as pd


def parse_raw_phenotype_file():
    pheno_file = 'combined.tsv'

    eids = []
    ages = []
    birth_ages = []

    session_code = 'session_id'

    age_code = 'scan_age'  #age at dMRI scan    #in weeks
    birth_code = 'birth_age' 

    alpha_file = '../../dhcp.csv'
    alphas = pd.read_csv(alpha_file)
    dmris = list(alphas.id)

    with open(pheno_file, 'r') as rfile:
        columns = next(rfile).split('\t')

        for line in rfile:
            words = line.split('\t')
            eid = words[0]
            session = words[columns.index(session_code)]
            full_id = '{0}_{1}'.format(eid, session)

            age = words[columns.index(age_code)]
            birth_age = words[columns.index(birth_code)]
            if full_id in dmris:
                eids.append(full_id)
                ages.append(age)
                birth_ages.append(birth_age)
    return eids, ages, birth_ages 



if __name__ == '__main__':
    eids, ages, birth_ages = parse_raw_phenotype_file()

    print('Number of individuals parsed: {0}'.format( len(eids)))

    df = pd.DataFrame(list(zip(*[eids, ages, birth_ages])))
    df.columns = ['id', 'age', 'birth_age']
    df.to_csv('phenotypes.csv', index=False)
