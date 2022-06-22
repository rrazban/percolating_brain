"""
Parse raw UK Biobank phenotype file and write out 
into a more easily readable file.

"""


import sys, os
import numpy as np
import pandas as pd


def parse_raw_phenotype_file():
    pheno_file = 'too_big_ukb25909.csv'

    eids = []
    ages = []
    diabetes = []
    bipolars = []
    genders = []

    age_code = '21003-2.0'  #age at dMRI scan
    diabetes_code = '2976-0.0'  #age at diabetes diagnosis
    bipolar_dep_code = '20126-0.0'  #categorical, indicates bipolar or depression severity
#   alzheimers_code = '42020-0.0'   #do not currently have in my raw ukb phenotype file
    gender_code = '31-0.0'


    alpha_file = '../../ukb.csv'
    alphas = pd.read_csv(alpha_file)
    dmris = list(alphas.id)

    with open(pheno_file, 'r') as rfile:
        columns = next(rfile).split('","')

        for line in rfile:
            words = line.split('","')
            eid = words[0][1:]  #remove first "

            age = words[columns.index(age_code)]
            diabete = words[columns.index(diabetes_code)]
            bipolar = words[columns.index(bipolar_dep_code)]
            gender = words[columns.index(gender_code)]
        
            if int(eid) in dmris:
                eids.append(eid)
                ages.append(age)
                diabetes.append(diabete)
                bipolars.append(bipolar)
                genders.append(gender)
    return eids, ages, diabetes, bipolars, genders



if __name__ == '__main__':
    eids, ages, diabetes, bipolars, genders = parse_raw_phenotype_file()

    print('Number of individuals parsed: {0}'.format( len(eids)))

    df = pd.DataFrame(list(zip(*[eids, ages, diabetes, bipolars, genders])))
    df.columns = ['id', 'age', 'diabetes', 'bipolar/depression', 'gender']
    df.to_csv('phenotypes.csv', index=False)
