"""
Parse raw ABCD phenotype file and write out 
into a more easily readable file.

"""


import sys, os
import numpy as np
import pandas as pd


def parse_raw_phenotype_file():
    pheno_file = './fmriresults01.txt' 
    phenotypes = pd.read_csv(pheno_file,  sep="\t")
    source = phenotypes.file_source
    fnames = []
    for x in source:
        almost = x.split('/')[-1]
        firstTest = almost.find('_')
        fname = almost[:almost.find('_', firstTest+1)]
        fnames.append(fname)
    return list(phenotypes.interview_age), fnames 



if __name__ == '__main__':
    ages, eids = parse_raw_phenotype_file()

    df = pd.DataFrame(list(zip(*[eids, ages])))
    df.columns = ['id', 'age']
    df.to_csv('phenotypes.csv', index=False)
