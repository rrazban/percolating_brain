"""
Generate the percolation probability curve based 
on theory that assumes no secondary cluster formation

"""


import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import lambertw

from theory import make_graph, plotout


n = 100	#number of nodes

def theory_no_alpha(ks):
    p_soln = 1 + lambertw(-(ks**2/2 + ks + 1) * np.exp(-(ks + 1)))
    p_soln[0] = 0
    return p_soln


if __name__ == '__main__':
    repeat = 10#00

    collection = np.arange(0, 14, 0.5)
    output = [[] for _ in collection]

    for r in range(repeat):
        ks, result = make_graph(rescale_p_ngc=1, max_k = 14)
        for i, a0 in enumerate(collection):

            indi = (np.abs(ks-a0).argmin())
            output[i].append(result[indi])

    pred_output = theory_no_alpha(collection)
    plotout(collection, output, pred_output, '', '')
