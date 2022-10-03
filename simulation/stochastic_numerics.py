"""
Calculate numerical solution for probability of 
a given graph with n nodes and E edges according
to our percolation theory.

"""


import sys
import matplotlib.pyplot as plt

import numpy as np



def transition(n, E, alpha, n_tot):
    p_sg = n*(n_tot - n)/alpha 
    p_gg = n*(n-1)/2 - E
    return p_sg/(p_sg + p_gg)

def set_probability(n, e, d_p):
    E_tot = int(n*(n-1)/2)

    if e>E_tot or n-e>=2:
        p = 0
    elif e in d_p[n]:
        p = d_p[n][e]
    else:
        print('huh?',n,e)
        sys.exit()
    return p


def numerics_theory(max_edges, alpha, n_tot):
    d_p = {}    #dictionary of p(n|E)
    for n in range(1, n_tot+1):
        d_p[n] = {}
    
    d_p[1][0]=1 #initial condition

    for n in range(2, n_tot+1):
        E_tot = int(n*(n-1)/2)
        E_tot = min(E_tot, max_edges)
        for E in range(n-2, E_tot):
 #           print('number of edges: {0}'.format(E))
            p_ne = set_probability(n, E, d_p)
            p_n1e = set_probability(n-1, E, d_p)

            d_p[n][E+1] = p_ne * (1 - transition(n, E, alpha, n_tot)) + p_n1e * transition(n-1,E, alpha, n_tot)
    return d_p
    
def transform_p_to_matrix(d_p, n_tot, E_tot):
    M_mean = np.zeros((n_tot+1, E_tot+1)) 
    M_var = np.zeros((n_tot+1, E_tot+1)) 
    for n in range(1, n_tot+1):
        for E,p in d_p[n].items():
            M_mean[n][E] += p * n/n_tot  
            M_var[n][E] += p * (n/n_tot)**2
    return M_mean, M_var

def get_mean_std(M_mean, M_var):
    mean = np.sum(M_mean, axis=0)
    std = (np.sum(M_var, axis=0) - (mean)**2)**(0.5)
    return mean, std

if __name__ == '__main__':
    repeat = 10#00

    n_tot = 100 
    max_k = 50  #maximum average degree of graph
    alpha = 11

    collection = np.arange(0, 20, 0.5)



    from theory import make_graph
    output = [[] for _ in collection]
    for r in range(repeat):
        print(r)	
        _, ks, result = make_graph(alpha, max_k, n_tot)
        for i, a0 in enumerate(collection): #save Ps closeset to target <k> for accurate error bars
            indi = (np.abs(ks-a0).argmin())
            output[i].append(result[indi])

    plt.violinplot(output, positions = collection, showmeans=True, showmedians=True)#, showextrema=False)


    E_tot = int(n_tot*max(collection)/2)#divide by 2 to avoid overcounting
    d_p = numerics_theory(E_tot, alpha, n_tot)
    M_mean, M_var = transform_p_to_matrix(d_p, n_tot, E_tot)
    mean, std = get_mean_std(M_mean, M_var)

    freq = int(E_tot/(len(collection)-1))
    plt.errorbar(collection, mean[::freq], yerr=std[::freq]) 

    plt.title('$\\alpha = ${0}'.format(alpha), fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.ylabel('probability in the giant cluster $P$', fontsize = 14 )
    plt.xlabel('average degree $\langle k \\rangle$', fontsize = 14 )
    plt.show()

