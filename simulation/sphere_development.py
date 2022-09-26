"""
Generate the percolation probability curve based 
on targeted attack of model based on fasciculation
and scaling. Theory that assumes no secondary 
cluster formation is also included

"""


import sys
import matplotlib.pyplot as plt
import numpy as np
import random


from theory import plotout, theory

sys.path.append('../analyze')
from Pcurve import preprocess, break_apart, sample_equidistant, curve_fit



d_k_alpha11 = {100: 50, 727:100}
d_radius_gr = {100: 1.001, 727: 1.0001}


def get_distances(M_density, coords):
    M_dist = np.zeros((n, n))
    pairs = np.nonzero(M_density)
    for i1, i2 in zip(pairs[0],pairs[1]):
        coord1 = coords[i1]
        coord2 = coords[i2]

        dist = np.linalg.norm(coord1-coord2)
        M_dist[i1,i2] = dist
        M_dist[i2,i1] = dist

    return M_dist


def random_coord(amt):
    dx = random.uniform(0,amt) * random.choice([-1,1])

    dy_lim = (amt**2 - dx**2)**(0.5)
    dy = random.uniform(0, dy_lim) * random.choice([-1,1])

    dz = (amt**2 - dx**2 - dy**2)**(0.5) * random.choice([-1,1])
    return np.array([dx, dy, dz]) 


def make_graph(n, rescale_p_ngc):

    ks = []
    Ps=[]
    M_density = np.zeros((n, n))
    coords = np.zeros((n, 3))

    M_adjacency = np.zeros((n,n))
    M_power = np.zeros((n,n))

    M_pg = np.zeros(n)
    M_radius_power = np.zeros(n)

    all_nodes = list(range(n)) 

    #parameters to grow graph
    max_k = d_k_alpha11[n] 
    p = max_k/n 

    initial_density = 1
    density_growth_rate = 1.001 #works for n=727(largest order of magnitude)    #also works for n=100 although can go larger including 1.1
    initial_radius = 1
    radius_growth_rate = d_radius_gr[n] 

    #make first edge
    ind1,ind2 = np.random.choice(all_nodes, size=2, replace=False)
    first_edge = ((ind1, ind2))

    avgdegree= 2/n
    P_one = 2* 1/n
    Gcc = [ind1, ind2]

    M_density[ind1][ind2] = initial_density 
    M_density[ind2][ind1] = M_density[ind1][ind2]

    coords[ind1] = random_coord(initial_radius) 
    coords[ind2] = random_coord(initial_radius) 

    M_adjacency[ind1][ind2] = 1
    M_adjacency[ind2][ind1] = 1
    M_pg[ind1]=1
    M_pg[ind2] = 1
    while avgdegree<max_k: #surprising that this looks worse for n=727, P curve is much more depressed
#    while P_one < 1:#(n-1)/n: #interesting for large n, it never can quite get to 1 proly cuz of rounding issues

        ind1 = random.choice(Gcc)	#one node must always be a part of the giant cluster
        ind2 = random.choice(all_nodes)	

        if ind2 in Gcc:
            if M_density[ind1][ind2]==0 and ind1!=ind2:    #avoid replicate edges reseting density back to zero
                if random.random() <  p/2:

                    M_density[ind1][ind2] = initial_density 
                    M_density[ind2][ind1] = M_density[ind1][ind2]#must have transpose so that dont overwrite

                    M_power += M_adjacency 
                    M_adjacency[ind1][ind2] = 1
 #                   M_adjacency[ind2][ind1] = 1    #no need, take the transpose later
                    M_radius_power += M_pg

                    avgdegree+=2/n
        else:
            if random.random() < p/rescale_p_ngc:

                    avgdegree+=2/n
                    P_one+=1/n
                    Gcc.append(ind2)

                    M_density[ind1][ind2] = initial_density
                    M_density[ind2][ind1] = M_density[ind1][ind2]

                    coords[ind2] =  random_coord(initial_radius)

                    M_power += M_adjacency 
                    M_adjacency[ind1][ind2] = 1

                    M_radius_power += M_pg
                    M_pg[ind2] = 1

                    ks.append(avgdegree)
                    Ps.append(P_one)

    M_density *= density_growth_rate**(M_power+M_power.T)

    coord_factor = radius_growth_rate ** M_radius_power
    coord_factor = np.column_stack((coord_factor, coord_factor, coord_factor))
    coords *= coord_factor #need to apply across row


    M_distance = get_distances(M_density, coords) #implicitly assume length and distance are interchangeable
    M_adjacency_dist = preprocess(M_distance)
#    plt.hist(np.log10(M_adjacency_dist[M_adjacency_dist!=0].flatten()))
 #   plt.show()
    min_ = np.min(np.log10(M_adjacency_dist[M_adjacency_dist!=0].flatten()))
    thresholds_dist = np.logspace(min_, np.log10(np.max(M_adjacency_dist)), 500) 
    k_distance, P_distance = break_apart(M_adjacency_dist, thresholds_dist)

    M_adjacency_density = preprocess(M_density)

    dmin_ = np.min(np.log10(M_adjacency_density[M_adjacency_density!=0].flatten()))
    thresholds_density = np.logspace(dmin_, np.log10(np.max(M_adjacency_density)), 500) 

    k_density, P_density = break_apart(M_adjacency_density, thresholds_density)

    print(avgdegree, P_one)  #make sure len(Gcc) is equal to n!!!!
    return ks, Ps, k_density, P_density, k_distance, P_distance


def get_output(collection, ks, Ps, output):
    for i, a0 in enumerate(collection):
        indi = (np.abs(ks-a0).argmin())
        output[i].append(Ps[indi])
    return output

def fit_alpha(ks, Ps):
    avgdegree, Pones = sample_equidistant(ks, Ps)
    popt, pcov = curve_fit(theory, avgdegree, Pones)

    return popt[0]

if __name__ == '__main__':
    repeat = 10#00

    n = 100
    alpha = 11

    if n==727:
        collection = np.arange(0, 30, 0.5)
    else:
        collection = np.arange(0, 20, 0.5)


    output_standard = [[] for _ in collection]
    output_density = [[] for _ in collection]
    output_distance = [[] for _ in collection]

    check_alpha_dist = False 
    alphas_distance = [] 
    alphas_density = [] 
    for r in range(repeat):
        print(r)

        k_standard, P_standard, k_density, P_density, k_distance, P_distance = make_graph(n, alpha)

        output_standard = get_output(collection, k_standard, P_standard, output_standard)
        output_density = get_output(collection, k_density, P_density, output_density)
        output_distance = get_output(collection, k_distance, P_distance, output_distance)
        if check_alpha_dist:
            alphas_distance.append(fit_alpha(k_distance, P_distance))
            alphas_density.append(fit_alpha(k_density, P_density))

    if check_alpha_dist:
        plt.hist(alphas_distance, label='distance')
        plt.hist(alphas_density, label='density')
        plt.legend()
        plt.show()

    pred_output = theory(collection, alpha)
    plotout(collection, output_standard, pred_output, alpha, 'standard')
    plotout(collection, output_density, pred_output,alpha , 'density')
    plotout(collection, output_distance, pred_output, alpha, 'distance')
