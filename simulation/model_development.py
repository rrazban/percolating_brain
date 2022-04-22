"""
Generate the percolation probability curve based 
on targeted attack of model based on fasciculation
and scaling. Theory that assumes no secondary 
cluster formation is also included

"""


import sys
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import random
from math import pi, cos, sin
from scipy.stats import spearmanr

from theory import plotout
from theory_no_alpha import theory_no_alpha

sys.path.append('../analyze')
from Pcurve import get_k_and_P, preprocess, break_apart


n = 100	#number of nodes


def initialize_edge(M, i1, i2, amt):
	M[i1,i2]=amt
	M[i2,i1]=amt
	return M

def grow_edge(M, amt):
	M[np.nonzero(M)]+=amt
	return M


def num_per_orbit(dr):
    max_ = 15
    d = {}
    for r in range(max_):
        num_per_circle = int(2*np.pi*r/dr)	#based on circumference of unit circle
        d[r] = num_per_circle
    return d

def get_available_orbits(d_orbits, total_orbits, M_orbit):
    orbits = []
    d_num = {}
    for i in range(1, total_orbits+1):
        num = len(np.where(M_orbit==i)[0])

        for _ in range(d_orbits[i]-num):
            orbits.append(i)
    return orbits	

def circle_coords_equidistant(r,dr):
#adapted from https://stackoverflow.com/questions/33510979/generator-of-evenly-spaced-points-in-a-circle-in-python

    h=0    #coordinates of origin
    k=0

    num_per_circle = int(2*np.pi*r/dr)	#based on circumference of unit circle

    coords = []	
    for j in range(num_per_circle):    
        theta = 2 * pi * j/num_per_circle
        coord = [h + cos(theta) * r, k + sin(theta) * r]
        coords.append(coord)
    return coords



def get_distances(M_density, M_orbit, dr):
    d_coords_orbit = {}	#all coords for a given orbit

    for r in range(1, int(np.max(M_orbit))+1):
        coords = circle_coords_equidistant(r, dr)
        d_coords_orbit[r] = coords

    d_coords = {}
    for i in range(n):
        if M_orbit[i]!=0:	#if equal to 0 then i guess this node has no edges
            coord = d_coords_orbit[M_orbit[i]][0]	#just choose the first item of list randomly
            d_coords[i] = coord
            d_coords_orbit[M_orbit[i]].remove(coord)

    M_dist = np.zeros((n, n))
    pairs = np.nonzero(M_density)
    for i1, i2 in zip(pairs[0],pairs[1]):
        coord1 = d_coords[i1]
        coord2 = d_coords[i2]
        dist = np.sqrt( (coord2[0] - coord1[0])**2 + (coord2[1] - coord1[1])**2 )   #makes it a distance, not a length
        M_dist[i1,i2] = dist

    return M_dist

def distance_v_density(pre_M_density, pre_M_dist):
	M_density = preprocess(pre_M_density)   #get rid of diagonals
	M_dist = preprocess(pre_M_dist)	

	conns = M_density[np.nonzero(M_density)]
	dists = M_dist[np.nonzero(M_dist)]
	rho, pval = spearmanr(conns, dists)
	plt.title('Model Simulation', fontsize=16)
	plt.scatter(dists, conns, label = '$\\rho=$ {0:.2f} ({1:.2E})\n$N=$ {2}'.format(rho, pval, len(conns)))
	plt.ylabel('tract density', fontsize=14)
	plt.xlabel('tract distance', fontsize=14)
	plt.xticks(fontsize=14)
	plt.yticks(fontsize=14)
	plt.legend(prop={'size':12})
	plt.yscale('log')
	plt.xscale('log')
	plt.tight_layout()
	plt.show()




def make_graph():
    M_density = np.zeros((n, n))
    M_orbit = np.zeros(n)	#have regions be randomly placed along circle	#randomly place at the end
                                #this should be called orbits!!!
    G=nx.Graph()
    G.add_nodes_from(range(n))

    max_k = 14
    p = max_k/n
    avgdegree = 0
    dr = 3
    d_orbits = num_per_orbit(dr)


    #make first edge
    ind1,ind2 = np.random.choice(range(n), size=2, replace=False)
    potential_edge = ((ind1, ind2))
    G.add_edge(*potential_edge)

    M_density = initialize_edge(M_density, ind1, ind2, 1)
    M_orbit[ind1] = 1
    M_orbit[ind2] = 1
    M_orbit = grow_edge(M_orbit, 1)  
    total_orbits = 2  #assumes that first orbit has radius of 1 (can only fit two regions)


    while avgdegree < max_k:
        Gcc = list(max(nx.connected_components(G), key=len))

        ind1 = random.choice(Gcc)	#one node must always be a part of the giant cluster
        ind2 = random.choice(list(range(n)))	#choose from any node
        potential_edge = ((ind1, ind2))

        if random.random() < p:

            G.add_edge(*potential_edge)
            avgdegree, P_one = get_k_and_P(G)

            if ind2 not in Gcc:

                available_orbits = get_available_orbits(d_orbits, total_orbits, M_orbit)

                if len(available_orbits)==0:  #current available orbits are full, shift +1
                    M_orbit = grow_edge(M_orbit, 1)
                    total_orbits += 1	
                    available_orbits = get_available_orbits(d_orbits, total_orbits, M_orbit)

#                   M_density = grow_edge(M_density, 1)
                    M_density = grow_edge(M_density, random.uniform(1,10)) #introduce some randomness so that do not have a lot of edges with the exact same density #becomes explosive percolation if no randomness

                M_orbit[ind2] = random.choice(available_orbits)
				
#            M_density = initialize_edge(M_density, ind1, ind2, 1)
            M_density = initialize_edge(M_density, ind1, ind2, random.uniform(1,10))


    M_distance = get_distances(M_density, M_orbit,dr) #implicitly assume length and distance are interchangeable

    distance_v_density(M_density, M_distance) #distance vs density correlation
    sys.exit()

    M_adjacency_dist = preprocess(M_distance)
    thresholds = np.arange(0, np.max(M_adjacency_dist), 1)
    k_distance, P_distance = break_apart(M_adjacency_dist, thresholds)

    M_adjacency_density = preprocess(M_density)
    thresholds = np.arange(0, np.max(M_adjacency_density), 1)
    k_density, P_density = break_apart(M_adjacency_density, thresholds)

    return k_distance, P_distance, k_density, P_density


if __name__ == '__main__':
    repeat = 10#00

    collection = np.arange(0, 14, 0.5)
    output = [[] for _ in collection]
    output1 = [[] for _ in collection]

    for r in range(repeat):
        print(r)

        ks, result, ks1, result1 = make_graph()
        for i, a0 in enumerate(collection):

            indi = (np.abs(ks-a0).argmin())
            output[i].append(result[indi])
            indi1 = (np.abs(ks1-a0).argmin())
            output1[i].append(result1[indi1])

    pred_output = theory_no_alpha(collection)
    plotout(collection, output1, pred_output, '', 'density')
    plotout(collection, output, pred_output, '', 'distance')

