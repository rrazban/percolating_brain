"""
Generate the percolation probability curve for
a random graph constrained in three dimensional
volume of the brain 

"""


import sys
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

from random_graph import plotout

sys.path.append('../../process/atlas')
import roi_coordinates

sys.path.append('../../analyze')
from Pcurve import break_apart


d_readin_coords = {'Harvard-Oxford': roi_coordinates.Harvard_Oxford, 'Desikan-Killiany': roi_coordinates.Desikan_Killiany, 'Talairach': roi_coordinates.Talairach}
#remove indices of background, white matter and orphan nodes
d_remove = {'Harvard-Oxford':[0, 1, 2, 8, 12, 13], 'Desikan-Killiany':[0, 1, 2, 40, 83], 'Talairach': [0, 1, 2, 3, 4, 5, 6, 14, 20, 24, 25, 26, 27, 28, 37, 39, 42, 46, 47, 49, 50, 51, 52, 53, 55, 56, 57, 58, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 80, 81, 82, 83, 84, 85, 86, 87, 91, 96, 98, 99, 101, 112, 113, 114, 115, 116, 117, 118, 119, 120, 128, 131, 133, 134, 137, 138, 148, 149, 150, 155, 156, 159, 160, 161, 162, 163, 164, 165, 167, 184, 189, 190, 193, 194, 195, 196, 200, 206, 207, 210, 211, 212, 220, 230, 231, 240, 243, 252, 253, 255, 256, 260, 261, 262, 267, 268, 272, 273, 289, 290, 291, 293, 304, 308, 309, 313, 317, 318, 319, 320, 321, 322, 326, 327, 331, 341, 344, 346, 347, 348, 353, 357, 359, 361, 370, 371, 374, 375, 387, 388, 389, 390, 391, 398, 400, 403, 404, 416, 424, 430, 434, 437, 438, 439, 441, 459, 461, 462, 467, 472, 485, 486, 489, 492, 493, 494, 496, 498, 502, 503, 504, 505, 506, 507, 508, 509, 510, 512, 513, 514, 515, 518, 519, 522, 523, 524, 526, 527, 530, 534, 536, 537, 538, 539, 540, 542, 543, 544, 545, 546, 547, 550, 551, 553, 554, 559, 562, 563, 564, 565, 569, 572, 573, 574, 576, 577, 580, 581, 582, 583, 588, 589, 590, 591, 593, 594, 595, 596, 597, 598, 601, 620, 621, 629, 630, 635, 642, 651, 652, 653, 657, 659, 662, 664, 666, 676, 677, 685, 687, 688, 694, 695, 696, 698, 701, 711, 712, 713, 714, 715, 716, 717, 718, 724, 731, 732, 736, 737, 743, 745, 746, 747, 753, 754, 765, 769, 775, 776, 777, 784, 786, 787, 788, 789, 790, 796, 798, 800, 801, 802, 811, 812, 817, 820, 821, 822, 831, 832, 833, 834, 842, 845, 850, 851, 854, 855, 860, 864, 865, 866, 867, 868, 870, 871, 874, 875, 876, 879, 880, 884, 885, 886, 888, 890, 905, 909, 917, 918, 919, 921, 922, 924, 925, 926, 929, 930, 931, 936, 937, 951, 954, 955, 956, 957, 960, 961, 962, 963, 964, 965, 966, 973, 976, 985, 988, 995, 996, 1001, 1002, 1020, 1041, 1042, 1047, 1051, 1052, 1054, 1055, 1056, 1057, 1058, 1059, 1060, 1063, 1064, 1073, 1074, 1084]}


def make_graph():
    k=14
    p=k/N
    G = nx.binomial_graph(N, p) #random graph
    edges = list(G.edges())

    distances = np.zeros((N, N))
    for i in range(N):
        for j in range(i+1, N):
            if (i, j) in edges:
                dist = np.linalg.norm(np.array(coords[i]) - np.array(coords[j]))	#probly faster to pick out terms rather than recalculate distance over and over	
                distances[i][j] = dist

#   plt.hist(distances.flatten())   #check out distribution
#   plt.show()

    thresholds = np.arange(0, np.max(distances), 0.02)
    avg_degrees, P_ones = break_apart(distances, thresholds)
    return avg_degrees, P_ones



if __name__ == '__main__':
    atlas = 'Harvard-Oxford'

    coords =  d_readin_coords[atlas]()

    remove = d_remove[atlas]
    #make sure goes from largest to smallest, or else indices change!
    for x in sorted(remove, reverse=True):
        del coords[x]
    N = len(coords)
	

    repeat = 10#0
    collection = np.arange(0, 14, 0.5)
    output = [[] for _ in collection]

    for r in range(repeat):
        print(r)
        ks, result = make_graph()
        for i, a0 in enumerate(collection):
            indi = (np.abs(ks-a0).argmin())
            output[i].append(result[indi])
        
    plt.title("{0} atlas".format(atlas), fontsize=16)
    plotout(collection, output, 'random ROI')
