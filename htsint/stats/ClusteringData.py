#!/usr/bin/env python


import numpy as np
import networkx as nx

###########################
### create scatter data ###
###########################   
n = 6                                                                 # specify the number of observations for each cluster  
np.random.seed(42)                                                    # provide a seed for the random number generator  
feature1 = np.array([np.random.normal(2,1,n)]).transpose()            # creates a colums vector from a gaussian for cluster 1
feature2 = np.array([np.random.normal(2,1,n)]).transpose()            # creates a colums vector from a gaussian for cluster 1 
cluster1 = np.hstack((feature1,feature2))                             # combines the col vectors into one matrix 
feature1 = np.array([np.random.normal(7,1,n)]).transpose()            # creates a colums vector from a gaussian for cluster 1 
feature2 = np.array([np.random.normal(7,1,n)]).transpose()            # creates a colums vector from a gaussian for cluster 1  
cluster2 = np.hstack((feature1,feature2))                             # combines the col vectors into one matrix  
dataScatter = np.vstack((cluster1,cluster2))                          # combines the clusters
dataScatterLabels = np.zeros(n * 2,dtype=int)                         # specifies the actual labels
dataScatterLabels[-n:] = 1

#############################  
### create circle data    ###    
#############################   
n = 12
feature1 = np.array([np.random.normal(4.8,0.2,n)]).transpose()
feature2 = np.array([np.random.normal(5.1,0.2,n)]).transpose()
cluster1 = np.hstack((feature1,feature2))
cluster2 = np.array([[3.5,4.0],[3.2,4.5],[3.2, 5.0], [3.5,5.5], [3.8,6.0],[4.1,6.5],[4.5,7.0],[4.9,7.0],[5.3,6.8],[5.7,6.4],
                     [6.1,6.0],[6.3,5.6],[6.3,5.2],[6.2,4.8],[6.1,4.2],[5.8,3.8],[5.4,3.4],[5.0,3.1],[4.6,3.1],[4.2,3.3],[3.8,3.6]])

n1,d1 = np.shape(cluster1)
n2,d2 = np.shape(cluster2)
dataCircle = np.vstack([cluster1,cluster2])
dataCircleLabels = np.hstack([np.zeros(n1,dtype=int),np.zeros(n2,dtype=int)+1])

##############################  
### create letters data    ### [3.35,3.15]
############################## 
cluster1 = np.array([[1.10,3.3],[1.18,3.3],[1.26,3.3],[1.26,3.2],[1.22,3.1],[1.14,3.1],[1.06,3.1],[0.98,3.1],
                     [0.90,3.1],[0.82,3.1],[0.75,3.15],[0.71,3.2],[0.7,3.3],[0.7,3.4],[0.7,3.5],[0.70,3.6],
                     [0.70,3.7],[0.73,3.8],[0.8,3.88],[0.9,3.9],[1.0,3.9],[1.07,3.88],[1.15,3.83],[1.21,3.75],[1.26,3.65]])

cluster2 = np.array([[2.0,3.9],[2.0,3.8],[2.0,3.7],[2.0,3.6],[2.0,3.5],[2.0,3.4],[2.0,3.3],[2.0,3.2],[2.0,3.1],
                     [2.06,3.8],[2.12,3.7],[2.18,3.6],[2.3,3.5],[2.3,3.4],[2.36,3.3],[2.42,3.2],[2.48,3.1],
                     [2.55,3.9],[2.55,3.8],[2.55,3.7],[2.55,3.6],[2.55,3.5],[2.55,3.4],[2.55,3.3],[2.55,3.2],[2.55,3.1]])

cluster3 = np.array([[3.3,3.9],[3.3,3.8],[3.3,3.7],[3.3,3.6],[3.3,3.5],[3.3,3.4],[3.3,3.3],[3.3,3.22],
                     [3.35,3.17],[3.45,3.1],[3.66,3.1],[3.75,3.17],
                     [3.8,3.9],[3.8,3.8],[3.8,3.7],[3.8,3.6],[3.8,3.5],[3.8,3.4],[3.8,3.3],[3.8,3.22]])

#cluster2[:,0] = cluster2[:,0] + 1.5
#cluster3[:,0] = cluster3[:,0] + 3.0
dataLetters = np.vstack([cluster1,cluster2,cluster3])
#dataLetters = np.vstack([cluster1,cluster2])
n1,d1 = np.shape(cluster1)
n2,d2 = np.shape(cluster2)
n3,d3 = np.shape(cluster3)
dataLettersLabels = np.hstack([np.zeros(n1,dtype=int),np.zeros(n2,dtype=int)+1,np.zeros(n3,dtype=int)+2])
#dataLettersLabels = np.hstack([np.zeros(n1,dtype=int),np.zeros(n2,dtype=int)+1])

#############################  
### create network data 1 ###  
#############################   
group1 = ['1','2','3','4','5','6','7',"J"]
group2 = ['9','10','11','12','13','14','15',"K"]

G = nx.Graph()
G.add_edge('1','2')
G.add_edge('1','3')
G.add_edge('1','5')
G.add_edge('1','6')
G.add_edge('2','3')
G.add_edge('3','4')
G.add_edge('4','5')
G.add_edge('7','J')
G.add_edge('J','K')
G.add_edge('K','9')
G.add_edge('7','1')
G.add_edge('9','10')
G.add_edge('10','11')
G.add_edge('10','12')
G.add_edge('12','15')
G.add_edge('12','14')
G.add_edge('11','15')
G.add_edge('11','13')
G.add_edge('13','15')

dataNetwork = G
dataNetworkLabels = np.zeros(len(G.nodes()),dtype=int)
for i in range(len(G.nodes())):
    n = G.nodes()[i]
    if n in group2:
        dataNetworkLabels[i] = 1
