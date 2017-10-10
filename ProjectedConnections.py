"""
Python 3

Created by L. Brochini for Neuromat on Feb 2017

This program is part of the Supplementary material to the paper:

 "Interaction graph estimation for the first olfactory relay of an insect"
L. Brochini, P.Hodara, C.Pouzat, A.Galves 
_______________________________
Generates "Projected Connections" Figure.

This procedure estimates the complete connectivity graph of a network of TOTALN neurons 
by computing the estimated graph only the for a subset of neurons at a time. Then, results are compiled to produce 
the final graph where the nature of connections is discriminated between direct 
connections (when j->i is identified for all subsets containing i and j)
absent connection (when j->i is identified as absent for all subsets containing i and j),
projected connection of first or second degree (explained below)  or inconclusives

defining:
 - I a the whole set of neurons
 - Delta_ij_F is Delta(i,n)(j) (see def. in paper) computed for the subset F where {i,j} in F
 (Note that, when computed for the whole graph, it is compared to the sensitivity measure epsilon
 i.e., when  Delta_ij_I >epsilon means that we accept j as a presynaptic neuron to i )
- D_ji_s it the set of all Delta_ij_F for all F in I where {i,j} belongs to F and |F|=s

The criterion for identifying a projection is:
if D_ij_s has at least one value below AND at least one value above epsilon, AND there 
is a path of s connections from j to i for the set where the j->i is identified as absent, then we say there is a projected connection 
of order s-2 from j to i

"""





#%%

import GEplots
import numpy as np
import genSimData
import GEallsubsets # computes GE for subsets of certain size
import time
import pickle
import pylab
from matplotlib import pyplot as plt
#%%
start = time.time()

TOTALN=10 #total number of neurons
Wmatrix,muNet2,PspontNet2=genSimData.setMatrixNet2() # generates simple connectivty matrix

Xdat=genSimData.GenSample(Wmatrix,200000,mu=muNet2,Pspont=PspontNet2) #generates sample process

genSimData.PlotRaster(Xdat['X'],"RasterN10",myaspect=20)
epsilon=0.1
xi=0.001
#%%
# GEsubs3 is a dictionary to be used in identifying first order projections
# Each key is a tuple with the index of 3 neurons and the correponding value is the 
# the estimated graph dictionary for the given sample for only the selected neuron in key.
# The keys correspond all 3 by 3 combinations of the entire set of 10 neurons.
# The 'maxDj' is subkey containing values for each link in the subset passed as key

GEsubs3=GEallsubsets.main(Xdat,TOTALN,3,epsilon, xi)


D_ji_3=[[np.nan]*(TOTALN) for i in range(TOTALN)]
DC=[list(np.ones(TOTALN)*False) for i in range(TOTALN)] # Direct connections map of conections j->i : lines correspond to presybaptic and cols to postsynaptic : attention inverted order from for loop
UC=[list(np.ones(TOTALN)*False) for i in range(TOTALN)] # Undirected connections
PC=[list(np.zeros(TOTALN,dtype=int)) for i in range(TOTALN)] # Projected Conections Graph 0: not a projection, 1: first order projection, 2: second order projection, nan: inconclusive
FirstProjCandidates=[] # connections candidate to projections
#%%
for i in range(TOTALN): # for each postsynaptic neuron i
    presyn=list(range(TOTALN))
    presyn.remove(i) # list of presynaptic candidates
    for j in presyn:
        klist=[k for k in range(TOTALN) if k not in [i,j]] #k= F \ {i,j} in crescent order for F in I

        #F3list is the list of F subsets of size 3 sets ordered by the value of k.
        #Each set is represented by an ordered tuple to be fed as keys to the GEsubs3 dictionary
        F3list=[tuple(sorted([i,j,k])) for k in range(TOTALN) if k not in [i,j]]


        D_ji_3[i][j]=[GEsubs3[key]['maxDj'][key.index(j),key.index(i)] for key in F3list ]# list of Delta_ij_F values for F in F3list
        if np.nanmax(D_ji_3[i][j])<=epsilon: # if maximum is less than epsilon connection is False
            DC[j][i]=False
        elif np.nanmin(D_ji_3[i][j])>epsilon: # if minimum is greater than epsilon connection is True
            DC[j][i]=True
        elif np.nanmax(D_ji_3[i][j])>epsilon and np.nanmin(D_ji_3[i][j])<epsilon:
        #    FirstProjCandidates.append([i,j]) # if it crosses epsilon line is a candidate to projection
            PC[j][i]=1 # if it crosses epsilon line is a candidate to projection
        elif False not in np.isnan(D_ji_3[i][j]): # all entries are nan values: completely inconclusive
            PC[j][i]=np.nan
        else:
            PC[j][i]=2 #other unforeseen case
            

#%%
end = time.time()
print('\nTotal time elapsed %.1f'%(end - start),'s') 

#%%
maxdelta=np.nanmax([np.nanmax(D_ji_3[i][j]) for i in range(10) for j in range(10)])
f, axarr = plt.subplots(TOTALN,TOTALN,figsize=(12,12))
for i in range(TOTALN):
     for j in range(TOTALN):
        axarr[i][j].plot(D_ji_3[j][i],'o-')
        axarr[i][j].plot([0,TOTALN-1],[epsilon,epsilon],'r')
        axarr[i][j].set_ylim([0,maxdelta])
        

pylab.savefig('MaxDeltas.pdf', format='pdf', dpi=1000)   

figname='Projections'
#titletext='True=Grey, False=White, 1st proj=Green, 2nd proj=Purple, Nan=Cyan'
GEplots.PlotOneMapProjections(DC,PC,figname,' ')



pickle.dump({'PC':PC,'DC':DC,'Dji3':D_ji_3,'epsilon':epsilon,'xi':xi},open('AnyProjConn.pkl','wb'))

#%%
