"""
Python 3

Created by L. Brochini for Neuromat on Feb 2017

This program is part of the Supplementary material to the paper:

 "Interaction graph estimation for the first olfactory relay of an insect"
L. Brochini, P.Hodara, C.Pouzat, A.Galves 
_______________________________
Generates "Projected Connections" Figure.

This procedure estimates the complete connectivity graph of a network of TOTALN neurons 
by computing only the GE for a subset of neurons at a time. Then, results are compiled to produce 
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

start = time.time()

TOTALN=10 #total number of neurons
Wmatrix=genSimData.setMatrixNet2() # generates simple connectivty matrix 
Xdat=genSimData.GenSample(Wmatrix,1000000,Pspont=0.04) #generates sample process

genSimData.PlotRaster(Xdat['X'],"RasterN10",myaspect=20)
epsilon=0.05
xi=0.001

# GEsubs3 is a dictionary to be used in identifying first order projections
# Each key is a tuple with the index of 3 neurons and the correponding value is the 
# the estimated graph dictionary for the given sample for only the selected neuron in key.
# The keys correspond all 3 by 3 combinations of the entire set of 10 neurons.
# The 'maxDj' is subkey containing values for each link in the subset passed as key
 
GEsubs3=GEallsubsets.main(Xdat,TOTALN,3,epsilon, xi)



#%%
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


        D_ji_3=[GEsubs3[key]['maxDj'][key.index(j),key.index(i)] for key in F3list ]# list of Delta_ij_F values for F in F3list
        if np.nanmax(D_ji_3)<=epsilon: # if maximum is less than epsilon connection is False
            DC[j][i]=False
        elif np.nanmin(D_ji_3)>epsilon: # if minimum is greater than epsilon connection is True
            DC[j][i]=True
        else:
            FirstProjCandidates.append([i,j]) # if it crosses epsilon line is a candidate to projection

#%%
for i in range(TOTALN): # Build undirected graph of identified direct connections.
    for j in range(TOTALN):
        UC[i][j]=(DC[i][j] or DC[j][i])

SecondProjCandidates=list(FirstProjCandidates)
for i,j in FirstProjCandidates:
    klist=[k for k in range(TOTALN) if k not in [i,j]]
    for k in klist:
        key=tuple(sorted([i,j,k]))
        dji=GEsubs3[key]['maxDj'][key.index(j),key.index(i)] # = Delta_ij_F for F={i,j,k}
        if dji<=epsilon:# if it is identified as absent connection for this subset ... 
            if (UC[k][j] and UC[k][i]) or (UC[j][k] and UC[k][i]) or (UC[i][k] and UC[k][j]):#... AND there are 2 direct connections path from j to i within this subset... 
                PC[j][i]=1 # ... then we consider it a first order projection
                SecondProjCandidates.remove([i,j])
                break

#%%            
OtherCases=list(SecondProjCandidates)

for i,j in SecondProjCandidates:
    GEsubs4=GEallsubsets.main(Xdat,TOTALN,4,epsilon, xi,SpecificPair=[i,j]) # runs GE for all 4 neurons subsets that contain i and j
    
    klist=[k for k in range(TOTALN) if k not in [i,j]]
    for k in klist:
        llist=[l for l in range(TOTALN) if l not in [i,j,k]]
        for l in llist:
            key=tuple(sorted([i,j,k,l]))
            dji=GEsubs4[key]['maxDj'][key.index(j),key.index(i)] # = Delta_ij_F for F={i,j,k,l}
            if dji<=epsilon: # if it is identified as absent connection for this subset ... 
                if (UC[j][k] and UC[k][l] and UC[l][i]) or (UC[j][l] and UC[l][k] and UC[k][i]): #... AND there are 3 direct connections path from j to i within this subset... 
                    PC[j][i]=2 # ... then we consider it a second order projection
                    OtherCases.remove([i,j])
                    break
        if PC[j][i]==2:
            break

for i,j in OtherCases:
    PC[j][i]=np.nan

figname='Projections'
titletext='True=Grey, False=White, 1st proj=Green, 2nd proj=Purple, Nan=Cyan'
GEplots.PlotOneMapProjections(DC,PC,figname,titletext)

end = time.time()
print('\nTotal time elapsed %.1f'%(end - start),'s') 


pickle.dump({'PC':PC,'DC':DC},open('ProjConn.pkl','wb'))

#%%
