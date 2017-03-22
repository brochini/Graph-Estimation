'''
Python 3

Created by L. Brochini for Neuromat on Feb 2017

This program is part of the Supplementary material to the paper:

"Interaction graph estimation for the first olfactory relay of an insect"
L. Brochini, P.Hodara, C.Pouzat, A.Galves 

Routine use by Projections program.

If SpecifiPair is set to None, then the program computes the graph estimation for all possible subsets of size N of a set of neurons [0,1,...,TOTALN-1].
If a SpecifiPair is given, then the program computes the graph estimation for all subsets of size N that contain {i,j}.

'''

import GE
import numpy as np
import itertools


#%%

def main(Xdat,TOTALN,N,epsilon, xi,SpecificPair=None): # Xdat is a given sample dictionary as provided by genSample function

    
    Wmatrix=Xdat['Wmatrix']
    if SpecificPair ==None: 
        comb=list(itertools.combinations(range(TOTALN),N)) #all possible subsets
    elif len(SpecificPair)==2:
        comb=[it for it in list(itertools.combinations(range(TOTALN),N)) if SpecificPair[0] in it if SpecificPair[1] in it] # all subsets that contain {i,j}

        
    V={} # dictionary of graph estimations for each subset in comb, where keys are tuples corresponding to the subset and values are dictionaries of GE routine output
    
    print("Now running Graph Estimation procedure to ", len(comb), " subsets of size ", N, " from a set of size", TOTALN,'neurons')

    
    for gind in range(len(comb)):
        
        GoodNeurons=comb[gind]    
        print('\n\n\n--------Subset #',gind,'corresponds to neurons',GoodNeurons,'\n\n')        
       
        X=Xdat['X'][GoodNeurons,:] # sample only for subset of neurons
        auxW=Wmatrix[GoodNeurons,:]
        subWmatrix=auxW[:,GoodNeurons]

        V[GoodNeurons]=GE.main([X],xi,epsilon,details=True, upperlim=20) #upperlim reduced in this case because it 
        V[GoodNeurons]['Wmatrix']=subWmatrix
            
#%%  includes maxDj into dictionary:
            
        deltamax=V[GoodNeurons]['maxDeltajPOSTPRE']
        maxDj=np.zeros([N,N])#inverting because recorded values are post-pre. I need pre-post to compare
        
        for postneuron in range(N):
            for preneuron in range(N):
                maxDj[preneuron,postneuron]=deltamax[postneuron][preneuron]
    
        V[GoodNeurons]['maxDj']=maxDj
                        
#%%Returns dictionary of Subset:GE    
    return(V)



