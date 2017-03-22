"""
Python 3 

Created by L. Brochini for Neuromat on Feb 2017

This program is part of the Supplementary material to the paper:

 "Interaction graph estimation for the first olfactory relay of an insect"
L. Brochini, P.Hodara, C.Pouzat, A.Galves 

__________________
Generates "Pruning" Figure.

Creates a pair of plots that correspond to "before pruning" and "after pruning" graph estimation.
The goal of pruning is to get rid of inconclusives.

Prunning procedure:
Run the GE for all neurons.
Then for each post synaptic neuron: if, for all presynaptic candidates there is at least one absent connection AND one 
inconclusive connection, then, add one of the nerons corresponding to the absent connection to a list of neurons to be pruned.
Next, rerun GE for the original set, minus the set of pruned neurons.
Then, do consecutive prunnings are done until there is nothing more to prune

"""

import GE
import GEplots
import genSimData
import pickle

# Set parameter used in graph estimation (GE)
epsilon=0.05 
xi=0.001 

# Create very simple network of 10 connected neurons
Wmat=genSimData.setMatrixNet2()
N=len(Wmat)
GElist=[] # List of graph estimations for no pruning, and consecutive prunings

# generate sample : More neurons and smaller sample size: pruning should be helpful in this case
D=genSimData.GenSample(Wmat,200000,Pspont=0.04) 
prunables=[[] for neuron in range(N)] # list (for each postsyn) of lists of presynaptic candidate neurons to be pruned at each step

GElist.append(GE.main([D['X']],xi,epsilon,details=False,upperlim=30)) # First, runing the GE without pruning
Titlearr=["No prune"]
prunecount=0

# Now using the same sample, we re-run the GE algorithm, excluding the list of pruned neurons at each pruning step. 
# The list increases at each step.
# When the list ceaces to increase, we obtain the final GE
while True: 
    needprune=False
    
    for post in range(N):
        GElist[-1][post,post]=-3 # ... just so there is no 0 in the diagonal
        col=list(GElist[-1][:,post])    
        if 0 in col and -2 in col: # if there is an absent (0) connection and an inconclusive (-2)
            inds= [i for i, j in enumerate(col) if j == 0 if i not in prunables[post]] # list of "no connections" in a column where there is an inconclusive
            if len(inds)>0: 
                needprune=True
                prunables[post].append(inds[0]) # take the first presyn candidate that has "no connection" to postsyn neuron and add it to the list of neurons to be pruned from GE estimation of this specific postsyn neuron
    
    if needprune:
        prunecount+=1
        print("\n\nPrune number:",prunecount,"\nNow prunning:",prunables)
        GElist.append(GE.main([D['X']],xi,epsilon,details=False,upperlim=30,PrunnedNeurons=prunables))   #recompute GE without pruned neurons     
        Titlearr.append("Prune "+str(prunecount))
    else:
        break
    
print("Number of GE to end pruning=",prunecount)
figname="PruningN10"   
GEplots.PlotLineDiffMaps(Wmat,[GElist[0],GElist[-1]],["Before pruning","After pruning"],figname,showax=True)

pickle.dump({'GElist':GElist,'Wmat':Wmat},open('Pruning.pkl','wb'))

     
