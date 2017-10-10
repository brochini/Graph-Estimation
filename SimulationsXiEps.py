"""
Python 3

Created by L. Brochini for Neuromat on Feb 2017

This program is part of the Supplementary material to the paper:

 "Interaction graph estimation for the first olfactory relay of an insect"
L. Brochini, P.Hodara, C.Pouzat, A.Galves 

___________________________________________________
Creates Figure 1.

The aim of this procedure is to evaluate the parameters epsilon and xi used for graph estimation for a specific sample size.
To do so, we create a network with known connections, generate a sample of the process of GL neurons connected with the 
chosen matrix and then evaluate if the estimator returns the correct graph with chosen parameters.

"""

import GE #Graph estimation module
import GEplots #Plotting module
import genSimData # Simulation of Galves-Locherbach neurons
import pickle

epsilonrange=[0.01, 0.05, 0.1, 0.15] # Sensitivity parameter range to be tested
xirange=[0.001, 0.01, 0.1, 0.15] # cutoff parameter range to be tested
V={}
nsteps=1000000
Wmat,muNet1,PspontNet1=genSimData.setMatrixNet1() # Preselected connectivity matrix. To change it, simply create an np.array of NxN with values in the range [-1,1]
D=genSimData.GenSample(Wmat,nsteps,mu=muNet1,Pspont=PspontNet1) # Provides sample of process of GL neurons connected with the given connectivity matrix and for n time steps with leakage parameter mu 

for epsilon in epsilonrange:
    for xi in xirange:
        print('\n\nParameters: xi=',xi,' epsilon=',epsilon)
        V[(epsilon,xi)]=GE.main([D['X']],xi,epsilon,details=False) #returns only the GE matrix
        
pickle.dump({'V':V,'Wmat':Wmat,'epsilonrange':epsilonrange,'xirange':xirange,'nsteps':nsteps},open('SimdataNet1.pkl','wb'))
GEplots.PlotMultDifs(Wmat,V,epsilonrange,xirange,"SimDataEpsXi"," ")  
