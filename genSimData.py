
"""
Python 3

Created by L. Brochini for Neuromat on Feb 2017

This program is part of the Supplementary material to the paper:

 "Interaction graph estimation for the first olfactory relay of an insect"
L. Brochini, P.Hodara, C.Pouzat, A.Galves 
_____________________________________

This program is used to simulate neuronal networks with stochastic dynamics as described in Galves and Locherbach 2013
This is the specific case where leakage is geometric (param. mu), leading to a markov chain with respect to potential (Brochini et al 2016)

Function SetMatrixNet1( and 2) can be used to determine connectivity matrices used to produce figure results in simulation section.
GenSample generates sample of GL neurons activity with given connection where the firing function (phi) is linear saturating of the type:
phi(U)=Pspont+U if U<1-Pspont and 
phi(U)=1 if U>= 1-Pspont
Returns the X (Nxn np.array of zeros and ones) series of spikes for the N neurons in n time steps

"""


import numpy as np
import matplotlib.pyplot as plt 
import pylab
import time


def setMatrixNet1():
    N=5
    Wmatrix=np.zeros([N,N])

    Wmatrix[0]=[0,0,0.1,0,0] 
    Wmatrix[1]=[0.1,0,0.3,0.4,0]
    Wmatrix[2]=[0,0.4,0,0.8,0]
    Wmatrix[3]=[0.3,0,0.1,0,0.5]
    Wmatrix[4]=[0.2,0,0.8,0,0]
    return Wmatrix

def setMatrixNet2():
    N=10
    w=0.5
    Connlist=['10','04','23','35','56','87','89'] # list of connections 1->0 etc..
    Wmat=np.zeros([N,N])
    for conn in Connlist:
        Wmat[int(conn[0]),int(conn[1])]=w
    return(Wmat)


def GenSample(Wmatrix,nsteps,mu=0.9,Pspont=0.02):
    start = time.time()      
    print("\n\nSimulating a network of GL neurons with the provided connectivity matrix.")
    print("Sample size will be ", nsteps)
     
    N=len(Wmatrix)
    Vn=np.zeros(N)
    Xn=np.zeros(N)
    
    Vn[0]=1
    Xn[0]=1
    
    X=np.zeros((N,nsteps), dtype=np.int)
    V=np.zeros((N,nsteps+1))
    
            
    def randcompare(x,r):
        if x+Pspont>r: # spontaneous firing rate is Pspont. added to the typical linear saturating: Phi(U)=Pspont+U
            return 1
        else:
            return 0

    
    for n in range(nsteps):   
         
        R=np.random.rand(N)
        Xn=np.array([randcompare(x,r) for (x,r) in zip(Vn,R)])   
        
        Vnplus1=np.squeeze(np.array(mu*Vn + np.dot(Xn,Wmatrix)))*(1-Xn)     
        
        Vn=Vnplus1
        
        X[0:N,n]=[int(elem) for elem in Xn]
        V[0:N,n+1]=Vnplus1    
        
    end = time.time()
    print('\n\n Time elapsed', end-start,'\n\n')        
    
    return {'X':X,'nsteps':nsteps,'mu':mu,'Pspont':Pspont,'Wmatrix':Wmatrix} 



def PlotRaster(X,figname,myaspect=20):
    
    print("\nPlase check a Raster plot representation in ",figname,".eps")    
    plt.matshow(X,cmap='Greys',  interpolation='nearest',aspect=myaspect)
    plt.xlim(100,min(1000,X.shape[1]))
    pylab.savefig(figname+'.eps', format='eps', dpi=1000)
    
    


