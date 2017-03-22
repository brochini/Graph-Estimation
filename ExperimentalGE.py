"""
Python 3

Created by L. Brochini for Neuromat on Feb 2017

This program is part of the Supplementary material to the paper:

 "Interaction graph estimation for the first olfactory relay of an insect"
L. Brochini, P.Hodara, C.Pouzat, A.Galves 
_____________________________________

This program generates the influence graph estimation figure for exeperimental data

Dataset used here corresponds to zenodo locust20010217 spontaneous tetrode D
data should be placed in subfolder:
LocustData/locust20010217_spont_tetD_uX.txt where X corresponds to well isolated neurons X=[1,2,3,4,7]

"""

import GE
import GEplots
import numpy as np 
from matplotlib import pyplot as plt
import pylab
import pickle

#%%
def SpkTimesOneTrial(GoodNeurons,fname):
    '''Gests spike times series from txt files and generates list of lists with each neurons spiketimes'''
    N=len(GoodNeurons)
    spktimes=[[[]] for neuron in range(N)] 
    
    for neuron in range(N):
        
        with open(fname+str(GoodNeurons[neuron])+".txt","r") as fp:        
            while True:
            
                line=fp.readline()
                if not line: break
                spktimes[neuron][0].append(float(line))

    return(spktimes)
                
#%%
                
def digitize(rate,N,totaltrials,spktimes):
    ''' Creates series of 0 and 1's based on spike times. rate= discretization window size'''
    X=[]
    nstepstotal=0
    superpos=[0 for neuron in range(N)]
    for trial in range(totaltrials):
        
        maxtime=max([spktimes[neuron][trial][-1] for neuron in range(N)])
        mintime=min([spktimes[neuron][trial][0] for neuron in range(N)])
        nsteps=int((maxtime-mintime)/rate)
        X0=np.zeros([N,nsteps+1],dtype=int)    
        nstepstotal+=nsteps    
        for neuron in range(N):
        
            timeint=[int((it-mintime)/rate) for it in spktimes[neuron][trial]];
            for el in timeint:
                if X0[neuron][el]==1:
                    superpos[neuron]+=1
                else:    
                    X0[neuron][el]=1

        X.append(X0)   
        
    return X,nstepstotal,superpos
#%% 

def PlotSuperposPerWindow(figname,myraterange,spktimes,myylim,ylab,xlab,tit):
    '''Plots Superposition ratio per discretization window size  '''
    
    N=len(spktimes) # number of neurons in set
    
    colors=['bo-','ro-','go-','mo-','co-','yo-','ko-','bo-','ro-','go-','mo-','co-','yo-','ko-'] 
    legend=[str(neuron) for neuron in range(N)]
    superposed=[[] for neuron in range(N)]
    for rate in myraterange:
        X,nstepstotal,superposrate=digitize(rate,N,totaltrials,spktimes)
        for neuron in range(N):        
            superposed[neuron].append(superposrate[neuron]/len(spktimes[neuron][0]))
     
    f1=plt.figure(figsize=(10,7))
    for neuron in range(N):
        plt.plot(myraterange,superposed[neuron],colors[neuron])
    
    plt.plot(myraterange,np.ones(len(myraterange))*0.01,'--')    
    plt.legend(legend,loc=2)
    plt.ylim(myylim)
    plt.xlim([min(myraterange),max(myraterange)])
    plt.ylabel(ylab)
    plt.xlabel(xlab)
        
    f1.suptitle(tit)    
    pylab.savefig(figname+'.eps', format='eps', dpi=1000)
    
    
#%%      
    
if __name__ == '__main__':
    
    
    fname="LocustData/locust20010217_spont_tetD_u"
    dataname='Loc217' # contains all trials from all experiments of spontaneous activity (same animal)
    
    GoodNeurons=[1,2,3,4,7] # for this dataset, these are the units which are best isolated. See nice ISI distributions below.
    N=len(GoodNeurons) # number of neurons in set
    totaltrials=1 # in this format all spktimes are saves as if there were only one trial
    xi=0.001 # xi parameter choice, see simulation results section
    epsilon=0.05 # epsilon parameter choice, see simulation results section
    acquisitionrate=15000 # Hz
    
    spktimes=SpkTimesOneTrial(GoodNeurons,fname) #create N lists with each neurons spike times

#%% Plotting ISI histograms for each neuron with useful informations
### Note distributions look nice for these neurons
    
    # create list of ISI of size up to xlim acquisition points 
    xlim=6000 
    ISIlim=[[] for neuron in range(N)] 
    for neuron in range(N):
        ISIlim[neuron]=[it for it in np.diff(spktimes[neuron][0]) if it <xlim]
    
    maxtime=max([spktimes[neuron][0][-1] for neuron in range(N)])
    mintime=min([spktimes[neuron][0][0] for neuron in range(N)])
    #Plot ISI histograms
    GEplots.PlotISIHist(ISIlim,200,xlim,dataname+'ISI','ISI histograms. Dataset 20010217 spont tetD. Total recording time='+str(int(maxtime/15000))+'s',GoodNeurons,acquisitionrate)

#%% Plotting superposition ratio as a function of the size of discretization window
# discretization window in acquisition time units. 

    figname=dataname+'_superpositions'    
    windowrange=range(1,200,1)
    ylab='($\#$ bins with superposed spikes)/(total $\#$ of spikes)'
    xlab='discretization window size (acquisition units= x6.67 e-5 s)'
    tit="Spikes superposition ratio per discretization window size"
    PlotSuperposPerWindow(figname,windowrange,spktimes,[0,0.015],ylab,xlab,tit)

    
    
#%% Choose discretization window size (in acquisition time units) as the one that roughly produces at most  1% superpositions    
    window=155 # in sampling time units 
    
    spkpart1=[] # first half of this dataset
    spkpart2=[] # second half of this dataset
    for neuron in range(N):
        aux1=[it for it in spktimes[neuron][0] if it< (maxtime-mintime)/2+mintime]
        aux2=[it for it in spktimes[neuron][0] if it>= (maxtime-mintime)/2+mintime]
        spkpart1.append([aux1])
        spkpart2.append([aux2])

    spktimes=[spktimes,spkpart1,spkpart2]
    Titlearr=['Whole Data','First half', 'second half']
    
    Vlist=[]
#%% Calling graph estimation procedure and saving figures
    for i in range(3):
        X,nstepstotal,superposrate=digitize(window,N,totaltrials,spktimes[i])
        Vlist.append(GE.main(X,xi,epsilon,details=False,upperlim=50))  #calls GE procedure. V is dictionary containing the estimated graph for dataset

    GEplots.PlotLineDiffMaps([0],Vlist,Titlearr,'ExperimentalGE',showax=True,experimental=True)      

    pickle.dump({'Vlist':Vlist,'window':window,'Titlearr':Titlearr},open('ExperimentalGE.pkl','wb'))
