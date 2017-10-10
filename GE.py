"""
Python 3

Created by Ludmila Brochini for Neuromat on Feb 2017

Statistical model selection of influence graph in a stochastic neural network.
This is the Influece graph estimation routine part of the Supplementary material to the paper:

 "Interaction graph estimation for the first olfactory relay of an insect"
Ludmila Brochini, Pierre Hodara, Christophe Pouzat, Antonio Galves 

____________________________

This program provides the interaction graph estimation (GE) from a given spike train

_____________________________

When called, main takes 4 sequential arguments:
    X               : list of np.arrays. Each element of the list corresponds to one trial data.
                      Each np array is a N x nsteps matrix of zeros and ones(1=spike,0=not spike) where 
                      N=number of neurons and nsteps=sample size
    xi              : xi parameter value
    epsilon         : epsilon parameter value
    details         : Default value is True. if set to False returns only Vest, the matrix of estimated connectivity. 
                      If True returns a dictionary with detailed outputs where Vest is a key.
    upper_lim       : maximum lenght of w to be verified. Can be set arbitrarly large with no consequence other than runtime increase. 
                      do not use small values that will cause to miss statistically important events. Default value is set to 50

main returns a dictionary where the key 'Vest' carries the estimated graph.
   Vest columns correspond to the postsynaptic neuron and lines to candidate presynaptic neurons.
   each element has one of the following values:
   0 : no connection
   1 : there is a connection
   -2: inconclusive due to lack of statistics. Events do not occur often enough to be able to exclude of accept j as presynaptic neuron


Procedure description:

For each pair of neurons (i,j), the method proposes to find if j is  presynaptic  to i
To do so the method infers if the history of j (since i fired for the last time up to time t)
 is relevant to predict if i will fire or not at t+1.

V is estimated through 

- X is an array of zeros and ones where 1 corresponds to the occurence of a spike 
- Wmatrix is the connectity matrix of the simulated data

W(i,l)= {w in {0,1}^({-l,...,-1} x Fi) } = set of all possible contexts in the matrix 
of zeros and ones where lines correspond to neurons in Fi and l columns, where l is the number of
 time steps in the past one has to take to find a spike of neuron i


"""

#%%
import numpy as np
import itertools
import time



class SpikeData: 
    """
    Data spike times, ISI etc
    
    usage:
    OBJ=SpikeData(X), where X i and N x nsteps 
    matrix of zeros and ones(1=spike,0=not spike)
    where N=number of neurons and nsteps=sample size

    Useful attibutes:
        spktime      : corresponds to [N]array of timestamps
        ISI          : [N]array of interspike intervals. The last value is a fake ISI, 
                       just records the time difference between the last spike and the end of the sample
        minstart     : first instant all neurons have fired 
        Fi           : [N]array of set of neurons excluding the postsynaptic neuron: Fi[i]= F\i 
        preneuron_map: zip of the index of a presynaptic neuron with the corresponding neuron in Fi
    Useful methods:    
        validISIlimitCumu: Based on cutoff parameter. Returns argmax(l,N(ISI>=l)>minfreq) given minfreq=nsteps**(0.5+xi)
    """

    def __init__(self, X,subsetpostneuron):# X is the np.array of {0,1} with lines=N and columns=nsteps
        self.X=X
        self.spktime=[]        #contains all instants that each neuron fired   
        self.minstart=0  # minimum time when each starts firing
        self.ISI=[]    #maximum ISI of each neuron determines the size of W(i,l)
        self.nsteps=self.X.shape[1] # sample size
        self.N=self.X.shape[0]      # number of neurons

        self.create()
        self.calc_minstart()
        self.buildISI()
        self.append_last_fake_ISI()
        self.F=tuple(range(self.N))
        self.Fi=[]
        self.Fialt=[]
        self.preneuron_map=[[] for i in subsetpostneuron]                    
        self.maxisi=[max(self.ISI[i]) for i in self.F]        
        
        for post_neuron in self.F:
            aux=list(self.F)
            aux.pop(post_neuron)
            self.Fi.append(tuple(aux))
            
        for post_neuron in range(self.N):
            aux=list(subsetpostneuron)
            aux.remove(subsetpostneuron[post_neuron])
            self.Fialt.append(tuple(aux))
            self.preneuron_map[post_neuron]=list(zip(range(len(self.Fialt[post_neuron])),self.Fialt[post_neuron]))
        
        
    def create(self):   # identifies spike times from X
        Laux=[]
        for i in range(self.N):
            t=0;
            l=[]
            while t<self.nsteps:
                if self.X[i][t] ==1:
                    l.append(t)
                t+=1
            Laux.append(np.array(l))
        self.spktime=np.array(Laux)
        
    def calc_minstart(self): #first instant all neurons have fired
        
        
        self.minstart=[min(self.spktime[i])+1 for i in range(self.N)]

    def buildISI(self): # creates attribute matrix of interspike intervals
        
        for i in range(self.N):
            isi=np.diff(self.spktime[i])
            self.ISI.append(isi)
            
    def append_last_fake_ISI(self): #appends a fake value at the end of ISI list to help procedure reach after last spike 
        for i in range(self.N):
            self.ISI[i]=np.append(self.ISI[i],len(self.X[i])-max(self.spktime[i])-1)


       
def validISIlimitCumu(ISI,minfreq,upperlim): 
                                    #returns limISI=argmax(l,N(ISI>=l)>minfreq).
                                    # in order to prevend dealing with huge l's, values are capped to l=50 or chosen value
                                    
    N=len(ISI)                                 #this is useful to pre-clean data since no N(w)_i can be greater than minfreq if
    count=[[] for neuron in range(N)]                        
    validISI=[]
    
    validISIcum=[]
    limISIcum=[]
    for neuron in range(N):
        validISI.append([])        
        D=ISI[neuron]
        count[neuron]=np.zeros(max(D)+1)
        for i in range(len(D)):
            count[neuron][D[i]]+=1

    cumul=[[] for neuron in range(N)]
    for neuron in range(N):
        validISIcum.append([])
        for l in range(len(count[neuron])):
           cumul[neuron].append(np.sum(count[neuron][l:]))
        for isi in range(1,len(cumul[neuron])):
            if cumul[neuron][isi]>=minfreq:
                validISIcum[neuron].append(isi)
    
    for neuron in range(N):
        if len(validISIcum[neuron])==0:
            limISIcum.append(1)
        else:
            limISIcum.append(min(max(validISIcum[neuron]),upperlim))
    
    return(limISIcum)    



class Wcollection:
    
    """
    Collection of all w for lines in Fi with l columns for a postsynaptic neuron i 
    and corresponding conditional probabilities, number of occurences, etc
    
    usage:
    OBJ=Wcollection(post_neuron,l) 
    
    Useful attibutes:
        count      : dictionary of {w:[N1w,N0w]}, where N1w (or N0w) is the number of times
                     a X[i]=1 (or 0) is observed after w in Fi x l
        
        prob       : [N-1]list of nested dictionries. Each element in prob corresponds to a presynaptic neuron
                     prob[index_preneuron]={wjc:{w:[p1w,Nw]}} where Nw is the total number of occurences 
                     of w in Fi x l, p1w is the conditional probability p(1|w) and wjc corrresponds to 
                     w with a deleted line corresponding to the candidate presynaptic neuron j
    Useful methods:
        count_all_lsize_sequences_after_i_fired : takes SpikeData object as argument.
        
        
    Obs:
    Wcollect= {w in {0,1}^({-l,...,-1} x Fi) } = set of all possible contexts in the matrix 
    zeros and ones where lines correspond to neurons in Fi and l columns, where l is the number of
    time steps in the past one has to take to find a spike of neuron i
    
        
    """

    def __init__(self,post_neuron,l):
        
        
        self.l=l                     #length in time of each w   
        self.post_neuron=post_neuron # ith neuron
        
        
        self.count={}
        self.prob=[]
        
         
    def add_event_to_w(self,w,a):     # registers N(a|w), a being 0 or 1  
                                      # count={w:[N(0|w),N(1|w))]}    
       
        wint=[]
        for k in range(len(w)):
            wint.append(np.sum(np.array([2**i*w[k][i] for i in range(len(w[k]))])))
        wtup=tuple(wint)
        
        if wtup not in self.count:
            self.count[wtup]=[0,0]              
        self.count[wtup][a]+=1 
       
       
    def count_all_lsize_sequences_after_i_fired(self,data,NeuronSubset): 
                
        l=self.l
        localpostneuron=NeuronSubset[self.post_neuron].index(self.post_neuron)

        lenspk=len(data.spktime[localpostneuron]) #number of times i fired 
        for k in range(lenspk):                  # for all spike times   
            t=data.spktime[localpostneuron][k]  #time reference when i fired. w beggins in t+1 and goes up to t+l
            if t>=data.minstart[localpostneuron]:                 #after all neurons fires at least once                              
                if data.ISI[localpostneuron][k] >l: # ISI(t)<l garantees neuron i does not fire between t and t+l
                    
                    w=data.X[data.Fi[localpostneuron],t+1:t+l+1] #=X(Fi,t+1:t+l)
                    a=int(data.X[localpostneuron][t+l+1]) # =0 or 1
                    self.add_event_to_w(w,a)

                    

    
    def OrderByWjc(self,N):

        for index_preneuron in range(N-1):
            self.prob.append({})  
                

        for W in self.count.keys():
            for index_preneuron in range(N-1): # prob is a list of lenght equal to the 
                                                        # number of presynaptic neuron
                 
                wjc=tuple([i for (i,j) in zip(W,range(N)) if j != index_preneuron])
                #wjc= wj complementary: Let j be a cadidate to presynaptic neuron.
                #and index_preneuron is the index corresponding to the jth neuron in the Fi set
                # then wjc is equal to w minus the line corresponding to index_preneuron
                
                if  wjc not in self.prob[index_preneuron].keys():
                    self.prob[index_preneuron][wjc]={}  # prob[index_preneuron] has wjc as keys. 
                                                        
                    
                self.prob[index_preneuron][wjc][W]=self.count[W]    #prob[index_preneuron][wjc] has w as key which value
                                                                         # is a list of a list with the conditional probability of 
                                                                         # having a spike in neuron i and the number of total occurences of w   
     

    def clean_collection_and_calc_probs(self,min_freq):
        wlist=[W for W in self.count.keys()] #need to generate list beause cant erase entry in a loop over keys
        for W in wlist:
            NW=self.count[W][1]+self.count[W][0]
            if NW < min_freq:
                del self.count[W]
            else:
                p1w=self.count[W][1]/NW
                self.count[W]=[p1w,NW]
             
#%%             

def main(X,xi,epsilon,details=False,upperlim=50,PrunnedNeurons=None):
    
    
    start = time.time()
#%%%%%%%%%%%%%%%%%%%%%%%
#    Taking data entry X and  transforming in manageble data objects of class SpikeData 
#
#
#    N is the number of neurons in the sample
#   totaltrials is the number of trials in the dataset
#   nstepstotal is the sample size which is the total number of steps for all trials.       
#   X  is a list of np.arrays of size. The list has size = totaltrials. 
#   Each trial yield a different np.array of size N x nsteps.
#   Each np.array may have different sample size, but the same number of neurons
#%%%%%%%%%%%%%%%%%%%%%%%%
    
    N=len(X[0]) #number of neurons in the sample
    dataset=[]

    totaltrials=len(X)
    
    print('N=',N)
    NeuronSubset=[list(range(N)) for neuron in range(N)]
    totalISI=[] # gathers ISI for all trials to latter compute limISI
    
    if PrunnedNeurons!=None: # if there is a list of elements to be removed, alter NeuronSubset to be analysed
        for post_neuron in range(N):
            for el in PrunnedNeurons[post_neuron]:
                NeuronSubset[post_neuron].remove(el)
    for post_neuron in range(N):
        totalISI.append([[] for neuron in NeuronSubset[post_neuron]])
    
 


    nstepstotal=0
    for trial in range(totaltrials):
        data=[]
        for post_neuron in range(N):
            data.append(SpikeData(X[0][NeuronSubset[post_neuron],:],NeuronSubset[post_neuron]))
            
        dataset.append(data)
        nstepstotal+=dataset[trial][0].nsteps
    
        for neuron in range(N):
            for post_neuron in range(N):
                for neuron in range(len(NeuronSubset[post_neuron])):
                    for it in dataset[trial][post_neuron].ISI[neuron]:
                        totalISI[post_neuron][neuron].append(it)                 
    
    min_freq=nstepstotal**(0.5+xi) # xi parameter here defines the cutoff value min_freq 
                                   # which is the minimum amount of times a pattern should occur
                                   # to provide a useful conditional probability estimation
    
    limISI=[]
    for post_neuron in range(N):
        limISI.append(validISIlimitCumu(totalISI[post_neuron],min_freq,upperlim))
    

#%% Building a list of lists Nw of objects of type Wcollection.  Basically it collects the number events N(w|a), a={0,1} for each neuron i
# and each w is a specific firing pattern of lenght l of all other neurons after i fired
# Counting w occurences...

    Nw=[]
    maxl=[]
    
    print('Building collection of patterns...')
    firstdata=True
    trial=1
    for data in dataset:
        trial+=1

        for post_neuron in range(N):
            if firstdata:
                Nw.append([])
            
            for l in range(limISI[post_neuron][NeuronSubset[post_neuron].index(post_neuron)]):
                if firstdata:
                    Nw[post_neuron].append(Wcollection(post_neuron,l))
                Nw[post_neuron][l].count_all_lsize_sequences_after_i_fired(data[post_neuron],NeuronSubset)   
                
        firstdata=False
    key= list(Nw[0][0].count.keys())
    print('First neuron has ',np.sum(Nw[0][0].count[key[0]]) ,' spikes')    
    

#%% ... and using the cutoff value min_freq to get rid of events that don't happen often enough and estimate conditional probabilities:

    for post_neuron in range(N):
        auxmaxl=0
        for l in range(len(Nw[post_neuron])):
            Nw[post_neuron][l].clean_collection_and_calc_probs(min_freq)  
                
            if len(Nw[post_neuron][l].count.keys())>0:
                auxmaxl=l
        maxl.append(auxmaxl)   # is the maximum valid value of l for each neuron after cuttoff. No reason to look at patterns longer than l. 
   
    for post_neuron in range(N):
        del Nw[post_neuron][(maxl[post_neuron]+1):]
        for l in range(maxl[post_neuron]):
            Nw[post_neuron][l].OrderByWjc(N)
            
#%% Now estimating delta values and comparing with the sensitivity measure epsilon to accecpt or reject a connection.
        
    Tau=Nw;
 
    Vdic={}# dictionary of estimated connectivities {post_neuron:[preneuronA,preneuronB,...]}                            
    Vest=np.zeros([N,N],dtype=int)
    Deltaj=[[] for neuron in range(N)]
    maxDeltaj=[[] for neuron in range(N)]

#%%

    for post_neuron in range(N):
        indpostneuron=NeuronSubset[post_neuron].index(post_neuron)
        
        Vdic[post_neuron]=[]
        Deltaj[post_neuron]=[[] for neuron in range(N)]
        maxDeltaj[post_neuron]=[np.nan for neuron in range(N)]


#%%

        if PrunnedNeurons!=None:
                
            for el in PrunnedNeurons[post_neuron]:
                    maxDeltaj[post_neuron][el]=0
#%%
        for index_preneuron in range(dataset[0][post_neuron].N-1):
            preneuronlist=[j for i,j in dataset[0][post_neuron].preneuron_map[indpostneuron] if i==index_preneuron ]
           
            preneuron=preneuronlist[0]
            finalpreneuron=preneuron
            for l in range(maxl[post_neuron]):
                
                for wjc in Tau[post_neuron][l].prob[index_preneuron].keys():
                    
                    v_in_Tauwjc=Tau[post_neuron][l].prob[index_preneuron][wjc]  # all v in Tau_i_wj 
                    if len(v_in_Tauwjc)>=2: # needs at least 2 w for comparison
                        for pair in itertools.combinations(v_in_Tauwjc,2):
                            p1w=Tau[post_neuron][l].prob[index_preneuron][wjc][pair[0]][0]
                            p1v=Tau[post_neuron][l].prob[index_preneuron][wjc][pair[1]][0]
                            Deltaj[post_neuron][preneuron].append(abs(p1w-p1v))
                            
#%%

            if Deltaj[post_neuron][finalpreneuron]:
                maxDeltaj[post_neuron][finalpreneuron]=max(Deltaj[post_neuron][finalpreneuron])
                
            if maxDeltaj[post_neuron][finalpreneuron]>epsilon: # accept j as presynaptic neuron           
                # now find out to which presynaptic neuron the index_preneuron corresponds:
                 
                Vdic[post_neuron].append(finalpreneuron)    
                Vest[finalpreneuron,post_neuron]=1
            elif maxDeltaj[post_neuron][finalpreneuron]<=epsilon:
                Vest[finalpreneuron,post_neuron]=0
            else:# means its NaN
                Vest[finalpreneuron,post_neuron]=-2
                
#%%                
                
    end = time.time()
    print('Time elapsed %.1f'%(end - start),'s')    
  
     
    if details==True:
        return ({'Vest':Vest,'Nw':Nw,'limISI': limISI,'xi':xi,'epsilon':epsilon,'maxl':maxl,'nsteps':nstepstotal,'totaltrials':totaltrials, 'DeltajPOSTPRE':Deltaj,'maxDeltajPOSTPRE':maxDeltaj,'PrunnedNeurons':PrunnedNeurons})      
    else:
        return (Vest)           
                                        
#%%