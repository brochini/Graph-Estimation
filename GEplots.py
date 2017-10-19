'''
Python 3

Created by L. Brochini for Neuromat on Feb 2017

This program is part of the Supplementary material to the paper:

 "Interaction graph estimation for the first olfactory relay of an insect"
L. Brochini, P.Hodara, C.Pouzat, A.Galves 
__________________________

Module of plotting functions used to represent connectivity matrices
obs: hatch and transparecny are not saved properly in eps and many other formats.
Formats that display correctly are pdf and svgz

Default is set to output black and white figures. To obtain color figures, pass default argument colors='yes' 
'''



import matplotlib.pyplot as plt
import pylab
from matplotlib import rc
import matplotlib.ticker as ticker

import numpy as np
mymat=np.array([[0,1,2,3],[3,0,1,2],[2,3,0,1],[1,2,3,0]])


#defining colormaps
def map_projections():
    
    from matplotlib.colors import LinearSegmentedColormap
    vmin=0
    vmax = 4.0
    cmap = LinearSegmentedColormap.from_list('mycmap', [(0 / vmax, 'white'), # False
                                                        (1 / vmax, 'gray'), # True
                                                        (2 / vmax, 'green'), # First order
                                                        (3 / vmax, [0.66,0,1]), # Second Order (Violet)
                                                        (4 / vmax, 'cyan'),   # Inconclusives due to only Nan values
                                                        ] # Inconclusives due to higher order projections or real False Positive/FalseNegative
                                            )
    
    return(cmap,vmin,vmax) 


def map_w_inconclusives():
    
    from matplotlib.colors import LinearSegmentedColormap
    vmin=0
    vmax = 5.0
    cmap = LinearSegmentedColormap.from_list('mycmap', [(0 / vmax, 'white'),
                                                        (1 / vmax, 'gray'),
                                                        (2 / vmax, 'blue'),
                                                        (3 / vmax, 'red'),
                                                        (4 / vmax, [0.6784313725490196, 0.8470588235294118, 0.9019607843137255]),
                                                        (5/vmax,   [ 0.98823529,  0.70588235,  0.63529412])]# light red
                                            )
    
    return(cmap,vmin,vmax) 


    
  
def maphatch():
    from matplotlib.colors import LinearSegmentedColormap
    vmin=0
    vmax=4
    cmap = LinearSegmentedColormap.from_list('mycmap', [(0 / (vmax-vmin), 'white'),
                                                        (1 /(vmax-vmin), 'black'),
                                                        (2 / (vmax-vmin), (0,0,0,0)),
                                                        (3 / (vmax-vmin), (0,0,0,0.3)), #total black with lots of transparency leads to light grey
                                                        (4 / (vmax-vmin),(0.7,0.7,0.7))]) 
    return(cmap,vmin,vmax) 


def maphatch2(): # one extra level of grey. Used to differentiate between inconclusives corresponding to true or absent connections
    from matplotlib.colors import LinearSegmentedColormap
    vmin=0
    vmax=5
    cmap = LinearSegmentedColormap.from_list('mycmap', [(0 / (vmax-vmin), 'white'),
                                                        (1 /(vmax-vmin), 'black'),
                                                        (2 / (vmax-vmin), (0,0,0,0)),
                                                        (3 / (vmax-vmin), (0,0,0,0.3)), #total black with lots of transparency leads to light grey
                                                        (4 / (vmax-vmin),(0.8,0.8,0.8)),
                                                        (5 / (vmax-vmin),(0.3,0.3,0.3))]) 
    return(cmap,vmin,vmax) 




def maphatchBase():
    from matplotlib.colors import LinearSegmentedColormap
    vmin=0
    vmax=4
    cmap = LinearSegmentedColormap.from_list('mycmap', [(0 / (vmax-vmin), (1,1,1)),
                                                        (1 /(vmax-vmin), (1,1,1)),
                                                        (2 / (vmax-vmin), (1,1,1)),
                                                        (3 / (vmax-vmin), (1,1,1)),
                                                        (4 / (vmax-vmin), (1,1,1))])
    return(cmap,vmin,vmax) 

def maphatchBase2():
    from matplotlib.colors import LinearSegmentedColormap
    vmin=0
    vmax=5
    cmap = LinearSegmentedColormap.from_list('mycmap', [(0 / (vmax-vmin), (1,1,1)),
                                                        (1 /(vmax-vmin), (1,1,1)),
                                                        (2 / (vmax-vmin), (1,1,1)),
                                                        (3 / (vmax-vmin), (1,1,1)),
                                                        (4 / (vmax-vmin), (1,1,1)),
                                                        (5 / (vmax-vmin), (1,1,1))])
    return(cmap,vmin,vmax) 


def PlotOneDiffMapsWInconclusives_ax(ax,Wmatrix,Vest,filename,titletext,showax=False,colors='no'):
    
    rc('text', usetex=True)
    rc('font', family='serif')
    
    N=len(Wmatrix)
    cmap1,vmin1,vmax1=maphatchBase2()
    cmap2,vmin2,vmax2=maphatch2()
    cmapc,vminc,vmaxc=map_w_inconclusives()
    
    Vtrue=np.zeros([N,N])
    for i in range(N):
        for j in range(N):
            if Wmatrix[i,j]!=0:
                Vtrue[i,j]=1
                
    
    diffmat=np.zeros([N,N])
    for i in range(N):
        for j in range(N):
            if Vest[i,j]==0:
                if Vtrue[i,j]==0:
                    diffmat[i,j]=0
                if Vtrue[i,j]==1:
                    diffmat[i,j]=2
            if Vest[i,j]==1:
                if Vtrue[i,j]==0:
                    diffmat[i,j]=3
                if Vtrue[i,j]==1:
                    diffmat[i,j]=1
            if Vest[i,j]==-2:
                if Vtrue[i,j]==0:
                    diffmat[i,j]=4
                if Vtrue[i,j]==1:
                    diffmat[i,j]=5
    
    if colors=='no':
        ax.pcolor(diffmat, cmap=cmap1, vmin=vmin1, vmax=vmax1, edgecolors='black' ,hatch='///')
        ax.pcolor(diffmat, cmap=cmap2, vmin=vmin2, vmax=vmax2, edgecolors='black')
    else:        
        ax.pcolor(diffmat, cmap=cmapc, vmin=0, vmax=vmaxc, edgecolors='black' )
        
    
    if showax:    
        ax.set_title(titletext,y=1.12)
        ax.set_ylabel("Presynaptic")
        ax.set_xlabel("Postsynaptic")      
        ax.set(frame_on=False, aspect=1, xticks=range(1,N+1), yticks=range(1,N+1))
    else:
        ax.set_title(titletext)
        
    ax.set(frame_on=False, aspect=1,xticks=[],yticks=[])
    ax.invert_yaxis() 
    if showax:
        ax.xaxis.set_label_position('top')
        ax.xaxis.tick_top()
        tickpos=list(np.array(range(N))+0.5)
        ticknames=[str(i) for i in range(1,N+1)]
        ax.tick_params(axis=u'both', which=u'both',length=0)    
        ax.xaxis.set_major_formatter(ticker.NullFormatter())
        ax.xaxis.set_minor_locator(ticker.FixedLocator(tickpos))
        ax.xaxis.set_minor_formatter(ticker.FixedFormatter(ticknames))
    
        ax.tick_params(axis=u'both', which=u'both',length=0)    
        ax.yaxis.set_major_formatter(ticker.NullFormatter())
        ax.yaxis.set_minor_locator(ticker.FixedLocator(tickpos))
        ax.yaxis.set_minor_formatter(ticker.FixedFormatter(ticknames))



def PlotMultDifs(Wmatrix,Vdic,epsilon_values,xi_values,filename,titletext,colors='no'):

    yplots=[i*len(xi_values)+1 for i in range(len(epsilon_values))] # subplot index that receive lateral title with epsilon value
    nkeys=len(Vdic.keys())    
    f, axarr = plt.subplots(len(epsilon_values),len(xi_values),figsize=(12,12))
    countsubplot=0
    for ind_epsilon in range(len(epsilon_values)):
        
        epsilon=epsilon_values[ind_epsilon]
        for ind_xi in range(len(xi_values)):
            countsubplot+=1
            xi=xi_values[ind_xi]
            PlotOneDiffMapsWInconclusives_ax( axarr[ind_epsilon,ind_xi],Wmatrix,Vdic[(epsilon,xi)],filename,'',showax=False,colors=colors)
            if countsubplot in range(1,int(nkeys/len(epsilon_values)+1)):
                uptitle=r"$\xi=$"+str(xi)
                axarr[ind_epsilon,ind_xi].set_title(uptitle)
            if countsubplot in yplots:
                ytitle=r"$\epsilon=$"+str(epsilon)
                axarr[ind_epsilon,ind_xi].set_ylabel(ytitle)              
        
    f.suptitle(titletext, fontsize=14)        
    pylab.savefig(filename+'.pdf', format='pdf', dpi=1000)        



def PlotLineDiffMaps(Wmatrix,Vlist,Titlearr,filename,showax=False,experimental=False,colors='no'):
    cols=len(Vlist)
    f, axarr = plt.subplots(1,cols,figsize=(5*cols,6))
    for i in range(cols):
        if experimental:
            PlotOneDiffMapsWInconclusives_ax( axarr[i],Vlist[i],Vlist[i],filename,Titlearr[i],showax=showax,colors=colors)
        else:
            PlotOneDiffMapsWInconclusives_ax( axarr[i],Wmatrix,Vlist[i],filename,Titlearr[i],showax=showax,colors=colors)
        
    pylab.savefig(filename+'.pdf', format='pdf', dpi=1000)        

    

def PlotOneMapProjections(DirectConn,ProjectConn,figname,titletext,colors='no'):
    
    rc('text', usetex=True)
    rc('font', family='serif')
    
    N=len(DirectConn)
    cmap1,vmin1,vmax1=maphatchBase()
    cmap2,vmin2,vmax2=maphatch()
    cmapc,vminc,vmaxc= map_projections()   
    
    MapProj=np.zeros([N,N])    
    for i in range(N):
        for j in range(N):
            if i !=j:
                if DirectConn[i][j]==True:
                    MapProj[i,j]=1
                else:
                    if ProjectConn[i][j]==1:
                        MapProj[i,j]=2
                    elif ProjectConn[i][j]==2:
                        MapProj[i,j]=3
                    elif DirectConn[i][j]==False and ProjectConn[i][j]==0:
                        MapProj[i,j]=0                    
                    else:
                        MapProj[i,j]=4

    fig=plt.figure()
    ax=fig.add_subplot(111)
    
    
    if colors=='no':
        ax.pcolor(MapProj, cmap=cmap1, vmin=0, vmax=vmax1, edgecolors='black',hatch='///')
        ax.pcolor(MapProj, cmap=cmap2, vmin=0, vmax=vmax2, edgecolors='black')
    else:
        ax.pcolor(MapProj, cmap=cmapc, vmin=0, vmax=vmaxc, edgecolors='black')
    
    ax.set_ylabel("Presynaptic")
    ax.set_xlabel("Postsynaptic")      
    ax.set(frame_on=False, aspect=1, xticks=range(1,N+1), yticks=range(1,N+1))
    ax.invert_yaxis() 
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()

    tickpos=list(np.array(range(N))+0.5)
    ticknames=[str(i) for i in range(1,N+1)]
    ax.tick_params(axis=u'both', which=u'both',length=0)    
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    ax.xaxis.set_minor_locator(ticker.FixedLocator(tickpos))
    ax.xaxis.set_minor_formatter(ticker.FixedFormatter(ticknames))

    ax.tick_params(axis=u'both', which=u'both',length=0)    
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
    ax.yaxis.set_minor_locator(ticker.FixedLocator(tickpos))
    ax.yaxis.set_minor_formatter(ticker.FixedFormatter(ticknames))
    ax.set_title(titletext, fontsize=14,y=1.08)        
    
    pylab.savefig(figname+'.pdf', format='pdf', dpi=1000)

        
def PlotISIHist(ISI,nbins,xlim,figname,title,GoodNeurons,acquisitionrate):
    N=len(ISI)
    f, axarr = plt.subplots(N,1,figsize=(10,10))
    for neuron in range(N):
        axarr[neuron].hist(np.array(ISI[neuron])/acquisitionrate,nbins)
        axarr[neuron].set_xlim([0,xlim/acquisitionrate])
        axarr[neuron].set_ylabel('neuron '+str(GoodNeurons[neuron]))
        annot='total spikes='+str(len(ISI[neuron])+1)
        posx=xlim/(2*acquisitionrate)
        posy=axarr[neuron].get_ylim()[1]*0.8
        axarr[neuron].annotate(annot,xy=(posx, posy))    
    axarr[N-1].set_xlabel('ISI (s)')
    f.suptitle(title, fontsize=14)  
    
    pylab.savefig(figname+'.pdf', format='pdf', dpi=1000)        
        

def PlotWmatrix(Wmatrix,figname,HighDiag=False,title=''):
    N=len(Wmatrix)
    myvmax=max(max([max(Wmatrix[i]) for i in range(len(Wmatrix))]),(-1)* min([min(Wmatrix[i]) for i in range(len(Wmatrix))]))
    myvmin=-myvmax
    if HighDiag:
        for neuron in range(N):
            Wmatrix[neuron,neuron]=myvmax
            
    fig=plt.figure();
    ax=fig.add_subplot(111);
    pplot=ax.pcolor(Wmatrix, cmap='jet', vmin=myvmin, vmax=myvmax);
    ax.set(aspect=1)
    
    ax.set_ylabel("Presynaptic")
    ax.set_xlabel("Postsynaptic")      
    ax.set(frame_on=False, aspect=1, xticks=range(1,N+1), yticks=range(1,N+1))
    ax.invert_yaxis() 
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()
    
    tickpos=list(np.array(range(N))+0.5)
    ticknames=[str(i) for i in range(1,N+1)]
    ax.tick_params(axis=u'both', which=u'both',length=0)    
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    ax.xaxis.set_minor_locator(ticker.FixedLocator(tickpos))
    ax.xaxis.set_minor_formatter(ticker.FixedFormatter(ticknames))
    
    ax.tick_params(axis=u'both', which=u'both',length=0)    
    ax.yaxis.set_major_formatter(ticker.NullFormatter())
    ax.yaxis.set_minor_locator(ticker.FixedLocator(tickpos))
    ax.yaxis.set_minor_formatter(ticker.FixedFormatter(ticknames))
    fig.suptitle(title)
    fig.colorbar(pplot)
    pylab.savefig(figname)
    