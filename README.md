# Graph-Estimation
Mar 22 2017
by L. Brochini for Neuromat

These programs are part of the Supplementary material, allowing you to reproduce figures in the paper:


 "Interaction graph estimation for the first olfactory relay of an insect"
 L. Brochini, P.Hodara, C.Pouzat, A.Galves 


Currently on arxiv
https://arxiv.org/abs/1612.05226

Requirements:
Python 3 with the following modules: matplotlib, numpy, pylab, time and itertools

To generate figures in order of appearance in the manuscript, run: SimulationsXiEps.py, Pruning.py, ProjectedConnections.py, ExperimentalGE.py

The following modules are called within these programs: GEallsubsets.py, genSimData.py, GE.py, GEplots.py

All these .py files should be in a folder containing a subfolder named
LocustData that must contain data files: locust20010217_spont_tetD_uX.txt where X=[1,2,3,4,7]

Acknoledgements: This content was produced as part of the activities of FAPESP  Research, Innovation and Dissemination Center for Neuromathematics (grant no 2013/07699-0 , S.Paulo Research Foundation). I also acknowledge Cnpq support (proc no 165828/2015-3) and FAPESP support (proc no 2016/24676-1).
