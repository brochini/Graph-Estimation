# Graph-Estimation

#### Oct 2017
#### by L. Brochini for Neuromat

These programs are part of the Supplementary material, allowing reproduction of figures in the paper:


 #### "Estimation of neuronal interaction graph from spike train data"
 #### Brochini L., Galves A., Hodara P., Ost G., Pouzat C. 


Currently on arxiv
https://arxiv.org/abs/1612.05226

Requirements:
Python 3 with the following modules: matplotlib, numpy, pylab, time and itertools

To generate figures in order of appearance in the manuscript, run: 
* SimulationsXiEps.py
* Pruning.py
* ProjectedConnections.py
* ExperimentalGE.py

The following modules are called within these programs: GEallsubsets.py, genSimData.py, GE.py, GEplots.py

All these .py files should be in a folder containing a subfolder named LocustData that must contain data files:

*locust20010217_spont_tetD_uX.txt where X=[1,2,3,4,7]*

which can be found in LocustData.tar.gz in this repository and contains the pre-processed data corresponding to the spike times for this set of 5 neurons during a trial of a experiment. The raw experimental data corresponds to multielectrode recordings of the spontaneous activity of neurons in the antennal lobe of the locust *Schistocerca americana* and can be found in zenodo
https://zenodo.org/record/21589#%20.WQoUCx1Jlz9

Christophe Pouzat provides a thorough description of the entire preprocessing and the implementation of the spike sorting procedure for this dataset:
https://christophe-pouzat.github.io/zenodo-locust-datasets-analysis/Locust_Analysis_with_R/locust20010217/Sorting_20010217_tetD.html

Acknoledgements: This content was produced as part of the activities of FAPESP  Research, Innovation and Dissemination Center for Neuromathematics (grant no 2013/07699-0 , S.Paulo Research Foundation). I also acknowledge Cnpq support (proc no 165828/2015-3) and FAPESP support (proc no 2016/24676-1).
