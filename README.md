# Information_processing_DBS_Cortex

Computational model from S. Valverde, M. Vandecasteele, C. Piette, W. Derousseaux, G. Gangarossa, A. Aristieta Arbelaiz, J. Touboul, B. Degos, L. Venance, "Deep brain stimulation-guided optogenetic rescue of parkinsonian symptoms", Nature Communications (in preparation). 

This simplified spiking network model of L5 motor cortex was developed in order to test its information processing capabilities in different conditions (control, Parkinson, Parkinson + DBS or opto-activation of cortical interneurons).  


1. System requirements
Simulations of network activity patterns were run in Matlab R2016b. More recent versions of Matlab can be used and downloaded from https://www.mathworks.com/downloads/.   

Simulations of classifier accuracy were run in Python 2.7 in a dedicated Anaconda environment. Scikit-learn (version 0.19.1 used in the simulations, easy download from https://scikit-learn.org/stable/install.html) is a necessary package for running the “Classifier_results.py” code. 

Simulations have been tested on OS.X version 10.9.5 and Ubuntu Linux 16.04.  
Typical install time of Matlab and Python on a “normal” desktop computer is about an hour long. 


2. Short demo
A short demo to visualise network activity patterns can be run from “Demo_code.m” in Matlab. The average simulation time on a “normal” desktop computer is about 2 min. 

The idea of this demo is to get familiar with the different network activity patterns obtained in the five conditions studied in the manuscript: control, parkinsonian, parkinsonian with DBS stimulation or opto-activation of PV or SST neurons. The network parameters are the same as those used in the manuscript. 

The precise parameters of DBS or opto-activation stimulation can easily be modified (and are by default those used in the manuscript for information processing measures).  Additional inputs injected to a subset of pyramidal cells can also be chosen by the user. 

By default, Figure 1 displays the stimulus injected. Then subsequent figures consist of raster plots similar to Figure 7.b for all the different conditions, as well as moving spike counts and binned raster plots of pyramidal cells (used in the manuscript for measuring information capacities). 


3. Detailed instructions
The .zip file contains different source codes (.m and .py files), including: 

1/ those simulating network activity patterns in different conditions: “Network_simulations_firing_rate.m” , ”Network_simulations_firing_rate_linear_scaling”, “Network_simulations_firing_rate_robustness_maps” simulate network activity in five different conditions (control, Parkinson, Parkinson+DBS or opto-stimulation), while varying the parameters of the stimulation (DBS frequency, DBS pulse duration and intensity or Opto-stimulation frequency) ;  “Network_simulations_additional_input.m” simulate network activity in the presence of additional inputs (constant, ramping inputs or Ornstein-Uhlenbeck noise) in those five conditions. 

2/ those displaying figures obtained from the simulations (Figure_.m files) and displayed in the manuscript (in order for them to run correctly, the user needs first to save data from the simulating source codes described in 1 and 3)

3/ those used at intermediate steps: “OU_noise_process.m” generates Ornstein-Uhlenbeck noise used as stimuli ; “Classifier_densification.m” is used to densify binned matrices of pyramidal cells responses; “Classifier_results.py” trains classical machine-learning algorithms to decode stimulus identity from densified matrices of pyramidal cells responses.  


Importantly, some of the simulations rely on specific files, found in the folder  “Simulations_init_files”. These concern simulations of network activity patterns in the presence of additional stimuli and are the following: 
- fixed connectivity matrices (“Connect_” files) for comparing results obtained for different stimuli and in different conditions (control/Parkinson/DBS/Opto-activation) 
- fixed choice of cells receiving each type of inputs (“receptor_cells_” files) for comparing results obtained in different conditions (control/Parkinson/DBS/Opto-activation)
- saved stochastic inputs (“OU_input_” files)
- weight_matrix_3 (used for densifying binned matrices of pyramidal cells responses in the “Classifier_densification.m” script)
One could easily decide to modify these files, and one should obtain very similar results, as those presented in the manuscript. 


Some simulations require a long processing time (in particular “Network_simulations_firing_rate_robustness_maps” and “Network_simulations_additional_input”). For the second source code, two matrices are normally saved and used for analyses on information capacities : 
- moving spike counts of all pyramidal cells across time (“Moving_proba_PSTH_” files), used for measuring for each trial the Pearson correlation coefficient. 
- binned raster plots of pyramidal cell activity for each trial (“V_PYR_binned_PSTH” files), which will then be densified (“Sampled_new_vector_” files). The classifiers will then be trained on these densified matrices to decode stimulus identity. 

The list of stimuli used for the manuscript is the following: 
Constant inputs:
1/ Constant input, 40 pA (n°1)

2/ Constant input, 50pA (n°2)

3/ Constant input, 70pA (n°3)

4/ Constant input, 60 pA (n°20)

Ramping stimuli
1/ Ramp de 0 à 160 pA (n°4)

2/ Ramp de 160 à 0 pA (n°5)

3/ Ramp de 0 à 120 (n°6) 

4/ Ramp de 120 à 0 (n°7)

5/ Ramp de 0 à 80 pA (n°8)
6/ Ramp de 80 à 0 pA (n°9)
7/ Ramp de 0 à 40 pA (n°10)
8/ Ramp de 40 à 0 pA (n°11)
OU noise
1/ OU_noise_1 (n°12)
2/ OU_noise_2 (n°13)
3/ OU_noise_3 (n°14)
4/ OU_noise_4 (n°15)
5/ OU_noise_5 (n°16)
6/ OU_noise_6 (n°17)
7/ OU_noise_7 (n°18)
8/ OU_noise_8 (n°19)


