#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 22:21:53 2018

@author: charlottepiette
"""

path0 = '/Volumes/LaCie 1/DBS_PAPER/DBS_NEW/Classifier_matrices' 

import os
import numpy as np
import random 

from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors.nearest_centroid import NearestCentroid
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn import svm


number_seeds = 5; 
number_splits = 5; 
max_i = 300; 

score_control_LR = np.zeros((number_seeds,number_splits))
score_noDBS_LR = np.zeros((number_seeds,number_splits))
score_DBS_LR = np.zeros((number_seeds,number_splits))
score_opto_SOM_LR = np.zeros((number_seeds,number_splits))
score_opto_PV_LR = np.zeros((number_seeds,number_splits))

score_control_LDA = np.zeros((number_seeds,number_splits))
score_noDBS_LDA = np.zeros((number_seeds,number_splits))
score_DBS_LDA = np.zeros((number_seeds,number_splits))
score_opto_SOM_LDA = np.zeros((number_seeds,number_splits))
score_opto_PV_LDA = np.zeros((number_seeds,number_splits))

score_control_NN = np.zeros((number_seeds,number_splits))
score_noDBS_NN = np.zeros((number_seeds,number_splits))
score_DBS_NN = np.zeros((number_seeds,number_splits))
score_opto_PV_NN = np.zeros((number_seeds,number_splits))
score_opto_SOM_NN = np.zeros((number_seeds,number_splits))

score_control_SVM = np.zeros((number_seeds,number_splits))
score_noDBS_SVM = np.zeros((number_seeds,number_splits))
score_DBS_SVM = np.zeros((number_seeds,number_splits))
score_opto_PV_SVM = np.zeros((number_seeds,number_splits))
score_opto_SOM_SVM = np.zeros((number_seeds,number_splits))

# Properties of the data to load
N_PYR = 800;
trial_number = 100;
stimulus_duration = 500;
number_repetitions = 100;  # after densification
Bin_Length = 50; 

# The numbers indicated here correspond to the "Sample_new_vector number" 
conditions = [1,2,20,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19] 

print(conditions)
data_control = []
data_DBS = []
data_noDBS = []
data_opto_PV = []
data_opto_SOM = []
target = []

m=0

for k in conditions : 
    os.chdir(path0)
    M1 = open('Sampled_new_vector_control_'+str(k))
    Spike_count_control = [ map(float,line.split(',')) for line in M1 ]        
    M1 = open('Sampled_new_vector_DBS_'+str(k)); 
    Spike_count_DBS = [ map(float,line.split(',')) for line in M1 ]      
    M1 = open('Sampled_new_vector_noDBS_'+str(k)); 
    Spike_count_noDBS = [ map(float,line.split(',')) for line in M1 ]
    data_control.append(Spike_count_control)
    data_DBS.append(Spike_count_DBS)
    data_noDBS.append(Spike_count_noDBS)
    target.append(m*np.ones(trial_number))
    M1 = open('Sampled_new_vector_opto_SOM_'+str(k)); 
    Spike_count_opto_SOM = [ map(float,line.split(',')) for line in M1 ]
    M1 = open('Sampled_new_vector_opto_PV_'+str(k)); 
    Spike_count_opto_PV = [ map(float,line.split(',')) for line in M1 ]
    data_opto_SOM.append(Spike_count_opto_SOM)
    data_opto_PV.append(Spike_count_opto_PV)

    m=m+1
    
data_control = np.vstack(data_control)
data_DBS = np.vstack(data_DBS)
data_noDBS = np.vstack(data_noDBS)
data_opto_SOM = np.vstack(data_opto_SOM)
data_opto_PV = np.vstack(data_opto_PV)
target = np.hstack(target)

for n in range(number_seeds): 

    seed = random.randint(1,100)
    print(seed)
    np.random.seed(seed)
    
    kfold = StratifiedKFold(n_splits=number_splits, shuffle=True, random_state=seed)
    
    f=0; 

    for train_index, test_index in kfold.split(data_opto_SOM,target): 
        
        ## CONTROL ##
        x_train, x_test = data_control[train_index,:],data_control[test_index,:]
        y_train,y_test = target[train_index],target[test_index]      
        mul_lr = LogisticRegression(multi_class='multinomial', solver='newton-cg',max_iter=300)      
        mul_lr.fit(x_train, y_train)  
        score_control_LR[n,f]=mul_lr.score(x_test, y_test)*100
        print(mul_lr.score(x_test,y_test))        
               
        clf = NearestCentroid(metric='euclidean',shrink_threshold=None)  
        clf.fit(x_train,y_train)
        score_control_NN[n,f]=clf.score(x_test,y_test)*100
    
        lda = LinearDiscriminantAnalysis(solver='svd')
        lda.fit(x_train,y_train)
        score_control_LDA[n,f]=lda.score(x_test,y_test)*100
        print(lda.score(x_test,y_test))
       
        svm_algo = svm.SVC(decision_function_shape='ovo',kernel='linear')
        svm_algo.fit(x_train,y_train)
        score_control_SVM[n,f]=svm_algo.score(x_test,y_test)*100

        ## DBS ##
        x_train, x_test = data_DBS[train_index,:],data_DBS[test_index,:]
        y_train,y_test = target[train_index],target[test_index]      

        mul_lr = LogisticRegression(multi_class='multinomial', solver='newton-cg',max_iter=max_i)
        x_train, x_test = data_DBS[train_index,:],data_DBS[test_index,:]
        y_train,y_test = target[train_index],target[test_index]      
        mul_lr.fit(x_train, y_train)
        score_DBS_LR[n,f] = mul_lr.score(x_test, y_test)*100
        print(mul_lr.score(x_test,y_test))        

        clf = NearestCentroid(metric='euclidean',shrink_threshold=None)  
        clf.fit(x_train,y_train) 
        score_DBS_NN[n,f] = clf.score(x_test,y_test)*100

        lda = LinearDiscriminantAnalysis(solver='svd')
        lda.fit(x_train,y_train)
        score_DBS_LDA[n,f]=lda.score(x_test,y_test)*100
        print(lda.score(x_test,y_test))

        svm_algo = svm.SVC(decision_function_shape='ovo',kernel='linear')
        svm_algo.fit(x_train,y_train)
        score_DBS_SVM[n,f]=svm_algo.score(x_test,y_test)*100

        ## NO DBS ##
        x_train, x_test = data_noDBS[train_index,:],data_noDBS[test_index,:]
        y_train,y_test = target[train_index],target[test_index]      
        mul_lr = LogisticRegression(multi_class='multinomial', solver='newton-cg',max_iter=max_i)
        mul_lr.fit(x_train, y_train)
        score_noDBS_LR[n,f] = mul_lr.score(x_test, y_test)*100
        print(mul_lr.score(x_test,y_test))        

        clf = NearestCentroid(metric='euclidean',shrink_threshold=None)  
        clf.fit(x_train,y_train) 
        score_noDBS_NN[n,f] = clf.score(x_test,y_test)*100

        lda = LinearDiscriminantAnalysis(solver='svd')
        lda.fit(x_train,y_train)
        score_noDBS_LDA[n,f]=lda.score(x_test,y_test)*100
        print(lda.score(x_test,y_test))

        svm_algo = svm.SVC(decision_function_shape='ovo',kernel='linear')
        svm_algo.fit(x_train,y_train)
        score_noDBS_SVM[n,f]=svm_algo.score(x_test,y_test)*100

        ## Opto SOM ##
        x_train, x_test = data_opto_SOM[train_index,:],data_opto_SOM[test_index,:]
        y_train,y_test = target[train_index],target[test_index]      
        mul_lr = LogisticRegression(multi_class='multinomial', solver='newton-cg',max_iter=max_i)
        mul_lr.fit(x_train, y_train)
        score_opto_SOM_LR[n,f] = mul_lr.score(x_test, y_test)*100
        print(mul_lr.score(x_test,y_test))        

        clf = NearestCentroid(metric='euclidean',shrink_threshold=None)  
        clf.fit(x_train,y_train)
        score_opto_SOM_NN[n,f] = clf.score(x_test,y_test)*100
     
        lda = LinearDiscriminantAnalysis(solver='svd')
        lda.fit(x_train,y_train)
        score_opto_SOM_LDA[n,f]=lda.score(x_test,y_test)*100
        print(lda.score(x_test,y_test))

        svm_algo = svm.SVC(decision_function_shape='ovo',kernel='linear')
        svm_algo.fit(x_train,y_train)
        score_opto_SOM_SVM[n,f]=svm_algo.score(x_test,y_test)*100
 
        ## Opto PV ##
        x_train, x_test = data_opto_PV[train_index,:],data_opto_PV[test_index,:]
        y_train,y_test = target[train_index],target[test_index]      
        mul_lr = LogisticRegression(multi_class='multinomial', solver='newton-cg',max_iter=max_i)
        mul_lr.fit(x_train, y_train)
        score_opto_PV_LR[n,f] = mul_lr.score(x_test, y_test)*100
        print(mul_lr.score(x_test,y_test))        
        
        clf = NearestCentroid(metric='euclidean',shrink_threshold=None)  
        clf.fit(x_train,y_train)
        score_opto_PV_NN[n,f] = clf.score(x_test,y_test)*100
  
        lda = LinearDiscriminantAnalysis(solver='svd')
        lda.fit(x_train,y_train)
        score_opto_PV_LDA[n,f]=lda.score(x_test,y_test)*100
        print(lda.score(x_test,y_test))

        svm_algo = svm.SVC(decision_function_shape='ovo',kernel='linear')
        svm_algo.fit(x_train,y_train)
        score_opto_PV_SVM[n,f]=svm_algo.score(x_test,y_test)*100
 
        f=f+1; 
                

## Save data 
np.savetxt('score_control_LR',score_control_LR,delimiter=',')
np.savetxt('score_DBS_LR',score_DBS_LR,delimiter=',')
np.savetxt('score_noDBS_LR',score_noDBS_LR,delimiter=',')
np.savetxt('score_opto_SOM_LR',score_opto_SOM_LR,delimiter=',')
np.savetxt('score_opto_PV_LR',score_opto_PV_LR,delimiter=',')
np.savetxt('score_control_NN',score_control_NN,delimiter=',')
np.savetxt('score_DBS_NN',score_DBS_NN,delimiter=',')
np.savetxt('score_noDBS_NN',score_noDBS_NN,delimiter=',')
np.savetxt('score_opto_SOM_NN',score_opto_SOM_NN,delimiter=',')
np.savetxt('score_opto_PV_NN',score_opto_PV_NN,delimiter=',')
np.savetxt('score_control_LDA',score_control_LDA,delimiter=',')
np.savetxt('score_DBS_LDA',score_DBS_LDA,delimiter=',')
np.savetxt('score_noDBS_LDA',score_noDBS_LDA,delimiter=',')
np.savetxt('score_opto_SOM_LDA',score_opto_SOM_LDA,delimiter=',')
np.savetxt('score_opto_PV_LDA',score_opto_PV_LDA,delimiter=',')
np.savetxt('score_control_SVM',score_control_SVM,delimiter=',')
np.savetxt('score_DBS_SVM',score_DBS_SVM,delimiter=',')
np.savetxt('score_noDBS_SVM',score_noDBS_SVM,delimiter=',')
np.savetxt('score_opto_SOM_SVM',score_opto_SOM_SVM,delimiter=',')
np.savetxt('score_opto_PV_SVM',score_opto_PV_SVM,delimiter=',')
