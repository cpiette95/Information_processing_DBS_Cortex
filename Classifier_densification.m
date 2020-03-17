% Code used for reducing the dimensionality of the binned raster plot of
% pyramidal cell activity (densification process)

% Properties of the loaded data:
number_trials=100;
number_cells=800; % number of pyramidal cells 
number_time_bins = 50 ; % (500 ms stimulus duration / 10ms, cf. bin_length = 10 )

% Properties of the output data (how much densification)
number_repetitions=100; 

% Generation of the weight matrix (uniform distribution, positive weights)
% weight_matrix=unifrnd(0,1,number_repetitions,number_cells); 
cd ('/Volumes/LaCie 1/DBS_PAPER/DBS_NEW/Simulations_init_files')
weight_matrix = load('weight_matrix_3'); % saved matrix used over all trials/all stimuli

input_type = 2; 
if input_type ==1
    % Constant input:
    path1 = '/Volumes/LaCie 1/DBS_PAPER/DBS_NEW/Simulations_init_files/Constant_stim';
    conditions = [1,2,3,4]; 
    saving_index = [1,2,3,20];
    
elseif input_type==2
    % Ramping input:
    path1 = '/Volumes/LaCie 1/DBS_PAPER/DBS_NEW/Simulations_init_files/Ramping_stim';
    conditions = 1:8; 
    saving_index = 4:11; 
    
elseif input_type==3
    % OU input:
   path1 = '/Volumes/LaCie 1/DBS_PAPER/DBS_NEW/Simulations_init_files/OU_stim';
    conditions = 1:8;
    saving_index = 12:19;
end

counter = 0; 
for k=conditions    
    counter = counter +1; 
    cd(path1)
    
    % Remember, for each condition, the first 200 rows correspond to the
    % cells being activated by the stimulus (cf. imagesc(reshape(M1(1,:),[number_cells,number_time_bins]));

    M1 = importdata(strcat('V_PYR_binned_PSTH_control_',num2str(k)));
    M2 = importdata(strcat('V_PYR_binned_PSTH_DBS_',num2str(k)));
    M3 = importdata(strcat('V_PYR_binned_PSTH_noDBS_',num2str(k)));   
    M4 = importdata(strcat('V_PYR_binned_PSTH_Opto_PV_',num2str(k)));
    M5 = importdata(strcat('V_PYR_binned_PSTH_Opto_SOM_',num2str(k)));
   
    M1_bis = zeros(number_trials*number_cells,number_time_bins); 
    M2_bis = zeros(number_trials*number_cells,number_time_bins); 
    M3_bis = zeros(number_trials*number_cells,number_time_bins); 
    M4_bis = zeros(number_trials*number_cells,number_time_bins); 
    M5_bis = zeros(number_trials*number_cells,number_time_bins); 
    
    for w=1:number_trials
        M1_small=reshape(M1(w,:),[number_cells,number_time_bins]); 
        M1_bis(1+(w-1)*number_cells:number_cells+(w-1)*number_cells,:) = M1_small;
        
        M1_small=reshape(M2(w,:),[number_cells,number_time_bins]);
        M2_bis(1+(w-1)*number_cells:number_cells+(w-1)*number_cells,:) = M1_small;
        
        M1_small=reshape(M3(w,:),[number_cells,number_time_bins]);
        M3_bis(1+(w-1)*number_cells:number_cells+(w-1)*number_cells,:) = M1_small;
        
        M1_small=reshape(M4(w,:),[number_cells,number_time_bins]);
        M4_bis(1+(w-1)*number_cells:number_cells+(w-1)*number_cells,:) = M1_small;
        
        M1_small=reshape(M5(w,:),[number_cells,number_time_bins]);
        M5_bis(1+(w-1)*number_cells:number_cells+(w-1)*number_cells,:) = M1_small;
        
    end
     
    % Ceate random summation of inputs (same weight matrix over all trials
    % and all stimuli)
    Sampled_new_vector_control = zeros(number_trials,number_repetitions,number_time_bins); 
    Sampled_new_vector_DBS = zeros(number_trials,number_repetitions,number_time_bins); 
    Sampled_new_vector_noDBS = zeros(number_trials,number_repetitions,number_time_bins); 
    Sampled_new_vector_opto_PV = zeros(number_trials,number_repetitions,number_time_bins); 
    Sampled_new_vector_opto_SOM = zeros(number_trials,number_repetitions,number_time_bins); 
       
    for w=1:number_trials
            trial_matrix = M1_bis(1+(w-1)*number_cells:number_cells+(w-1)*number_cells,:); 
            weighted_trial_matrix = weight_matrix*trial_matrix;   
            Sampled_new_vector_control(w,:,:) = weighted_trial_matrix ;
            
            trial_matrix = M2_bis(1+(w-1)*number_cells:number_cells+(w-1)*number_cells,:); 
            weighted_trial_matrix = weight_matrix*trial_matrix; 
            Sampled_new_vector_DBS(w,:,:)= weighted_trial_matrix; 
            
            trial_matrix = M3_bis(1+(w-1)*number_cells:number_cells+(w-1)*number_cells,:);
            weighted_trial_matrix = weight_matrix*trial_matrix ; 
            Sampled_new_vector_noDBS(w,:,:)=weighted_trial_matrix; 
            
            trial_matrix = M4_bis(1+(w-1)*number_cells:number_cells+(w-1)*number_cells,:); 
            weighted_trial_matrix = weight_matrix*trial_matrix; 
            Sampled_new_vector_opto_PV(w,:,:)=weighted_trial_matrix; 
            
            trial_matrix = M5_bis(1+(w-1)*number_cells:number_cells+(w-1)*number_cells,:); 
            weighted_trial_matrix = weight_matrix*trial_matrix; 
            Sampled_new_vector_opto_SOM(w,:,:)= weighted_trial_matrix; 
    end

    cd '/Volumes/LaCie 1/DBS_PAPER/DBS_NEW/Classifier_matrices_basic'
    dlmwrite(strcat('Sampled_new_vector_control_',num2str(saving_index(counter))),Sampled_new_vector_control);
    dlmwrite(strcat('Sampled_new_vector_DBS_',num2str(saving_index(counter))),Sampled_new_vector_DBS);
    dlmwrite(strcat('Sampled_new_vector_noDBS_',num2str(saving_index(counter))),Sampled_new_vector_noDBS);
    dlmwrite(strcat('Sampled_new_vector_opto_PV_',num2str(saving_index(counter))),Sampled_new_vector_opto_PV);
    dlmwrite(strcat('Sampled_new_vector_opto_SOM_',num2str(saving_index(counter))),Sampled_new_vector_opto_SOM);

end
