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

intensity={'Opto_PV_200','Opto_SOM_200','Opto_PV_400','Opto_SOM_400','Opto_PV_800','Opto_SOM_800'};
input_type = 1; 
if input_type ==1
    % Constant input:
    path1 =  '/Volumes/LaCie 1/DBS_PAPER/DBS_NEW/Simulations_init_files/Constant_stim_opto';
    conditions = [1,2,3,4]; 
    saving_index = [1,2,3,20];
    
elseif input_type==2
    % Ramping input:
    path1 = '/Volumes/LaCie 1/DBS_PAPER/DBS_NEW/Simulations_init_files/Ramping_stim_opto';
    conditions = 1:8; 
    saving_index = 4:11; 
    
elseif input_type==3
    % OU input:
    path1 = '/Volumes/LaCie 1/DBS_PAPER/DBS_NEW/Simulations_init_files/OU_stim_opto'; 
    conditions = 1:8;
    saving_index = 12:19;
end


for Z=1:length(intensity)
counter = 0; 
for k=conditions    
    counter = counter +1; 
    cd(path1)
    
    % Remember, for each condition, the first 200 rows correspond to the
    % cells being activated by the stimulus (cf. imagesc(reshape(M1(1,:),[number_cells,number_time_bins]));
    M4 = importdata(strcat('V_PYR_binned_PSTH_',intensity{Z},'_',num2str(k)));
    M4_bis = zeros(number_trials*number_cells,number_time_bins); 
    
    for w=1:number_trials

        M1_small=reshape(M4(w,:),[number_cells,number_time_bins]);
        M4_bis(1+(w-1)*number_cells:number_cells+(w-1)*number_cells,:) = M1_small;        
    end
     
    % Ceate random summation of inputs (same weight matrix over all trials
    % and all stimuli)
    Sampled_new_vector_opto = zeros(number_trials,number_repetitions,number_time_bins); 
       
    for w=1:number_trials
            trial_matrix = M4_bis(1+(w-1)*number_cells:number_cells+(w-1)*number_cells,:); 
            weighted_trial_matrix = weight_matrix*trial_matrix; 
            Sampled_new_vector_opto(w,:,:)=weighted_trial_matrix; 
    end

    cd ('/Volumes/LaCie 1/DBS_PAPER/DBS_NEW/Classifier_matrices_opto')
    dlmwrite(strcat('Sampled_new_vector_',intensity{Z},'_',num2str(saving_index(counter))),Sampled_new_vector_opto);

end
end
