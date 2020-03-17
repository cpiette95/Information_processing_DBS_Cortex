%% CODE FOR DISPLAYING FIGURE S7b (after running Network_simulations_additional_input, or using Saved_simulation_results)

clear
number_conditions = 11;
label={'control','noDBS','DBS','Opto_PV_200', 'Opto_PV_400','Opto_PV_600','Opto_PV_800','Opto_SOM_200','Opto_SOM_400','Opto_SOM_600','Opto_SOM_800'};

trial_number = 100;
N_PYR  = 800;
dt = 0.05;
stimulus_duration = 500;

% Figure S7b left: type = 1 (ramping stimuli)
% Figure S7b right: type = 2 (OU noise)
type = 2; %[1,2];

Corr_1 = {};
Stats = {}; 
av_corr_1_median = {}; 

for U=1:length(type)
    
    if type(U)==1
        path1 = '/Volumes/LaCie 1/DBS_PAPER/DBS_NEW/Simulations_init_files/Ramping_stim_opto';
        path2 = '/Volumes/LaCie 1/DBS_PAPER/DBS_NEW/Simulations_init_files/Ramping_stim';
        interesting_stimuli_index = 1:8;
        order_stimuli =1:8;
        % Ramping stimuli
        stimulus(1,:) = linspace(0,160,stimulus_duration/dt);
        stimulus(2,:) = linspace(160,0,stimulus_duration/dt);
        stimulus(3,:) = linspace(0,120,stimulus_duration/dt);
        stimulus(4,:) = linspace(120,0,stimulus_duration/dt);
        stimulus(5,:) = linspace(0,80,stimulus_duration/dt);
        stimulus(6,:) = linspace(80,0,stimulus_duration/dt);
        stimulus(7,:) = linspace(0,40,stimulus_duration/dt);
        stimulus(8,:) = linspace(40,0,stimulus_duration/dt);
        
    elseif type(U)==2
        path1 = '/Volumes/LaCie 1/DBS_PAPER/DBS_NEW/Simulations_init_files/OU_stim_opto';
        path2 = '/Volumes/LaCie 1/DBS_PAPER/DBS_NEW/Simulations_init_files/OU_stim';
       cd('/Volumes/LaCie 1/DBS_PAPER/DBS_NEW/Simulations_init_files/');
        interesting_stimuli_index = 1:8;
        order_stimuli = [6,1,4,5,8,7,3,2];% order of the stimulus similar to the display of the manuscript
        
        % OU stimuli
        clear stimulus
        stimulus(1,:) = load('OU_input_1');
        stimulus(2,:) = load('OU_input_2');
        stimulus(3,:) = load('OU_input_3');
        stimulus(4,:) = load('OU_input_4');
        stimulus(5,:) = load('OU_input_5');
        stimulus(6,:) = load('OU_input_6');
        stimulus(7,:) = load('OU_input_7');
        stimulus(8,:) = load('OU_input_8');
        stimulus = stimulus(:,1:end-1);
    end
    
    m=1;
    for w=order_stimuli ; 
        
        stim = stimulus(w,200:end-200);
        cd(path1)
        
        for p=1:number_conditions
            
            if p==1 || p==2 || p==3
                cd(path2)
                M1 = importdata(strcat('Moving_proba_PSTH_',label{p},'_',num2str(w)));
            elseif p==4 || p==5 || p==7 || p==8 || p==9 || p==11
                cd(path1)
                M1 = importdata(strcat('Moving_proba_PSTH_',label{p},'_',num2str(w)));
            elseif p==6
                cd(path2)
                M1 = importdata(strcat('Moving_proba_PSTH_Opto_PV_',num2str(w)));
            elseif p==10
                cd(path2)
                M1 = importdata(strcat('Moving_proba_PSTH_Opto_SOM_',num2str(w)));
            end
            
            for o=1:100
                Proba_moving_PSTH=M1(o,200:end-200); % stimulus
                signal = corrcoef(Proba_moving_PSTH,stim);
                Corr_1{U,p}(m,o) = signal(1,2);
            end
            
            av_corr_1_median{U}(m,p) = nanmedian(Corr_1{U,p}(m,:));
        end
        
        [pval,tbl,stats] = kruskalwallis([Corr_1{U,1}(m,:);Corr_1{U,2}(m,:);Corr_1{U,3}(m,:);Corr_1{U,4}(m,:);Corr_1{U,5}(m,:);Corr_1{U,6}(m,:);Corr_1{U,7}(m,:);Corr_1{U,8}(m,:);Corr_1{U,9}(m,:);Corr_1{U,10}(m,:);Corr_1{U,11}(m,:)]',[],'off');
        c=multcompare(stats);
        Stats{U,m} = c(:,[1,2,6]);
        
        m=m+1;
    end
    
end


% Figure S7 (left)
figure();
imagesc(av_corr_1_median{1}(:,[1,2,3,4,5,6,7,8,9,10,11]))
colormap(hot)
colorbar
caxis([0 0.65])
title('Ramping stimuli')

% Figure S7 (right)
figure();
imagesc(av_corr_1_median{2}(:,[1,2,3,4,5,6,7,8,9,10,11]))
colormap(hot)
caxis([0 0.65])
title('OU stimuli')

