%% CODE FOR DISPLAYING FIGURE 7C or S6 (after running Network_simulations_additional_input)

clear
number_conditions = 5;
label={'Control','DBS off','DBS on','Opto PV','Opto SOM'};

trial_number = 100;
N_PYR  = 800;
dt = 0.05;
stimulus_duration = 500;

% Figure 7c left: type = 1 (ramping stimuli)
% Figure 7c right: type = 2 (OU noise)
type = 1; % [1,2];
more_DBS = 1; % check opto PV/SST


Corr_1 = {};
Stats = {}; 
av_corr_1_median = {}; 

for U=1:length(type)
    
    if type(U)==1
        path1 = '/Volumes/LaCie 1/DBS_PAPER/DBS_NEW/Simulations_init_files/Ramping_stim';
        interesting_stimuli_index = 1:12;
        order_stimuli =1:12; % order of the stimulus similar to the display of the manuscript
        
        % Ramping stimuli
        stimulus(1,:) = linspace(0,160,stimulus_duration/dt);
        stimulus(2,:) = linspace(160,0,stimulus_duration/dt);
        stimulus(3,:) = linspace(0,120,stimulus_duration/dt);
        stimulus(4,:) = linspace(120,0,stimulus_duration/dt);
        stimulus(5,:) = linspace(0,80,stimulus_duration/dt);
        stimulus(6,:) = linspace(80,0,stimulus_duration/dt);
        stimulus(7,:) = linspace(0,40,stimulus_duration/dt);
        stimulus(8,:) = linspace(40,0,stimulus_duration/dt);
        stimulus(9,:) = linspace(0,20,stimulus_duration/dt);
        stimulus(10,:) = linspace(20,0,stimulus_duration/dt);
        stimulus(11,:) = linspace(0,10,stimulus_duration/dt);
        stimulus(12,:) = linspace(10,0,stimulus_duration/dt);

    elseif type(U)==2
        path1 = '/Volumes/LaCie 1/DBS_PAPER/DBS_NEW/Simulations_init_files/OU_stim';
        cd '/Volumes/LaCie 1/DBS_PAPER/DBS_NEW/Simulations_init_files';
        interesting_stimuli_index = 1:12;
        order_stimuli = [6,1,4,5,8,7,3,2,12,9,11,10]; % order of the stimulus similar to the display of the manuscript
        
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
        stimulus(9,:) = load('OU_input_9');
        stimulus(10,:) = load('OU_input_10');
        stimulus(11,:) = load('OU_input_11');
        stimulus(12,:) = load('OU_input_12');
        stimulus = stimulus(:,1:end-1);
    end
    
    m=1;
    for w=order_stimuli ; 
        
        stim = stimulus(w,200:end-200);
        cd(path1)
        
        for p=1:number_conditions
            
            if p==1
                M1 = importdata(strcat('Moving_proba_PSTH_control_',num2str(w)));
            elseif p==2
                M1 = importdata(strcat('Moving_proba_PSTH_noDBS_',num2str(w)));
            elseif p==3
                M1 = importdata(strcat('Moving_proba_PSTH_DBS_',num2str(w)));
            end
            
            if more_DBS==1
                if p==4
                    M1 = importdata(strcat('Moving_proba_PSTH_Opto_PV_',num2str(w)));
                elseif p==5
                    M1 = importdata(strcat('Moving_proba_PSTH_Opto_SOM_',num2str(w)));
                end
                
            end
            
            for o=1:100
                Proba_moving_PSTH=M1(o,200:end-200); % stimulus
                signal = corrcoef(Proba_moving_PSTH,stim);
                Corr_1{U,p}(m,o) = signal(1,2);
            end
            
            av_corr_1_median{U}(m,p) = nanmedian(Corr_1{U,p}(m,:));
        end
        
        [pval,tbl,stats] = kruskalwallis([Corr_1{U,1}(m,:);Corr_1{U,2}(m,:);Corr_1{U,3}(m,:);Corr_1{U,4}(m,:);Corr_1{U,5}(m,:)]',[],'off');
        c=multcompare(stats);
        Stats{U,m} = c(:,[1,2,6]);
        
        m=m+1;
    end
    
end



% Figures 6 - Heatmaps
% 6E (left)
figure();
imagesc(av_corr_1_median{1}(:,[1,2,3,4,5]))
colormap(hot)
colorbar
caxis([0 0.65])
title('Ramping stimuli')
% 6E (right)
figure();
imagesc(av_corr_1_median{1}(:,[1,2,3,4,5]))
colormap(hot)
caxis([0 0.65])
title('OU stimuli')



% Figure S6 - Boxplot
func_tris = @(x) nanmedian(x); 
x=[]; 
group = []; 
counter = 1; 
counter_bis = 1; 
order_ramp = 1:12; 
order_noise = [6,1,4,5,8,7,3,2,12,9,11,10];
for w=order_ramp; 
    for p=1:5
        x = [x ; Corr_1{p}(w,:)] ;
        group = [group ; repmat(counter,trial_number,1)];
        counter = counter +1;
    end
end
positions = [1 1.25 1.5 1.75 2 2.75 3 3.25 3.5 3.75 4.5 4.75 5 5.25 5.5 6.25 6.5 6.75 7 7.25 8 8.25 8.5 8.75 9 9.75 10 10.25 10.5 10.75 11.5 11.75 12 12.25 12.5 13.25 13.5 13.75 14 14.25 15 15.25 15.5 15.75 16 16.75 17 17.25 17.5 17.75 18.5 18.75 19 19.25 19.5 20.25 20.5 20.75 21 21.25] ; % 22 22.25 22.5 22.75 23];  
figure();
xlabel('Stimuli for different conditions')
ylabel('Correlation coefficient (on stim only)')
line([0 22], [0 0])
boxplot(x',group,'PlotStyle','compact','positions',positions,'MedianStyle','target','BoxStyle','filled','Colors','rgbmc','OutlierSize',0); 
ylim([-0.4 0.8])



