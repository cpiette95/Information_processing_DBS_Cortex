%% CODE FOR DISPLAYING FIGURE 6B, S6A/B (after running Network_simulations_firing_rate)

% DBS stimulation
DBS_intensity = [80, 100, 120, 160,200,240];
DBS_pulse_duration = [1, 1.5, 2, 2.5, 3];
DBS_frequency = [18, 67, 130, 180];

for m=3% 1:3
    
    clearvars -except m DBS_intensity DBS_pulse_duration DBS_frequency
    if m==1
        
        % Figure 4B
        conditions = 'DBS';
        parameter = DBS_frequency;
        parameter_label = 'DBS_frequency';
        Labels = {'Control','PD','18Hz','67Hz','130Hz','180Hz'};
        
    elseif m==2
        % Figure S5A
        conditions = 'DBS';
        parameter = DBS_intensity;
        parameter_label = 'DBS_intensity';
        Labels = {'Control','PD','80 pA','100 pA','120 pA', '160 pA','200 pA','240 pA'} ;
        
    elseif m==3
        % Figure S5B
        conditions = 'DBS';
        parameter = DBS_pulse_duration;
        parameter_label = 'DBS_pulse_duration';
        Labels = {'Control','PD','1 ms', '1.5 ms','2 ms','2.5 ms','3 ms'};
    end
    
    firing_rate_control(:,1) = mean(importdata('new_firing_rate_control_1'),2);
    firing_rate_control(:,2) = std(importdata('new_firing_rate_control_1'),0,2);
    firing_rate_noDBS(:,1) = mean(importdata('new_firing_rate_noDBS_1'),2);
    firing_rate_noDBS(:,2) = std(importdata('new_firing_rate_noDBS_1'),0,2);
    rate = importdata(strcat('new_firing_rate_',conditions,'_',parameter_label)) ; 
    firing_rate_PYR(1:length(parameter),1)= mean(rate(1:length(parameter),:),2);
    firing_rate_PYR(1:length(parameter),2)= std(rate(1:length(parameter),:),0,2);
    firing_rate_PV(1:length(parameter),1)= mean(rate(1+length(parameter):2*length(parameter),:),2);
    firing_rate_PV(1:length(parameter),2)= std(rate(1+length(parameter):2*length(parameter),:),0,2);
    firing_rate_SOM(1:length(parameter),1)= mean(rate(1+2*length(parameter):3*length(parameter),:),2);
    firing_rate_SOM(1:length(parameter),2)= std(rate(1+2*length(parameter):3*length(parameter),:),0,2);
    
    figure()
    suptitle(parameter_label); hold on ; 
    %subplot(1,3,1); hold on 
    mean_rate = [firing_rate_control(1,1),firing_rate_noDBS(1,1),firing_rate_PYR(:,1)'];
    std_rate = [firing_rate_control(1,2),firing_rate_noDBS(1,2),firing_rate_PYR(:,2)'];
    bar(1:length(parameter)+2,mean_rate)
    errorbar(1:length(parameter)+2,mean_rate,std_rate,'.')
    set(gca,'XTick',1:length(parameter)+2);
    set(gca,'XTickLabel',Labels)
    ylabel('Pyramidal cells')
    ylim([0 max(mean_rate)*1.2])
    
    subplot(1,3,2);
    hold on
    mean_rate = [firing_rate_control(2,1),firing_rate_noDBS(2,1),firing_rate_PV(:,1)'];
    std_rate = [firing_rate_control(2,2),firing_rate_noDBS(2,2),firing_rate_PV(:,2)'];
    bar(1:length(parameter)+2,mean_rate)
    errorbar(1:length(parameter)+2,mean_rate,std_rate,'.')
    set(gca,'XTick',1:length(parameter)+2);
    set(gca,'XTickLabel',Labels)
    ylabel('PV cells')
    ylim([0 max(mean_rate)*1.2])
    
    subplot(1,3,3);
    hold on
    mean_rate = [firing_rate_control(3,1),firing_rate_noDBS(3,1),firing_rate_SOM(:,1)'];
    std_rate = [firing_rate_control(3,2),firing_rate_noDBS(3,2),firing_rate_SOM(:,2)'];
    bar(1:length(parameter)+2,mean_rate)
    errorbar(1:length(parameter)+2,mean_rate,std_rate,'.')
    set(gca,'XTick',1:length(parameter)+2);
    set(gca,'XTickLabel',Labels)
    ylabel('SST cells')
    ylim([0 max(mean_rate)*1.2])
    
end

clear
% Opto stimulation
Opto_frequency = [13, 67, 130];
Opto_intensity = [200, 400, 600,800];
for m=3:4
    
    clearvars -except m Opto_frequency Opto_intensity
    
    if m==1
        % Figure S6A
        conditions = 'Opto_SOM';
        parameter = Opto_frequency;
        parameter_label = 'Opto_frequency';
        Labels = {'Control','PD','13Hz','67Hz','130Hz'};
    elseif m==2
        % Figure S6B
        conditions = 'Opto_PV';
        parameter = Opto_frequency;
        parameter_label = 'Opto_frequency';
        Labels = {'Control','PD','13Hz','67Hz','130Hz'};
        
    elseif m==3
        % Figure S7A
        conditions = 'Opto_SOM';
        parameter = Opto_intensity;
        parameter_label = 'Opto_intensity';
        Labels = {'Control','PD','200 pA','400 pA','600 pA','800 pA'};
    elseif m==4
        % Figure S7B
        conditions = 'Opto_PV';
        parameter = Opto_intensity;
        parameter_label = 'Opto_intensity';
        Labels = {'Control','PD','200 pA','400 pA','600 pA','800 pA'};
    end
    
    firing_rate_control(:,1) = mean(importdata('new_firing_rate_control_1'),2);
    firing_rate_control(:,2) = std(importdata('new_firing_rate_control_1'),0,2);
    firing_rate_noDBS(:,1) = mean(importdata('new_firing_rate_noDBS_1'),2);
    firing_rate_noDBS(:,2) = std(importdata('new_firing_rate_noDBS_1'),0,2);
    rate = importdata(strcat('new_firing_rate_',conditions,'_',parameter_label));
    firing_rate_PYR(1:length(parameter),1)= mean(rate(1:length(parameter),:),2);
    firing_rate_PYR(1:length(parameter),2)= std(rate(1:length(parameter),:),0,2);
    if m==1 || m== 3
        firing_rate_PV(1:length(parameter),1)= mean(rate(1+length(parameter):2*length(parameter),:),2);
        firing_rate_PV(1:length(parameter),2)= std(rate(1+length(parameter):2*length(parameter),:),0,2);
        firing_rate_SOM(1:length(parameter),1)= mean(rate(1+3*length(parameter):4*length(parameter),:),2);
        firing_rate_SOM(1:length(parameter),2)= std(rate(1+3*length(parameter):4*length(parameter),:),0,2);
    elseif m==2 || m== 4
        firing_rate_PV(1:length(parameter),1)= mean(rate(1+2*length(parameter):3*length(parameter),:),2);
        firing_rate_PV(1:length(parameter),2)= std(rate(1+2*length(parameter):3*length(parameter),:),0,2);
        firing_rate_SOM(1:length(parameter),1)= mean(rate(1+3*length(parameter):4*length(parameter),:),2);
        firing_rate_SOM(1:length(parameter),2)= std(rate(1+3*length(parameter):4*length(parameter),:),0,2);
    end
    
    
    figure()
    suptitle(conditions)
    subplot(1,3,1); hold on
    mean_rate = [firing_rate_control(1,1),firing_rate_noDBS(1,1),firing_rate_PYR(:,1)'];
    std_rate = [firing_rate_control(1,2),firing_rate_noDBS(1,2),firing_rate_PYR(:,2)'];
    bar(1:length(parameter)+2,mean_rate)
    errorbar(1:length(parameter)+2,mean_rate,std_rate,'.')
    set(gca,'XTick',1:length(parameter)+2);
    set(gca,'XTickLabel',Labels)
    ylabel('Pyramidal cells')
    ylim([0 max(mean_rate)*1.2])
    
    subplot(1,3,2);
    hold on
    mean_rate = [firing_rate_control(2,1),firing_rate_noDBS(2,1),firing_rate_PV(:,1)'];
    std_rate = [firing_rate_control(2,2),firing_rate_noDBS(2,2),firing_rate_PV(:,2)'];
    bar(1:length(parameter)+2,mean_rate)
    errorbar(1:length(parameter)+2,mean_rate,std_rate,'.')
    set(gca,'XTick',1:length(parameter)+2);
    set(gca,'XTickLabel',Labels)
    ylabel('PV cells')
    ylim([0 max(mean_rate)*1.2])
    
    subplot(1,3,3);
    hold on
    mean_rate = [firing_rate_control(3,1),firing_rate_noDBS(3,1),firing_rate_SOM(:,1)'];
    std_rate = [firing_rate_control(3,2),firing_rate_noDBS(3,2),firing_rate_SOM(:,2)'];
    bar(1:length(parameter)+2,mean_rate)
    errorbar(1:length(parameter)+2,mean_rate,std_rate,'.')
    set(gca,'XTick',1:length(parameter)+2);
    set(gca,'XTickLabel',Labels)
    ylabel('SST cells')
    ylim([0 max(mean_rate)*1.2])
end


