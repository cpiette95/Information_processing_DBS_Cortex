%% CODE FOR DISPLAYING FIGURE S5.B (after running Network_simulations_firing_rate_linear_scaling)

frequency_range = [18,67,130,200];
total_pulse_duration = [1, 2, 3];
total_DBS_intensity = [40,80, 120, 160,200,240];
trial_number = 20;

%cd('/Saved_simulation_results/Firing_rates/')
no_DBS = load('new_firing_rate_noDBS_1'); 
firing_rate_Pyr_mean_std(1,1) = 0; 
firing_rate_Pyr_mean_std(1,2) = mean(no_DBS(1,:)); 
firing_rate_Pyr_mean_std(1,3) = std(no_DBS(1,:)); 

%cd('/Saved_simulation_results/Firing_rates/Linear_scaling')
counter = 2;
for w=1:length(total_DBS_intensity)
    for m=1:length(total_pulse_duration)
        av_firing_rate_PYR_DBS = load(strcat('new_firing_rate_DBS_amp=',num2str(total_DBS_intensity(w)),'pA_pulse_dur=',num2str(total_pulse_duration(m)),'ms_freq_range=all'));      
        for k=1:length(frequency_range)
            firing_rate_Pyr_mean_std(counter,1) = total_DBS_intensity(w)*total_pulse_duration(m)*1e-3*frequency_range(k); 
            firing_rate_Pyr_mean_std(counter,2) = mean(av_firing_rate_PYR_DBS(k,:)); 
            firing_rate_Pyr_mean_std(counter,3) = std(av_firing_rate_PYR_DBS(k,:)); 
            counter = counter +1; 
        end
    end
end

firing_rate_Pyr_mean_std = sortrows(firing_rate_Pyr_mean_std,1); 

figure(); 
xlabel('Total amount of current injected during 1 second')
ylabel('Pyramidal cell firing rate (Hz)')
P = polyfit(firing_rate_Pyr_mean_std(1:end-12,1),firing_rate_Pyr_mean_std(1:end-12,2),1);
yfit = P(1)*firing_rate_Pyr_mean_std(1:end-12,1)+P(2);
hold on;
plot(firing_rate_Pyr_mean_std(1:end-12,1),yfit,'r-.','LineWidth',3);
errorbar(firing_rate_Pyr_mean_std(:,1),firing_rate_Pyr_mean_std(:,2),firing_rate_Pyr_mean_std(:,3),'ko')
errorbar(firing_rate_Pyr_mean_std(52,1),firing_rate_Pyr_mean_std(52,2),firing_rate_Pyr_mean_std(52,3),'bo')
xlim([0 130])
fitlm(firing_rate_Pyr_mean_std(1:end-12,1),firing_rate_Pyr_mean_std(1:end-12,2))
