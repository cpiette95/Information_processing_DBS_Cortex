%% CODE FOR REPRODUCING FIGURE S5.B

dt = 0.05;
T = 1200;
N_times=ceil(T/dt);
Time_vector = 0:dt:T;
incr=1;

N_PYR_REF = 800;
N_PV_REF = 120;
N_SOM_REF = 80;

N_PYR = ceil(N_PYR_REF*incr);
N_PV = ceil(N_PV_REF*incr);
N_SOM = ceil(N_SOM_REF*incr);
N_total = N_PYR + N_SOM + N_PV;

% DBS inputs
frequency_range = [18,67,130,200];
total_pulse_duration = [1, 2, 3];
total_DBS_intensity = [40,80,120,160,200,240];
trial_number = 20;  

counter = 0 ; 
current_injected = zeros(length(frequency_range)*length(total_DBS_intensity)*length(total_pulse_duration),1); 
firing_rate_Pyr_mean_std = zeros(length(frequency_range)*length(total_DBS_intensity)*length(total_pulse_duration),2); 

av_firing_rate_PYR_DBS = zeros(length(total_pulse_duration),trial_number);
av_firing_rate_PV_DBS = zeros(length(total_pulse_duration),trial_number);
av_firing_rate_SOM_DBS = zeros(length(total_pulse_duration),trial_number);

for w=1:length(total_DBS_intensity)
    
    for m=1:length(total_pulse_duration)
        
        for k=1:length(frequency_range)
             
            disp(['Counter=',num2str(counter)]);
            % DBS parameters (pA)
            I0_PYR = total_DBS_intensity(w);
            I0_PV = total_DBS_intensity(w);
            I0_SOM = total_DBS_intensity(w);
            disp(['Intensity=',num2str(total_DBS_intensity(w))]);
            frequency = 1e-3*frequency_range(k);
            disp(['Frequency=',num2str(frequency_range(k))]);
            pulse_duration = total_pulse_duration(m);
            disp(['Pulse duration=',num2str(total_pulse_duration(m))]);
            delay_PV = 0;
            delay_SOM = 2;
            delay_PYR = 2;

            counter = counter +1; 
            current_injected(counter) = total_DBS_intensity(w)*frequency_range(k)*total_pulse_duration(m); 
            
            stim_on_PYR = ones(round(pulse_duration/dt),400);
            stim_off_PYR = zeros((round(1/(frequency*dt))-size(stim_on_PYR,1)),400);
            stimulus_template_PYR = vertcat(stim_on_PYR,stim_off_PYR);
            stimulus_total_PYR = repmat(I0_PYR*stimulus_template_PYR',1,ceil(ceil(T/dt)/size(stimulus_template_PYR,1)));
            
            stim_off_pre_PYR_bis = zeros(round(delay_PYR/dt),400);
            stim_on_PYR_bis = ones(round(pulse_duration/dt),400);
            stim_off_PYR_bis = zeros(round(1/(frequency*dt) - size(stim_on_PYR_bis,1) -size(stim_off_pre_PYR_bis,1)),400);
            stimulus_template_PYR_bis = vertcat(stim_off_pre_PYR_bis,stim_on_PYR_bis,stim_off_PYR_bis);
            stimulus_total_PYR_bis= repmat(I0_PYR*stimulus_template_PYR_bis',1,ceil(ceil(T/dt)/size(stimulus_template_PYR_bis,1)));
            
            stimulus_total_PYR = [stimulus_total_PYR;stimulus_total_PYR_bis];
            
            stim_on_PV = ones(round(pulse_duration/dt),N_PV);
            stim_off_PV = zeros(round(1/(frequency*dt)-size(stim_on_PYR,1)),N_PV);
            stimulus_template_PV = vertcat(stim_on_PV,stim_off_PV);
            stimulus_total_PV = repmat(I0_PV*stimulus_template_PV',1,ceil(ceil(T/dt)/size(stimulus_template_PV,1)));
            
            stim_off_pre_SOM = zeros(round(delay_SOM/dt),N_SOM);
            stim_on_SOM = ones(round(pulse_duration/dt),N_SOM);
            stim_off_SOM = zeros(round(1/(frequency*dt) - size(stim_on_SOM,1) -size(stim_off_pre_SOM,1)),N_SOM);
            stimulus_template_SOM = vertcat(stim_off_pre_SOM,stim_on_SOM,stim_off_SOM);
            stimulus_total_SOM = repmat(I0_SOM*stimulus_template_SOM',1,ceil(ceil(T/dt)/size(stimulus_template_PV,1)));
            
            
            for p=1:trial_number
                
                EL = -60;
                Ee = 0;
                Ei = -80;
                
                sigma = 5;
                
                DeltaT_PYR = 1;
                DeltaT_PV = 1;
                DeltaT_SOM = 5;
                
                taue = 3;
                taui = 5;
                synaptic_delay = 1;
                
                gL_PYR = 6;
                gL_SOM = 5;
                gL_PV = 5;
                gL_Int = 5;
                C_PYR =  180;
                C_Int = 80;
                taum_PYR = C_PYR/gL_PYR;
                taum_Int = C_Int/gL_Int;
                
                % Specific parameters of AdEx model
                VT_PYR = -49; 
                VT_PV = -52;
                VT_SOM = -53;
                
                VR_PYR = -60;
                VR_PV = -60;
                VR_SOM = -60;
                
                Vcut = min([VT_PYR,VT_PV,VT_SOM])+40;
                
                a_PYR = 4;
                a_PV = 0;
                a_SOM = 4;
                
                b_PYR = 100;
                b_PV = 0;
                b_SOM = 90;
                
                tauw_PYR = 100;
                tauw_PV = 15;
                tauw_SOM = 40;
                
                % Connectivity weights
                we_PYR_PYR = 0.5;
                we_PYR_PV = 1;
                wi_PV_PYR = 2.8;
                wi_PV_PV = 2.5;
                wi_SOM_PYR = 2.2;
                wi_SOM_PV = 2.4;
                wi_SOM_SOM = 1.6;
                wi_PV_SOM = 1.6;
                
                % External constant inputs
                Iext_PYR = 100;
                Iext_PV = 50;
                Iext_SOM = 25;
                
                Amp_spike = 30;
                
                % Connectivity matrices (rows:presynaptic, columns:post-synaptic) 
                p_PYR_PYR = 0.5;
                p_PV_PV = 0.6;
                p_PYR_PV = 0.4;
                p_PV_PYR = 0.4;
                p_SOM_PYR = 0.4;
                p_SOM_PV = 0.3;
                p_SOM_SOM = 0.1;
                p_PV_SOM = 0.1;
                
                Connect_PYR_PYR = double((rand(N_PYR)>(1-p_PYR_PYR))) ; 
                Connect_PYR_PYR(logical(eye(size(Connect_PYR_PYR)))) = 0;
                Connect_PV_PV = double((rand(N_PV)>(1-p_PV_PV))) ;
                Connect_PV_PV(logical(eye(size(Connect_PV_PV,1)))) =0;
                Connect_SOM_SOM = double((rand(N_SOM)>(1-p_SOM_SOM))) ;
                Connect_SOM_SOM(logical(eye(size(Connect_SOM_SOM,1)))) =0;
                
                Connect_PYR_PV = double((rand(N_PYR,N_PV)>(1-p_PYR_PV)));
                Connect_PV_PYR = double((rand(N_PV,N_PYR)>(1-p_PV_PYR)));
                Connect_SOM_PYR = double((rand(N_SOM,N_PYR)>(1-p_SOM_PYR)));
                Connect_SOM_PV = double((rand(N_SOM,N_PV)>(1-p_SOM_PV)));
                Connect_PV_SOM = double((rand(N_PV,N_SOM)>(1-p_PV_SOM)));
                
                V_PYR = zeros(N_PYR,N_times);
                w_PYR = zeros(N_PYR,N_times);
                V_PV = zeros(N_PV,N_times);
                w_PV = zeros(N_PV,N_times);
                V_SOM = zeros(N_SOM,N_times);
                w_SOM = zeros(N_SOM,N_times);
                               
                % Initial conditions
                sigma_V_init=5;
                sigma_W_init=5;
                
                V_PYR(:,1) = EL + sigma_V_init*rand(N_PYR,1);
                V_PV(:,1) = EL + sigma_V_init*rand(N_PV,1);
                V_SOM(:,1) = EL + sigma_V_init*rand(N_SOM,1);
                
                w_PYR(:,1) =  a_PYR*(V_PYR(:,1)-EL)+ sigma_W_init*rand(N_PYR,1);
                w_PV(:,1) =  a_PV*(V_PV(:,1)-EL)+ sigma_W_init*rand(N_PV,1);
                w_SOM(:,1) =  a_SOM*(V_SOM(:,1)-EL)+ sigma_W_init*rand(N_SOM,1);
                
                n_spikes_init=1;
                
                ge_PYR = zeros(N_PYR,N_times);
                gi_PYR = zeros(N_PYR,N_times);
                ge_PV = zeros(N_PV,N_times);
                gi_PV = zeros(N_PV,N_times);
                gi_SOM = zeros(N_SOM,N_times);
                                
                ge_PYR(:,1) = we_PYR_PYR*N_PYR_REF/N_PYR * rand(N_PYR,1)*n_spikes_init;
                gi_PYR(:,1) = (wi_PV_PYR*N_PV_REF/N_PV + wi_SOM_PYR*N_SOM_REF/N_SOM)/2*rand(N_PYR,1)*n_spikes_init;
                ge_PV(:,1) = we_PYR_PV*N_PYR_REF/N_PYR * rand(N_PV,1)*n_spikes_init;
                gi_PV(:,1) = (wi_PV_PV*N_PV_REF/N_PV + wi_SOM_PV*N_SOM_REF/N_SOM)/2*rand(N_PV,1)*n_spikes_init;
                gi_SOM(:,1) = (wi_SOM_SOM*N_SOM_REF/N_SOM+wi_PV_SOM*N_PV_REF/N_PV)/2*rand(N_SOM,1)*n_spikes_init;
                
                result_exc_spikes_PYR = zeros(N_PYR,N_times);
                result_inh_spikes_PYR = zeros(N_PYR,N_times);
                result_exc_spikes_PV = zeros(N_PV,N_times);
                result_inh_spikes_PV = zeros(N_PV,N_times);
                result_inh_spikes_SOM = zeros(N_SOM,N_times);
                
                i=1;
                for t=1:dt:(T-dt)
                    %% Update and save the voltage + adaptation variables
                    
                    dzeta_PYR = randn(1,N_PYR);  % random input to each neuron
                    dzeta_PV = randn(1,N_PV);
                    dzeta_SOM = randn(1,N_SOM);
                    
                    index_PYR = find(V_PYR(:,i)<Vcut); % distinction between spiking neurons (-> update to Vreset) and non-spiking neurons (-> differential equation)
                    index_PV = find(V_PV(:,i)<Vcut);
                    index_SOM = find(V_SOM(:,i)<Vcut);
                                       
                    V_PYR(index_PYR,i+1)=V_PYR(index_PYR,i)+dt*(gL_PYR*(EL-V_PYR(index_PYR,i))+gL_PYR*DeltaT_PYR*exp((V_PYR(index_PYR,i)-VT_PYR)/DeltaT_PYR) + ge_PYR(index_PYR,i).*(Ee-V_PYR(index_PYR,i)) + gi_PYR(index_PYR,i).*(Ei-V_PYR(index_PYR,i)) + Iext_PYR + stimulus_total_PYR(index_PYR,i) - w_PYR(index_PYR,i))/C_PYR + (sigma*sqrt(dt)/sqrt(taum_PYR))*dzeta_PYR(index_PYR)';
                    w_PYR(index_PYR,i+1)=w_PYR(index_PYR,i)+dt*(a_PYR*(V_PYR(index_PYR,i)-EL) - w_PYR(index_PYR,i))/tauw_PYR ;

                    V_PV(index_PV,i+1)=V_PV(index_PV,i)+dt*(gL_PV*(EL-V_PV(index_PV,i))+gL_PV*DeltaT_PV*exp((V_PV(index_PV,i)-VT_PV)/DeltaT_PV) + ge_PV(index_PV,i).*(Ee-V_PV(index_PV,i)) + gi_PV(index_PV,i).*(Ei-V_PV(index_PV,i)) + Iext_PV + stimulus_total_PV(index_PV,i) - w_PV(index_PV,i))/C_Int + (sigma*sqrt(dt)/sqrt(taum_Int))*dzeta_PV(index_PV)';
                    w_PV(index_PV,i+1)=w_PV(index_PV,i)+dt*(a_PV*(V_PV(index_PV,i)-EL) - w_PV(index_PV,i))/tauw_PV ;
                    
                    V_SOM(index_SOM,i+1)=V_SOM(index_SOM,i)+dt*(gL_SOM*(EL-V_SOM(index_SOM,i))+gL_SOM*DeltaT_SOM*exp((V_SOM(index_SOM,i)-VT_SOM)/DeltaT_SOM) + gi_SOM(index_SOM,i).*(Ei-V_SOM(index_SOM,i))+ Iext_SOM + stimulus_total_SOM(index_SOM,i) - w_SOM(index_SOM,i))/C_Int + (sigma*sqrt(dt)/sqrt(taum_Int))*dzeta_SOM(index_SOM)';
                    w_SOM(index_SOM,i+1)=w_SOM(index_SOM,i)+dt*(a_SOM*(V_SOM(index_SOM,i)-EL) - w_SOM(index_SOM,i))/tauw_SOM ;
                    
                    
                    %% Update the inputs received at each synapses at every timestep
                    
                    Spikes_PYR = double(V_PYR(:,i)>Vcut);  % Boolean: spiking or not spiking
                    Spikes_PV = double(V_PV(:,i)>Vcut);
                    Spikes_SOM = double(V_SOM(:,i)>Vcut);
                    
                    
                    %output vector: number of inputs (0,1,2, etc.) received by each post-synaptic neurons, then multiplied by the weight of the connection
                    % how to normalize the weights? (here used the average number of
                    % pre-synaptic neurons targeting each postsynaptic neuron...
                    result_exc_spikes_PYR(:,i) = we_PYR_PYR*N_PYR_REF/N_PYR.*(Connect_PYR_PYR'*Spikes_PYR(:)); % attention, utilisation de la transpose de Connect
                    result_inh_spikes_PYR(:,i) = wi_PV_PYR*N_PV_REF/N_PV.*(Connect_PV_PYR'*Spikes_PV(:))+wi_SOM_PYR*N_SOM_REF/N_SOM.*(Connect_SOM_PYR'*Spikes_SOM(:));
                    result_exc_spikes_PV(:,i) = we_PYR_PV*N_PYR_REF/N_PYR.*(Connect_PYR_PV'*Spikes_PYR(:));
                    result_inh_spikes_PV(:,i) = wi_PV_PV*N_PV_REF/N_PV.*(Connect_PV_PV'*Spikes_PV(:))+wi_SOM_PV*N_SOM_REF/N_SOM.*(Connect_SOM_PV'*Spikes_SOM(:));
                    result_inh_spikes_SOM(:,i) = wi_SOM_SOM*N_SOM_REF/N_SOM.*(Connect_SOM_SOM'*Spikes_SOM(:))+wi_PV_SOM*N_PV_REF/N_PV.*(Connect_PV_SOM'*Spikes_PV(:));
                    
                    
                    if i>synaptic_delay/dt
                        ge_PYR(:,i+1) = ge_PYR(:,i) + dt*(-ge_PYR(:,i)/taue) + result_exc_spikes_PYR(:,i-synaptic_delay/dt);                        
                        gi_PYR(:,i+1) = gi_PYR(:,i) + dt*(-gi_PYR(:,i)/taui) + result_inh_spikes_PYR(:,i-synaptic_delay/dt) ;
                        ge_PV(:,i+1) = ge_PV(:,i) + dt*(-ge_PV(:,i)/taue) + result_exc_spikes_PV(:,i-synaptic_delay/dt);
                        gi_PV(:,i+1) = gi_PV(:,i) + dt*(-gi_PV(:,i)/taui) + result_inh_spikes_PV(:,i-synaptic_delay/dt);
                        gi_SOM(:,i+1) = gi_SOM(:,i) + dt*(-gi_SOM(:,i)/taui) + result_inh_spikes_SOM(:,i-synaptic_delay/dt);
                    else
                        ge_PYR(:,i+1) = ge_PYR(:,i) + dt*(-ge_PYR(:,i)/taue) ;
                        gi_PYR(:,i+1) = gi_PYR(:,i) + dt*(-gi_PYR(:,i)/taui);
                        ge_PV(:,i+1) = ge_PV(:,i) + dt*(-ge_PV(:,i)/taue) ;
                        gi_PV(:,i+1) = gi_PV(:,i) + dt*(-gi_PV(:,i)/taui) ;
                        gi_SOM(:,i+1) = gi_SOM(:,i) + dt*(-gi_SOM(:,i)/taui) ;                           
                    end
                    
                    index_PYR_s = find(Spikes_PYR);
                    index_PV_s = find(Spikes_PV);
                    index_SOM_s = find(Spikes_SOM);
                    
                    V_PYR(index_PYR_s,i)= Amp_spike;
                    V_PYR(index_PYR_s,i+1)=  VR_PYR;
                    w_PYR(index_PYR_s,i+1)= w_PYR(index_PYR_s,i)+b_PYR;
                    
                    V_PV(index_PV_s,i)= Amp_spike;
                    V_PV(index_PV_s,i+1)=  VR_PV;
                    w_PV(index_PV_s,i+1)= w_PV(index_PV_s,i)+b_PV;
                    
                    V_SOM(index_SOM_s,i)= Amp_spike;
                    V_SOM(index_SOM_s,i+1)=  VR_SOM;
                    w_SOM(index_SOM_s,i+1)= w_SOM(index_SOM_s,i)+b_SOM;                    
                    
                    i=i+1;
                    
                end
                
                %% Raster plots
                
                % With DBS
                PYR_spike_number_DBS = length(find(V_PYR(:,4000:end)==Amp_spike)); % retrait des 200 premi?res ms
                PV_spike_number_DBS = length(find(V_PV(:,4000:end)==Amp_spike));
                SOM_spike_number_DBS = length(find(V_SOM(:,4000:end)==Amp_spike));
                                
                mean_firing_rate_PYR_DBS = PYR_spike_number_DBS/((T-200)*1e-3*N_PYR); % retrait des 200 premi?res ms
                mean_firing_rate_PV_DBS = PV_spike_number_DBS/((T-200)*1e-3*N_PV);
                mean_firing_rate_SOM_DBS = SOM_spike_number_DBS/((T-200)*1e-3*N_SOM);
                
                av_firing_rate_PYR_DBS(k,p) = mean_firing_rate_PYR_DBS;
                av_firing_rate_PV_DBS(k,p) = mean_firing_rate_PV_DBS;
                av_firing_rate_SOM_DBS(k,p) = mean_firing_rate_SOM_DBS;

            end
                        
        end
        M1 = [av_firing_rate_PYR_DBS ; av_firing_rate_PV_DBS ; av_firing_rate_SOM_DBS];
        dlmwrite(strcat('new_firing_rate_DBS_amp=',num2str(total_DBS_intensity(w)),'pA_pulse_dur=',num2str(total_pulse_duration(m)),'ms_freq_range=all'),M1);
    end
end


