%% Code used for reproducing figures showing changes in firing rates in different conditions while different parameters of the stimulation are changed (frequency, intensity, pulse duration)

% The choice of conditions can be: 'control', 'noDBS', 'DBS', 'Opto_PV' or
% 'Opto_SOM'
% The choice of parameters to modify : for DBS stimulation: DBS_intensity,
% DBS_pulse duration and DBS_frequency (for recapitulating Figure 6B, Figure
% S5A, Figure S5B) ; for opto-stimulation: Opto_frequency (for
% recapitulating Figure 6D and S7A)


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

trial_number = 20;
conditions = 'DBS' ; %Opto_SOM'; % or noDBS or DBS or Opto_PV or Opto_SOM


%% DBS stimulation
% DBS stimulation (default: 120 pA, 130 Hz, 2 ms)
DBS_intensity = [80, 120, 160, 200,240];  %pA
DBS_pulse_duration = [1, 1.5, 2, 2.5, 3]; %ms
DBS_frequency = [18, 67, 130, 180]; %Hz


%% Opto stimulation
% Opto stimulation (default: 600 pA, 67 Hz, 3 ms)
Opto_intensity = [200, 400, 600,800];  %pA
Opto_frequency = [13, 67, 130]; %Hz

% % Choice of parameter
parameter = DBS_pulse_duration; % or DBS_pulse_duration or DBS_frequency or [1] (for control and PD)
str_parameter = 'DBS_pulse_duration'; %DBS_intensity'; %'Opto_frequency'; %Opto_frequency';

av_firing_rate_PYR = zeros(length(parameter),trial_number);
av_firing_rate_PV = zeros(length(parameter),trial_number);
av_firing_rate_SOM = zeros(length(parameter),trial_number);


for k=1:length(parameter)
    
    stimulus_total_PYR = zeros(N_PYR,N_times);
    stimulus_total_PV = zeros(N_PV,N_times);
    stimulus_total_SOM = zeros(N_SOM,N_times);
    
    % Default values for DBS
    I0_PYR = 120;
    I0_PYR_bis = 120;
    I0_PV = 120;
    I0_SOM = 120;
    delay_PV = 0;
    delay_SOM = 2;
    delay_PYR = 2;
    frequency = 130*1e-3;
    pulse_duration = 2;
    proportion_activated = 1 ;
    
    if (strcmp(conditions,'noDBS') || strcmp(conditions,'control'))  
        I0_PYR = 0;
        I0_PYR_bis = 0;
        I0_PV = 0;
        I0_SOM = 0;
        
    %%    DBS conditions
    elseif strcmp(conditions,'DBS') && strcmp(str_parameter,'DBS_intensity')
        I0_PYR = DBS_intensity(k);
        I0_PYR_bis = DBS_intensity(k);
        I0_PV = DBS_intensity(k);
        I0_SOM = DBS_intensity(k);
    elseif strcmp(conditions,'DBS') && strcmp(str_parameter,'DBS_frequency')
        frequency = DBS_frequency(k)*1e-3;
    elseif strcmp(conditions,'DBS') && strcmp(str_parameter,'DBS_pulse_duration')
        pulse_duration = DBS_pulse_duration(k);
        
        % Opto conditions
    elseif strcmp(conditions,'Opto_PV')
        I0_PYR = 0;
        I0_PYR_bis = 0;
        I0_PV = 600;
        I0_SOM = 0;
        frequency = 67*1e-3;
        delay_PV = 0;
        delay_SOM = 0;
        delay_PYR = 0;
        pulse_duration = 3;
        proportion_activated = 0.5 ;
        
        if strcmp(str_parameter,'Opto_intensity')
            I0_PV = Opto_intensity(k);
        elseif strcmp(str_parameter,'Opto_frequency')
            frequency = Opto_frequency(k)*1e-3;
        end
        
    elseif strcmp(conditions,'Opto_SOM')
        I0_PYR = 0;
        I0_PYR_bis = 0;
        I0_PV = 0;
        I0_SOM = 600;
        delay_PV = 0;
        delay_SOM = 0;
        delay_PYR = 0;
        frequency = 67*1e-3;
        pulse_duration = 3;
        proportion_activated = 0.5 ;
        
        if strcmp(str_parameter,'Opto_intensity')
            I0_SOM = Opto_intensity(k);
        elseif strcmp(str_parameter,'Opto_frequency')
            frequency = Opto_frequency(k)*1e-3;
        end
        
    else
        disp('Error')
    end
    
    
    pre_stim_PYR_template = vertcat(ones(round(pulse_duration/dt),1),zeros(round(1/(frequency*dt)-round(pulse_duration/dt)),1));
    stim_PYR_template = repmat(I0_PYR*pre_stim_PYR_template',1,ceil(ceil(T/dt)/length(pre_stim_PYR_template)));
    stim_PV_template = repmat(I0_PV*pre_stim_PYR_template',1,ceil(ceil(T/dt)/length(pre_stim_PYR_template)));
    
    pre_stim_PYR_delayed_template = vertcat(zeros(round(delay_PYR/dt),1),ones(round(pulse_duration/dt),1),zeros(round(1/(frequency*dt) - round(delay_PYR/dt) - round(pulse_duration/dt)),1));
    stim_PYR_delayed_template = repmat(I0_PYR_bis*pre_stim_PYR_delayed_template',1,ceil(ceil(T/dt)/length(pre_stim_PYR_delayed_template)));
    
    pre_stim_SOM_template = vertcat(zeros(round(delay_SOM/dt),1),ones(round(pulse_duration/dt),1),zeros(round(1/(frequency*dt) - round(delay_SOM/dt) - round(pulse_duration/dt)),1));
    stim_SOM_template = repmat(I0_SOM*pre_stim_SOM_template',1,ceil(ceil(T/dt)/length(pre_stim_SOM_template)));
    
    stimulus_total_PYR = [repmat(stim_PYR_template,proportion_activated*ceil(N_PYR/2),1) ; repmat(stim_PYR_delayed_template,proportion_activated*ceil(N_PYR/2),1) ; zeros(N_PYR - proportion_activated*N_PYR,length(stim_PYR_template))];
    stimulus_total_PV = [repmat(stim_PV_template,proportion_activated*N_PV,1) ; zeros(N_PV - proportion_activated*N_PV,length(stim_PV_template))];
    stimulus_total_SOM = [repmat(stim_SOM_template,proportion_activated*N_SOM,1); zeros(N_SOM - proportion_activated*N_SOM,length(stim_SOM_template))];
    
    
    
    for p=1:trial_number
        
        % Model parameters
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
        
        % Specific parameters
        
        if strcmp(conditions,'control')
            VT_PYR = -48;
        else
            VT_PYR = -49;
        end
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
        
        % External inputs
        Iext_PYR = 100;
        Iext_PV = 50;
        Iext_SOM = 25;
        
        Amp_spike = 30;
        
        % Connectivity matrices (les lignes: pr?-synaptiques, les colonnes:
        % post-synaptiques)
        p_PYR_PYR = 0.5;
        p_PV_PV = 0.6; 
        p_PYR_PV = 0.4;
        p_PV_PYR = 0.4;
        p_SOM_PYR = 0.4;
        p_SOM_PV = 0.3;
        p_SOM_SOM = 0.1;
        p_PV_SOM = 0.1;
        
        Connect_PYR_PYR = double((rand(N_PYR)>(1-p_PYR_PYR))) ; % matrix of 0s and 1s (no need for it to be symetrical)
        Connect_PYR_PYR(logical(eye(size(Connect_PYR_PYR)))) = 0;% fill with 0: no i and i connectivity
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
        
        % Same initial conditions for DBS on
        
        V_PYR(:,1) = EL + sigma_V_init*rand(N_PYR,1);
        V_PV(:,1) = EL + sigma_V_init*rand(N_PV,1);
        V_SOM(:,1) = EL + sigma_V_init*rand(N_SOM,1);
        
        w_PYR(:,1) =  a_PYR*(V_PYR(:,1)-EL)+ sigma_W_init*rand(N_PYR,1);
        w_PV(:,1) =  a_PV*(V_PV(:,1)-EL)+ sigma_W_init*rand(N_PV,1);
        w_SOM(:,1) =a_SOM*(V_SOM(:,1)-EL)+ sigma_W_init*rand(N_SOM,1);
        
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
            
            % Update the voltage/recovery variables of non-spiking neurons
            dzeta_PYR = randn(1,N_PYR);
            dzeta_PV = randn(1,N_PV);
            dzeta_SOM = randn(1,N_SOM);
            
            index_PYR = find(V_PYR(:,i)<Vcut);
            index_PV = find(V_PV(:,i)<Vcut);
            index_SOM = find(V_SOM(:,i)<Vcut);
            
            V_PYR(index_PYR,i+1)=V_PYR(index_PYR,i)+dt*(gL_PYR*(EL-V_PYR(index_PYR,i))+gL_PYR*DeltaT_PYR*exp((V_PYR(index_PYR,i)-VT_PYR)/DeltaT_PYR) + ge_PYR(index_PYR,i).*(Ee-V_PYR(index_PYR,i)) + gi_PYR(index_PYR,i).*(Ei-V_PYR(index_PYR,i)) + Iext_PYR + stimulus_total_PYR(index_PYR,i) - w_PYR(index_PYR,i))/C_PYR + (sigma*sqrt(dt)/sqrt(taum_PYR))*dzeta_PYR(index_PYR)';
            w_PYR(index_PYR,i+1)=w_PYR(index_PYR,i)+dt*(a_PYR*(V_PYR(index_PYR,i)-EL) - w_PYR(index_PYR,i))/tauw_PYR ;
            
            V_PV(index_PV,i+1)=V_PV(index_PV,i)+dt*(gL_PV*(EL-V_PV(index_PV,i))+gL_PV*DeltaT_PV*exp((V_PV(index_PV,i)-VT_PV)/DeltaT_PV) + ge_PV(index_PV,i).*(Ee-V_PV(index_PV,i)) + gi_PV(index_PV,i).*(Ei-V_PV(index_PV,i)) + Iext_PV + stimulus_total_PV(index_PV,i) - w_PV(index_PV,i))/C_Int + (sigma*sqrt(dt)/sqrt(taum_Int))*dzeta_PV(index_PV)';
            w_PV(index_PV,i+1)=w_PV(index_PV,i)+dt*(a_PV*(V_PV(index_PV,i)-EL) - w_PV(index_PV,i))/tauw_PV ;
            
            V_SOM(index_SOM,i+1)=V_SOM(index_SOM,i)+dt*(gL_SOM*(EL-V_SOM(index_SOM,i))+gL_SOM*DeltaT_SOM*exp((V_SOM(index_SOM,i)-VT_SOM)/DeltaT_SOM) + gi_SOM(index_SOM,i).*(Ei-V_SOM(index_SOM,i))+ Iext_SOM + stimulus_total_SOM(index_SOM,i) - w_SOM(index_SOM,i))/C_Int + (sigma*sqrt(dt)/sqrt(taum_Int))*dzeta_SOM(index_SOM)';
            w_SOM(index_SOM,i+1)=w_SOM(index_SOM,i)+dt*(a_SOM*(V_SOM(index_SOM,i)-EL) - w_SOM(index_SOM,i))/tauw_SOM ;
            
            
            % Update the conductance-based synapses
            Spikes_PYR = double(V_PYR(:,i)>Vcut);
            Spikes_PV = double(V_PV(:,i)>Vcut);
            Spikes_SOM = double(V_SOM(:,i)>Vcut);
            
            result_exc_spikes_PYR(:,i) = we_PYR_PYR*N_PYR_REF/N_PYR.*(Connect_PYR_PYR'*Spikes_PYR(:));
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
                gi_PYR(:,i+1) = gi_PYR(:,i) + dt*(-gi_PYR(:,i)/taui) ;
                ge_PV(:,i+1) = ge_PV(:,i) + dt*(-ge_PV(:,i)/taue) ;
                gi_PV(:,i+1) = gi_PV(:,i) + dt*(-gi_PV(:,i)/taui) ;
                gi_SOM(:,i+1) = gi_SOM(:,i) + dt*(-gi_SOM(:,i)/taui) ;
                
            end
            
            % Reset of spiking neurons
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
        
        %% Firing rate
        if strcmp(conditions,'control') || strcmp(conditions,'noDBS')  || strcmp(conditions,'DBS')
            
            PYR_spike_number = length(find(V_PYR(:,ceil(200/dt):end)==Amp_spike)); % start from 200ms to end of stimulation
            PV_spike_number = length(find(V_PV(:,ceil(200/dt):end)==Amp_spike));
            SOM_spike_number = length(find(V_SOM(:,ceil(200/dt):end)==Amp_spike));
            
            mean_firing_rate_PYR = PYR_spike_number/((T-200)*1e-3*N_PYR);
            disp(['Firing rate PYR =', num2str(mean_firing_rate_PYR)])
            mean_firing_rate_PV = PV_spike_number/((T-200)*1e-3*N_PV);
            disp(['Firing rate PV =', num2str(mean_firing_rate_PV)])
            mean_firing_rate_SOM = SOM_spike_number/((T-200)*1e-3*N_SOM);
            disp(['Firing rate SOM =', num2str(mean_firing_rate_SOM)])
            
            av_firing_rate_PYR(k,p) = mean_firing_rate_PYR;
            av_firing_rate_PV(k,p) = mean_firing_rate_PV;
            av_firing_rate_SOM(k,p) = mean_firing_rate_SOM;
            
        elseif strcmp(conditions,'Opto_PV')
            
            PYR_spike_number = length(find(V_PYR(:,ceil(200/dt):end)==Amp_spike)); % start from 200ms to end of stimulation
            PV_spike_number = length(find(V_PV(:,ceil(200/dt):end)==Amp_spike));
            PV_activated_spike_number = length(find(V_PV(1:proportion_activated*N_PV,ceil(200/dt):end)==Amp_spike));
            SOM_spike_number = length(find(V_SOM(:,ceil(200/dt):end)==Amp_spike));
            
            mean_firing_rate_PYR = PYR_spike_number/((T-200)*1e-3*N_PYR);
            mean_firing_rate_PV = PV_spike_number/((T-200)*1e-3*N_PV);
            mean_firing_rate_PV_activated = PV_activated_spike_number/((T-200)*1e-3*proportion_activated*N_PV);
            mean_firing_rate_SOM = SOM_spike_number/((T-200)*1e-3*N_SOM);
            
            av_firing_rate_PYR(k,p) = mean_firing_rate_PYR;
            av_firing_rate_PV(k,p) = mean_firing_rate_PV;
            av_firing_rate_PV_activated(k,p) = mean_firing_rate_PV_activated;
            av_firing_rate_SOM(k,p) = mean_firing_rate_SOM;
            
        elseif strcmp(conditions,'Opto_SOM')
            
            PYR_spike_number = length(find(V_PYR(:,ceil(200/dt):end)==Amp_spike)); % start from 200ms to end of stimulation
            PV_spike_number = length(find(V_PV(:,ceil(200/dt):end)==Amp_spike));
            SOM_spike_number = length(find(V_SOM(:,ceil(200/dt):end)==Amp_spike));
            SOM_activated_spike_number = length(find(V_SOM(1:proportion_activated*N_SOM,4000:end)==Amp_spike));
            
            mean_firing_rate_PYR = PYR_spike_number/((T-200)*1e-3*N_PYR);
            mean_firing_rate_PV = PV_spike_number/((T-200)*1e-3*N_PV);
            mean_firing_rate_SOM = SOM_spike_number/((T-200)*1e-3*N_SOM);
            mean_firing_rate_SOM_activated = SOM_activated_spike_number/((T-200)*1e-3*proportion_activated*N_SOM);
            
            av_firing_rate_PYR(k,p) = mean_firing_rate_PYR;
            av_firing_rate_PV(k,p) = mean_firing_rate_PV;
            av_firing_rate_SOM(k,p) = mean_firing_rate_SOM;
            av_firing_rate_SOM_activated(k,p) = mean_firing_rate_SOM_activated;
        end
                
        
    end
    
end

if strcmp(conditions,'control') || strcmp(conditions,'noDBS')  || strcmp(conditions,'DBS')
    M = [av_firing_rate_PYR ; av_firing_rate_PV ; av_firing_rate_SOM];
elseif strcmp(conditions,'Opto_PV');
    M = [av_firing_rate_PYR ; av_firing_rate_PV ; av_firing_rate_PV_activated ;av_firing_rate_SOM ];
elseif strcmp(conditions,'Opto_SOM');
    M = [av_firing_rate_PYR ; av_firing_rate_PV ; av_firing_rate_SOM ; av_firing_rate_SOM_activated];
end
dlmwrite(strcat('new_firing_rate_',conditions,'_',str_parameter),M);





