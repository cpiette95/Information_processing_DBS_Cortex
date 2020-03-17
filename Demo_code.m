% Demo-code for visualizing network activity patterns in the different
% conditions (control, Parkinson, Parkinson+DBS, Parkinson+Opto-activation
% of interneurons)

% Choice of conditions:
conditions = {'control','noDBS','DBS','Opto_PV','Opto_SOM'};
% by default: DBS stimulation consists of 120 pA pulses of 2ms duration,
% repeated at 130Hz, targeting all neurons (pyramidal cells, SST and PV
% interneurons)
% Opto stimulation consists of 400 pA pulses of 3ms duration, repeated at
% 67 Hz and targeting only half of PV or SST interneurons

% Simulation parameters

T = 900; % duration of the simulation
dt = 0.05;
N_times=ceil(T/dt);
Time_vector = 0:dt:T;
Bin_Length =10 ; % ms bins (-> 50 bins for 500 ms stimulus duration)
incr=1;

% Number of cells in the network
N_PYR_REF = 800;
N_PV_REF = 120;
N_SOM_REF = 80;
N_PYR = ceil(N_PYR_REF*incr);
N_PV = ceil(N_PV_REF*incr);
N_SOM = ceil(N_SOM_REF*incr);
N_total = N_PYR + N_SOM + N_PV;


% Choice of additional stimuli:
% The stimulus can be a constant stimulus (type 1), a ramping input (type 2), a Ornstein-Uhlenbeck noise
% (type 3). However, the user can easily add any other inputs to the
% network. It is injected in a subset of pyramidal cells (but again, this
% parameters can be easily adjusted).

type = 3;
number_cell = 200; % number of pyramidal cells receiving the stimulus
receptor_PYR = randsample(N_PYR,number_cell) ; % index of pyramidal cells receiving the additional stimulus
stimulus_duration = 500; % duration during which the stimulus is presented to the network
Time_position = 200; % start of the stimulus injection

if type==1
    stimulus_intensity = input('What is the intensity of the constant stimulus being injected?');  % expected input: eg. 60 (will result in 60 pA injection)
elseif type==2
    start_ramp = input('What is the intensity for the start of the ramp?');
    end_ramp = input('What is the intensity for the end of the ramp?');
    stimulus_intensity = linspace(start_ramp,end_ramp,stimulus_duration/dt+1);
elseif type==3
    th = 0.1; % time constant
    mu = 120; % mean of OU noise
    sig = 14; % variance of OU noise
    t = 0:dt:stimulus_duration;
    x = zeros(1,length(t)-1);
    for i = 1:length(t)-1
        x(i+1) = x(i)+th*(mu-x(i))*dt+sig*sqrt(dt)*randn;
    end
    stimulus_intensity = x ;
end


for U=1:length(conditions)
    
    % Initialisation
    stimulus = zeros(N_PYR,N_times); % Additional stimulus to Pyramidal cells
    stimulus_total_PYR = zeros(N_PYR,N_times);
    stimulus_total_PV = zeros(N_PV,N_times);
    stimulus_total_SOM = zeros(N_SOM,N_times);
        
    i=1;
    for w=1:N_PYR
        if ismember(w,receptor_PYR)==1
            stimulus(w,Time_position/dt:Time_position/dt + stimulus_duration/dt) = stimulus_intensity; % additional stimuli only for a subset of pyramidal cells
        else
            non_receptor_PYR(i)=w ;
            i=i+1;
        end
    end
    figure(); plot(stimulus(receptor_PYR(1),:))
    
    % Default values
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
    
    if (strcmp(conditions{U},'noDBS') || strcmp(conditions{U},'control'))
        I0_PYR = 0;
        I0_PYR_bis = 0;
        I0_PV = 0;
        I0_SOM = 0;
    elseif strcmp(conditions{U},'DBS')
        I0_PYR = 120;
        I0_PYR_bis = 120;
        I0_PV = 120;
        I0_SOM = 120;
        frequency = 130*1e-3;
        pulse_duration = 2;
        proportion_activated = 1 ;
    elseif strcmp(conditions{U},'Opto_PV')
        I0_PYR = 0;
        I0_PYR_bis = 0;
        I0_PV = 600;
        I0_SOM = 0;
        delay_PV = 0;
        delay_SOM = 0;
        delay_PYR = 0;
        frequency = 67*1e-3;
        pulse_duration = 3;
        proportion_activated = 0.5 ;
    elseif strcmp(conditions{U},'Opto_SOM')
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
    end
    
    % DBS or Opto-activation inputs  (stimulus_total)
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
    
    % Parameters (detailed description found in Supplementary information of the Manuscript)
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
    
    % Connectivity matrices (rows: pre-synaptic, columns:post-synaptic)
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
    V_PYR(:,1) = EL + sigma_V_init*rand(N_PYR,1);
    V_PV(:,1) = EL + sigma_V_init*rand(N_PV,1);
    V_SOM(:,1) = EL + sigma_V_init*rand(N_SOM,1);
    w_PYR(:,1) =  a_PYR*(V_PYR(:,1)-EL)+ sigma_W_init*rand(N_PYR,1);
    w_PV(:,1) =  a_PV*(V_PV(:,1)-EL)+ sigma_W_init*rand(N_PV,1);
    w_SOM(:,1) = a_SOM*(V_SOM(:,1)-EL)+ sigma_W_init*rand(N_SOM,1);
    
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
        
        V_PYR(index_PYR,i+1)=V_PYR(index_PYR,i)+dt*(gL_PYR*(EL-V_PYR(index_PYR,i))+gL_PYR*DeltaT_PYR*exp((V_PYR(index_PYR,i)-VT_PYR)/DeltaT_PYR) + ge_PYR(index_PYR,i).*(Ee-V_PYR(index_PYR,i)) + gi_PYR(index_PYR,i).*(Ei-V_PYR(index_PYR,i)) + Iext_PYR + stimulus(index_PYR,i) + stimulus_total_PYR(index_PYR,i) - w_PYR(index_PYR,i))/C_PYR + (sigma*sqrt(dt)/sqrt(taum_PYR))*dzeta_PYR(index_PYR)';
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
    V_PYR_receptor_DBS = V_PYR(receptor_PYR,:);
    idx = find(V_PYR_receptor_DBS==Amp_spike);
    [row_PYR_receptor_DBS,col]=ind2sub(size(V_PYR_receptor_DBS),idx);
    PYR_spike_time_receptor_DBS = col*dt;
    V_PYR_non_receptor_DBS = V_PYR(non_receptor_PYR,:);
    idx = find(V_PYR_non_receptor_DBS==Amp_spike);
    [row_PYR_DBS,col]=ind2sub(size(V_PYR_non_receptor_DBS),idx);
    PYR_spike_time_DBS  = col*dt ;
    idx = find(V_PV==Amp_spike);
    [row_PV_DBS,col]=ind2sub(size(V_PV),idx);
    PV_spike_time_DBS = col*dt;
    idx = find(V_SOM==Amp_spike);
    [row_SOM_DBS,col]=ind2sub(size(V_SOM),idx);
    SOM_spike_time_DBS = col*dt;
    
    figure(); title(strcat('Raster plot for ',conditions(U)))
    hold on
    plot(PYR_spike_time_receptor_DBS,row_PYR_receptor_DBS,'k.','MarkerSize',10)
    plot(PYR_spike_time_DBS,number_cell+row_PYR_DBS,'k.','MarkerSize',10)
    plot(PV_spike_time_DBS,N_PYR+row_PV_DBS,'r.','MarkerSize',10)
    plot(SOM_spike_time_DBS,N_PYR+N_PV+row_SOM_DBS,'g.','MarkerSize',10)
    xlabel('Time (ms)')
    ylabel('Neuron index')
    %xlim([500 1400])
    %ylim([0 1010])
    hold off;
    
    
    
    %% Spike binning during stimulus
    V_PYR_binary=zeros(N_PYR,N_times);
    V_PYR_binary(V_PYR==Amp_spike)=1;
    Spike_binned_PSTH=zeros(N_PYR,ceil(stimulus_duration/Bin_Length)-1);
    v=1;
    for z=0:((stimulus_duration)/Bin_Length-1)
        Spike_binned_PSTH(:,v)=sum(V_PYR_binary(:,(Time_position)/dt+z*Bin_Length/dt:(Time_position)/dt+(z+1)*Bin_Length/dt),2);
        v=v+1;
    end
    Spike_binned_PSTH=vertcat(Spike_binned_PSTH(receptor_PYR,:),Spike_binned_PSTH(non_receptor_PYR,:));
    Spike_binning_PSTH =Spike_binned_PSTH;
    
    % Moving Proba during stimulus
    Count_Spike_PSTH = movsum(V_PYR_binary(:,(Time_position)/dt:(Time_position+stimulus_duration)/dt-1),Bin_Length/dt,2);
    Moving_Proba_spike_PSTH = sum(Count_Spike_PSTH,1);
    
    
    figure();
    suptitle(strcat(conditions{U}))
    subplot(2,1,1)
    plot(Moving_Proba_spike_PSTH)
    title('Moving Spike count of all Pyramidal cells')
    subplot(2,1,2)
    imagesc(Spike_binning_PSTH)
    title('Binned raster plot of all Pyramidal cells')
end


