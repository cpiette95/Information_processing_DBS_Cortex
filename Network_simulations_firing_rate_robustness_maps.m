%% CODE FOR REPRODUCING FIGURES 6.C, S5.C

dt = 0.05;
T = 1200;
N_times=ceil(T/dt);
Time_vector = 0:dt:T;
incr = 1;

N_PYR_REF = 800;
N_PV_REF = 120;
N_SOM_REF = 80;

N_PYR = ceil(N_PYR_REF*incr);
N_PV = ceil(N_PV_REF*incr);
N_SOM = ceil(N_SOM_REF*incr);
N_total = N_PYR + N_SOM + N_PV;

trial_number = 20;

% We test threee different conditions displayed in three heatmpas of the
% paper
% 1: corresponds to Figure S5.C (right) - varying input to PV cells and
% half of Pyramidal cells
% 2: corresponds to Figure S5.C (left) - varying input to SST cells and
% half of Pyramidal cells
% 3: corresponds to Figure 6.B - varying input to PV/SST cells, half
% Pyramidal cells receive a given input


for m=1:3
    
    if m==1 && m==2
        proportion_factor_1 = 0:0.1:1;
        proportion_factor_2 = -1:0.1:1;
    elseif m==3
        proportion_factor_1 = 0:0.1:1;
        proportion_factor_2 = 0:0.1:1;
    end
 
    av_firing_rate_PYR_DBS = {} ;
    av_firing_rate_PV_DBS = {} ;
    av_firing_rate_SOM_DBS = {};
    
    for k=1:length(proportion_factor_2)
        
        for j=1:length(proportion_factor_1)
            
            disp(['Proportion of DBS current to P1 =',num2str(proportion_factor_2(k))]);
            disp(['Proportion of DBS current to P2 =',num2str(proportion_factor_1(j))]);
            
            % DBS parameters (current intensity, pA)
            if m==1
                I0_PYR = 120;
                I0_PYR_bis = I0_PYR*proportion_factor_1(j);
                I0_PV = I0_PYR*proportion_factor_2(k);
                I0_SOM = 120;
            elseif m==2
                I0_PYR = 120;
                I0_PYR_bis = I0_PYR*proportion_factor_1(j);
                I0_PV = 120;
                I0_SOM = I0_PYR*proportion_factor_2(k);
                
            elseif m==3
                I0_PYR = 120;
                I0_PYR_bis = 0;
                I0_PV = I0_PYR*proportion_factor_1(j);
                I0_SOM = I0_PYR*proportion_factor_2(k);
            end
            
            % DBS parameters
            delay_PV = 0;
            delay_SOM = 2;
            delay_PYR = 2;
            frequency = 130*1e-3; % 130Hz stimulation
            pulse_duration = 2; % 2ms pulse duration
            
            stim_on_PYR = ones(round(pulse_duration/dt),ceil(N_PYR/2));
            stim_off_PYR = zeros((round(1/(frequency*dt))-size(stim_on_PYR,1)),ceil(N_PYR/2));
            stimulus_template_PYR = vertcat(stim_on_PYR,stim_off_PYR);
            stimulus_total_PYR = repmat(I0_PYR*stimulus_template_PYR',1,ceil(ceil(T/dt)/size(stimulus_template_PYR,1)));
            
            stim_off_pre_PYR_bis = zeros(round(delay_PYR/dt),ceil(N_PYR/2));
            stim_on_PYR_bis = ones(round(pulse_duration/dt),ceil(N_PYR/2));
            stim_off_PYR_bis = zeros(round(1/(frequency*dt) - size(stim_on_PYR_bis,1) -size(stim_off_pre_PYR_bis,1)),ceil(N_PYR/2));
            stimulus_template_PYR_bis = vertcat(stim_off_pre_PYR_bis,stim_on_PYR_bis,stim_off_PYR_bis);
            stimulus_total_PYR_bis= repmat(I0_PYR_bis*stimulus_template_PYR_bis',1,ceil(ceil(T/dt)/size(stimulus_template_PYR_bis,1)));
            
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
                
                % Specific parameters
                VT_PYR = -49; %-48 for control
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
                
                V_PYR_DBS = zeros(N_PYR,N_times);
                w_PYR_DBS = zeros(N_PYR,N_times);
                V_PV_DBS = zeros(N_PV,N_times);
                w_PV_DBS = zeros(N_PV,N_times);
                V_SOM_DBS = zeros(N_SOM,N_times);
                w_SOM_DBS = zeros(N_SOM,N_times);
                
                % Initial conditions
                sigma_V_init=5;
                sigma_W_init=5;
                
                V_PYR_DBS(:,1) = EL + sigma_V_init*rand(N_PYR,1);
                V_PV_DBS(:,1) = EL + sigma_V_init*rand(N_PV,1);
                V_SOM_DBS(:,1) = EL + sigma_V_init*rand(N_SOM,1);
                
                w_PYR_DBS(:,1) =  a_PYR*(V_PYR_DBS(:,1)-EL)+ sigma_W_init*rand(N_PYR,1);
                w_PV_DBS(:,1) =  a_PV*(V_PV_DBS(:,1)-EL)+ sigma_W_init*rand(N_PV,1);
                w_SOM_DBS(:,1) =a_SOM*(V_SOM_DBS(:,1)-EL)+ sigma_W_init*rand(N_SOM,1);
                
                n_spikes_init=1;
                
                ge_PYR_DBS = zeros(N_PYR,N_times);
                gi_PYR_DBS = zeros(N_PYR,N_times);
                ge_PV_DBS = zeros(N_PV,N_times);
                gi_PV_DBS = zeros(N_PV,N_times);
                gi_SOM_DBS = zeros(N_SOM,N_times);
                
                ge_PYR_DBS(:,1) = we_PYR_PYR*N_PYR_REF/N_PYR * rand(N_PYR,1)*n_spikes_init;
                gi_PYR_DBS(:,1) = (wi_PV_PYR*N_PV_REF/N_PV + wi_SOM_PYR*N_SOM_REF/N_SOM)/2*rand(N_PYR,1)*n_spikes_init;
                ge_PV_DBS(:,1) = we_PYR_PV*N_PYR_REF/N_PYR * rand(N_PV,1)*n_spikes_init;
                gi_PV_DBS(:,1) = (wi_PV_PV*N_PV_REF/N_PV + wi_SOM_PV*N_SOM_REF/N_SOM)/2*rand(N_PV,1)*n_spikes_init;
                gi_SOM_DBS(:,1) = (wi_SOM_SOM*N_SOM_REF/N_SOM+wi_PV_SOM*N_PV_REF/N_PV)/2*rand(N_SOM,1)*n_spikes_init;
                
                result_exc_spikes_PYR_DBS = zeros(N_PYR,N_times);
                result_inh_spikes_PYR_DBS = zeros(N_PYR,N_times);
                result_exc_spikes_PV_DBS = zeros(N_PV,N_times);
                result_inh_spikes_PV_DBS = zeros(N_PV,N_times);
                result_inh_spikes_SOM_DBS = zeros(N_SOM,N_times);
                
                i=1;
                
                for t=1:dt:(T-dt)
                    
                    % Generation of random and indendepent noisy inputs for all neurons
                    dzeta_PYR = randn(1,N_PYR);
                    dzeta_PV = randn(1,N_PV);
                    dzeta_SOM = randn(1,N_SOM);
                    
                    % Update of the voltage variables of non-spiking neurons
                    index_PYR_DBS = find(V_PYR_DBS(:,i)<Vcut);
                    index_PV_DBS = find(V_PV_DBS(:,i)<Vcut);
                    index_SOM_DBS = find(V_SOM_DBS(:,i)<Vcut);
                    
                    V_PYR_DBS(index_PYR_DBS,i+1)=V_PYR_DBS(index_PYR_DBS,i)+dt*(gL_PYR*(EL-V_PYR_DBS(index_PYR_DBS,i))+gL_PYR*DeltaT_PYR*exp((V_PYR_DBS(index_PYR_DBS,i)-VT_PYR)/DeltaT_PYR) + ge_PYR_DBS(index_PYR_DBS,i).*(Ee-V_PYR_DBS(index_PYR_DBS,i)) + gi_PYR_DBS(index_PYR_DBS,i).*(Ei-V_PYR_DBS(index_PYR_DBS,i)) + Iext_PYR + stimulus_total_PYR(index_PYR_DBS,i) - w_PYR_DBS(index_PYR_DBS,i))/C_PYR + (sigma*sqrt(dt)/sqrt(taum_PYR))*dzeta_PYR(index_PYR_DBS)';
                    w_PYR_DBS(index_PYR_DBS,i+1)=w_PYR_DBS(index_PYR_DBS,i)+dt*(a_PYR*(V_PYR_DBS(index_PYR_DBS,i)-EL) - w_PYR_DBS(index_PYR_DBS,i))/tauw_PYR ;
                    
                    V_PV_DBS(index_PV_DBS,i+1)=V_PV_DBS(index_PV_DBS,i)+dt*(gL_PV*(EL-V_PV_DBS(index_PV_DBS,i))+gL_PV*DeltaT_PV*exp((V_PV_DBS(index_PV_DBS,i)-VT_PV)/DeltaT_PV) + ge_PV_DBS(index_PV_DBS,i).*(Ee-V_PV_DBS(index_PV_DBS,i)) + gi_PV_DBS(index_PV_DBS,i).*(Ei-V_PV_DBS(index_PV_DBS,i)) + Iext_PV + stimulus_total_PV(index_PV_DBS,i) - w_PV_DBS(index_PV_DBS,i))/C_Int + (sigma*sqrt(dt)/sqrt(taum_Int))*dzeta_PV(index_PV_DBS)';
                    w_PV_DBS(index_PV_DBS,i+1)=w_PV_DBS(index_PV_DBS,i)+dt*(a_PV*(V_PV_DBS(index_PV_DBS,i)-EL) - w_PV_DBS(index_PV_DBS,i))/tauw_PV ;
                    
                    V_SOM_DBS(index_SOM_DBS,i+1)=V_SOM_DBS(index_SOM_DBS,i)+dt*(gL_SOM*(EL-V_SOM_DBS(index_SOM_DBS,i))+gL_SOM*DeltaT_SOM*exp((V_SOM_DBS(index_SOM_DBS,i)-VT_SOM)/DeltaT_SOM) + gi_SOM_DBS(index_SOM_DBS,i).*(Ei-V_SOM_DBS(index_SOM_DBS,i))+ Iext_SOM + stimulus_total_SOM(index_SOM_DBS,i) - w_SOM_DBS(index_SOM_DBS,i))/C_Int + (sigma*sqrt(dt)/sqrt(taum_Int))*dzeta_SOM(index_SOM_DBS)';
                    w_SOM_DBS(index_SOM_DBS,i+1)=w_SOM_DBS(index_SOM_DBS,i)+dt*(a_SOM*(V_SOM_DBS(index_SOM_DBS,i)-EL) - w_SOM_DBS(index_SOM_DBS,i))/tauw_SOM ;
                    
                    % Update of the conductance-based synapses
                    Spikes_PYR_DBS = double(V_PYR_DBS(:,i)>Vcut);
                    Spikes_PV_DBS = double(V_PV_DBS(:,i)>Vcut);
                    Spikes_SOM_DBS = double(V_SOM_DBS(:,i)>Vcut);
                    
                    result_exc_spikes_PYR_DBS(:,i) = we_PYR_PYR*N_PYR_REF/N_PYR.*(Connect_PYR_PYR'*Spikes_PYR_DBS(:));
                    result_inh_spikes_PYR_DBS(:,i) = wi_PV_PYR*N_PV_REF/N_PV.*(Connect_PV_PYR'*Spikes_PV_DBS(:))+wi_SOM_PYR*N_SOM_REF/N_SOM.*(Connect_SOM_PYR'*Spikes_SOM_DBS(:));
                    result_exc_spikes_PV_DBS(:,i) = we_PYR_PV*N_PYR_REF/N_PYR.*(Connect_PYR_PV'*Spikes_PYR_DBS(:));
                    result_inh_spikes_PV_DBS(:,i) = wi_PV_PV*N_PV_REF/N_PV.*(Connect_PV_PV'*Spikes_PV_DBS(:))+wi_SOM_PV*N_SOM_REF/N_SOM.*(Connect_SOM_PV'*Spikes_SOM_DBS(:));
                    result_inh_spikes_SOM_DBS(:,i) = wi_SOM_SOM*N_SOM_REF/N_SOM.*(Connect_SOM_SOM'*Spikes_SOM_DBS(:))+wi_PV_SOM*N_PV_REF/N_PV.*(Connect_PV_SOM'*Spikes_PV_DBS(:));
                    
                    if i>synaptic_delay/dt
                        ge_PYR_DBS(:,i+1) = ge_PYR_DBS(:,i) + dt*(-ge_PYR_DBS(:,i)/taue) + result_exc_spikes_PYR_DBS(:,i-synaptic_delay/dt);
                        gi_PYR_DBS(:,i+1) = gi_PYR_DBS(:,i) + dt*(-gi_PYR_DBS(:,i)/taui) + result_inh_spikes_PYR_DBS(:,i-synaptic_delay/dt) ;
                        ge_PV_DBS(:,i+1) = ge_PV_DBS(:,i) + dt*(-ge_PV_DBS(:,i)/taue) + result_exc_spikes_PV_DBS(:,i-synaptic_delay/dt);
                        gi_PV_DBS(:,i+1) = gi_PV_DBS(:,i) + dt*(-gi_PV_DBS(:,i)/taui) + result_inh_spikes_PV_DBS(:,i-synaptic_delay/dt);
                        gi_SOM_DBS(:,i+1) = gi_SOM_DBS(:,i) + dt*(-gi_SOM_DBS(:,i)/taui) + result_inh_spikes_SOM_DBS(:,i-synaptic_delay/dt);
                    else
                        ge_PYR_DBS(:,i+1) = ge_PYR_DBS(:,i) + dt*(-ge_PYR_DBS(:,i)/taue) ;
                        gi_PYR_DBS(:,i+1) = gi_PYR_DBS(:,i) + dt*(-gi_PYR_DBS(:,i)/taui) ;
                        ge_PV_DBS(:,i+1) = ge_PV_DBS(:,i) + dt*(-ge_PV_DBS(:,i)/taue) ;
                        gi_PV_DBS(:,i+1) = gi_PV_DBS(:,i) + dt*(-gi_PV_DBS(:,i)/taui) ;
                        gi_SOM_DBS(:,i+1) = gi_SOM_DBS(:,i) + dt*(-gi_SOM_DBS(:,i)/taui) ;
                    end
                    
                    % Reset after spike emission
                    index_PYR_s_DBS = find(Spikes_PYR_DBS);
                    index_PV_s_DBS = find(Spikes_PV_DBS);
                    index_SOM_s_DBS = find(Spikes_SOM_DBS);
                    
                    V_PYR_DBS(index_PYR_s_DBS,i)= Amp_spike;
                    V_PYR_DBS(index_PYR_s_DBS,i+1)= VR_PYR;
                    w_PYR_DBS(index_PYR_s_DBS,i+1)= w_PYR_DBS(index_PYR_s_DBS,i)+b_PYR;
                    
                    V_PV_DBS(index_PV_s_DBS,i)= Amp_spike;
                    V_PV_DBS(index_PV_s_DBS,i+1)= VR_PV;
                    w_PV_DBS(index_PV_s_DBS,i+1)= w_PV_DBS(index_PV_s_DBS,i)+b_PV;
                    
                    V_SOM_DBS(index_SOM_s_DBS,i)= Amp_spike;
                    V_SOM_DBS(index_SOM_s_DBS,i+1)= VR_SOM;
                    w_SOM_DBS(index_SOM_s_DBS,i+1)= w_SOM_DBS(index_SOM_s_DBS,i)+b_SOM;
                    
                    i=i+1;
                    
                end
                
                
                %% Average firing rate
                PYR_spike_number_DBS = length(find(V_PYR_DBS(:,ceil(200/dt):end)==Amp_spike)); % start from 200ms to end of stimulation
                PV_spike_number_DBS = length(find(V_PV_DBS(:,ceil(200/dt):end)==Amp_spike));
                SOM_spike_number_DBS = length(find(V_SOM_DBS(:,ceil(200/dt):end)==Amp_spike));
                
                mean_firing_rate_PYR_DBS = PYR_spike_number_DBS/((T-200)*1e-3*N_PYR);
                mean_firing_rate_PV_DBS = PV_spike_number_DBS/((T-200)*1e-3*N_PV);
                mean_firing_rate_SOM_DBS = SOM_spike_number_DBS/((T-200)*1e-3*N_SOM);
                
                av_firing_rate_PYR_DBS{k,j}(p) = mean_firing_rate_PYR_DBS;
                av_firing_rate_PV_DBS{k,j}(p) = mean_firing_rate_PV_DBS;
                av_firing_rate_SOM_DBS{k,j}(p) = mean_firing_rate_SOM_DBS;
                
            end
            
        end
        
    end
    
    save(strcat('firing_rate_DBS_130Hz_120pA_proportion_current_m=',num2str(m)), 'av_firing_rate_PYR_DBS', 'av_firing_rate_PV_DBS','av_firing_rate_SOM_DBS');

    A=load('new_firing_rate_noDBS_1'); 
    for r=1:length(proportion_factor_2)
        for s=1:length(proportion_factor_1)
            firing_rate_Pyr{m}(r,s) = mean(av_firing_rate_PYR_DBS{r,s})/mean(A(1,:));
            firing_rate_Pyr_bis{m}(r,s) = std(av_firing_rate_PYR_DBS{r,s}) ;
            firing_rate_PV{m}(r,s) = mean(av_firing_rate_PV_DBS{r,s})/mean(A(2,:));
            firing_rate_SST{m}(r,s) = mean(av_firing_rate_SOM_DBS{r,s})/mean(A(3,:));
        end
    end
    
    if m==1 % PYR and PV (Figure S5 right)
        figure() ; hold on
        imagesc(firing_rate_Pyr{m}')
        ylabel('% of I(DBS) applied to n=400 PYR')
        xlabel('% of I(DBS) applied to all PV')
        ylabels = {'0','20','40','60','80','100'};
        xlabels = {'-100','-80','-60','-40','-20','0','20','40','60','80','100'};
        set(gca,'XTick',1:2:21,'XTickLabel',xlabels);
        set(gca,'YTick',1:2:11,'YTickLabel',ylabels);
        caxis manual
        
    elseif m==2 %PYR and SST (Figure S5 left)
        figure() ; hold on
        imagesc(firing_rate_Pyr{m}')
        ylabel('Fraction of I(DBS) applied to n=400 PYR')
        xlabel('Fraction of I(DBS) applied to all SOM')
        ylabels = {'0','20','40','60','80','100'};
        xlabels = {'-100','-80','-60','-40','-20','0','20','40','60','80','100'};
        set(gca,'XTick',1:2:21,'XTickLabel',xlabels);
        set(gca,'YTick',1:2:11,'YTickLabel',ylabels);
        caxis manual
    elseif m==3 %PV and SST (Figure 6B)
        figure() ; hold on
        imagesc(firing_rate_Pyr{m})
        xlabel('Fraction of I(DBS) applied to all PV')
        ylabel('Fraction of I(DBS) applied to all SST')
        xlabels = {'0','20','40','60','80','100'};
        ylabels = {'0','20','40','60','80','100'};
        set(gca,'XTick',1:2:11,'XTickLabel',xlabels);
        set(gca,'YTick',1:2:11,'YTickLabel',ylabels);
        caxis manual
        caxis([0.4 1.6])
    end
       
end
