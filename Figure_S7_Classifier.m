%% CODE FOR DISPLAYING FIGURE S7c(after running Classifier_results_opto)


%cd('Saved_simulation_results/Information_processing/Classifier_scores')
number_classifiers = 1; % NN, LR, LDA and SVM
number_conditions = 11 ; % control, PD, PD+DBS, PD+Opto_SST, PD+Opto_PV
total_repetitions = 25; % 5*10

classifier = 'SVM'
cd('/Volumes/LaCie 1/DBS_PAPER/DBS_NEW/Classifier_matrices_basic')
M1={};    
M1{1} = importdata(strcat('score_control_',classifier),',');
M1{2} = importdata(strcat('score_noDBS_',classifier),',');
M1{3} = importdata(strcat('score_DBS_',classifier),',');
M1{6} = importdata(strcat('score_opto_PV_',classifier),',');
M1{10} = importdata(strcat('score_opto_SOM_',classifier),',');
cd('/Volumes/LaCie 1/DBS_PAPER/DBS_NEW/Classifier_matrices_opto')
M1{4} = importdata(strcat('score_opto_PV_200pA_0.5_',classifier),',');
M1{5} = importdata(strcat('score_opto_PV_400pA_0.5_',classifier),',');
M1{7} = importdata(strcat('score_opto_PV_800pA_0.5_',classifier),',');
M1{8} = importdata(strcat('score_opto_SOM_200pA_0.5_',classifier),',');
M1{9} = importdata(strcat('score_opto_SOM_400pA_0.5_',classifier),',');
M1{11} = importdata(strcat('score_opto_SOM_800pA_0.5_',classifier),',');

score_final = zeros(2*number_classifiers,number_conditions);
for p=1:number_conditions
    score_final(1,p) = nanmean(M1{p}(:)); 
    score_final(2,p) = nanstd(M1{p}(:)); 
end

% Statistical tests
M=M1;
[p,tbl,stats] = anova1([M{1}(:) M{2}(:) M{3}(:) M{4}(:) M{5}(:) M{6}(:) M{7}(:) M{8}(:) M{9}(:) M{10}(:) M{11}(:)],[],'on');
c=multcompare(stats);
Summary_Stats = c(:,[1,2,6]);


% Figure S7 - Classifier accuracy
figure() ; hold on ;
size_b=30;
scatter(ones(total_repetitions,1),M1{1}(:),size_b,'r','o')
scatter(1.1*ones(total_repetitions,1),M1{2}(:),size_b,'b','o')
scatter(1.2*ones(total_repetitions,1),M1{3}(:),size_b,'g','o')
scatter(1.6*ones(total_repetitions,1),M1{4}(:),size_b,'c','o')
scatter(1.7*ones(total_repetitions,1),M1{5}(:),size_b,'c','o')
scatter(1.8*ones(total_repetitions,1),M1{6}(:),size_b,'c','o')
scatter(1.9*ones(total_repetitions,1),M1{7}(:),size_b,'c','o')
scatter(2.4*ones(total_repetitions,1),M1{8}(:),size_b,'m','o')
scatter(2.5*ones(total_repetitions,1),M1{9}(:),size_b,'m','o')
scatter(2.6*ones(total_repetitions,1),M1{10}(:),size_b,'m','o')
scatter(2.7*ones(total_repetitions,1),M1{11}(:),size_b,'m','o')
e=errorbar([1,1.1,1.2,1.6,1.7,1.8,1.9,2.4,2.5,2.6,2.7],score_final(1,:),score_final(2,:),'s','MarkerSize',10,...
    'MarkerEdgeColor','black','MarkerFaceColor','none')
e.Color='black';
e.LineWidth=1;
ylabel('Classifier accuracy (%)')
ylim([0 100])
xlim([0.8 3])




