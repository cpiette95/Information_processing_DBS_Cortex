%% CODE FOR DISPLAYING FIGURE 7D (after running Classifier_results)


%cd('Saved_simulation_results/Information_processing/Classifier_scores')
number_classifiers = 4; % NN, LR, LDA and SVM
number_conditions = 5 ; % control, PD, PD+DBS, PD+Opto_SST, PD+Opto_PV
total_repetitions = 25; % 5*10

% M1: NN ; M2: LR ; M3 : LDA ; M4 : SVM
M1={};    
M1{1} = importdata(strcat('score_control_NN'),',');
M1{2} = importdata(strcat('score_noDBS_NN'),',');
M1{3} = importdata(strcat('score_DBS_NN'),',');
M1{4} = importdata(strcat('score_opto_PV_NN'),',');
M1{5} = importdata(strcat('score_opto_SOM_NN'),',');
M2={};    
M2{1} = importdata(strcat('score_control_LR'),',');
M2{2} = importdata(strcat('score_noDBS_LR'),',');
M2{3} = importdata(strcat('score_DBS_LR'),',');
M2{4} = importdata(strcat('score_opto_PV_LR'),',');
M2{5} = importdata(strcat('score_opto_SOM_LR'),','); 
M3={};    
M3{1} = importdata(strcat('score_control_LDA'),',');
M3{2} = importdata(strcat('score_noDBS_LDA'),',');
M3{3} = importdata(strcat('score_DBS_LDA'),',');
M3{4} = importdata(strcat('score_opto_PV_LDA'),',');
M3{5} = importdata(strcat('score_opto_SOM_LDA'),','); 
M4={};    
M4{1} = importdata(strcat('score_control_SVM'),',');
M4{2} = importdata(strcat('score_noDBS_SVM'),',');
M4{3} = importdata(strcat('score_DBS_SVM'),',');
M4{4} = importdata(strcat('score_opto_PV_SVM'),',');
M4{5} = importdata(strcat('score_opto_SOM_SVM'),','); 

score_final = zeros(2*number_classifiers,number_conditions);
for p=1:number_conditions
    score_final(1,p) = nanmean(M1{p}(:)); 
    score_final(2,p) = nanstd(M1{p}(:)); 
    score_final(3,p) = nanmean(M2{p}(:)); 
    score_final(4,p) = nanstd(M2{p}(:));
    score_final(5,p) = nanmean(M3{p}(:)); 
    score_final(6,p) = nanstd(M3{p}(:));
    score_final(7,p) = nanmean(M4{p}(:)); 
    score_final(8,p) = nanstd(M4{p}(:)); 
end


% Statistical tests
for U=1:4
    if U==1
        M=M1;
    elseif U==2
        M=M2;
    elseif U==3
        M=M3;
    elseif U==4
        M=M4;
    end
    [p,tbl,stats] = anova1([M{1}(:) M{2}(:) M{3}(:) M{4}(:) M{5}(:)],[],'on');
    c=multcompare(stats);
    Summary_Stats{U} = c(:,[1,2,6]);
end


% Figure 4E - Classifier accuracy
figure() ; hold on ;
size_b=10;
scatter(ones(total_repetitions,1),M1{1}(:),size_b,'r','o')
scatter(1.1*ones(total_repetitions,1),M1{2}(:),size_b,'b','o')
scatter(1.2*ones(total_repetitions,1),M1{3}(:),size_b,'g','o')
scatter(1.3*ones(total_repetitions,1),M1{4}(:),size_b,'c','o')
scatter(1.4*ones(total_repetitions,1),M1{5}(:),size_b,'m','o')

scatter(2*ones(total_repetitions,1),M2{1}(:),size_b,'r','o')
scatter(2.1*ones(total_repetitions,1),M2{2}(:),size_b,'b','o')
scatter(2.2*ones(total_repetitions,1),M2{3}(:),size_b,'g','o')
scatter(2.3*ones(total_repetitions,1),M2{4}(:),size_b,'c','o')
scatter(2.4*ones(total_repetitions,1),M2{5}(:),size_b,'m','o')

scatter(3*ones(total_repetitions,1),M3{1}(:),size_b,'r','o')
scatter(3.1*ones(total_repetitions,1),M3{2}(:),size_b,'b','o')
scatter(3.2*ones(total_repetitions,1),M3{3}(:),size_b,'g','o')
scatter(3.3*ones(total_repetitions,1),M3{4}(:),size_b,'c','o')
scatter(3.4*ones(total_repetitions,1),M3{5}(:),size_b,'m','o')

scatter(4*ones(total_repetitions,1),M4{1}(:),size_b,'r','o')
scatter(4.1*ones(total_repetitions,1),M4{2}(:),size_b,'b','o')
scatter(4.2*ones(total_repetitions,1),M4{3}(:),size_b,'g','o')
scatter(4.3*ones(total_repetitions,1),M4{4}(:),size_b,'c','o')
scatter(4.4*ones(total_repetitions,1),M4{5}(:),size_b,'m','o')

e=errorbar(1:4,score_final([1,3,5,7],1),score_final([2,4,6,8],1),'s','MarkerSize',10,...
    'MarkerEdgeColor','black','MarkerFaceColor','none')
e.Color='black';
e.LineWidth=1;
e=errorbar([1.1,2.1,3.1,4.1],score_final([1,3,5,7],2),score_final([2,4,6,8],2),'s','MarkerSize',10,...
    'MarkerEdgeColor','black','MarkerFaceColor','none')
e.Color='black';
e.LineWidth=1;
e=errorbar([1.2,2.2,3.2,4.2],score_final([1,3,5,7],3),score_final([2,4,6,8],3),'s','MarkerSize',10,...
    'MarkerEdgeColor','black','MarkerFaceColor','none')
e.Color='black';
e.LineWidth=1;
e=errorbar([1.3,2.3,3.3,4.3],score_final([1,3,5,7],4),score_final([2,4,6,8],4),'s','MarkerSize',10,...
    'MarkerEdgeColor','black','MarkerFaceColor','none')
e.Color='black'
e.LineWidth=1;
e=errorbar([1.4,2.4,3.4,4.4],score_final([1,3,5,7],5),score_final([2,4,6,8],5),'s','MarkerSize',10,...
    'MarkerEdgeColor','black','MarkerFaceColor','none')
e.Color='black'
e.LineWidth=1;
set(gca,'Xtick',[1 2 3 4 5])
set(gca,'XtickLabel',{'Algo NN', 'Algo MLR','Algo LDA','Algo SVM'})
ylabel('Classifier accuracy (%)')
legend('Control','PD','DBS 130 Hz','Opto PV','Opto SOM'); 
ylim([50 100])




