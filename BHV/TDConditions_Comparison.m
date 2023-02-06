% DeltaT task
% fMRI experiment
% BHV Analysis: Compare TD conditions (TD0, TD1, TD2-20)
% Ji Sun Kim

%% Directory Set-up
clear all; clc;

sbj_dir = 'E:\Delta_T_Analysis\Behavioral\Subject_data_final';
sbj_dir_RT= 'E:\Delta_T_Analysis\Behavioral\Subject_data_final_RT';
ori_dir = 'D:\Delta_T_Analysis_Past\GLM_Nov\Onset_input';
% ori_dir = 'E:\Delta_T_Analysis\GLM_Nov\Onset_input';

% Defining subjects
% load('E:\Delta_T_Analysis\GLM_Nov\Onset_input\sbj_run_list_April.mat');
load('E:\Delta_T_Analysis\Retrieval_only\sbj_run_list_Jan.mat');

save_dir = 'E:\Delta_T_Analysis\Behavioral\Sim_Consec_Seq';
% save_dir = 'E:\Delta_T_Analysis\Behavioral\Sim_Seq';
% ~mkdir(save_dir);

%% Accuracy: 0 vs. 1 vs. 2~20

for sbj_ct=1:length(sbj_run_list)
    
    % load onset.mat
    load([ori_dir '\' num2str(sbj_run_list{1,sbj_ct}) '_onset.mat']);
    
    % Accuracy: Temporal Distance
    TD_0_amount = 0; TD_1_amount = 0; TD_seq_amount = 0;
    TD_0_corr = 0; TD_1_corr = 0; TD_seq_corr = 0;
    
    for run_ct=1:length(sbj_onset.Ret_data)
        for trial_ct=1:20
            if sbj_onset.Ret_data{run_ct}{trial_ct+1,14} == 0
                TD_0_amount = TD_0_amount+1;
                if sbj_onset.Ret_data{run_ct}{trial_ct+1,22} == 1
                    TD_0_corr = TD_0_corr+1;
                else
                    continue
                end
            elseif sbj_onset.Ret_data{run_ct}{trial_ct+1,14} == 1
                TD_1_amount = TD_1_amount+1;
                if sbj_onset.Ret_data{run_ct}{trial_ct+1,22} == 1
                    TD_1_corr = TD_1_corr+1;
                else
                    continue
                end
            elseif sbj_onset.Ret_data{run_ct}{trial_ct+1,14} > 1
                TD_seq_amount = TD_seq_amount+1;
                if sbj_onset.Ret_data{run_ct}{trial_ct+1,22} == 1
                    TD_seq_corr = TD_seq_corr+1;
                else
                    continue
                end
            end
        end  %-- end of for trial_ct
    end  % end of for run_ct
    
    TD_0_acc(sbj_ct) = TD_0_corr/TD_0_amount;
    TD_1_acc(sbj_ct) = TD_1_corr/TD_1_amount;
    TD_seq_acc(sbj_ct) = TD_seq_corr/TD_seq_amount;
    
end  %-- end of for sbj_ct

save_dir = 'E:\Delta_T_Analysis\Behavioral\Sim_Consec_Seq';
save([save_dir '\TD_0_acc.mat'], 'TD_0_acc');
save([save_dir '\TD_1_acc.mat'], 'TD_1_acc');
save([save_dir '\TD_seq_acc.mat'], 'TD_seq_acc');

          
%% Response Time: 0 vs. 1 vs. 2~20 

for sbj_ct=1:length(sbj_run_list)
    
    % load onset.mat
    load([ori_dir '\' num2str(sbj_run_list{1,sbj_ct}) '_onset.mat']);
    
    % RT: Temporal Distance
    TD_0_rt = []; TD_1_rt = []; TD_seq_rt= [];
    
    for run_ct=1:length(sbj_onset.Ret_data)
        for trial_ct=1:20
            if sbj_onset.Ret_data{run_ct}{trial_ct,24} == 0
                continue
            else
                if sbj_onset.Ret_data{run_ct}{trial_ct+1,14} == 0
                    TD_0_rt = [TD_0_rt sbj_onset.Ret_data{run_ct}{trial_ct+1,24}];
                elseif sbj_onset.Ret_data{run_ct}{trial_ct+1,14} == 1
                    TD_1_rt = [TD_1_rt sbj_onset.Ret_data{run_ct}{trial_ct+1,24}];
                elseif sbj_onset.Ret_data{run_ct}{trial_ct+1,14} > 1
                    TD_seq_rt = [TD_seq_rt sbj_onset.Ret_data{run_ct}{trial_ct+1,24}];
                end
            end
        end  %-- end of for trial_ct
    end  % end of for run_ct
    
    mean_TD_0_rt(sbj_ct) = mean(TD_0_rt);
    mean_TD_1_rt(sbj_ct) = mean(TD_1_rt);
    mean_TD_seq_rt(sbj_ct) = mean(TD_seq_rt);
    
end  %-- end of for sbj_ct

save_dir = 'E:\Delta_T_Analysis\Behavioral\Sim_Consec_Seq';
save([save_dir '\mean_TD_0_rt.mat'], 'mean_TD_0_rt');
save([save_dir '\mean_TD_1_rt.mat'], 'mean_TD_1_rt');
save([save_dir '\mean_TD_seq_rt.mat'], 'mean_TD_seq_rt');

%% Statistical Analysis (Bonferroni correction)

% Accuracy
load([save_dir '\TD_0_acc.mat']);
load([save_dir '\TD_1_acc.mat']);
load([save_dir '\TD_seq_acc.mat']);

data = table(TD_0_acc', TD_1_acc', TD_seq_acc', ...
    'VariableNames', {'TD0', 'TD1', 'TDseq'});
w = table(categorical([1 2 3].'), 'VariableNames', {'temporal'});
rm = fitrm(data, 'TD0-TDseq ~ 1', 'WithinDesign', w);
AT = ranova(rm, 'WithinModel', 'temporal');
disp(anovaTable(AT,'Value'));

pre_data = [TD_0_acc' TD_1_acc' TD_seq_acc'];
% pairs = [1 2; 1 3; 2 3]; alpha = 0.05; tail = 0;
% [h,p,sigPairs] = ttest_bonf(pre_data, pairs, alpha, tail)

[p,tbl,stats] = anova1(pre_data);
fprintf(['Accuracy: F(2,21) = ' num2str(tbl{2,5}) ', p = %0.3f \n'], p);

[results, means,h, gnames] = multcompare(stats, "CType", "bonferroni");
% tbl = array2table(means, "RowNames", {'TD=0', 'TD=1', 'TD=2'}, ...
%     "VariableNames", ["Mean", "Standard Error"])
tbl = array2table(results, "VariableNames", ...
     ["Group1", "Group2", "Lower Limit", "Difference","Upper Limit", "P-value"]);
[h1,p1,ci,stats] = ttest(pre_data(:,2), pre_data(:,1));
fprintf(['Accuracy: TD0 vs. TD1: t(23) = %0.3f p = %0.3f \n'], stats.tstat, tbl{1,6});
[h2,p2,ci,stats] = ttest(pre_data(:,1), pre_data(:,3));
fprintf(['Accuracy: TD0 vs. TD1: t(23) = %0.3f p = %0.3f \n'], stats.tstat, tbl{2,6});
[h3,p3,ci,stats] = ttest(pre_data(:,2), pre_data(:,3));
fprintf(['Accuracy: TD0 vs. TD1: t(23) = %0.3f p = %0.3f \n'], stats.tstat, tbl{3,6});


% Response time
load([save_dir '\mean_TD_0_rt.mat']);
load([save_dir '\mean_TD_1_rt.mat']);
load([save_dir '\mean_TD_seq_rt.mat']);

data = table(mean_TD_0_rt', mean_TD_1_rt', mean_TD_seq_rt', ...
    'VariableNames', {'TD0', 'TD1', 'TDseq'});
w = table(categorical([1 2 3].'), 'VariableNames', {'temporal'});
rm = fitrm(data, 'TD0-TDseq ~ 1', 'WithinDesign', w);
AT = ranova(rm, 'WithinModel', 'temporal');
disp(anovaTable(AT,'Value'));

clear pre_data
pre_data = [mean_TD_0_rt' mean_TD_1_rt' mean_TD_seq_rt'];

[p,tbl,stats] = anova1(pre_data);
fprintf(['RT: F(2,21) = ' num2str(tbl{2,5}) ', p = %0.3f \n'], p);

[results, means,h, gnames] = multcompare(stats, "CType", "bonferroni");

tbl = array2table(results, "VariableNames", ...
     ["Group1", "Group2", "Lower Limit", "Difference","Upper Limit", "P-value"]);

[h1,p1,ci,stats] = ttest(pre_data(:,1), pre_data(:,2));
fprintf(['RT: TD0 vs. TD1: t(23) = %0.3f p = %0.3f \n'], stats.tstat, p1);
[h2,p2,ci,stats] = ttest(pre_data(:,1), pre_data(:,3));
fprintf(['RT: TD0 vs. TD1: t(23) = %0.3f p = %0.3f \n'], stats.tstat, p2);
[h3,p3,ci,stats] = ttest(pre_data(:,2), pre_data(:,3));
fprintf(['RT: TD0 vs. TD1: t(23) = %0.3f p = %0.3f \n'], stats.tstat, p3);
            
