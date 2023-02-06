% DeltaT task
% fMRI experiment
% BHV Analysis: Compare TD conditions (TD0, TD1, TD2-20)
% Ji Sun Kim

%% Directory Set-up
clear all; clc;

% Set-up directory for:
% sbj_dir
% sbj_dir_rt
% ori_dir
% save_dir

% Define subjects by loading sbj_list file

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
            
%% Plot results 

%% Accuracy

% load accuracy
load([save_dir '\TD_0_acc.mat']);
load([save_dir '\TD_1_acc.mat']);
load([save_dir '\TD_seq_acc.mat']);

% plot data
data = {TD_0_acc TD_1_acc TD_seq_acc};
[avg,err] = jh_mean_err(2,data);

labels = categorical({'TD=0', 'TD=1', 'TD=2~20'});
labels = reordercats(labels,string(labels));
y = avg;
figure;
ylim([0 1]);
resultsbar = bar(labels,y,'FaceColor','flat', 'BarWidth', 0.7);
ax = gca;

% ylabel('Accuracy', 'FontWeight', 'bold');

set(gca,'FontSize',18, 'linewidth',2, 'FontWeight','bold', 'box','off');
ax.XAxis.FontSize=20;
ax.XAxis.FontWeight='bold';

ax.XAxis.TickLength = [0 0];

%     ax.Children(2).CData(3,:) = [0.7 0.7 0.7];
resultsbar.CData(1,:) = [1 0.8 0.8];
resultsbar.CData(2,:) = [0.8 0.9 0.8];
resultsbar.CData(3,:) = [0.8 0.9 1];
hold on
bar(y(1:3), 'FaceAlpha',0,'EdgeColor',[0.3 0.5 0.6],'LineWidth',4, 'BarWidth', 0.7);
bar(y(1:2),'FaceAlpha',0,'EdgeColor',[0.3 0.6 0.3],'LineWidth',4, 'BarWidth', 0.7);
bar(y(1),'FaceAlpha',0,'EdgeColor',[0.8 0.1 0],'LineWidth',4, 'BarWidth', 0.7);


for i=1:3
    hold on
    plot([i i],[avg(i)-err(i) avg(i)+err(i)], 'Color', 'k', 'linewidth', 2);
end


hold on
yline(0.33,'k--','linewidth',2);

saveas(gcf, ['E:\Delta_T_Analysis\Behavioral\Figs_Jan2023\TDcondition_acc.png']);


%% RT

% load response time
load([save_dir '\mean_TD_0_rt.mat']);
load([save_dir '\mean_TD_1_rt.mat']);
load([save_dir '\mean_TD_seq_rt.mat']);

% plot data
data = {mean_TD_0_rt mean_TD_1_rt mean_TD_seq_rt};
[avg,err] = jh_mean_err(2,data);

labels = categorical({'TD=0', 'TD=1', 'TD=2~20'});
labels = reordercats(labels,string(labels));
y = avg;
figure;
resultsbar = bar(labels,y,'FaceColor','flat', 'BarWidth', 0.7);
ax = gca;
ax.XAxis.FontSize=15;
ax.XAxis.FontWeight='bold';
% ylabel('Response Time (ms)', 'FontWeight', 'bold');
ax = gca;
set(gca,'FontSize',18, 'linewidth',2, 'FontWeight','bold', 'box','off');
ax.XAxis.FontSize=20; ax.XAxis.FontWeight='bold';
ax.XAxis.TickLength = [0 0];

resultsbar.CData(1,:) = [1 0.8 0.8];
resultsbar.CData(2,:) = [0.8 0.9 0.8];
resultsbar.CData(3,:) = [0.8 0.9 1];
hold on
bar(y(1:3), 'FaceAlpha',0,'EdgeColor',[0.3 0.5 0.6],'LineWidth',4, 'BarWidth', 0.7);
bar(y(1:2),'FaceAlpha',0,'EdgeColor',[0.3 0.6 0.3],'LineWidth',4, 'BarWidth', 0.7);
bar(y(1),'FaceAlpha',0,'EdgeColor',[0.8 0.1 0],'LineWidth',4, 'BarWidth', 0.7);
hold on

ylim([0 3000]);


for i=1:3
    hold on
    plot([i i],[avg(i)-err(i) avg(i)+err(i)], 'Color', 'k', 'linewidth', 2);
end

saveas(gcf, ['E:\Delta_T_Analysis\Behavioral\Figs_Jan2023\TDcondition_RT.png']);
