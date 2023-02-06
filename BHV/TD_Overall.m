% DeltaT task
% fMRI experiment
% BHV Analysis: Calculate mean Acc & RT for each TD (0~20)
% Ji Sun Kim

%% Directory setup

clear all; clc;

% Set-up directory for:
% sbj_dir
% sbj_dir_rt
% ori_dir
% save_dir

% Define subjects by loading sbj_list file (sbj_run_list.mat)
  
  
%% Calculate Accuracy for each TD
for i =1:length(sbj_run_list)
    
    clear sbj; clear Enc_data_xls; clear Ret_data_xls;
    clear accuracy; clear accuracy_run; clear accuracy_overall;
    
    % Load data
    for sess=1:2
        for run=1:4
            Enc_data_xls{sess,run} = readcell([ori_Enc '\xls\Exp_results_Subj_' num2str(sbj_run_list{1,i}) '_Sess_' num2str(sess) '_Run_' num2str(run) '.xls']);
            Ret_data_xls{sess,run} = readcell([ori_Ret '\xls\Exp_results_Subj_' num2str(sbj_run_list{1,i}) '_Sess_' num2str(sess) '_Run_' num2str(run) '.xls']);
        end
    end  
    Enc_data = {Enc_data_xls{1,1:4} Enc_data_xls{2,1:4}};
    Ret_data = {Ret_data_xls{1,1:4} Ret_data_xls{2,1:4}};
    
    for sbj_run = 1: length(sbj_run_list{2,i})
        sbj.Enc_data{1,sbj_run} = Enc_data{sbj_run_list{2,i}{sbj_run}};
        sbj.Ret_data{1,sbj_run} = Ret_data{sbj_run_list{2,i}{sbj_run}};
    end
    
    % Bins
    bin_1_amount = 0;  bin_2_amount = 0;  bin_3_amount = 0;
    bin_4_amount = 0;  bin_5_amount = 0;  bin_6_amount = 0;
    bin_7_amount = 0;  bin_8_amount = 0;  bin_9_amount = 0;
    bin_10_amount = 0;  bin_11_amount = 0;  bin_12_amount = 0;
    bin_13_amount = 0;  bin_14_amount = 0;  bin_15_amount = 0;
    bin_16_amount = 0;  bin_17_amount = 0;  bin_18_amount = 0;
    bin_19_amount = 0;  bin_20_amount = 0;
    bin_0_amount= 0;
    
    bin_1_correct = 0;  bin_2_correct = 0;  bin_3_correct = 0;
    bin_4_correct = 0;  bin_5_correct = 0;  bin_6_correct = 0;
    bin_7_correct = 0;  bin_8_correct = 0;  bin_9_correct = 0;
    bin_10_correct = 0;  bin_11_correct = 0;  bin_12_correct = 0;
    bin_13_correct = 0;  bin_14_correct = 0;  bin_15_correct = 0;
    bin_16_correct = 0;  bin_17_correct = 0;  bin_18_correct = 0;
    bin_19_correct = 0;  bin_20_correct = 0;
    bin_0_correct = 0;
    
    for sbj_run=1:length(sbj_run_list{2,i})
        TempDiff_acc{sbj_run}= [cell2mat(sbj.Ret_data{sbj_run}(2:21,14)) sbj.accuracy{sbj_run}];
    end
    sbj.TempDiff_acc = TempDiff_acc;
    
    for run = 1:length(sbj_run_list{2,i})
        for trial = 1:20
            if sbj.TempDiff_acc{run}(trial,1) == 0
                bin_0_amount = bin_0_amount+1;
                bin_0_correct = bin_0_correct+sbj.TempDiff_acc{run}(trial,2);
            elseif sbj.TempDiff_acc{run}(trial,1) == 1
                bin_1_amount = bin_1_amount+1;
                bin_1_correct = bin_1_correct+sbj.TempDiff_acc{run}(trial,2);
            elseif sbj.TempDiff_acc{run}(trial,1) == 2
                bin_2_amount = bin_2_amount+1;
                bin_2_correct = bin_2_correct+sbj.TempDiff_acc{run}(trial,2);
            elseif sbj.TempDiff_acc{run}(trial,1) == 3
                bin_3_amount = bin_3_amount+1;
                bin_3_correct = bin_3_correct+sbj.TempDiff_acc{run}(trial,2);
            elseif sbj.TempDiff_acc{run}(trial,1) == 4
                bin_4_amount = bin_4_amount+1;
                bin_4_correct = bin_4_correct+sbj.TempDiff_acc{run}(trial,2);
            elseif sbj.TempDiff_acc{run}(trial,1) == 5
                bin_5_amount = bin_5_amount+1;
                bin_5_correct = bin_5_correct+sbj.TempDiff_acc{run}(trial,2);
            elseif sbj.TempDiff_acc{run}(trial,1) == 6
                bin_6_amount = bin_6_amount+1;
                bin_6_correct = bin_6_correct+sbj.TempDiff_acc{run}(trial,2);
            elseif sbj.TempDiff_acc{run}(trial,1) == 7
                bin_7_amount = bin_7_amount+1;
                bin_7_correct = bin_7_correct+sbj.TempDiff_acc{run}(trial,2);
            elseif sbj.TempDiff_acc{run}(trial,1) == 8
                bin_8_amount = bin_8_amount+1;
                bin_8_correct = bin_8_correct+sbj.TempDiff_acc{run}(trial,2);
            elseif sbj.TempDiff_acc{run}(trial,1) == 9
                bin_9_amount = bin_9_amount+1;
                bin_9_correct = bin_9_correct+sbj.TempDiff_acc{run}(trial,2);
            elseif sbj.TempDiff_acc{run}(trial,1) == 10
                bin_10_amount = bin_10_amount+1;
                bin_10_correct = bin_10_correct+sbj.TempDiff_acc{run}(trial,2);
            elseif sbj.TempDiff_acc{run}(trial,1) == 11
                bin_11_amount = bin_11_amount+1;
                bin_11_correct = bin_11_correct+sbj.TempDiff_acc{run}(trial,2);
            elseif sbj.TempDiff_acc{run}(trial,1) == 12
                bin_12_amount = bin_12_amount+1;
                bin_12_correct = bin_12_correct+sbj.TempDiff_acc{run}(trial,2);
            elseif sbj.TempDiff_acc{run}(trial,1) == 13
                bin_13_amount = bin_13_amount+1;
                bin_13_correct = bin_13_correct+sbj.TempDiff_acc{run}(trial,2);
            elseif sbj.TempDiff_acc{run}(trial,1) == 14
                bin_14_amount = bin_14_amount+1;
                bin_14_correct = bin_14_correct+sbj.TempDiff_acc{run}(trial,2);
            elseif sbj.TempDiff_acc{run}(trial,1) == 15
                bin_15_amount = bin_15_amount+1;
                bin_15_correct = bin_15_correct+sbj.TempDiff_acc{run}(trial,2);
            elseif sbj.TempDiff_acc{run}(trial,1) == 16
                bin_16_amount = bin_16_amount+1;
                bin_16_correct = bin_16_correct+sbj.TempDiff_acc{run}(trial,2);
            elseif sbj.TempDiff_acc{run}(trial,1) == 17
                bin_17_amount = bin_17_amount+1;
                bin_17_correct = bin_17_correct+sbj.TempDiff_acc{run}(trial,2);
            elseif sbj.TempDiff_acc{run}(trial,1) == 18
                bin_18_amount = bin_18_amount+1;
                bin_18_correct = bin_18_correct+sbj.TempDiff_acc{run}(trial,2);
            elseif sbj.TempDiff_acc{run}(trial,1) == 19
                bin_19_amount = bin_19_amount+1;
                bin_19_correct = bin_19_correct+sbj.TempDiff_acc{run}(trial,2);
            elseif sbj.TempDiff_acc{run}(trial,1) == 20
                bin_20_amount = bin_20_amount+1;
                bin_20_correct = bin_20_correct+sbj.TempDiff_acc{run}(trial,2);
            end
        end
    end
    
    sbj.TempDiff.bin_0_amount = bin_0_amount;   sbj.TempDiff.bin_0_correct = bin_0_correct;
    sbj.TempDiff.bin_1_amount = bin_1_amount;   sbj.TempDiff.bin_1_correct = bin_1_correct;
    sbj.TempDiff.bin_2_amount = bin_2_amount;   sbj.TempDiff.bin_2_correct = bin_2_correct;
    sbj.TempDiff.bin_3_amount = bin_3_amount;   sbj.TempDiff.bin_3_correct = bin_3_correct;
    sbj.TempDiff.bin_4_amount = bin_4_amount;   sbj.TempDiff.bin_4_correct = bin_4_correct;
    sbj.TempDiff.bin_5_amount = bin_5_amount;   sbj.TempDiff.bin_5_correct = bin_5_correct;
    sbj.TempDiff.bin_6_amount = bin_6_amount;   sbj.TempDiff.bin_6_correct = bin_6_correct;
    sbj.TempDiff.bin_7_amount = bin_7_amount;   sbj.TempDiff.bin_7_correct = bin_7_correct;
    sbj.TempDiff.bin_8_amount = bin_8_amount;   sbj.TempDiff.bin_8_correct = bin_8_correct;
    sbj.TempDiff.bin_9_amount = bin_9_amount;   sbj.TempDiff.bin_9_correct = bin_9_correct;
    sbj.TempDiff.bin_10_amount = bin_10_amount;   sbj.TempDiff.bin_10_correct = bin_10_correct;
    sbj.TempDiff.bin_11_amount = bin_11_amount;   sbj.TempDiff.bin_11_correct = bin_11_correct;
    sbj.TempDiff.bin_12_amount = bin_12_amount;   sbj.TempDiff.bin_12_correct = bin_12_correct;
    sbj.TempDiff.bin_13_amount = bin_13_amount;   sbj.TempDiff.bin_13_correct = bin_13_correct;
    sbj.TempDiff.bin_14_amount = bin_14_amount;   sbj.TempDiff.bin_14_correct = bin_14_correct;
    sbj.TempDiff.bin_15_amount = bin_15_amount;   sbj.TempDiff.bin_15_correct = bin_15_correct;
    sbj.TempDiff.bin_16_amount = bin_16_amount;   sbj.TempDiff.bin_16_correct = bin_16_correct;
    sbj.TempDiff.bin_17_amount = bin_17_amount;   sbj.TempDiff.bin_17_correct = bin_17_correct;
    sbj.TempDiff.bin_18_amount = bin_18_amount;   sbj.TempDiff.bin_18_correct = bin_18_correct;
    sbj.TempDiff.bin_19_amount = bin_19_amount;   sbj.TempDiff.bin_19_correct = bin_19_correct;
    sbj.TempDiff.bin_20_amount = bin_20_amount;   sbj.TempDiff.bin_20_correct = bin_20_correct;
    
    % Temporal Difference & accuracy
    Tempdiff_avg = [bin_0_correct/bin_0_amount bin_1_correct/bin_1_amount bin_2_correct/bin_2_amount ...
        bin_3_correct/bin_3_amount bin_4_correct/bin_4_amount bin_5_correct/bin_5_amount ...
        bin_6_correct/bin_6_amount bin_7_correct/bin_7_amount bin_8_correct/bin_8_amount ...
        bin_9_correct/bin_9_amount bin_10_correct/bin_10_amount bin_11_correct/bin_11_amount ...
        bin_12_correct/bin_12_amount bin_13_correct/bin_13_amount bin_14_correct/bin_14_amount ...
        bin_15_correct/bin_15_amount bin_16_correct/bin_16_amount bin_17_correct/bin_17_amount ...
        bin_18_correct/bin_18_amount bin_19_correct/bin_19_amount bin_20_correct/bin_20_amount];
    
    sbj.Tempdiff_avg = Tempdiff_avg;
    
    % save structure
    save([sbj_dir '\' num2str(sbj_run_list{1,i}) '_final.mat'], 'sbj');
    fprintf('saved\n');
end


%% Calculate RT for each TD

for sbj_ct =1:length(sbj_run_list)

    % load onset.mat
    load([ori_dir '\' num2str(sbj_run_list{1,sbj_ct}) '_onset.mat']);

    % Accuracy: Temporal Difference
    % Bins
    bin_1_amount = 0;  bin_2_amount = 0;  bin_3_amount = 0;
    bin_4_amount = 0;  bin_5_amount = 0;  bin_6_amount = 0;
    bin_7_amount = 0;  bin_8_amount = 0;  bin_9_amount = 0;
    bin_10_amount = 0;  bin_11_amount = 0;  bin_12_amount = 0;
    bin_13_amount = 0;  bin_14_amount = 0;  bin_15_amount = 0;
    bin_16_amount = 0;  bin_17_amount = 0;  bin_18_amount = 0;
    bin_19_amount = 0;  bin_20_amount = 0;
    bin_0_amount= 0;
    
    bin_1_rt = [];  bin_2_rt = [];  bin_3_rt = [];
    bin_4_rt = [];  bin_5_rt = [];  bin_6_rt = [];
    bin_7_rt = [];  bin_8_rt = [];  bin_9_rt = [];
    bin_10_rt = [];  bin_11_rt = [];  bin_12_rt = [];
    bin_13_rt = [];  bin_14_rt = [];  bin_15_rt = [];
    bin_16_rt = [];  bin_17_rt = [];  bin_18_rt = [];
    bin_19_rt = [];  bin_20_rt = [];
    bin_0_rt = [];
    
    
    for run_ct = 1:length(sbj_onset.Ret_data)
        for trial_ct = 1:20
            if sbj_onset.Ret_data{run_ct}{trial_ct+1,14} == 0 && sbj_onset.Ret_data{run_ct}{trial_ct+1,22} == 1
                bin_0_amount = bin_0_amount+1;
                bin_0_rt = [bin_0_rt sbj_onset.Ret_data{run_ct}{trial_ct+1,24}];
            elseif sbj_onset.Ret_data{run_ct}{trial_ct+1,14} == 1 && sbj_onset.Ret_data{run_ct}{trial_ct+1,22} == 1
                bin_1_amount = bin_1_amount+1;
                bin_1_rt = [bin_1_rt sbj_onset.Ret_data{run_ct}{trial_ct+1,24}];
            elseif sbj_onset.Ret_data{run_ct}{trial_ct+1,14} == 2 && sbj_onset.Ret_data{run_ct}{trial_ct+1,22} == 1
                bin_2_amount = bin_2_amount+1;
                bin_2_rt = [bin_2_rt sbj_onset.Ret_data{run_ct}{trial_ct+1,24}];
            elseif sbj_onset.Ret_data{run_ct}{trial_ct+1,14} == 3 && sbj_onset.Ret_data{run_ct}{trial_ct+1,22} == 1
                bin_3_amount = bin_3_amount+1;
                bin_3_rt = [bin_3_rt sbj_onset.Ret_data{run_ct}{trial_ct+1,24}];
            elseif sbj_onset.Ret_data{run_ct}{trial_ct+1,14} == 4 && sbj_onset.Ret_data{run_ct}{trial_ct+1,22} == 1
                bin_4_amount = bin_4_amount+1;
                bin_4_rt = [bin_4_rt sbj_onset.Ret_data{run_ct}{trial_ct+1,24}];
            elseif sbj_onset.Ret_data{run_ct}{trial_ct+1,14} == 5 && sbj_onset.Ret_data{run_ct}{trial_ct+1,22} == 1
                bin_5_amount = bin_5_amount+1;
                bin_5_rt = [bin_5_rt sbj_onset.Ret_data{run_ct}{trial_ct+1,24}];
            elseif sbj_onset.Ret_data{run_ct}{trial_ct+1,14} == 6 && sbj_onset.Ret_data{run_ct}{trial_ct+1,22} == 1
                bin_6_amount = bin_6_amount+1;
                bin_6_rt = [bin_6_rt sbj_onset.Ret_data{run_ct}{trial_ct+1,24}];
            elseif sbj_onset.Ret_data{run_ct}{trial_ct+1,14} == 7 && sbj_onset.Ret_data{run_ct}{trial_ct+1,22} == 1
                bin_7_amount = bin_7_amount+1;
                bin_7_rt = [bin_7_rt sbj_onset.Ret_data{run_ct}{trial_ct+1,24}];
            elseif sbj_onset.Ret_data{run_ct}{trial_ct+1,14} == 8 && sbj_onset.Ret_data{run_ct}{trial_ct+1,22} == 1
                bin_8_amount = bin_8_amount+1;
                bin_8_rt = [bin_8_rt sbj_onset.Ret_data{run_ct}{trial_ct+1,24}];
            elseif sbj_onset.Ret_data{run_ct}{trial_ct+1,14} == 9 && sbj_onset.Ret_data{run_ct}{trial_ct+1,22} == 1
                bin_9_amount = bin_9_amount+1;
                bin_9_rt = [bin_9_rt sbj_onset.Ret_data{run_ct}{trial_ct+1,24}];
            elseif sbj_onset.Ret_data{run_ct}{trial_ct+1,14} == 10 && sbj_onset.Ret_data{run_ct}{trial_ct+1,22} == 1
                bin_10_amount = bin_10_amount+1;
                bin_10_rt = [bin_10_rt sbj_onset.Ret_data{run_ct}{trial_ct+1,24}];
            elseif sbj_onset.Ret_data{run_ct}{trial_ct+1,14} == 11 && sbj_onset.Ret_data{run_ct}{trial_ct+1,22} == 1
                bin_11_amount = bin_11_amount+1;
                bin_11_rt = [bin_11_rt sbj_onset.Ret_data{run_ct}{trial_ct+1,24}];
            elseif sbj_onset.Ret_data{run_ct}{trial_ct+1,14} == 12 && sbj_onset.Ret_data{run_ct}{trial_ct+1,22} == 1
                bin_12_amount = bin_12_amount+1;
                bin_12_rt = [bin_12_rt sbj_onset.Ret_data{run_ct}{trial_ct+1,24}];
            elseif sbj_onset.Ret_data{run_ct}{trial_ct+1,14} == 13 && sbj_onset.Ret_data{run_ct}{trial_ct+1,22} == 1
                bin_13_amount = bin_13_amount+1;
                bin_13_rt = [bin_13_rt sbj_onset.Ret_data{run_ct}{trial_ct+1,24}];
            elseif sbj_onset.Ret_data{run_ct}{trial_ct+1,14} == 14 && sbj_onset.Ret_data{run_ct}{trial_ct+1,22} == 1
                bin_14_amount = bin_14_amount+1;
                bin_14_rt = [bin_14_rt sbj_onset.Ret_data{run_ct}{trial_ct+1,24}];
            elseif sbj_onset.Ret_data{run_ct}{trial_ct+1,14} == 15 && sbj_onset.Ret_data{run_ct}{trial_ct+1,22} == 1
                bin_15_amount = bin_15_amount+1;
                bin_15_rt = [bin_15_rt sbj_onset.Ret_data{run_ct}{trial_ct+1,24}];
            elseif sbj_onset.Ret_data{run_ct}{trial_ct+1,14} == 16 && sbj_onset.Ret_data{run_ct}{trial_ct+1,22} == 1
                bin_16_amount = bin_16_amount+1;
                bin_16_rt = [bin_16_rt sbj_onset.Ret_data{run_ct}{trial_ct+1,24}];
            elseif sbj_onset.Ret_data{run_ct}{trial_ct+1,14} == 17 && sbj_onset.Ret_data{run_ct}{trial_ct+1,22} == 1
                bin_17_amount = bin_17_amount+1;
                bin_17_rt = [bin_17_rt sbj_onset.Ret_data{run_ct}{trial_ct+1,24}];
            elseif sbj_onset.Ret_data{run_ct}{trial_ct+1,14} == 18 && sbj_onset.Ret_data{run_ct}{trial_ct+1,22} == 1
                bin_18_amount = bin_18_amount+1;
                bin_18_rt = [bin_18_rt sbj_onset.Ret_data{run_ct}{trial_ct+1,24}];
            elseif sbj_onset.Ret_data{run_ct}{trial_ct+1,14} == 19 && sbj_onset.Ret_data{run_ct}{trial_ct+1,22} == 1
                bin_19_amount = bin_19_amount+1;
                bin_19_rt = [bin_19_rt sbj_onset.Ret_data{run_ct}{trial_ct+1,24}];
            elseif sbj_onset.Ret_data{run_ct}{trial_ct+1,14} == 20 && sbj_onset.Ret_data{run_ct}{trial_ct+1,22} == 1
                bin_20_amount = bin_20_amount+1;
                bin_20_rt = [bin_20_rt sbj_onset.Ret_data{run_ct}{trial_ct+1,24}];
            end
        end
    end

    if bin_2_amount == 0
        bin_2_rt = NaN;
    end
    if bin_3_amount == 0
        bin_3_rt = NaN;
    end
    if bin_4_amount == 0
        bin_4_rt = NaN;
    end
    if bin_5_amount == 0
        bin_5_rt = NaN;
    end
    if bin_6_amount == 0
        bin_6_rt = NaN;
    end
    if bin_7_amount == 0
        bin_7_rt = NaN;
    end
    if bin_8_amount == 0
        bin_8_rt = NaN;
    end
    if bin_9_amount == 0
        bin_9_rt = NaN;
    end
    if bin_10_amount == 0
        bin_10_rt = NaN;
    end
    if bin_11_amount == 0
        bin_11_rt = NaN;
    end
    if bin_12_amount == 0
        bin_12_rt = NaN;
    end
    if bin_13_amount == 0
        bin_13_rt = NaN;
    end
    if bin_14_amount == 0
        bin_14_rt = NaN;
    end
    if bin_15_amount == 0
        bin_15_rt = NaN;
    end
    if bin_16_amount == 0
        bin_16_rt = NaN;
    end
    if bin_17_amount == 0
        bin_17_rt = NaN;
    end
    if bin_18_amount == 0
        bin_18_rt = NaN;
    end
    if bin_19_amount == 0
        bin_19_rt = NaN;
    end
    if bin_20_amount == 0
        bin_20_rt = NaN;
    end

    % Temporal Difference & accuracy
    Tempdiff_avg = [mean(bin_0_rt) mean(bin_1_rt) mean(bin_2_rt) ...
        mean(bin_3_rt) mean(bin_4_rt) mean(bin_5_rt) ...
        mean(bin_6_rt) mean(bin_7_rt) mean(bin_8_rt) ...
        mean(bin_9_rt) mean(bin_10_rt) mean(bin_11_rt) ...
        mean(bin_12_rt) mean(bin_13_rt) mean(bin_14_rt) ...
        mean(bin_15_rt) mean(bin_16_rt) mean(bin_17_rt) ...
        mean(bin_18_rt) mean(bin_19_rt) mean(bin_20_rt)];
    
    sbj.Tempdiff_avg = Tempdiff_avg;

    rt_Tempdiff_sbj = [rt_Tempdiff_sbj ; Tempdiff_avg];

    
    % save structure
    save([sbj_dir '\' num2str(sbj_run_list{1,sbj_ct}) '_TD_rt.mat'], 'sbj');
    fprintf([num2str(sbj_run_list{1,sbj_ct}) ' saved\n']);

end

save([ori_dir '\rt_Tempdiff_sbj.mat'], 'rt_Tempdiff_sbj');

%% Plot figures

%% Accuracy
TD_acc = [];

for sbj_ct=1:length(sbj_run_list)

    load([sbj_dir '\' num2str(sbj_run_list{1,i}) '_final.mat']);
    TD_acc(sbj_ct,:) = sbj.Tempdiff_avg;

end  %-- end of for sbj_ct

for td_ct=1:21
    data{1,td_ct} = TD_acc(:,td_ct);
end  %-- end of for td_ct
[avg,err] = jh_mean_err(2,data);

[p1,s1] = polyfit(1:20, avg(2:21)', 1);
y1 = polyval(p1,1:20);

figure;
% scatter([1:20], avg(2:21),'filled');
hold on; sim_scat = scatter([0], mean(TD_acc(:,1)), 'filled', 'red', 'h'); sim_scat.SizeData = 200;
hold on; er = errorbar([0], avg(1),err(1),'.', 'LineWidth',2); er.Color = 'red';
hold on; jbfill([1:20],(avg(2:21)+err(2:21)),(avg(2:21)-err(2:21)),[0.5 0.7 0.9],[0.5 0.7 0.9], 0, 0.2);
hold on; plot([1:20], avg(2:21), 'k--', 'linewidth', 3);
hold on; plot([1:20], y1, 'b-.', 'linewidth', 2.5);
xlim([-0.5 20]);
% ylabel('Accuracy'); xlabel('Temporal Distance');
set(gca,'FontName','Arial','FontSize',18,'fontweight','bold','linewidth',2, 'box','off');

%% RT
TD_rt= [];

for sbj_ct=1:length(sbj_run_list)

    load([sbj_dir '\' num2str(sbj_run_list{1,sbj_ct}) '_TD_rt.mat']);
    TD_rt(sbj_ct,:) = sbj.Tempdiff_avg;

end  %-- end of for sbj_ct

for td_ct=1:21
    data{1,td_ct} = TD_rt(:,td_ct);
end  %-- end of for td_ct
[avg,err] = jh_mean_err(2,data);

[p1,s1] = polyfit(1:20, avg(2:21)', 1);
y1 = polyval(p1,1:20);

figure;
% scatter([1:20], avg(2:21),'filled');
hold on; sim_scat = scatter([0], avg(1), 'filled', 'red', 'h'); sim_scat.SizeData = 200;
hold on; er = errorbar([0], avg(1),err(1),'.', 'LineWidth',2); er.Color = 'red';
hold on; jbfill([1:20],(avg(2:21)+err(2:21)),(avg(2:21)-err(2:21)),[0.5 0.7 0.9],[0.5 0.7 0.9], 0, 0.2);
hold on; plot([1:20], avg(2:21), 'k--', 'linewidth', 3);
hold on; plot([1:20], y1, 'b-.', 'linewidth', 2.5);
xlim([-0.5 20]);
ylim([1600 3000]);
% ylabel('Response Time (ms)'); xlabel('Temporal Distance');
set(gca,'FontName','Arial','FontSize',18,'fontweight','bold','linewidth',2, 'box','off');     
