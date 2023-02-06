% DeltaT task
% fMRI experiment
% Create regressor/onset files
% Ji Sun Kim

%% Directory set-up

clear all; clc

% Set-up directory for:
% ori_dir
% save_dir
% rp_dir

% load sbj_run_list.mat

%% Create SOT

save_save_dir = [save_dir '\corr_incorr_0_1_2to20_RT'];

for sbj_ct=1:length(sbj_run_list)

    % load onset.mat
    load([ori_dir '\' num2str(sbj_run_list{1,sbj_ct}) '_onset.mat']);

    reg_vector = {};
    run_reg = {};
    run_reg_ct = 0;

    for run_ct=1:length(sbj_run_list{2,sbj_ct})

        corr_0 = []; corr_1 = []; corr_2to20 = [];
        incorr_0 = []; incorr_1 = []; incorr_2to20 = []; miss = [];
        corr_0_rt = []; corr_1_rt = []; corr_2to20_rt = [];
        incorr_0_rt = []; incorr_1_rt = []; incorr_2to20_rt = []; miss_rt = [];

        for trial_ct=1:20

            if strcmp(sbj_onset.Ret_data{run_ct}{trial_ct+1,23},'Correct') == 1
                
                if sbj_onset.Ret_data{run_ct}{trial_ct+1,14}==0
                    corr_0 = [corr_0 sbj_onset.Ret_data{run_ct}{trial_ct+1,3} - 206];
                    corr_0_rt = [corr_0_rt sbj_onset.Ret_data{run_ct}{trial_ct+1,24}];
                elseif sbj_onset.Ret_data{run_ct}{trial_ct+1,14}==1
                    corr_1 = [corr_1 sbj_onset.Ret_data{run_ct}{trial_ct+1,3} - 206];
                    corr_1_rt = [corr_1_rt sbj_onset.Ret_data{run_ct}{trial_ct+1,24}];
                elseif sbj_onset.Ret_data{run_ct}{trial_ct+1,14}>1
                    corr_2to20 = [corr_2to20 sbj_onset.Ret_data{run_ct}{trial_ct+1,3} - 206];
                    corr_2to20_rt = [corr_2to20_rt sbj_onset.Ret_data{run_ct}{trial_ct+1,24}];
                end

            elseif strcmp(sbj_onset.Ret_data{run_ct}{trial_ct+1,23},'Incorret') == 1

                if sbj_onset.Ret_data{run_ct}{trial_ct+1,14}==0
                    incorr_0 = [incorr_0 sbj_onset.Ret_data{run_ct}{trial_ct+1,3} - 206];
                    incorr_0_rt = [incorr_0_rt sbj_onset.Ret_data{run_ct}{trial_ct+1,24}];
                elseif sbj_onset.Ret_data{run_ct}{trial_ct+1,14}==1
                    incorr_1 = [incorr_1 sbj_onset.Ret_data{run_ct}{trial_ct+1,3} - 206];
                    incorr_1_rt = [incorr_1_rt sbj_onset.Ret_data{run_ct}{trial_ct+1,24}];
                elseif sbj_onset.Ret_data{run_ct}{trial_ct+1,14}>1
                    incorr_2to20 = [incorr_2to20 sbj_onset.Ret_data{run_ct}{trial_ct+1,3} - 206];
                    incorr_2to20_rt = [incorr_2to20_rt sbj_onset.Ret_data{run_ct}{trial_ct+1,24}];
                end

            elseif strcmp(sbj_onset.Ret_data{run_ct}{trial_ct+1,23},'Miss') == 1

                miss = [miss sbj_onset.Ret_data{run_ct}{trial_ct+1,3} - 206];
                miss_rt = [miss_rt 4000];
%                 miss_rt = [miss_rt sbj_onset.Ret_data{run_ct}{trial_ct+1,24}];

            end  %-- end of if Correct

        end  %-- end of for trial_ct

        % create regressor vectors
        reg_ct=0;
        if isempty(corr_0) == 0
            reg_ct=reg_ct+1;
            run_reg_ct=run_reg_ct+1;
            run_reg{1,run_ct}{1,reg_ct} = 1;
            run_reg{1,run_ct}{2,reg_ct} = 'Corr_0';
            run_reg{1,run_ct}{3,reg_ct} = corr_0;
            run_reg{1,run_ct}{4,reg_ct} = corr_0_rt;
            reg_vector{1,run_reg_ct} = run_ct;
            reg_vector{2,run_reg_ct} = 1;
            reg_vector{3,run_reg_ct} = 'Corr_0';
            reg_vector{4,run_reg_ct} = corr_0;
            run_reg_ct=run_reg_ct+1;
            reg_vector{1,run_reg_ct} = run_ct;
            reg_vector{2,run_reg_ct} = 11;
            reg_vector{3,run_reg_ct} = 'Corr_0_RT_PM';
            reg_vector{4,run_reg_ct} = corr_0_rt;
        end
        if isempty(corr_1) == 0
            reg_ct=reg_ct+1;
            run_reg_ct=run_reg_ct+1;
            run_reg{1,run_ct}{1,reg_ct} = 2;
            run_reg{1,run_ct}{2,reg_ct} = 'Corr_1';
            run_reg{1,run_ct}{3,reg_ct} = corr_1;
            run_reg{1,run_ct}{4,reg_ct} = corr_1_rt;
            reg_vector{1,run_reg_ct} = run_ct;
            reg_vector{2,run_reg_ct} = 2;
            reg_vector{3,run_reg_ct} = 'Corr_1';
            reg_vector{4,run_reg_ct} = corr_1;
            run_reg_ct=run_reg_ct+1;
            reg_vector{1,run_reg_ct} = run_ct;
            reg_vector{2,run_reg_ct} = 12;
            reg_vector{3,run_reg_ct} = 'Corr_1_RT_PM';
            reg_vector{4,run_reg_ct} = corr_1_rt;
        end
        if isempty(corr_2to20) == 0
            reg_ct=reg_ct+1;
            run_reg_ct=run_reg_ct+1;
            run_reg{1,run_ct}{1,reg_ct} = 3;
            run_reg{1,run_ct}{2,reg_ct} = 'Corr_2to20';
            run_reg{1,run_ct}{3,reg_ct} = corr_2to20;
            run_reg{1,run_ct}{4,reg_ct} = corr_2to20_rt;
            reg_vector{1,run_reg_ct} = run_ct;
            reg_vector{2,run_reg_ct} = 3;
            reg_vector{3,run_reg_ct} = 'Corr_2to20';
            reg_vector{4,run_reg_ct} = corr_2to20;
            run_reg_ct=run_reg_ct+1;
            reg_vector{1,run_reg_ct} = run_ct;
            reg_vector{2,run_reg_ct} = 13;
            reg_vector{3,run_reg_ct} = 'Corr_2to20_RT_PM';
            reg_vector{4,run_reg_ct} = corr_2to20_rt;
        end
        if isempty(incorr_0) == 0
            reg_ct=reg_ct+1;
            run_reg_ct=run_reg_ct+1;
            run_reg{1,run_ct}{1,reg_ct} = 4;
            run_reg{1,run_ct}{2,reg_ct} = 'Incorr_0';
            run_reg{1,run_ct}{3,reg_ct} = incorr_0;
            run_reg{1,run_ct}{4,reg_ct} = incorr_0_rt;
            reg_vector{1,run_reg_ct} = run_ct;
            reg_vector{2,run_reg_ct} = 4;
            reg_vector{3,run_reg_ct} = 'Incorr_0';
            reg_vector{4,run_reg_ct} = incorr_0;
            run_reg_ct=run_reg_ct+1;
            reg_vector{1,run_reg_ct} = run_ct;
            reg_vector{2,run_reg_ct} = 14;
            reg_vector{3,run_reg_ct} = 'Incorr_0_RT_PM';
            reg_vector{4,run_reg_ct} = incorr_0_rt;
        end
        if isempty(incorr_1) == 0
            reg_ct=reg_ct+1;
            run_reg_ct=run_reg_ct+1;
            run_reg{1,run_ct}{1,reg_ct} = 5;
            run_reg{1,run_ct}{2,reg_ct} = 'Incorr_1';
            run_reg{1,run_ct}{3,reg_ct} = incorr_1;
            run_reg{1,run_ct}{4,reg_ct} = incorr_1_rt;
            reg_vector{1,run_reg_ct} = run_ct;
            reg_vector{2,run_reg_ct} = 5;
            reg_vector{3,run_reg_ct} = 'Incorr_1';
            reg_vector{4,run_reg_ct} = incorr_1;
            run_reg_ct=run_reg_ct+1;
            reg_vector{1,run_reg_ct} = run_ct;
            reg_vector{2,run_reg_ct} = 15;
            reg_vector{3,run_reg_ct} = 'Incorr_1_RT_PM';
            reg_vector{4,run_reg_ct} = incorr_1_rt;
        end
        if isempty(incorr_2to20) == 0
            reg_ct=reg_ct+1;
            run_reg_ct=run_reg_ct+1;
            run_reg{1,run_ct}{1,reg_ct} = 6;
            run_reg{1,run_ct}{2,reg_ct} = 'Incorr_2to20';
            run_reg{1,run_ct}{3,reg_ct} = incorr_2to20;
            run_reg{1,run_ct}{4,reg_ct} = incorr_2to20_rt;
            reg_vector{1,run_reg_ct} = run_ct;
            reg_vector{2,run_reg_ct} = 6;
            reg_vector{3,run_reg_ct} = 'Incorr_2to20';
            reg_vector{4,run_reg_ct} = incorr_2to20;
            run_reg_ct=run_reg_ct+1;
            reg_vector{1,run_reg_ct} = run_ct;
            reg_vector{2,run_reg_ct} = 16;
            reg_vector{3,run_reg_ct} = 'Incorr_2to20_RT_PM';
            reg_vector{4,run_reg_ct} = incorr_2to20_rt;
        end
        if isempty(miss) == 0
            reg_ct=reg_ct+1;
            run_reg_ct=run_reg_ct+1;
            run_reg{1,run_ct}{1,reg_ct} = 7;
            run_reg{1,run_ct}{2,reg_ct} = 'Miss';
            run_reg{1,run_ct}{3,reg_ct} = miss;
            run_reg{1,run_ct}{4,reg_ct} = miss_rt;
            reg_vector{1,run_reg_ct} = run_ct;
            reg_vector{2,run_reg_ct} = 7;
            reg_vector{3,run_reg_ct} = 'Miss';
            reg_vector{4,run_reg_ct} = miss;
            run_reg_ct=run_reg_ct+1;
            reg_vector{1,run_reg_ct} = run_ct;
            reg_vector{2,run_reg_ct} = 17;
            reg_vector{3,run_reg_ct} = 'Miss_RT_PM';
            reg_vector{4,run_reg_ct} = miss_rt;
        end

        % add motion regressors
        for mcreg_ct=1:6
            reg_vector{1,run_reg_ct+mcreg_ct} = run_ct;
            reg_vector{2,run_reg_ct+mcreg_ct} = 9;
            reg_vector{3,run_reg_ct+mcreg_ct} = 'Motion Reg';
            reg_vector{4,run_reg_ct+mcreg_ct} = 0;
        end
        run_reg_ct=run_reg_ct+6;

        rpfile = load([rp_dir '\' num2str(sbj_run_list{1,sbj_ct}) '\spike_regs_rp_modified_run-0' num2str(sbj_run_list{2,sbj_ct}{run_ct}) '.txt']);
        if size(rpfile,2) > 6
            for spikereg_ct=1:size(rpfile,2)-6
                reg_vector{1,run_reg_ct+spikereg_ct} = run_ct;
                reg_vector{2,run_reg_ct+spikereg_ct} = 10;
                reg_vector{3,run_reg_ct+spikereg_ct} = 'Spike Reg';
                reg_vector{4,run_reg_ct+spikereg_ct} = 0;
            end
        end
        run_reg_ct=run_reg_ct+size(rpfile,2)-6;

    end  %-- end of for run_ct

     % add run regressors
    for runreg_ct=1:length(sbj_run_list{2,sbj_ct})
        reg_vector{1,run_reg_ct+runreg_ct} = 0;
        reg_vector{2,run_reg_ct+runreg_ct} = 21;
        reg_vector{3,run_reg_ct+runreg_ct} = 'Run Reg';
        reg_vector{4,run_reg_ct+runreg_ct} = 0;
    end

    % save .mat
    save([save_save_dir '\' num2str(sbj_run_list{1,sbj_ct}) '_run_reg.mat'], 'run_reg');
    save([save_save_dir '\' num2str(sbj_run_list{1,sbj_ct}) '_reg_vector.mat'], 'reg_vector');


end  %-- end of for sbj_ct
