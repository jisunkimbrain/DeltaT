% DeltaT task
% fMRI experiment
% GLM Analysis
% Ji Sun Kim

% Modeling (1st-level Model for Within-Subject Analysis)

%%%%%%%%%%%%%%%%%%% Delta T Batch Script %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% SPM12 version %%%%%%%%%%%%%%%%%%%%%%%%

% To regress: all 3 TD categories
% Association (TD=0), Consecutive (TD=1), Longer TD (TD=2~20)

% With RT modeled as Parametric Modulators


%% Defining Parameters
clear all;
clc;
dur = 4;  %-- 0 or 4 (mvpa => 0)

% Define subjects by loading sbj_run_list.mat

% Set-up directory for:
% raw_dir
% first_level_dir
% onset_dir
% rp_dir
% batch_save_dir

%% Run Modeling & Estimation session

for sbj_num=1:length(sbj_run_list)
    
    clear all_run; clear run_dir;
    clear matlabbatch;
    
    %% Defining pathway
    Rawdata_dir = [raw_dir num2str(sbj_run_list{1,sbj_num})];

    dir_first_level = [first_level_dir '\' num2str(sbj_run_list{1,sbj_num})];

    if ~exist(dir_first_level)
        mkdir(dir_first_level)
    end
    
    % Defining directories for scans of each run
    cd(Rawdata_dir);
    dir_sub_struct = dir;
    filter = '\w*3x3x3\w*';
    ct_run = 1;
    for ct_dir_sub = 3:length(dir_sub_struct)
        str = dir_sub_struct(ct_dir_sub).name;
        matchStr = regexp(str,filter,'match');
        if isempty(matchStr) == 1
            continue
        else
            all_run{ct_run} = str;
            ct_run = ct_run+1;
        end
    end
    
    for run_ct=1:length(sbj_run_list{2,sbj_num})
        run_dir{run_ct} = all_run{sbj_run_list{2,sbj_num}{run_ct}};
    end
    
    rpfiles = {};
    for run_num=1:length(sbj_run_list{2,sbj_num})
        rpfiles{run_num,1} = [rp_dir '\' num2str(sbj_run_list{1,sbj_num}) '\spike_regs_rp_modified_run-0' num2str(sbj_run_list{2,sbj_num}{run_num}) '.txt'];
    end
    
    
    %% Conditions for Modeling: Corr vs Incorr (each TD)

    %     load([onset_dir '\' num2str(sbj_run_list{1,sbj_num}) '_reg_vector.mat']);
    load([onset_dir '\' num2str(sbj_run_list{1,sbj_num}) '_run_reg.mat']);
    
    %% Start Batch
    Batch_Start_Time = datestr(now)
    
    %% Model Specification
    matlabbatch{1}.spm.stats.fmri_spec.dir = {dir_first_level};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
    
    % Functional images
    filter_s = 's*.img'; %-- 'ua*.img' for MVPA 
    for run_num=1:length(sbj_run_list{2,sbj_num})
        Rawdata_dir_run = [Rawdata_dir '\' run_dir{run_num}];
        dir_raw_struct = dir(fullfile(Rawdata_dir_run, filter_s));
        cell_dir = struct2cell(dir_raw_struct);
        sfiles = strcat(cell_dir(2,:)', {filesep}, cell_dir(1,:)');
        % only include scans from the beginning of retrieval phase (#104)
        sfiles(1:103) = [];

        % Data & Design
        % Subject/Session
        % Scanes
        matlabbatch{1}.spm.stats.fmri_spec.sess(run_num).scans = sfiles;
        
        % Conditions
        % Retrieval: TD_0
        for reg_ct=1:size(run_reg{1,run_num},2)
            matlabbatch{1}.spm.stats.fmri_spec.sess(run_num).cond(reg_ct).name = run_reg{1,run_num}{2,reg_ct};
            matlabbatch{1}.spm.stats.fmri_spec.sess(run_num).cond(reg_ct).onset = run_reg{1,run_num}{3,reg_ct};
            matlabbatch{1}.spm.stats.fmri_spec.sess(run_num).cond(reg_ct).duration = dur;
            matlabbatch{1}.spm.stats.fmri_spec.sess(run_num).cond(reg_ct).tmod = 0;
            matlabbatch{1}.spm.stats.fmri_spec.sess(run_num).cond(reg_ct).pmod.name = 'RT';
            matlabbatch{1}.spm.stats.fmri_spec.sess(run_num).cond(reg_ct).pmod.param = run_reg{1,run_num}{4,reg_ct};
            matlabbatch{1}.spm.stats.fmri_spec.sess(run_num).cond(reg_ct).pmod.poly = 1;
            matlabbatch{1}.spm.stats.fmri_spec.sess(run_num).cond(reg_ct).orth = 1;
        end  %-- end of for reg_ct
        
        % Multiple Conditions
        matlabbatch{1}.spm.stats.fmri_spec.sess(run_num).multi = {''};
        % Regressors
        matlabbatch{1}.spm.stats.fmri_spec.sess(run_num).regress = struct('name', {}, 'val', {});
        % Multiple regressors
        matlabbatch{1}.spm.stats.fmri_spec.sess(run_num).multi_reg = {rpfiles{run_num,1}};
        % High-pass filter
        matlabbatch{1}.spm.stats.fmri_spec.sess(run_num).hpf = 128;
    end
    
    % Basic Functions
    % Canonical HRF
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    % Model derivatives (no derivatives)
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    % Model Interactions (Volterra)
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    % Global normalisation
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    % Masking threshold
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
    % Explicit mask
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    % Serial correlations
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    
    
    %% Model Estimation
    matlabbatch{2}.spm.stats.fmri_est.spmmat = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    
    %% Saving directory/batch script
    %     save([batch_save_dir '\Modeling_Spec_Est_Sec_Dur_' num2str(dur) '_Ret_only_Corr_0_1_2to20_PM_RT_' num2str(sbj_run_list{1,sbj_num}) '.mat'], 'matlabbatch');
    save([batch_save_dir '\Modeling_Spec_Est_Sec_Dur_' num2str(dur) '_Ret_only_Corr_0_1_2to20_PM_RT_ua_' num2str(sbj_run_list{1,sbj_num}) '.mat'], 'matlabbatch');

    spm_jobman('run', matlabbatch);
    
end %-- end of for sbj_num

% Batch_End_Time = datestr(now)
