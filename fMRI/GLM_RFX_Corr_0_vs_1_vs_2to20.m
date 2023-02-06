% DeltaT task
% fMRI experiment
% GLM Analysis
% Ji Sun Kim

% Random Effect Analysis (2nd-level Model for Between-Subject Analysis

%%%%%%%%%%%%%%%%%%%% MemRew Batch Script %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% SPM12 version %%%%%%%%%%%%%%%%%%%%%%%%

% Association (TD=0), Consecutive (TD=1), Longer TD (TD=2~20)

clear all;
clc;
dur = 4;  %-- 0 or 4
                          
Batch_Start_Time = datestr(now)
                          
%% Defining Parameters
clear all;
clc;
dur = 4;  %-- 0 or 4 (mvpa => 0)

% Define subjects by loading sbj_run_list.mat

% Set-up directory for:
% dir_first_level
% batch_save_dir
% ori_dir

% Same as Contrast batch code
conname{1} = 'Corr_0';
conname{2} = 'Corr_1';
conname{3} = 'Corr_2to20';
conname{4} = 'Corr_0_vs_Corr_1';
conname{5} = 'Corr_1_vs_Corr_0';
conname{6} = 'Corr_0_vs_Corr_2to20';
conname{7} = 'Corr_2to20_vs_Corr_0';
conname{8} = 'Corr_1_vs_Corr_2to20';
conname{9} = 'Corr_2to20_vs_Corr_1';
conname{10} = 'Corr_0_vs_Corr_Seq';
conname{11} = 'Corr_Seq_vs_Corr_0';

% directory that image files will be saved
dir_rfx = [pri_dir '\2nd_Level_RFX'];
if ~exist(dir_rfx)
    mkdir(dir_rfx);
end


%% Initialize SPM
spm('Defaults', 'fMRI');
spm_jobman('initcfg');

clear matlabbatch;

for ct_con=1:length(conname)
    ct_spec =  4 * ct_con - 3;
    ct_est =   ct_spec + 1;
    ct_ct =   ct_est + 1;
    ct_rep = ct_ct + 1;

    % make output directory
    dir_img = fullfile(dir_rfx, conname{ct_con});
    if ~exist(dir_img)
        mkdir(dir_img);
    end
    
    %% Factorial design specification
    % Directory
    matlabbatch{ct_spec}.spm.stats.factorial_design.dir = {dir_img};
    
    % Design
    % Scans
    % change directory when changing dur
    for ct_sub=1:length(sbj_run_list)
        conimgs(ct_sub,1) = cellstr(fullfile(dir_first_level, num2str(sbj_run_list{1,ct_sub}), ['\con_' num2str(ct_con, '%0.4d') '.nii,1']));
    end

    matlabbatch{ct_spec}.spm.stats.factorial_design.des.t1.scans = conimgs;
    
    % Covariates
    matlabbatch{ct_spec}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    % Multiple covariates
    matlabbatch{ct_spec}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    % Masking
    matlabbatch{ct_spec}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{ct_spec}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{ct_spec}.spm.stats.factorial_design.masking.em = {''};
    % Global calculation
    matlabbatch{ct_spec}.spm.stats.factorial_design.globalc.g_omit = 1;
    % Global normalisation
    matlabbatch{ct_spec}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    % Normalisation
    matlabbatch{ct_spec}.spm.stats.factorial_design.globalm.glonorm = 1;
    
    %% Model Estimation (2nd-level RFX)
    matlabbatch{ct_est}.spm.stats.fmri_est.spmmat = {[dir_img filesep 'SPM.mat']};
    matlabbatch{ct_est}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{ct_est}.spm.stats.fmri_est.method.Classical = 1;
    
    %% Contrast Manager (2nd-level RFX)
    matlabbatch{ct_ct}.spm.stats.con.spmmat = {[dir_img filesep 'SPM.mat']};
    % Name
    matlabbatch{ct_ct}.spm.stats.con.consess{1}.tcon.name = conname{ct_con};
    % Weights vector
    matlabbatch{ct_ct}.spm.stats.con.consess{1}.tcon.weights = 1;
    % Replicate over sessions (default; Don't replicate)
    matlabbatch{ct_ct}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    % Delete existing contrasts (Yes)
    matlabbatch{ct_ct}.spm.stats.con.delete = 1;
    
    %% Results Report
    matlabbatch{ct_rep}.spm.stats.results.spmmat = {[dir_img filesep 'SPM.mat']};
    matlabbatch{ct_rep}.spm.stats.results.conspec.titlestr = conname{ct_con};
    matlabbatch{ct_rep}.spm.stats.results.conspec.contrasts = 1;
    matlabbatch{ct_rep}.spm.stats.results.conspec.threshdesc = 'FWE';
    matlabbatch{ct_rep}.spm.stats.results.conspec.thresh = 0.05;   % GUI로는 0.005로 했음
    matlabbatch{ct_rep}.spm.stats.results.conspec.extent = 10;    % Nah.m에서는 05로 했음
    matlabbatch{ct_rep}.spm.stats.results.conspec.mask = struct('contrasts', {}, 'thresh', {}, 'mtype', {});
    matlabbatch{ct_rep}.spm.stats.results.units = 1;
    matlabbatch{ct_rep}.spm.stats.results.print = true;
    
end

%% Run & save batch

save([batch_save_dir '\RFX_Sec_Dur_', num2str(dur), '_Ret_only_Corr_0_1_2to20.mat'], 'matlabbatch');

spm_jobman('run', matlabbatch);

Batch_End_Time = datestr(now)
