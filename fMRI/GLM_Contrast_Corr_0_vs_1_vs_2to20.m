% DeltaT task
% fMRI experiment
% GLM Analysis
% Ji Sun Kim

% Contrast Manager (1st-level Model for Within-Subject Analysis)

%%%%%%%%%%%%%%%%%%% Delta T Batch Script %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% SPM12 version %%%%%%%%%%%%%%%%%%%%%%%%

% To contrast among: all 3 TD categories
% Association (TD=0), Consecutive (TD=1), Longer TD (TD=2~20)

%% Defining Parameters
clear all;
clc;
dur = 4;  %-- 0 or 4 (mvpa => 0)

% Define subjects by loading sbj_run_list.mat

% Set-up directory for:
% ori_dir
% onset_dir
% batch_save_dir

%% Initialize SPM
spm('Defaults', 'fMRI');
spm_jobman('initcfg');

clear matlabbatch;

%% Contrast Manager: 1st-level Analysis
for sbj_num=1:length(sbj_run_list)
    
    % 1st-level subject directory
    sbj_dir = [ori_dir '\' num2str(sbj_run_list{1,sbj_num})];
    matlabbatch{sbj_num}.spm.stats.con.spmmat = {[sbj_dir '\SPM.mat']};
    
    % load sbj onset data
    load([onset_dir '\' num2str(sbj_run_list{1,sbj_num}) '_reg_vector.mat']);
    load([onset_dir '\' num2str(sbj_run_list{1,sbj_num}) '_run_reg.mat']);
    
    % Set up contrast weight
    % make sure contrast weights add up to 0
    condigit = zeros(11,length(reg_vector));
    conname = {};
    %     f_condigit = [];
    
    Corr_TD_0_idx = find(cell2mat(reg_vector(2,:))==1);
    Corr_TD_1_idx = find(cell2mat(reg_vector(2,:))==2);
    Corr_TD_2to20_idx = find(cell2mat(reg_vector(2,:))==3);
    
    for con_ct=1:11
        for reg_ct=1:length(reg_vector)
            
            if con_ct==1  %-- Corr_TD_0
                conname{con_ct} = 'Corr_0';
                for idx_ct=1:length(Corr_TD_0_idx)
                    condigit(con_ct,Corr_TD_0_idx(idx_ct)) = 1/length(Corr_TD_0_idx);
                end
            end  %-- end of if con_ct==1
            
            if con_ct==2  %-- Corr_TD_1
                conname{con_ct} = 'Corr_1';
                for idx_ct=1:length(Corr_TD_1_idx)
                    condigit(con_ct,Corr_TD_1_idx(idx_ct)) = 1/length(Corr_TD_1_idx);
                end
            end  %-- end of if con_ct==1
            
            if con_ct==3  %-- Corr_TD_2to20
                conname{con_ct} = 'Corr_2to20';
                for idx_ct=1:length(Corr_TD_2to20_idx)
                    condigit(con_ct,Corr_TD_2to20_idx(idx_ct)) = 1/length(Corr_TD_2to20_idx);
                end
            end  %-- end of if con_ct==1
            
            if con_ct==4  %-- Corr_TD_0 - Corr_TD_1
                conname{con_ct} = 'Corr_0_vs_Corr_1';
                for idx_ct=1:length(Corr_TD_0_idx)
                    condigit(con_ct,Corr_TD_0_idx(idx_ct)) = length(Corr_TD_1_idx);
                end
                for idx_ct=1:length(Corr_TD_1_idx)
                    condigit(con_ct,Corr_TD_1_idx(idx_ct)) = -length(Corr_TD_0_idx);
                end
            end  %-- end of if con_ct==1
            
            if con_ct==5  %-- Corr_TD_1 - Corr_TD_0
                conname{con_ct} = 'Corr_1_vs_Corr_0';
                for idx_ct=1:length(Corr_TD_0_idx)
                    condigit(con_ct,Corr_TD_0_idx(idx_ct)) = -length(Corr_TD_1_idx);
                end
                for idx_ct=1:length(Corr_TD_1_idx)
                    condigit(con_ct,Corr_TD_1_idx(idx_ct)) = length(Corr_TD_0_idx);
                end
            end
            
            if con_ct==6  %-- Corr_TD_0 - Corr_TD_2to20
                conname{con_ct} = 'Corr_0_vs_Corr_2to20';
                for idx_ct=1:length(Corr_TD_0_idx)
                    condigit(con_ct,Corr_TD_0_idx(idx_ct)) = length(Corr_TD_2to20_idx);
                end
                for idx_ct=1:length(Corr_TD_2to20_idx)
                    condigit(con_ct,Corr_TD_2to20_idx(idx_ct)) = - length(Corr_TD_0_idx);
                end
            end
            
            if con_ct==7  %-- Corr_TD_2to20 - Corr_TD_0
                conname{con_ct} = 'Corr_2to20_vs_Corr_0';
                for idx_ct=1:length(Corr_TD_0_idx)
                    condigit(con_ct,Corr_TD_0_idx(idx_ct)) = -length(Corr_TD_2to20_idx);
                end
                for idx_ct=1:length(Corr_TD_2to20_idx)
                    condigit(con_ct,Corr_TD_2to20_idx(idx_ct)) = length(Corr_TD_0_idx);
                end
            end
            
            if con_ct==8  %-- Corr_TD_1 - Corr_TD_2to20
                conname{con_ct} = 'Corr_1_vs_Corr_2to20';
                for idx_ct=1:length(Corr_TD_1_idx)
                    condigit(con_ct,Corr_TD_1_idx(idx_ct)) = length(Corr_TD_2to20_idx);
                end
                for idx_ct=1:length(Corr_TD_2to20_idx)
                    condigit(con_ct,Corr_TD_2to20_idx(idx_ct)) = -length(Corr_TD_1_idx);
                end
            end
            
            if con_ct==9  %-- Corr_TD_2to20 - Corr_TD_1
                conname{con_ct} = 'Corr_2to20_vs_Corr_1';
                for idx_ct=1:length(Corr_TD_1_idx)
                    condigit(con_ct,Corr_TD_1_idx(idx_ct)) = -length(Corr_TD_2to20_idx);
                end
                for idx_ct=1:length(Corr_TD_2to20_idx)
                    condigit(con_ct,Corr_TD_2to20_idx(idx_ct)) = length(Corr_TD_1_idx);
                end
            end
            
            if con_ct==10
                conname{con_ct} = 'Corr_0_vs_Corr_Seq';
                for idx_ct=1:length(Corr_TD_0_idx)
                    condigit(con_ct,Corr_TD_0_idx(idx_ct)) = length(Corr_TD_1_idx)+length(Corr_TD_2to20_idx);
                end
                for idx_ct=1:length(Corr_TD_1_idx)
                    condigit(con_ct,Corr_TD_1_idx(idx_ct)) = -length(Corr_TD_0_idx);
                end
                for idx_ct=1:length(Corr_TD_2to20_idx)
                    condigit(con_ct,Corr_TD_2to20_idx(idx_ct)) = -length(Corr_TD_0_idx);
                end
            end
            
            if con_ct==11
                conname{con_ct} = 'Corr_Seq_vs_Corr_0';
                for idx_ct=1:length(Corr_TD_0_idx)
                    condigit(con_ct,Corr_TD_0_idx(idx_ct)) = -length(Corr_TD_1_idx)-length(Corr_TD_2to20_idx);
                end
                for idx_ct=1:length(Corr_TD_1_idx)
                    condigit(con_ct,Corr_TD_1_idx(idx_ct)) = length(Corr_TD_0_idx);
                end
                for idx_ct=1:length(Corr_TD_2to20_idx)
                    condigit(con_ct,Corr_TD_2to20_idx(idx_ct)) = length(Corr_TD_0_idx);
                end
            end
            
        end  %-- end of for reg_ct
    end  %-- end of for con_ct
    
    
    % Contrast Sessions
    for consess_ct=1:length(conname)
        % Name
        matlabbatch{sbj_num}.spm.stats.con.consess{consess_ct}.tcon.name = conname{1,consess_ct};
        % Weights vector
        matlabbatch{sbj_num}.spm.stats.con.consess{consess_ct}.tcon.weights = condigit(consess_ct,:);
        % Replicate over sessions (default; Don't replicate)
        matlabbatch{sbj_num}.spm.stats.con.consess{consess_ct}.tcon.sessrep = 'none';
        % Delete existing contrasts (Yes)
        matlabbatch{sbj_num}.spm.stats.con.delete = 1;
    end  %-- end of consess_ct
 
end

%% Run & save batch

save([batch_save_dir '\Contrast_Sec_Dur_' num2str(dur) '_Ret_only_Corr_0_1_2to20_PM_RT.mat'], 'matlabbatch');

spm_jobman('run', matlabbatch);

% Batch_End_Time = datestr(now)
