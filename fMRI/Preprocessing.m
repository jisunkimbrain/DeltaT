% DeltaT task
% fMRI experiment
% Preprocessing using SPM
% Ji Sun Kim

%% Defining pathway and subjects
clear all; clc;

Batch_Start_Time = datestr(now)

% start at directory: RawData;
dir_struct = dir;
dir_root = cd;

% Define subjects by loading sbj_list file

% Defining subject, run, volume
for ct_sub = 1:length(sbj_run_list)
    dir_sub = [dir_root '\' num2str(sbj_run_list{1,ct_sub})];
    cd(dir_sub);
    dir_sub_struct = dir;
    filter = '\w*3x3x3\w*';
    run_ct = 1;
    for ct_dir_sub = 3:length(dir_sub_struct)
        str = dir_sub_struct(ct_dir_sub).name;
        matchStr = regexp(str,filter,'match');
        if isempty(matchStr) == 1
            continue
        else
            run_dir{run_ct} = str;
            run_ct = run_ct+1;
        end
    end
    
    % functional image
     for ct_run=1:length(run_dir)
         dir_sub_run = [dir_sub '\' run_dir{ct_run}];
         cd(dir_sub_run);
         dir_sub_run_struct = dir(dir_sub_run);
         nvol(ct_run) = length(cellstr(ls('K*.img')));
         dir_sub_run_fimg_struct = dir('K*.img');
         for ct_vol=1:nvol
             scans{ct_run}{ct_vol,1} = [dir_sub_run '\' dir_sub_run_fimg_struct(ct_vol).name];
         end
     end
    
    % fieldmap image
    cd(dir_sub);
    filter ='\w*gre_field\w*';
    ct=1;
    for ct_dir_sub = 3:length(dir_sub_struct)
        str = dir_sub_struct(ct_dir_sub).name;
        matchStr = regexp(str,filter,'match');
        if isempty(matchStr) == 1
            continue
        else
            fieldmap_dir{ct} = str;
            ct = ct + 1;
        end
    end
    
    fieldmap_AP_dir = [dir_sub '\' fieldmap_dir{2}];
    cd(fieldmap_AP_dir);
    fm_vdm_img_dir = dir('vdm5_fpm*.img');
    fm_vdm_img = [fieldmap_AP_dir '\' fm_vdm_img_dir.name];
    
    % structural image (T1)
    cd(dir_sub);
    filter = '\w*T1\w*';
    ct=1;
    for ct_dir_sub = 3:length(dir_sub_struct)
        str = dir_sub_struct(ct_dir_sub).name;
        matchStr = regexp(str,filter,'match');
        if isempty(matchStr) == 1
            continue
        else
            structural_dir{ct} = str;
            ct = ct + 1;
        end
    end
    
    T1_dir = [dir_sub '\' structural_dir{1}];
    cd(T1_dir);
    T1_img_dir = dir('*.img');
    T1_img = [T1_dir '\' T1_img_dir.name];
    
    cd(dir_root);
    
    %% Initialize SPM
    
    spm('Defaults', 'fMRI');
    spm_jobman('initcfg');
    
    %% 1. Slice Timing
    
    matlabbatch{1}.spm.temporal.st.scans = scans;
    
    % Number of slices
    matlabbatch{1}.spm.temporal.st.nslices = 44;
    % TR
    matlabbatch{1}.spm.temporal.st.tr = 2;
    % TA = TR-(TR/Number of slices)
    matlabbatch{1}.spm.temporal.st.ta = 1.95454545454545;
    % Slice Order
    matlabbatch{1}.spm.temporal.st.so = [2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42 44 1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39 41 43];
    % Reference Slice
    matlabbatch{1}.spm.temporal.st.refslice = 22;
    % Filename Prefix
    matlabbatch{1}.spm.temporal.st.prefix = 'a';
    
    %% 2. Realignment & Unwarp
    
    for ct_run=1:length(run_dir)
        matlabbatch{2}.spm.spatial.realignunwarp.data(ct_run).scans(1) = cfg_dep(['Slice Timing: Slice Timing Corr. Images (Sess ' num2str(ct_run) ')'], substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{ct_run}, '.','files'));
        matlabbatch{2}.spm.spatial.realignunwarp.data(ct_run).pmscan = {[fm_vdm_img ',1']};
    end
    
    % Estimation Options (Default)
    % Quality
    matlabbatch{2}.spm.spatial.realignunwarp.eoptions.quality = 0.9;
    % Separation
    matlabbatch{2}.spm.spatial.realignunwarp.eoptions.sep = 4;
    % Smoothing (FWHM)
    matlabbatch{2}.spm.spatial.realignunwarp.eoptions.fwhm = 5;
    % Num Passes
    matlabbatch{2}.spm.spatial.realignunwarp.eoptions.rtm = 0;
    % Interpolation
    matlabbatch{2}.spm.spatial.realignunwarp.eoptions.einterp = 2;
    % Wrapping
    matlabbatch{2}.spm.spatial.realignunwarp.eoptions.ewrap = [0 0 0];
    % Weighting
    matlabbatch{2}.spm.spatial.realignunwarp.eoptions.weight = '';
    
    % Unwarp Estimation Options
    % Basis Functions
    matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
    % Regularisation
    matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
    % Reg. Factor
    matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
    % Jacobian deformations
    matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.jm = 0;
    % First-order effects
    matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
    % Second-order effects
    matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.sot = [];
    % Smoothing for unwarp (FWHM)
    matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
    % Re-estimate movement params
    matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.rem = 1;
    % Number of Iterations
    matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.noi = 5;
    % Taylor expansion point
    matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
    
    % Unwarp Reslicing Options
    % Resliced images (unwarp)?
    matlabbatch{2}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
    % Interpolation
    matlabbatch{2}.spm.spatial.realignunwarp.uwroptions.rinterp = 4;
    % Wrapping
    matlabbatch{2}.spm.spatial.realignunwarp.uwroptions.wrap = [0 0 0];
    % Masking
    matlabbatch{2}.spm.spatial.realignunwarp.uwroptions.mask = 1;
    % Filename Prefix
    matlabbatch{2}.spm.spatial.realignunwarp.uwroptions.prefix = 'u';
    
    %% 3. Coregister: Estimate
    % Reference Image
    matlabbatch{3}.spm.spatial.coreg.estimate.ref = {[T1_img ',1']};
    % Source Image
    matlabbatch{3}.spm.spatial.coreg.estimate.source(1) = cfg_dep('Realign & Unwarp: Unwarped Mean Image', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','meanuwr'));
    % Other Images
    matlabbatch{3}.spm.spatial.coreg.estimate.other = {''};
    % Estimation Options
    % Objective Function
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    % Separation
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    % Tolerances
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    % Histogram Smoothing
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    
    %% 4. Segment
    
    % Data
    % Channel
    % Volumes
    matlabbatch{4}.spm.spatial.preproc.channel.vols(1) = cfg_dep('Coregister: Estimate: Coregistered Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
    % Bias regularisation
    matlabbatch{4}.spm.spatial.preproc.channel.biasreg = 0.001;
    % Bias FWHM
    matlabbatch{4}.spm.spatial.preproc.channel.biasfwhm = 60;
    % Save Bias Corrected
    matlabbatch{4}.spm.spatial.preproc.channel.write = [0 0];
    
    % Tissues
    % Tissue
    % Tissue probability map
    matlabbatch{4}.spm.spatial.preproc.tissue(1).tpm = {'C:\Users\Jisun_PC\Documents\MATLAB\MATLAB_Toolboxes\spm12\tpm\TPM.nii,1'};
    % Num Gaussians
    matlabbatch{4}.spm.spatial.preproc.tissue(1).ngaus = 1;
    % Native Tissue
    matlabbatch{4}.spm.spatial.preproc.tissue(1).native = [1 0];
    % Warped Tissue
    matlabbatch{4}.spm.spatial.preproc.tissue(1).warped = [0 0];
    % Tissue
    matlabbatch{4}.spm.spatial.preproc.tissue(2).tpm = {'C:\Users\Jisun_PC\Documents\MATLAB\MATLAB_Toolboxes\spm12\tpm\TPM.nii,2'};
    matlabbatch{4}.spm.spatial.preproc.tissue(2).ngaus = 1;
    matlabbatch{4}.spm.spatial.preproc.tissue(2).native = [1 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(2).warped = [0 0];
    % Tissue
    matlabbatch{4}.spm.spatial.preproc.tissue(3).tpm = {'C:\Users\Jisun_PC\Documents\MATLAB\MATLAB_Toolboxes\spm12\tpm\TPM.nii,3'};
    matlabbatch{4}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{4}.spm.spatial.preproc.tissue(3).native = [1 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(3).warped = [0 0];
    % Tissue
    matlabbatch{4}.spm.spatial.preproc.tissue(4).tpm = {'C:\Users\Jisun_PC\Documents\MATLAB\MATLAB_Toolboxes\spm12\tpm\TPM.nii,4'};
    matlabbatch{4}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{4}.spm.spatial.preproc.tissue(4).native = [1 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(4).warped = [0 0];
    % Tissue
    matlabbatch{4}.spm.spatial.preproc.tissue(5).tpm = {'C:\Users\Jisun_PC\Documents\MATLAB\MATLAB_Toolboxes\spm12\tpm\TPM.nii,5'};
    matlabbatch{4}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{4}.spm.spatial.preproc.tissue(5).native = [1 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(5).warped = [0 0];
    % Tissue
    matlabbatch{4}.spm.spatial.preproc.tissue(6).tpm = {'C:\Users\Jisun_PC\Documents\MATLAB\MATLAB_Toolboxes\spm12\tpm\TPM.nii,6'};
    matlabbatch{4}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{4}.spm.spatial.preproc.tissue(6).native = [0 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(6).warped = [0 0];
    
    % Warping & MRF
    % MRF Parameter
    matlabbatch{4}.spm.spatial.preproc.warp.mrf = 1;
    % Clean Up
    matlabbatch{4}.spm.spatial.preproc.warp.cleanup = 1;
    % Warping Regularisation
    matlabbatch{4}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    % Affine Regularisation
    matlabbatch{4}.spm.spatial.preproc.warp.affreg = 'mni';
    % Smoothness
    matlabbatch{4}.spm.spatial.preproc.warp.fwhm = 0;
    % Sampling distance
    matlabbatch{4}.spm.spatial.preproc.warp.samp = 3;
    % Deformation Fields
    matlabbatch{4}.spm.spatial.preproc.warp.write = [0 1];
    
    %% 5. Normalize: Write
    
    % Deformation Field
    matlabbatch{5}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
    % Images to Write
    for ct_run=1:length(run_dir)
        matlabbatch{5}.spm.spatial.normalise.write.subj.resample(ct_run) = cfg_dep(['Realign & Unwarp: Unwarped Images (Sess ' num2str(ct_run) ')'], substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{ct_run}, '.','uwrfiles'));
    end
    matlabbatch{5}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Realign & Unwarp: Unwarped Images (Sess 1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','uwrfiles'));
    matlabbatch{5}.spm.spatial.normalise.write.subj.resample(2) = cfg_dep('Realign & Unwarp: Unwarped Images (Sess 2)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{2}, '.','uwrfiles'));
    matlabbatch{5}.spm.spatial.normalise.write.subj.resample(3) = cfg_dep('Realign & Unwarp: Unwarped Images (Sess 3)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{3}, '.','uwrfiles'));
    matlabbatch{5}.spm.spatial.normalise.write.subj.resample(4) = cfg_dep('Realign & Unwarp: Unwarped Images (Sess 4)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{4}, '.','uwrfiles'));
    matlabbatch{5}.spm.spatial.normalise.write.subj.resample(5) = cfg_dep('Realign & Unwarp: Unwarped Images (Sess 5)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{5}, '.','uwrfiles'));
    matlabbatch{5}.spm.spatial.normalise.write.subj.resample(6) = cfg_dep('Realign & Unwarp: Unwarped Images (Sess 6)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{6}, '.','uwrfiles'));
    matlabbatch{5}.spm.spatial.normalise.write.subj.resample(7) = cfg_dep('Realign & Unwarp: Unwarped Images (Sess 7)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{7}, '.','uwrfiles'));
    matlabbatch{5}.spm.spatial.normalise.write.subj.resample(8) = cfg_dep('Realign & Unwarp: Unwarped Images (Sess 8)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{8}, '.','uwrfiles'));
    
    % Writing Options
    % Bounding box
    matlabbatch{5}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
        78 76 85];
    % Voxel sizes
    matlabbatch{5}.spm.spatial.normalise.write.woptions.vox = [3 3 3.3];
    % Interpolation
    matlabbatch{5}.spm.spatial.normalise.write.woptions.interp = 4;
    % Filename Prefix
    matlabbatch{5}.spm.spatial.normalise.write.woptions.prefix = 'w';
    
    %% 6. Smooth
    % Images to smooth
    matlabbatch{6}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
    % FWHM
    matlabbatch{6}.spm.spatial.smooth.fwhm = [6 6 6];
    % Data Type
    matlabbatch{6}.spm.spatial.smooth.dtype = 0;
    % Implicit masking
    matlabbatch{6}.spm.spatial.smooth.im = 0;
    % Filename prefix
    matlabbatch{6}.spm.spatial.smooth.prefix = 's';
    
    %% Saving directory/batch script
    
    save(['prep_batch_' num2str(sbj_run_list{1,ct_sub}) '.mat'], 'matlabbatch');
    
    cd(dir_root);
    
    spm_jobman('run', matlabbatch);
    
end

Batch_End_Time = datestr(now)

end
