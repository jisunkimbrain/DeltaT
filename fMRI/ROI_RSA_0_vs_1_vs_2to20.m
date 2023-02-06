% DeltaT task
% fMRI experiment
% ROI Representational Similarity Analysis
% Ji Sun Kim

%% Basic Set-up

clear all; clc

% Define subjects by loading sbj_run_list.mat

% Set-up directory for:
% ori_dir
% results_dir
% roi_dir
% first_level_dir

mask_name = {'mri\raparc_aseg_lh_hippocampus.nii', 'mri\raparc_aseg_rh_hippocampus.nii', ...
    'mri\rlh-anterior-hippocampus.nii', 'mri\rrh-anterior-hippocampus.nii', ...
    'mri\rlh-posterior-hippocampus.nii', 'mri\rrh-posterior-hippocampus.nii', ...
    'mri\rlh_hp_ca1.nii', 'mri\rrh_hp_ca1.nii', ...
    'mri\rlh_hp_ca23dg.nii', 'mri\rrh_hp_ca23dg.nii', ...
    'mri\rlh_hp_sub.nii', 'mri\rrh_hp_sub.nii', ...
    'mri\rlh_hp_wholesub.nii', 'mri\rrh_hp_wholesub.nii', ...
    'mri\raparc_aseg_lh_entorhinal.nii', 'mri\raparc_aseg_rh_entorhinal.nii', ...
    'mri\raparc_aseg_lh_parahippocampal.nii', 'mri\raparc_aseg_rh_parahippocampal.nii'};
roi_name = {'Left Hippocampus', 'Right Hippocampus', ...
    'Left aHippocampus', 'Right aHippocampus', ...
    'Left pHippocampus', 'Right pHippocampus', ...
    'Left CA1', 'Right CA1', ...
    'Left CA23DG', 'Right CA23DG', ...
    'Left Subiculum', 'Right Subiculum', ...
    'Left Subiculum (whole)', 'Right Subiculum (whole)', ...
    'Left Entorhinal cortex', 'Right Entorhinal cortex', ...
    'Left Parahippocampal gyrus', 'Right Parahippocampal gyrus', ...
    'Bilateral Hippocampus', 'Bilateral aHippocampus', 'Bilateral pHippocampus', ...
    'Biltaeral CA1', 'Bilateral CA23DG', 'Bilateral Subiculum', 'Bilateral Subiculum (whole)' ...
    'Bilateral Entorhinal cortex', 'Bilateral Parahippocampal gyrus'};

conditions = {'TD0', 'TD1', 'TD2to20'};

%% Set up ROI: Read data for ROI masks

for ct_sub=1:length(sbj_run_list)

    for ct_roi=1:length(roi_name)

        if ct_roi < 19

            roi_file = [roi_dir '\' num2str(sbj_run_list{1,ct_sub}) '\' mask_name{ct_roi}];
            Y = spm_read_vols(spm_vol(roi_file),1);
            indx = find(Y>0);
            [x,y,z] = ind2sub(size(Y),indx);
            XYZ = [x y z]';

        elseif ct_roi == 19

            L_roi = [roi_dir '\' num2str(sbj_run_list{1,ct_sub}) '\' mask_name{1}];
            R_roi = [roi_dir '\' num2str(sbj_run_list{1,ct_sub}) '\' mask_name{2}];
            XYZ = Combine_Bilateral_ROI_coordinates(L_roi, R_roi);

        elseif ct_roi == 20

            L_roi = [roi_dir '\' num2str(sbj_run_list{1,ct_sub}) '\' mask_name{3}];
            R_roi = [roi_dir '\' num2str(sbj_run_list{1,ct_sub}) '\' mask_name{4}];
            XYZ = Combine_Bilateral_ROI_coordinates(L_roi, R_roi);

        elseif ct_roi == 21

            L_roi = [roi_dir '\' num2str(sbj_run_list{1,ct_sub}) '\' mask_name{5}];
            R_roi = [roi_dir '\' num2str(sbj_run_list{1,ct_sub}) '\' mask_name{6}];
            XYZ = Combine_Bilateral_ROI_coordinates(L_roi, R_roi);

        elseif ct_roi==22

            L_roi = [roi_dir '\' num2str(sbj_run_list{1,ct_sub}) '\' mask_name{7}];
            R_roi = [roi_dir '\' num2str(sbj_run_list{1,ct_sub}) '\' mask_name{8}];
            XYZ = Combine_Bilateral_ROI_coordinates(L_roi, R_roi);

        elseif ct_roi==23

            L_roi = [roi_dir '\' num2str(sbj_run_list{1,ct_sub}) '\' mask_name{9}];
            R_roi = [roi_dir '\' num2str(sbj_run_list{1,ct_sub}) '\' mask_name{10}];
            XYZ = Combine_Bilateral_ROI_coordinates(L_roi, R_roi);

        elseif ct_roi==24

            L_roi = [roi_dir '\' num2str(sbj_run_list{1,ct_sub}) '\' mask_name{11}];
            R_roi = [roi_dir '\' num2str(sbj_run_list{1,ct_sub}) '\' mask_name{12}];
            XYZ = Combine_Bilateral_ROI_coordinates(L_roi, R_roi);

        elseif ct_roi==25

            L_roi = [roi_dir '\' num2str(sbj_run_list{1,ct_sub}) '\' mask_name{13}];
            R_roi = [roi_dir '\' num2str(sbj_run_list{1,ct_sub}) '\' mask_name{14}];
            XYZ = Combine_Bilateral_ROI_coordinates(L_roi, R_roi);

        elseif ct_roi==26

            L_roi = [roi_dir '\' num2str(sbj_run_list{1,ct_sub}) '\' mask_name{15}];
            R_roi = [roi_dir '\' num2str(sbj_run_list{1,ct_sub}) '\' mask_name{16}];
            XYZ = Combine_Bilateral_ROI_coordinates(L_roi, R_roi);

        elseif ct_roi==27

            L_roi = [roi_dir '\' num2str(sbj_run_list{1,ct_sub}) '\' mask_name{17}];
            R_roi = [roi_dir '\' num2str(sbj_run_list{1,ct_sub}) '\' mask_name{18}];
            XYZ = Combine_Bilateral_ROI_coordinates(L_roi, R_roi);

        end
        roi_xyz{ct_sub,ct_roi} = XYZ;

    end

end  %-- end of for roi_name

%% Step 1. Computing and visualizing RDMs


%% Compute RDM for each region and for each individual

% Extract activation pattern for each condition, each region, each individual
for ct_sub=1:length(sbj_run_list)

    % load functional data of each condition of each subject
    sbj_first_dir = [first_level_dir '\' num2str(sbj_run_list{1,ct_sub})];
    TD_0_img = [sbj_first_dir '\con_0001.nii'];
    TD_1_img = [sbj_first_dir '\con_0002.nii'];
    TD_2to20_img = [sbj_first_dir '\con_0003.nii'];

    for ct_roi=1:length(roi_name)

        % TD=0
        cond_pattern{ct_sub,ct_roi}(:,1) = spm_get_data(TD_0_img,roi_xyz{ct_sub,ct_roi});
        % TD=1
        cond_pattern{ct_sub,ct_roi}(:,2) = spm_get_data(TD_1_img,roi_xyz{ct_sub,ct_roi});
        % TD=2~20
        cond_pattern{ct_sub,ct_roi}(:,3) = spm_get_data(TD_2to20_img,roi_xyz{ct_sub,ct_roi});

    end  %-- end of for ct_roi

end  %-- end of for ct_sub
save([results_dir '\cond_pattern.mat'], 'cond_pattern');


% Compute RDM for each region and each individual
for roi_ct=1:length(roi_name)

    for sbj_ct=1:length(sbj_run_list)

        for cond_i=1:numel(conditions)
            for cond_j=1:numel(conditions)
                if cond_i ~= cond_j
                roi{roi_ct}.rdms(cond_i,cond_j,sbj_ct) = 1 - corr( cond_pattern{sbj_ct,roi_ct}(:,cond_i), ...
                    cond_pattern{sbj_ct,roi_ct}(:,cond_j), 'Rows', 'complete');

                % to make RDMs symmetric
                roi{roi_ct}.rdms(cond_j,cond_i,sbj_i) = roi{roi_i}.rdms(cond_i,cond_j,sbj_i);

                end  %-- construct RDM
            end
        end  %-- end of for cond_i, cond_j

    end  %-- end of for sbj_ct

%     % print progress
%     fprintf(['Calculation done for ROI(' num2str(roi_ct) '): ' roi_name{roi_ct} '\n']);

end  %-- end of for roi_ct

save([results_dir '\roi.mat'], 'roi');

% % print progress
% fprintf(['Calculation done for subject(' num2str(sbj_ct) '): ' num2str(sbj_run_list{1,sbj_ct}) '\n']);


%% Visualize mean RDMs

my_colormap = customcolormap([0 0.4 1], {'#FA7A22', '#F5EDCB', '#9CABA2'});
colorbar;
colormap(my_colormap);

figure;
set(gcf, 'position', [1         800        1500         250], 'color', 'w')

for roi_ct = 1:numel(roi_foranalysis)
    subplot(1,4,roi_ct)
    imagesc(mean(roi{roi_foranalysis(roi_ct)}.rdms,3));   
    caxis([0.5 0.7]);
    set(gca, 'xtick', 1:3, 'ytick', 1:3, 'LineWidth', 2);
    set(gca, 'xticklabel', conditions, 'XTickLabelRotation', 60, ...
        'yticklabel', conditions);
    set(gca, 'FontWeight', 'bold', 'FontSize', 15);
    colorbar;
end

% Create model RDMs

% 1. TD=1 closer to TD=2~20 (sequence)
model{1}.rdms = [0 1 1; 1 0 0; 1 0 0];
model{1}.name = 'Seq';
                     
% 2. TD=1 closer to TD=0 (simultaneous association)
model{2}.rdms = [0 0 1; 0 0 1; 1 1 0];
model{2}.name = 'Prox';

% 3. TD=2~20 closer to TD=0 (conditions with similar accuracy)
model{3}.rdms = [0 1 0; 1 0 1; 0 1 0];
model{3}.name = 'Acc';

figure;
set(gcf, 'position', [1         755        1050         250], 'color', 'w')

for i = 1:numel(model)
    if i == 1, disp('Model Representational Dissimilarity Matrices'); end
    subplot(1,3,i);
    imagesc(model{i}.rdms);
    colorbar;
    caxis([0 1]);
    colormap(flipud(pink))
    title(model{i}.name, 'FontWeight', 'bold', 'FontSize', 15);
    set(gca, 'xtick', 1:3, 'ytick', 1:3, 'LineWidth', 2);
    set(gca, 'xticklabel', conditions, 'XTickLabelRotation', 90, ...
        'yticklabel', conditions);
end

%% Compare the ROI neural RDMs with model RDMs

% ROIs
for roi_i = 1:numel(roi)
    rdms.dat(:,:,roi_i) = mean(roi{roi_i}.rdms,3); % average across subjects
    rdms.names{roi_i} = roi_name{roi_i};           % ROI names
end

% models
for model_i = 1:numel(model)
    rdms.dat(:,:,model_i+numel(roi)) = model{model_i}.rdms;
    rdms.names{model_i+numel(roi)} = model{model_i}.name;
end

% Preallocation
r_models = NaN(size(rdms.dat,3),size(rdms.dat,3));
upper_triang_idx = triu(true(3,3),1);

% Calculate Kendall's tau a (recommendation by Kriegeskorte)
for i = 1:size(rdms.dat,3)
    xx = squeeze(rdms.dat(:,:,i));
    for j = 1:size(rdms.dat,3)
        if i ~= j
            yy = squeeze(rdms.dat(:,:,j));
            r_models(i,j) = rankCorr_Kendall_taua(xx(upper_triang_idx), ...
                yy(upper_triang_idx));
        else
            r_models(i,j) = 1;
        end
    end
    rdms.flatten_dat(:,i) = xx(upper_triang_idx);
end


%% Step 3. Statistical Inference

%% Prep for bootstrap tests

% flatten the roi matrix
for roi_i = 1:numel(roi)
    roi{roi_i}.rdms_flatten = [];
    for subj_i = 1:size(roi{roi_i}.rdms,3)
        a = roi{roi_i}.rdms(:,:,subj_i);
        roi{roi_i}.rdms_flatten(:,subj_i) = a(upper_triang_idx);
        % rdms.flatten_dat
    end
end


%% Run bootstrap tests

model_names = {'Seq', 'Prox', 'Acc'}; % use a short model name

for roi_i = 1:numel(roi)
    clear boot_vals;
    for model_i = 1:numel(model)
        boot_rdmsmean = bootstrp(10000, @mean, roi{roi_i}.rdms_flatten')';
        for iter_i = 1:size(boot_rdmsmean,2)
            boot_vals(iter_i,1) = rankCorr_Kendall_taua(...
                model{model_i}.rdms(upper_triang_idx), ...
                boot_rdmsmean(:,iter_i));
        end
        eval(['roi{roi_i}.bootmean_' model_names{model_i} ' = mean(boot_vals);']);
        eval(['roi{roi_i}.bootste_' model_names{model_i} ' = std(boot_vals);']);
        eval(['roi{roi_i}.bootZ_' model_names{model_i} ' = roi{roi_i}.bootmean_' ...
            model_names{model_i} './roi{roi_i}.bootste_' model_names{model_i} ';']);
        eval(['roi{roi_i}.bootP_' model_names{model_i} ...
            '= 2 * (1 - normcdf(abs(roi{roi_i}.bootZ_' ...
            model_names{model_i} ')));']);
        eval(['roi{roi_i}.ci95_' model_names{model_i} ...
            ' = [prctile(boot_vals, 2.5); prctile(boot_vals, 97.5)];']);
    end
    
    roi{roi_i}.noise_ceiling(1,1) = mean(corr(mean(roi{roi_i}.rdms_flatten,2), ...
        roi{roi_i}.rdms_flatten));

    for ii = 1:size(roi{roi_i}.rdms_flatten,2)
        temp = roi{roi_i}.rdms_flatten;
        temp(:,ii) = [];
        lower_bound_i(ii,1) = corr(mean(temp,2), roi{roi_i}.rdms_flatten(:,ii));
    end
    
    roi{roi_i}.noise_ceiling(1,2) = mean(lower_bound_i);
end

%% Plot results

for roi_i = 1:numel(roi_name)

    figure;
     
    y = [roi{roi_i}.bootmean_0_1 roi{roi_i}.bootmean_1_2to20 roi{roi_i}.bootmean_0_2to20];
    e = [roi{roi_i}.bootste_0_1 roi{roi_i}.bootste_1_2to20 roi{roi_i}.bootste_0_2to20];
    p = [roi{roi_i}.bootP_0_1 roi{roi_i}.bootP_1_2to20 roi{roi_i}.bootP_0_2to20];

    bar_wani_2016(y, e, .8, 'errbar_width', 0, 'ast', p, 'use_samefig', ...
        'ylim', [-.6 1]);
    patch(get(gca, 'xlim'), roi{roi_i}.noise_ceiling, [.7 .7 .7]);
    title(roi_name{roi_i}, 'FontSize', 15);
    ylabel('correlation (Kendall''s tau a)', 'FontSize', 15);
    xticklabels(model_name);
    ylim([-1 1]);

    saveas(gcf, [results_dir '\Stats\' roi_name{roi_i} '.bmp']);

    close all;

end
