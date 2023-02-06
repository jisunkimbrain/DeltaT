% DeltaT task
% fMRI experiment
% ROI Univariate analyses
% Ji Sun Kim

%% Define pathway

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

%% How many voxels?

ROI_voxel_number = [];

for ct_sub=1:length(sbj_run_list)
    for ct_roi=1:length(mask_name)

        ROI = [roi_dir '\' num2str(sbj_run_list{1,ct_sub}) '\' mask_name{ct_roi}];
        Y = spm_read_vols(spm_vol(ROI),1);
        indx = find(Y>0);
        ROI_voxel_number(ct_sub,ct_roi) = length(indx);

    end  %-- end of for ct_roi
end  %-- end of for ct_sub

for ct_sub=1:length(sbj_run_list)
    for ct_roi=19:length(roi_name)

        a = ct_roi-18;
        b = 1 + 2*(a-1);

        ROI_voxel_number(ct_sub,ct_roi) = ROI_voxel_number(ct_sub,b) + ROI_voxel_number(ct_sub,b+1);

    end  %-- end of for ct_roi
end  %-- end of for ct_sub

% calculate mean & sd
for ct_roi=1:length(roi_name)

    roi_mean = mean(ROI_voxel_number(:,ct_roi));
    roi_std = std(ROI_voxel_number(:,ct_roi));
    
    fprintf([roi_name{ct_roi} ': mean ' num2str(roi_mean) ' std ' num2str(roi_std) '\n']);

end

%% Extract mean parameter estimate

clear TD_meanbeta_0_1to20

for ct_sub=1:length(sbj_run_list)

    sbj_first_dir =[first_level_dir '\' num2str(sbj_run_list{1,ct_sub})];

    for ct_roi=1:length(roi_name)

        if ct_roi < 19

            roi_file = [roi_dir '\' num2str(sbj_run_list{1,ct_sub}) '\' mask_name{ct_roi}];

            for td_ct=1:2
                a = sprintf('%04d', td_ct);
                beta_file = [sbj_first_dir '\con_' a '.nii'];
                TD_meanbeta_0_1to20{ct_roi}(ct_sub, td_ct) = Extract_ROI_Data(roi_file, beta_file);
            end  %-- end of for td_ct

        elseif ct_roi == 19

            L_roi = [roi_dir '\' num2str(sbj_run_list{1,ct_sub}) '\' mask_name{1}];
            R_roi = [roi_dir '\' num2str(sbj_run_list{1,ct_sub}) '\' mask_name{2}];
            XYZ = Combine_Bilateral_ROI_coordinates(L_roi, R_roi);

            for td_ct=1:2
                a = sprintf('%04d', td_ct);
                beta_file = [sbj_first_dir '\con_' a '.nii'];
                TD_meanbeta_0_1to20{ct_roi}(ct_sub, td_ct) = nanmean(spm_get_data(beta_file, XYZ),2);
            end


        elseif ct_roi == 20

            L_roi = [roi_dir '\' num2str(sbj_run_list{1,ct_sub}) '\' mask_name{3}];
            R_roi = [roi_dir '\' num2str(sbj_run_list{1,ct_sub}) '\' mask_name{4}];
            XYZ = Combine_Bilateral_ROI_coordinates(L_roi, R_roi);

            for td_ct=1:2
                a = sprintf('%04d', td_ct);
                beta_file = [sbj_first_dir '\con_' a '.nii'];
                TD_meanbeta_0_1to20{ct_roi}(ct_sub, td_ct) = nanmean(spm_get_data(beta_file, XYZ),2);
            end

        elseif ct_roi == 21

            L_roi = [roi_dir '\' num2str(sbj_run_list{1,ct_sub}) '\' mask_name{5}];
            R_roi = [roi_dir '\' num2str(sbj_run_list{1,ct_sub}) '\' mask_name{6}];
            XYZ = Combine_Bilateral_ROI_coordinates(L_roi, R_roi);

            for td_ct=1:2
                a = sprintf('%04d', td_ct);
                beta_file = [sbj_first_dir '\con_' a '.nii'];
                TD_meanbeta_0_1to20{ct_roi}(ct_sub, td_ct) = nanmean(spm_get_data(beta_file, XYZ),2);
            end

        elseif ct_roi==22

            L_roi = [roi_dir '\' num2str(sbj_run_list{1,ct_sub}) '\' mask_name{7}];
            R_roi = [roi_dir '\' num2str(sbj_run_list{1,ct_sub}) '\' mask_name{8}];
            XYZ = Combine_Bilateral_ROI_coordinates(L_roi, R_roi);

            for td_ct=1:2
                a = sprintf('%04d', td_ct);
                beta_file = [sbj_first_dir '\con_' a '.nii'];
                TD_meanbeta_0_1to20{ct_roi}(ct_sub, td_ct) = nanmean(spm_get_data(beta_file, XYZ),2);
            end

        elseif ct_roi==23

            L_roi = [roi_dir '\' num2str(sbj_run_list{1,ct_sub}) '\' mask_name{9}];
            R_roi = [roi_dir '\' num2str(sbj_run_list{1,ct_sub}) '\' mask_name{10}];
            XYZ = Combine_Bilateral_ROI_coordinates(L_roi, R_roi);

            for td_ct=1:2
                a = sprintf('%04d', td_ct);
                beta_file = [sbj_first_dir '\con_' a '.nii'];
                TD_meanbeta_0_1to20{ct_roi}(ct_sub, td_ct) = nanmean(spm_get_data(beta_file, XYZ),2);
            end

        elseif ct_roi==24

            L_roi = [roi_dir '\' num2str(sbj_run_list{1,ct_sub}) '\' mask_name{11}];
            R_roi = [roi_dir '\' num2str(sbj_run_list{1,ct_sub}) '\' mask_name{12}];
            XYZ = Combine_Bilateral_ROI_coordinates(L_roi, R_roi);

            for td_ct=1:2
                a = sprintf('%04d', td_ct);
                beta_file = [sbj_first_dir '\con_' a '.nii'];
                TD_meanbeta_0_1to20{ct_roi}(ct_sub, td_ct) = nanmean(spm_get_data(beta_file, XYZ),2);
            end

        elseif ct_roi==25

            L_roi = [roi_dir '\' num2str(sbj_run_list{1,ct_sub}) '\' mask_name{13}];
            R_roi = [roi_dir '\' num2str(sbj_run_list{1,ct_sub}) '\' mask_name{14}];
            XYZ = Combine_Bilateral_ROI_coordinates(L_roi, R_roi);

            for td_ct=1:2
                a = sprintf('%04d', td_ct);
                beta_file = [sbj_first_dir '\con_' a '.nii'];
                TD_meanbeta_0_1to20{ct_roi}(ct_sub, td_ct) = nanmean(spm_get_data(beta_file, XYZ),2);
            end

        elseif ct_roi==26

            L_roi = [roi_dir '\' num2str(sbj_run_list{1,ct_sub}) '\' mask_name{15}];
            R_roi = [roi_dir '\' num2str(sbj_run_list{1,ct_sub}) '\' mask_name{16}];
            XYZ = Combine_Bilateral_ROI_coordinates(L_roi, R_roi);

            for td_ct=1:2
                a = sprintf('%04d', td_ct);
                beta_file = [sbj_first_dir '\con_' a '.nii'];
                TD_meanbeta_0_1to20{ct_roi}(ct_sub, td_ct) = nanmean(spm_get_data(beta_file, XYZ),2);
            end

        elseif ct_roi==27

            L_roi = [roi_dir '\' num2str(sbj_run_list{1,ct_sub}) '\' mask_name{17}];
            R_roi = [roi_dir '\' num2str(sbj_run_list{1,ct_sub}) '\' mask_name{18}];
            XYZ = Combine_Bilateral_ROI_coordinates(L_roi, R_roi);

            for td_ct=1:2
                a = sprintf('%04d', td_ct);
                beta_file = [sbj_first_dir '\con_' a '.nii'];
                TD_meanbeta_0_1to20{ct_roi}(ct_sub, td_ct) = nanmean(spm_get_data(beta_file, XYZ),2);
            end
        
        end


    end  %-- end of for ct_roi

    fprintf(['Done for sbj: ' num2str(sbj_run_list{1,ct_sub}) '\n']);

end  %-- end of for ct_sub

%% Statistical testing

for ct_roi=1:length(roi_name)

    pre_data = TD_meanbeta_0_1to20{1,ct_roi};

    [h1,p1,ci1,stats1] = ttest(pre_data(:,1), pre_data(:,2));
    fprintf([roi_name{ct_roi} '(' num2str(ct_roi) ') : significant difference of activation between Corr 0 & Corr 1to20: t(23): %0.3f , p = %0.3f \n'], stats1.tstat, p1);


end  %-- end of for ct_roi

%% Plot figures

for ct_roi=1:length(roi_name)

    if ismember(ct_roi,HPC_roi_indx) == 1

    pre_data = TD_meanbeta_0_1to20{1,ct_roi};

    % plot results
    data = {pre_data(:,1) pre_data(:,2)};
    [avg,err] = jh_mean_err(2,data);

    labels = categorical({'TD0', 'TDseq'});
    labels = reordercats(labels,string(labels));
    y = avg;
    figure;
    resultsbar = bar(labels,y,'FaceColor','flat', 'BarWidth', 0.7);
    hold on
    set(gca,'FontSize',20, 'linewidth',2, 'FontName', 'Arial', 'FontWeight', 'bold', 'box','off');

    ax = gca;
    ax.XAxis.FontSize=20;
    ax.XAxis.FontWeight='bold';
    ax.XAxis.FontName = 'arial';
    ax = gca;
    
    resultsbar.CData(1,:) = [1 0.8 0.8];
    resultsbar.CData(2,:) = [0.8 0.9 1];
    hold on
    bar(y(1:2), 'FaceAlpha',0,'EdgeColor',[0.3 0.5 0.6],'LineWidth',4, 'BarWidth', 0.7);
    bar(y(1),'FaceAlpha',0,'EdgeColor',[0.8 0.1 0],'LineWidth',4, 'BarWidth', 0.7);

    ylim([-0.02 0.3]);

    hold on
    yline(0,'k','linewidth',2);

    for i=1:2
        hold on
        plot([i i],[avg(i)-err(i) avg(i)+err(i)], 'Color', 'k', 'linewidth', 2.5);
    end

    else
        continue
    end

end  %-- end of for ct_roi

               
