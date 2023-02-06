% DeltaT task
% fMRI experiment
% BHV Analysis: Calculate correlation & partial correlation among TD conditions(TD0, TD1, TD2-20)
% Ji Sun Kim

%% Directory Set-up
clear all; clc;

% Set-up directory for:
% sbj_dir
% sbj_dir_rt
% ori_dir
% save_dir

% Define subjects by loading sbj_list file

%% Stats: Correlation & Partial correlation

% Accuracy
load([save_dir '\TD_0_acc.mat']);
load([save_dir '\TD_1_acc.mat']);
load([save_dir '\TD_seq_acc.mat']);
pre_data = [TD_0_acc' TD_1_acc' TD_seq_acc'];

load([save_dir '\mean_TD_0_rt.mat']);
load([save_dir '\mean_TD_1_rt.mat']);
load([save_dir '\mean_TD_seq_rt.mat']);
pre_data = [mean_TD_0_rt' mean_TD_1_rt' mean_TD_seq_rt'];

[rho_corr,pval_corr] = corr(pre_data);
[rho_partial,pval_partial] = partialcorr(pre_data);

fprintf(['Correlation: TD0 & TD1: r = %0.3f , p = %0.3f \n'], rho_corr(2,1), pval_corr(2,1));
fprintf(['Correlation: TD0 & TD2~20: r = %0.3f , p = %0.3f \n'], rho_corr(3,1), pval_corr(3,1));
fprintf(['Correlation: TD1 & TD2~20: r = %0.3f , p = %0.3f \n'], rho_corr(2,3), pval_corr(2,3));

fprintf(['Partial Correlation: TD0 & TD1: r = %0.3f , p = %0.3f \n'], rho_partial(2,1), pval_partial(2,1));
fprintf(['Partial Correlation: TD0 & TD2~20: r = %0.3f , p = %0.3f \n'], rho_partial(3,1), pval_partial(3,1));
fprintf(['Partial Correlation: TD1 & TD2~20: r = %0.3f , p = %0.3f \n'], rho_partial(2,3), pval_partial(2,3));


%% Plot results
% Upper triangle: correlation, lower triangle: partial correlation

%% Accuracy

% load accuracy
load([save_dir '\TD_0_acc.mat']);
load([save_dir '\TD_1_acc.mat']);
load([save_dir '\TD_seq_acc.mat']);

pre_data = [TD_0_acc' TD_1_acc' TD_seq_acc'];

[rho_partial,pval_partial] = partialcorr(pre_data);
[rho_corr,pval_corr] = corr(pre_data);
rho = [1 rho_partial(1,2) rho_partial(1,3); ...
    rho_corr(2,1) 1 rho_partial(2,3); ...
    rho_corr(3,1) rho_corr(3,2) 1];
pval = [1 pval_partial(1,2) pval_partial(1,3); ...
    pval_corr(2,1) 1 pval_partial(2,3); ...
    pval_corr(3,1) pval_corr(3,2) 1];

figure;
h = imagesc(rho);
h.AlphaData = 0.7;
colormap(flipud(hot));
set(gca,'xtick',[]);
set(gca,'ytick',[]);
c = colorbar();
drawnow
alphaVal = 0.7;
cdata = c.Face.Texture.CData;
% Change the 4th channel (alpha channel) to 10% of it's initial value (255)
cdata(end,:) = uint8(alphaVal * cdata(end,:));
% Ensure that the display respects the alpha channel
c.Face.Texture.ColorType = 'truecoloralpha';
% Update the color data with the new transparency information
c.Face.Texture.CData = cdata;
% c.FaceColor.FaceAlpha = 0.7;
c.LineWidth = 1.3;
c.FontName = 'arial'; c.FontWeight = 'bold'; c.FontSize = 10;


%% Response Time

load([save_dir '\mean_TD_0_rt.mat']);
load([save_dir '\mean_TD_1_rt.mat']);
load([save_dir '\mean_TD_seq_rt.mat']);

pre_data = [mean_TD_0_rt' mean_TD_1_rt' mean_TD_seq_rt'];

% rest is the same


