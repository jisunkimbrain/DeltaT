# DeltaT

# Repository details
Hello!
This repository contains codes used to collect data & perform analyses for the following paper:

[Kim, J.S., Lee, S.A. (2022) Hippocampal orchestration of associative and sequential memory networks for episodic retrieval](https://drive.google.com/file/d/1eP4APrUKljEK1VQNjqVDwKhwZt8aEMkF/view?usp=sharing)

All scripts are written by J.S. Kim

Please contact me through kimjisun@snu.ac.kr if you have any questions about the codes

The scripts here are based on Matlab r2021b.

# Script Folders
## TASK
*These scripts were used to run the "DeltaT" Task. The task runs based on Cogent 2000 v125*

- `fMRI_DeltaT_CondGeneration.m`
- `fMRI_DeltaT_GenerateStimSets.m`
- `fMRI_DeltaT_TaskPlayer.m`
- `fMRI_DeltaT_Initializer.m`

## BHV

- `TDConditions_Comparison.m`
- `TDConditions_Correlation.m`
- `TD_Overall.m`

## fMRI
*Preprocessing & General Linear Model (GLM) analyses were run using SPM12*

### Whole-brain Univariate Analyses
- `Preprocessing.m`
- `GLM_Modeling_Corr_0_vs_1_vs_2to20.m`
- `GLM_Contrast_Corr_0_vs_1_vs_2to20.m`
- `GLM_RFX_Corr_0_vs_1_vs_2to20.m`

### ROI: Univariate Analyses

- `Extract_ROI_Data.m`
- `ROI_Univariate_Asso_vs_Seq.m`

### ROI: Multivariate Analyses

*RSA Analyses were performed using (and codes are based on) functions from [CanlabCore](https://www.mathworks.com/matlabcentral/fileexchange/72750-canlabcore?s_tid=FX_rc3_behav)*


*Colorbars were created using the function [Custom Colormap](https://www.mathworks.com/matlabcentral/fileexchange/69470-custom-colormap)*

- `ROI_RSA_0_vs_1_vs_2to20.m`

### gPPI Analyses
*gPPI analysis was performed by using the CONN GUI*
