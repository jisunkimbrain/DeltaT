% DeltaT task
% fMRI experiment
% Initializer
% Ji Sun Kim

clear
clc

SID = input('Subjects Number: ');

[Parm] = fMRI_T0_CondGeneration(SID);
[Parm] = fMRI_T0_GenerateStimSets(SID, Parm);

if ~exist('Exp_Parm')
    mkdir('Exp_Parm')
end

movefile('Parm_*.mat', 'Exp_Parm');

if ~exist('BHV_Results_ENC')
    mkdir('BHV_Results_ENC\xls')
    mkdir('BHV_Results_ENC\mat')
    mkdir('BHV_Results_RET\xls')
    mkdir('BHV_Results_RET\mat')
    mkdir('BHV_Results_Del\xls')
    mkdir('BHV_Results_DEL\mat')
end
