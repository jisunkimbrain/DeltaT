% DeltaT task
% fMRI experiment
% Condition generator
% Ji Sun Kim

function [Parm] = fMRI_T0_CondGeneration(subject)

% Subject number
Parm.subj = subject;

% Parameter set-up
Parm.nsess = 2;
Parm.nrun = 4;
Parm.ntrial_Enc = 28; % 42 stims
Parm.ntrial_Del = 10;
Parm.ntrial_Ret = 20;

% Encoding stimuli duration
Parm.Stimdur_Enc = 2;

% Encoding ISI
ISI_Enc = 3.5;
Parm.ISI_Enc = ISI_Enc;

% Retrieval stimuli duration
Parm.Stimdur_Ret = 4;

% Retrieval ISI (jittering)
ISI_Ret_init = num2cell([ones(1,10), ones(1,10)*2]);
for sess=1:Parm.nsess
    for run=1:Parm.nrun
        ISI_Ret{sess,run} = ISI_Ret_init(randperm(20));
    end
end
Parm.ISI_Ret = ISI_Ret;

% TR (Enc, Ret)
for sess = 1:Parm.nsess
    for run=1:Parm.nrun
        for trial = 1:Parm.ntrial_Enc
            if trial==1
                TR_Enc{sess,run}{trial} = 1;
            else
                TR_Enc{sess,run}{trial} = TR_Enc{sess,run}{trial-1} + (Parm.ISI_Enc + Parm.Stimdur_Enc +0.5)/2;
            end
        end
    end
end
Parm.TR_Enc = TR_Enc;

for sess = 1:Parm.nsess
    for run=1:Parm.nrun
        for trial = 1:Parm.ntrial_Ret
            if trial==1
                TR_Ret{sess,run}{trial} = 1;
            else
                TR_Ret{sess,run}{trial} = TR_Ret{sess,run}{trial-1} + Parm.ISI_Ret{sess,run}{trial} + 2;
            end
        end
    end
end
Parm.TR_Ret = TR_Ret;

% Encoding Condition: Location (1; up-down / 2; down-up)
EncCond_location = {1 1 1 1 2 2 2 2};
Parm.EncCond_location = reshape(EncCond_location(randperm(length(EncCond_location))), 2,4);

% Encoding Condition: Item Number (1; one item left / 2; two items / 3; one item right)
% Retrieval Condition: Type of Question (1; T=1 / 2; Assocation T=0 / 3; TOJ)

% Pre-determined Sets
% Set 1
Set1_EncCond_item = {1 2 3 2 2 2 2 3 3 2 1 2 2 2 1 2 2 3 1 1 2 3 2 3 1 1 2 1};
Set1_RetCond_Q = {3 3 2 2 3 3 1 1 1 3 3 3 1 3 2 3 2 1 3 2};
% Set 2
Set2_EncCond_item = {1 3 2 3 3 2 2 2 3 2 1 2 2 1 3 2 2 3 1 1 2 1 1 2 2 2 2 3};
Set2_RetCond_Q = {2 1 3 2 2 3 3 1 3 1 3 1 3 1 3 3 2 3 2 3};
% Set 3
Set3_EncCond_item = {3 3 3 2 2 2 1 2 1 2 2 1 2 3 2 3 2 1 2 1 3 2 3 2 1 2 2 1};
Set3_RetCond_Q = {3 3 1 1 2 3 3 3 2 2 1 3 3 3 1 3 2 1 3 2};
% Set 4
Set4_EncCond_item = {3 3 2 1 2 2 2 3 3 1 2 2 1 3 2 1 2 2 2 1 1 3 2 2 2 2 3 3};
Set4_RetCond_Q = {2 3 1 3 3 1 2 3 3 3 1 1 3 2 2 1 3 3 3 2};
% Set 5
Set5_EncCond_item = {3 1 1 1 2 2 3 3 2 2 3 2 2 2 2 1 2 2 3 3 1 2 2 2 1 3 2 3};
Set5_RetCond_Q = {1 3 1 1 3 3 3 1 3 2 2 1 3 3 3 3 2 2 2 3};
% Set 6
Set6_EncCond_item = {3 2 2 2 1 3 2 2 3 3 1 3 2 2 2 2 2 2 2 1 3 1 2 1 2 3 1 1};
Set6_RetCond_Q = {2 2 2 2 1 3 3 3 1 3 3 3 3 3 1 2 3 1 1 3};
% Set 7
Set7_EncCond_item = {1 2 1 2 1 2 2 1 2 1 2 3 2 1 3 2 2 2 3 2 3 2 1 2 3 2 3 3};
Set7_RetCond_Q = {3 2 2 3 1 1 3 3 1 3 1 3 3 3 3 1 2 2 3 2};
% Set 8
Set8_EncCond_item = {1 3 3 1 2 2 2 2 2 1 2 2 3 2 3 2 2 1 2 1 2 3 1 2 3 1 2 1};
Set8_RetCond_Q = {3 1 3 3 3 2 3 3 2 1 3 3 1 3 1 2 3 1 2 2};
% Set 9
Set9_EncCond_item = {1 3 2 3 2 3 3 2 2 2 3 2 1 1 2 2 2 1 1 2 2 2 2 2 1 3 1 3};
Set9_RetCond_Q = {2 3 1 1 1 2 3 1 3 3 2 3 3 3 2 3 3 3 1 2};

% Set Condition
Set_pre_pre = {1 2 3 4 5 6 7 8 9};
Set_pre = Set_pre_pre(randperm(length(Set_pre_pre)));
Parm.SetCond = reshape(Set_pre(1:8),2,4);
for i=1:8
    randomize_trialorder{i} = randperm(20);
end

% Encoding & Retrieval set up
for sess=1:Parm.nsess
    for run=1:Parm.nrun
        if Parm.SetCond{sess,run} == 1
            Parm.EncCond_item{sess,run} = Set1_EncCond_item(randomize_trialorder{run});
            Parm.RetCond_Q{sess,run} =  Set1_RetCond_Q(randomize_trialorder{run});
        elseif Parm.SetCond{sess,run} == 2
            Parm.EncCond_item{sess,run} = Set2_EncCond_item(randomize_trialorder{run});
            Parm.RetCond_Q{sess,run} =  Set2_RetCond_Q(randomize_trialorder{run});
        elseif Parm.SetCond{sess,run} == 3
            Parm.EncCond_item{sess,run} = Set3_EncCond_item(randomize_trialorder{run});
            Parm.RetCond_Q{sess,run} =  Set3_RetCond_Q(randomize_trialorder{run});
        elseif Parm.SetCond{sess,run} == 4
            Parm.EncCond_item{sess,run} = Set4_EncCond_item(randomize_trialorder{run});
            Parm.RetCond_Q{sess,run} =  Set4_RetCond_Q(randomize_trialorder{run});
        elseif Parm.SetCond{sess,run} == 5
            Parm.EncCond_item{sess,run} = Set5_EncCond_item(randomize_trialorder{run});
            Parm.RetCond_Q{sess,run} =  Set5_RetCond_Q(randomize_trialorder{run});
        elseif Parm.SetCond{sess,run} == 6
            Parm.EncCond_item{sess,run} = Set6_EncCond_item(randomize_trialorder{run});
            Parm.RetCond_Q{sess,run} =  Set6_RetCond_Q(randomize_trialorder{run});
        elseif Parm.SetCond{sess,run} == 7
            Parm.EncCond_item{sess,run} = Set7_EncCond_item(randomize_trialorder{run});
            Parm.RetCond_Q{sess,run} =  Set7_RetCond_Q(randomize_trialorder{run});
        elseif Parm.SetCond{sess,run} == 8
            Parm.EncCond_item{sess,run} = Set8_EncCond_item(randomize_trialorder{run});
            Parm.RetCond_Q{sess,run} =  Set8_RetCond_Q(randomize_trialorder{run});
        elseif Parm.SetCond{sess,run} == 9
            Parm.EncCond_item{sess,run} = Set9_EncCond_item(randomize_trialorder{run});
            Parm.RetCond_Q{sess,run} =  Set9_RetCond_Q(randomize_trialorder{run});
        elseif Parm.SetCond{sess,run} == 10
            Parm.EncCond_item{sess,run} = Set10_EncCond_item(randomize_trialorder{run});
            Parm.RetCond_Q{sess,run} =  Set10_RetCond_Q(randomize_trialorder{run});
        end
    end
end

end
