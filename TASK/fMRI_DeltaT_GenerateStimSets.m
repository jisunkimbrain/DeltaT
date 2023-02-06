% DeltaT task
% fMRI experiment
% Stimuli list generator
% Ji Sun Kim


function [Parm] = fMRI_T0_GenerateStimSets(subject, Parm)

% Load picture list
ori_dir = cd;  %-- start from Temporal_difference_task folder
Stim_dir = [ori_dir '\Stimuli_set\JObjects\JObjects_jpgresized'];
dir_struct = dir(Stim_dir);

% Stimulation set
for i = 1:size(dir_struct,1)-2
    Stim_name{i,1} = dir_struct(i+2).name;
end

% Randomize picture list
random_order = randperm(length(Stim_name));
random_stim_all = Stim_name(random_order);
Parm.random_order = random_order;
Parm.random_stim_all = random_stim_all;

% Stimulation set
Parm.Enc_stim{1,1} = random_stim_all(1:42);
Parm.Enc_stim{1,2} = random_stim_all(43:84);
Parm.Enc_stim{1,3} = random_stim_all(85:126);
Parm.Enc_stim{1,4} = random_stim_all(253:294);
Parm.Enc_stim{2,1} = random_stim_all(127:168);
Parm.Enc_stim{2,2} = random_stim_all(169:210);
Parm.Enc_stim{2,3} = random_stim_all(211:252);
Parm.Enc_stim{2,4} = random_stim_all(295:336);

% Encoding: set up stimuli (left & right)
for sess=1:Parm.nsess
    for run=1:Parm.nrun
        i=1;
        for trial=1:Parm.ntrial_Enc
            if Parm.EncCond_item{sess,run}{trial} == 1
                Enc_stim_left{sess,run}{trial} = Parm.Enc_stim{sess,run}{i};
                Enc_stim_right{sess,run}{trial} = nan;
                i=i+1;
            elseif Parm.EncCond_item{sess,run}{trial} == 2
                Enc_stim_left{sess,run}{trial} = Parm.Enc_stim{sess,run}{i};
                Enc_stim_right{sess,run}{trial} = Parm.Enc_stim{sess,run}{i+1};
                i=i+2;
            elseif Parm.EncCond_item{sess,run}{trial} == 3
                Enc_stim_left{sess,run}{trial} = nan;
                Enc_stim_right{sess,run}{trial} = Parm.Enc_stim{sess,run}{i};
                i=i+1;
            end
        end
    end
end
Parm.Enc_stim_left = Enc_stim_left;
Parm.Enc_stim_right = Enc_stim_right;

% Retrieval: set up stimuli
for sess=1:Parm.nsess
    for run=1:Parm.nrun
        
        % Different combinations for different Sets
        if Parm.SetCond{sess,run} == 1
            Ret_stim_list_1 = {4 13 12 27 5 4 7 15 20 5 17 3 9 10 6 18 2 16 11 14};
            Ret_stim_list_2 = {21 24 12 27 7 22 8 16 21 25 23 13 10 26 6 23 2 17 19 14};
            Ret_stim_1_Enc_location = {1 1 1 1 2 2 1 1 1 1 1 2 2 2 1 2 1 2 1 1};
            Ret_stim_2_Enc_location = {1 2 2 2 2 2 2 1 2 1 1 2 1 1 2 2 2 2 1 2};
        elseif Parm.SetCond{sess,run} == 2
            Ret_stim_list_1 = {6 26 12 3 24 5 4 7 8 18 2 12 20 10 13 17 21 10 16 7};
            Ret_stim_list_2 = {6 27 15 3 24 14 23 8 22 19 9 13 26 11 17 25 21 25 16 27};
            Ret_stim_1_Enc_location = {1 1 2 1 1 2 2 1 2 2 2 1 1 2 2 1 1 1 1 2};
            Ret_stim_2_Enc_location = {2 2 2 2 2 1 1 1 1 1 2 1 2 1 2 1 2 2 2 1};
        elseif Parm.SetCond{sess,run} == 3
            Ret_stim_list_1 = {9 16 24 17 26 2 12 7 10 6 21 5 5 17 13 22 8 3 4 11};
            Ret_stim_list_2 = {15 19 25 18 26 20 24 27 10 6 22 15 19 23 14 27 8 4 13 11};
            Ret_stim_1_Enc_location = {1 2 1 1 1 2 1 1 1 1 2 2 1 2 2 2 1 2 2 1};
            Ret_stim_2_Enc_location = {2 2 1 1 2 1 2 2 2 2 1 1 1 2 2 1 2 1 1 2};
        elseif Parm.SetCond{sess,run} == 4
            Ret_stim_list_1 = {24 2 25 10 4 17 5 3 7 16 12 26 8 6 23 17 3 7 11 15};
            Ret_stim_list_2 = {24 19 26 12 9 18 5 19 22 20 13 27 21 6 23 18 14 11 25 15};
            Ret_stim_1_Enc_location = {1 2 2 1 1 2 1 2 1 1 1 2 2 1 1 1 1 2 1 1};
            Ret_stim_2_Enc_location = {2 1 1 2 2 2 2 2 2 1 1 2 1 2 2 1 2 2 1 2};
        elseif Parm.SetCond{sess,run} == 5
            Ret_stim_list_1 = {2 10 6 14 27 4 9 18 8 23 12 26 6 14 9 10 5 13 16 22};
            Ret_stim_list_2 = {3 21 7 15 17 20 24 19 15 23 12 27 18 17 22 16 5 13 11 24};
            Ret_stim_1_Enc_location = {1 1 1 2 1 1 1 1 2 1 1 2 2 1 2 2 1 1 1 1};
            Ret_stim_2_Enc_location = {1 1 2 2 2 2 1 2 1 2 2 2 2 1 2 1 2 2 2 2};
        elseif Parm.SetCond{sess,run} == 6
            Ret_stim_list_1 = {2 19 3 23 7 25 25 5 26 10 11 14 7 4 4 16 13 14 17 8};
            Ret_stim_list_2 = {2 19 3 23 8 6 15 9 27 22 20 17 13 24 18 16 21 15 18 12};
            Ret_stim_1_Enc_location = {1 1 1 1 1 1 2 1 2 2 1 1 2 1 2 1 1 2 1 1};
            Ret_stim_2_Enc_location = {2 2 2 2 2 2 1 2 1 1 1 2 2 1 2 2 2 2 1 2};
        elseif Parm.SetCond{sess,run} == 7
            Ret_stim_list_1 = {5 13 9  15 6 17 8 7 20 6 25 12 11 3 11 19 2 16 10 4};
            Ret_stim_list_2 = {23 13 9 22 7 18 27 24 21 17 26 14 24 18 22 20 2 16 26 4};
            Ret_stim_1_Enc_location = {1 1 1 2 2 2 1 1 1 1 2 2 1 1 2 2 1 1 1 1};
            Ret_stim_2_Enc_location = {1 2 2 2 2 1 2 1 2 1 1 1 2 2 1 2 2 2 2 2};
        elseif Parm.SetCond{sess,run} == 8
            Ret_stim_list_1 = {9 5 5 12 17 19 12 8 16 2 27 13 17 4 14 11 9 26 6 7};
            Ret_stim_list_2 = {20 10 21 24 24 19 14 25 16 3 8 23 18 22 15 11 21 27 6 7};
            Ret_stim_1_Enc_location = {2 2 1 1 2 1 2 1 1 2 2 2 1 1 1 1 1 1 1 1};
            Ret_stim_2_Enc_location = {1 1 1 1 2 2 2 2 2 2 2 1 1 2 2 2 2 2 2 2};
        elseif Parm.SetCond{sess,run} == 9
            Ret_stim_list_1 = {3 8 11 20 13 9 7 4 5 22 15 10 6 17 23 2 8 12 17 16};
            Ret_stim_list_2 = {3 21 12 21 14 9 27 5 24 24 15 20 19 26 23 10 25 22 18 16};
            Ret_stim_1_Enc_location = {1 2 2 2 1 1 2 2 2 1 1 2 2 1 1 2 1 1 2 1};
            Ret_stim_2_Enc_location = {2 2 2 1 1 2 1 1 2 1 2 1 1 2 2 1 1 2 1 2};
        end
        
        % Retrieval: trials
        for trial=1:Parm.ntrial_Ret
            
            if Ret_stim_1_Enc_location{trial} == 1
                Ret_stim_1{trial} = Parm.Enc_stim_left{sess,run}{Ret_stim_list_1{trial}};
            elseif Ret_stim_1_Enc_location{trial} == 2
                Ret_stim_1{trial} = Parm.Enc_stim_right{sess,run}{Ret_stim_list_1{trial}};
            end
            
            if Ret_stim_2_Enc_location{trial} == 1
                Ret_stim_2{trial} = Parm.Enc_stim_left{sess,run}{Ret_stim_list_2{trial}};
            elseif Ret_stim_2_Enc_location{trial} == 2
                Ret_stim_2{trial} = Parm.Enc_stim_right{sess,run}{Ret_stim_list_2{trial}};
            end
            
            answer_shuffle = randi([1 2], 1);
            if answer_shuffle == 1
                Ret_stim_left{sess,run}{trial} = Ret_stim_1{trial};
                Ret_stim_right{sess,run}{trial} = Ret_stim_2{trial};
            elseif answer_shuffle == 2
                Ret_stim_left{sess,run}{trial} = Ret_stim_2{trial};
                Ret_stim_right{sess,run}{trial} = Ret_stim_1{trial};
            end
        end
        
    end
end
Parm.Ret_stim_left = Ret_stim_left;
Parm.Ret_stim_right = Ret_stim_right;

% Delay Arithmetics set
% Equation sets
numbersread = readcell([ori_dir '\Stimuli_set\Delay_arithmetics.xlsx']);
numbersread = numbersread(2:81,:);
[m,n] = size(numbersread);
idx = randperm(m);

for i=1:m
    numbers(i,:) = numbersread(idx(i),:);
end

Del_numbers{1,1} = numbers(1:10,1:6);
Del_numbers{1,2} = numbers(11:20,1:6);
Del_numbers{1,3} = numbers(21:30,1:6);
Del_numbers{1,4} = numbers(61:70,1:6);
Del_numbers{2,1} = numbers(31:40,1:6);
Del_numbers{2,2} = numbers(41:50,1:6);
Del_numbers{2,3} = numbers(51:60,1:6);
Del_numbers{2,4} = numbers(71:80,1:6);
Parm.Del_numbers = Del_numbers;

% Retrieval Answer (1; TOJ(before) / 2; Association / 3; TOJ(after))
% Retrieval Spatial Condition (1; location constant / 2; location shuffle)
% bring this to Ret.m
for sess=1:Parm.nsess
    for run=1:Parm.nrun
        for trial=1:Parm.ntrial_Ret
            index_left_left = find(strcmp(Parm.Enc_stim_left{sess,run}, Parm.Ret_stim_left{sess,run}{trial}));
            index_left_right = find(strcmp(Parm.Enc_stim_right{sess,run}, Parm.Ret_stim_left{sess,run}{trial}));
            TF_left_left = isempty(index_left_left);
            TF_left_right = isempty(index_left_right);
            index_right_left = find(strcmp(Parm.Enc_stim_left{sess,run}, Parm.Ret_stim_right{sess,run}{trial}));
            index_right_right = find(strcmp(Parm.Enc_stim_right{sess,run}, Parm.Ret_stim_right{sess,run}{trial}));
            TF_right_left = isempty(index_right_left);
            TF_right_right = isempty(index_right_right);
            if TF_left_left == 1
                index_left = index_left_right;
                RetCond_location{sess,run}{trial,1} = 2;
                RetCond_location{sess,run}{trial,2} = 'left_shuffle';
            elseif TF_left_right == 1
                index_left = index_left_left;
                RetCond_location{sess,run}{trial,1} = 1;
                RetCond_location{sess,run}{trial,2} = 'left_constant';
            end
            
            if TF_right_left == 1
                index_right = index_right_right;
                RetCond_location{sess,run}{trial,3} = 1;
                RetCond_location{sess,run}{trial,4} = 'right_constant';
            elseif TF_right_right == 1
                index_right = index_right_left;
                RetCond_location{sess,run}{trial,3} = 2;
                RetCond_location{sess,run}{trial,4} = 'right_shuffle';
            end
                        
            if TF_left_right == 1
                if TF_right_left == 1
                    RetCond_location{sess,run}{trial,5} = 1;
                    RetCond_location{sess,run}{trial,6} = 'constant';
                else
                    RetCond_location{sess,run}{trial,5} = 2;
                    RetCond_location{sess,run}{trial,6} = 'shuffle';
                end
            else
                RetCond_location{sess,run}{trial,5} = 2;
                RetCond_location{sess,run}{trial,6} = 'shuffle';
            end
            
            if index_left == index_right
                RetCond_location{sess,run}{trial,7} = 2;
                RetCond_location{sess,run}{trial,8} = 'Asso';
                RetCond_distance{sess,run}{trial} = 0;
            elseif index_left < index_right
                RetCond_location{sess,run}{trial,7} = 1;
                RetCond_location{sess,run}{trial,8} = 'TOJ(before)';
                RetCond_distance{sess,run}{trial} = index_right - index_left;
            elseif index_left > index_right
                RetCond_location{sess,run}{trial,7} = 3;
                RetCond_location{sess,run}{trial,8} = 'TOJ(after}';
                RetCond_distance{sess,run}{trial} = index_left - index_right;
            end
            
            RetCond_location{sess,run}{trial,9} = index_left;
            RetCond_location{sess,run}{trial,10} = index_right;
            
        end
    end
end

Parm.RetCond_location = RetCond_location;
Parm.RetCond_distance = RetCond_distance;

save(['Parm_' num2str(subject)], 'Parm')
