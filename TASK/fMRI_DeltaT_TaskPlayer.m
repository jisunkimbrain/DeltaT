% DeltaT task
% fMRI experiment
% Task Player
% Ji Sun Kim

% Task instructions are in Korean

%% Exp setup
SID = input('Subject ID: ');
load(['Exp_Parm\Parm_' num2str(SID) '.mat']);

% emergency starter for interpution
start_sess = 1;   % change this for each session
start_run = 1;
start_trial = 1;

% Results header info
Headers_Enc = {'Trial', 'TR', 'Time', 'ItemNum', 'Stim_name_left', 'Stim_name_right'};
Headers_Del = {'Trial', 'Time', 'Question', 'Option1', 'Option2', 'Answer', 'Response', ...
    'Correctness', 'Correctness(str)', 'RT'};
Headers_Ret = {'Trial', 'TR', 'Time', 'Stim_left_RET', 'Stim_right_RET', ...
    'Stim_left_RET_idx', 'Spatial_Condition_left', 'Spatial_Condition_left(str)', ...
    'Enc_Condition_left', ...
    'Stim_right_RET_idx', 'Spatial_Condition_right', 'Spatial_Conditon_right(str)', ...
    'Enc_Condition_right', ...
    'Temp_Diff', ...
    'Condition_Spatial', 'Condition_Spatial(str)', ...
    'Condition_Question', 'Condition_Question(str)', 'Correct_Answer(str)', ...
    'Response', 'Response(str)', 'Correctness', 'Correctness(str)', 'RT'};

% display setup
config_display(1,2,[1,1,1],[0,0,0]); % (window, resolution, [background], [foreground]) [204,204,204]
config_display(0);
% config_display(0,2,[1,1,1],[0,0,0]);
config_keyboard;

% directory
ori_dir = cd;  %-- start from Temporal_difference_task folder
Stim_dir = [ori_dir '\Stimuli_set\JObjects\JObjects_jpgresized'];
Box_dir = [ori_dir '\Stimuli_set\JObjects'];

start_cogent

%% Task

for sess = start_sess
%     if sess == start_sess
%         start_cogent
%     end
    
    for run = start_run:Parm.nrun
        
        %% Encoding session
        % directory
%         Stim_dir = '\\kaistleelab.synology.me\jisun1006\Tasks\Temporal_difference_task\Object_database\JObjects\JObjects_jpgresized';
%         Box_dir = '\\kaistleelab.synology.me\jisun1006\Tasks\Temporal_difference_task\Object_database\JObjects';
%         
        % message before task
        settextstyle('HY울릉도B',35);
        preparestring('나오는 물체들을 잘 기억해주세요',1,0,70); % temporary
        preparestring('시작',1,0,-30);
        drawpict(1);
        
        % response key setup
        map = getkeymap; readkeys;
        mapp=[map.Space, map.E, map.W, map.N, map.Escape]; % KAIST fMRI 4-button setting
        
        waitkeydown(inf,19); % scanner sync (s)
        tic
        clearpict(1);
        
        for trial = start_trial:Parm.ntrial_Enc
            
            % emergency code initializer
            if run > start_run
                start_trial = 1;
            end
            
            % stimuli set
            box = [Box_dir, '\box1.jpg'];
            
            % fixation
            settextstyle('HY울릉도B', 30);
            setforecolour(0, 0, 0);
            preparestring('+',1,0,50);
            
            % location condition (box) & spatial cue
            if Parm.EncCond_location{sess,run} == 1      % (1; up-down)
                loadpict(box,1,-220,180);
                loadpict(box,1,220,-60);
                if Parm.EncCond_item{sess,run}{trial} == 1
                    settextstyle('HY울릉도B', 20);
                    setforecolour(0, 0, 0);
                    preparestring('+',1,-220,180);
                elseif Parm.EncCond_item{sess,run}{trial} == 2
                    settextstyle('HY울릉도B', 20);
                    setforecolour(0, 0, 0);
                    preparestring('+',1,-220,180);
                    preparestring('+',1,220,-60);
                elseif Parm.EncCond_item{sess,run}{trial} == 3
                    settextstyle('HY울릉도B', 20);
                    setforecolour(0, 0, 0);
                    preparestring('+',1,220,-60);
                end
            elseif Parm.EncCond_location{sess,run} == 2      % (2; down-up)
                loadpict(box,1,-220,-60);
                loadpict(box,1,220,180);
                if Parm.EncCond_item{sess,run}{trial} == 1
                    settextstyle('HY울릉도B', 20);
                    setforecolour(0, 0, 0);
                    preparestring('+',1,-220,-60);
                elseif Parm.EncCond_item{sess,run}{trial} == 2
                    settextstyle('HY울릉도B', 20);
                    setforecolour(0, 0, 0);
                    preparestring('+',1,-220,-60);
                    preparestring('+',1,220,180);
                elseif Parm.EncCond_item{sess,run}{trial} == 3
                    settextstyle('HY울릉도B', 20);
                    setforecolour(0, 0, 0);
                    preparestring('+',1,220,180);
                end
            end
            
            % present spatial cue
            start_point = time;   % task start
            drawpict(1);
            timepoint = toc
            wait(500);
            
            % stimuli set
            if Parm.EncCond_location{sess,run} == 1      % (1; up-down)
                if Parm.EncCond_item{sess,run}{trial} == 1
                    loadpict([Stim_dir '\' Parm.Enc_stim_left{sess,run}{trial}], 1, -220, 180);
                elseif Parm.EncCond_item{sess,run}{trial} == 2
                    loadpict([Stim_dir '\' Parm.Enc_stim_left{sess,run}{trial}], 1, -220, 180);
                    loadpict([Stim_dir '\' Parm.Enc_stim_right{sess,run}{trial}], 1, 220, -60);
                elseif Parm.EncCond_item{sess,run}{trial} == 3
                    loadpict([Stim_dir '\' Parm.Enc_stim_right{sess,run}{trial}], 1, 220, -60);
                end
            elseif Parm.EncCond_location{sess,run} == 2      % (2; down-up)
                if Parm.EncCond_item{sess,run}{trial} == 1
                    loadpict([Stim_dir '\' Parm.Enc_stim_left{sess,run}{trial}], 1, -220, -60);
                elseif Parm.EncCond_item{sess,run}{trial} == 2
                    loadpict([Stim_dir '\' Parm.Enc_stim_left{sess,run}{trial}], 1, -220, -60);
                    loadpict([Stim_dir '\' Parm.Enc_stim_right{sess,run}{trial}], 1, 220, 180);
                elseif Parm.EncCond_item{sess,run}{trial} == 3
                    loadpict([Stim_dir '\' Parm.Enc_stim_right{sess,run}{trial}], 1, 220, 180);
                end
            end
            
            % present stimuli
            setforecolour(0, 0, 0);
            % s key wait, inf
            % waitkeydown(inf,mapp.S);
            drawpict(1);
            
            % ready for response key press
            [key, t, n] = waitkeydown(Parm.Stimdur_Enc*1000, mapp);
            
            % prevent multiple key press
            if n==1
                response = key(1);
                response_time = t(1) - start_point;
                setforecolour(1, 0, 0);
                preparestring(' 버튼을 잘못 누르셨습니다 ',1,0,-220);
                drawpict(1);
                wait(3000 - response_time);
                clearpict(1);
            elseif n>1
                response = key(1);
                response_time = t(1) - start_point;
                setforecolour(1, 0, 0);
                preparestring(' 버튼을 잘못 누르셨습니다 ',1,0,-220);
                drawpict(1);
                wait(3000 - response_time);
                clearpict(1);
            else
                response = 0;
                response_time = 0;
            end
            
            % Stopper
            if n~=0
                if key == map.Escape
                    stop_cogent;
                end
            end
            
            % save results (Results cell) Headers = {'Trial', 'TR', 'ItemNum', 'Stim_name_left', 'Stim_name_right'};
            Results{trial, 1} = trial;
            Results{trial, 2} = Parm.TR_Enc{sess,run}{trial};
            Results{trial, 3} = timepoint;
            Results{trial, 4} = Parm.EncCond_item{sess,run}{trial};
            Results{trial, 5} = char(Parm.Enc_stim_left{sess,run}{trial});
            Results{trial, 6} = char(Parm.Enc_stim_right{sess,run}{trial});
            
            % fixation drawing
            clearpict(1);
            settextstyle('HY울릉도B', 30);
            setforecolour(0, 0, 0);
            preparestring('+',1,0,50);
            drawpict(1);
            wait(Parm.ISI_Enc * 1000);
            clearpict(1);
            
        end
        
        % saving result with mat file
        Results = [Headers_Enc; Results];
        % prevent overwriting
        if ~exist(['BHV_Results_ENC\mat\Exp_results_Subj_' num2str(SID) '_Sess_' num2str(sess) '_Run_' num2str(run) '.mat'])
            save(['BHV_Results_ENC\mat\Exp_results_Subj_' num2str(SID) '_Sess_' num2str(sess) '_Run_' num2str(run) '.mat'], 'Results');
            xlswrite(['BHV_Results_ENC\xls\Exp_results_Subj_' num2str(SID) '_Sess_' num2str(sess) '_Run_' num2str(run) '.xls'], Results);
        else
            save(['BHV_Results_ENC\mat\Exp_results_Subj_' num2str(SID) '_Sess_' num2str(sess) '_Run_' num2str(run) '_2.mat'], 'Results');
            xlswrite(['BHV_Results_ENC\xls\Exp_results_Subj_' num2str(SID) '_Sess_' num2str(sess) '_Run_' num2str(run) '_2.xls'], Results);
        end
        clear Results
        
        % message after task
        settextstyle('HY울릉도B', 40);
        setforecolour(0, 0, 0);
        preparestring('잠시 후 다음 실험이 진행됩니다', 1, 0, 0);
        drawpict(1);
        toc
        wait(3000);
        clearpict(1);
        
        %% Delay task
        % message before task
        settextstyle('HY울릉도B',35);
        preparestring('계산이 맞는지 답변해주세요',1,0,70);
        preparestring('시작',1,0,-30);
        drawpict(1);
        wait(5000); % get rid of this when scanning
        toc
        
        % response key setup
        map = getkeymap; readkeys;
        mapp=[map.Space, map.E, map.W, map.N, map.Escape]; % KAIST fMRI 4-button setting
        
%         waitkeydown(inf,19); % scanner sync
        clearpict(1);
        
        for trial = start_trial:Parm.ntrial_Del
            
            % emergency code initializer
            if run > start_run
                start_trial = 1;
            end
            
            % present stimuli
            start_point = time;   % task start
            setforecolour(0, 0, 0);
            settextstyle('HY울릉도B',40);
            preparestring([num2str(Parm.Del_numbers{sess,run}{trial,1}),' ', Parm.Del_numbers{sess,run}{trial,2},' ', num2str(Parm.Del_numbers{sess,run}{trial,3}),' ='],1,0,10);
            preparestring([num2str(Parm.Del_numbers{sess,run}{trial,4}), '                 ', num2str(Parm.Del_numbers{sess,run}{trial,5})],1,0,-240);
            drawpict(1);
            timepoint = toc
            
            % ready for response key press
            [key, t, n] = waitkeydown(3000, mapp);
            
            % prevent multiple key press
            if n==1
                response = key(1);
                response_time = t(1) - start_point ;
            elseif n>1
                response = key(1);
                response_time = t(1) - start_point ;
            else
                response = 0;
                response_time = 0;
            end
            
            % get response
            if response == map.E
                setforecolour(1, 0, 0);
                preparestring('★                     ',1,0,-240);
                drawpict(1);
                wait(3000 - response_time);          % change this back to 3000
                clearpict(1);
            elseif response == map.N
                setforecolour(1, 0, 0);
                preparestring('                    ★',1,0,-240);
                drawpict(1);
                wait(3000 - response_time);
                clearpict(1);
            elseif response == map.W
                setforecolour(1, 0, 0);
                preparestring(' 버튼을 잘못 누르셨습니다 ',1,0,-180);
                drawpict(1);
                wait(3000 - response_time);
                clearpict(1);
            elseif response == 0
                clearpict(1);
            end
            
            % Stopper
            if n~=0
                if key == map.Escape
                    stop_cogent;
                end
            end
            
            % Headers_Del = {'Trial', 'Question', 'Option1', 'Option2', 'Answer', 'Response', ... 'Correctness', 'Correctness(str)', 'RT'};
    
            % save results (Results cell)
            Results{trial, 1} = trial;
            Results{trial, 2} = timepoint;
            Results{trial, 3} = [num2str(Parm.Del_numbers{sess,run}{trial,1}) Parm.Del_numbers{sess,run}{trial,2} num2str(Parm.Del_numbers{sess,run}{trial,3}), '='];
            Results{trial, 4} = Parm.Del_numbers{sess,run}{trial,4};
            Results{trial, 5} = Parm.Del_numbers{sess,run}{trial,5};

            Results{trial, 6} = Parm.Del_numbers{sess,run}{trial,6};
            
            if response ~= 0
                if response == map.E
                    Results{trial, 7} = 1;     %-- Response (1; left, 2; right)
                    if Parm.Del_numbers{sess,run}{trial,6} == 1
                        Results{trial,8} = 1;
                        Results{trial,9} = 'Corr';
                    elseif Parm.Del_numbers{sess,run}{trial,6} == 2
                        Results{trial,8} = 2;
                        Results{trial, 9} = 'Incorr';
                    end
                elseif response == map.N
                    Results{trial, 7} = 2;
                    if Parm.Del_numbers{sess,run}{trial,6} == 1
                        Results{trial,8} = 2;
                        Results{trial,9} = 'Incorr';
                    elseif Parm.Del_numbers{sess,run}{trial,6} == 2
                        Results{trial,8} = 1;
                        Results{trial,9} = 'Corr';
                    end
                elseif response == map.W
                    Results{trial, 7} = 3;    %-- Response (0; Miss)
                    Results{trial, 8} = 2;    %-- 2; Incorr & Miss
                    Results{trial, 9} = 'Wrong button';
                    Results{trial, 10} = 0;
                end
                Results{trial, 10} = response_time;
            else
                Results{trial, 7} = 0;    %-- Response (0; Miss)
                Results{trial, 8} = 2;    %-- 2; Incorr & Miss
                Results{trial, 9} = 'Miss';
                Results{trial, 10} = 0;
            end
        end
        
        % saving result with mat file
        Results = [Headers_Del; Results];
        T = table(Results);
        
        % prevent overwriting
        if ~exist(['BHV_Results_Del\mat\Exp_results_Subj_' num2str(SID) '_Sess_' num2str(sess) '_Run_' num2str(run) '.mat'])
            save(['BHV_Results_Del\mat\Exp_results_Subj_' num2str(SID) '_Sess_' num2str(sess) '_Run_' num2str(run) '.mat'], 'Results');
            writetable(T, ['BHV_Results_Del\xls\Exp_results_Subj_' num2str(SID) '_Sess_' num2str(sess) '_Run_' num2str(run) '.xls']);
        else
            save(['BHV_Results_Del\mat\Exp_results_Subj_' num2str(SID) '_Sess_' num2str(sess) '_Run_' num2str(run) '_2.mat'], 'Results');
            writetable(T,['BHV_Results_Del\xls\Exp_results_Subj_' num2str(SID) '_Sess_' num2str(sess) '_Run_' num2str(run) '_2.xls']);
        end
        clear Results
        
        % message after task
        settextstyle('HY울릉도B', 40);
        setforecolour(0, 0, 0);
        preparestring('잠시 후 다음 실험이 진행됩니다', 1, 0, 0);
        drawpict(1);
        toc
        wait(3000);
        clearpict(1);
        
        %% Retrieval task
%         Stim_dir = '\\kaistleelab.synology.me\jisun1006\Tasks\Temporal_difference_task\Object_database\JObjects\JObjects_jpgresized';
%         Box_dir = '\\kaistleelab.synology.me\jisun1006\Tasks\Temporal_difference_task\Object_database\JObjects';

        % message before task
        clearpict(1);
        settextstyle('HY울릉도B',30);
        preparestring('두 개의 사물 중 어떤게 먼저 등장했는지 답변해주세요',1,0,60);
        preparestring('왼쪽               동시에               오른쪽',1,0,-60);
        drawpict(1);
        wait(5000); % get rid of this when scanning
        toc
        
        % response key setup
        map = getkeymap; readkeys;
        mapp=[map.Space, map.E, map.W, map.N, map.Escape]; % KAIST fMRI 4-button setting
        
%         waitkeydown(inf,19); % s
        clearpict(1);
        
        for trial = start_trial:Parm.ntrial_Ret
            
            % emergency code initializer
            if run > start_run
                start_trial = 1;
            end
            
            % stimuli set
            fname_left = fullfile(Stim_dir, Parm.Ret_stim_left{sess,run}{trial});
            box = fullfile(Box_dir, 'box1.jpg');
            fname_right = fullfile(Stim_dir, Parm.Ret_stim_right{sess,run}{trial});
            
            start_point = time;   % task start
            
            % stimulati location set
            if Parm.EncCond_location{sess,run} == 1      % (1; up-down)
                loadpict(box,1,-220,180);
                loadpict(box,1,220,-60);
                loadpict(fname_left,1,-220,180);
                loadpict(fname_right,1,220,-60);
            elseif Parm.EncCond_location{sess, run} == 2      % (2; down-up)
                loadpict(box,1,-220,-60);
                loadpict(box,1,220,180);
                loadpict(fname_left,1,-220,-60);
                loadpict(fname_right,1,220,180);
            end
            
            setforecolour(0, 0, 0);
            settextstyle('HY울릉도B',30);
            preparestring('어떤 물체가 먼저 등장했나요?',1,0,-205);
            preparestring('왼쪽               동시에               오른쪽',1,0,-240);
            
            % present stimuli
            drawpict(1);
            timepoint = toc
            
            % ready for response key press
            [key, t, n] = waitkeydown(4000, mapp);
            
            % prevent multiple key press
            if n==1
                response = key(1);
                response_time = t(1) - start_point ;
            elseif n>1
                response = key(1);
                response_time = t(1) - start_point ;
            else
                response = 0;
                response_time = 0;
            end
            
            % get response
            if response == map.E
                setforecolour(1, 0, 0);
                preparestring('★                                                ',1,0,-240);
                drawpict(1);
                wait(Parm.Stimdur_Ret*1000 - response_time);
                clearpict(1);
            elseif response == map.W
                setforecolour(1, 0, 0);
                preparestring('★   ',1,0,-240);
                drawpict(1);
                wait(Parm.Stimdur_Ret*1000 - response_time);
                clearpict(1);
            elseif response == map.N
                setforecolour(1, 0, 0);
                preparestring('                                             ★',1,0,-240);
                drawpict(1);
                wait(Parm.Stimdur_Ret*1000 - response_time);
                clearpict(1);
            elseif response == 0
                clearpict(1);
            end
            
            % Stopper
            if n~=0
                if key == map.Escape
                    stop_cogent;
                end
            end
            
            % save results (Results cell)
            Results{trial, 1} = trial;
            Results{trial, 2} = Parm.TR_Ret{sess,run}{trial};
            Results{trial, 3} = timepoint;
            Results{trial, 4} = char(Parm.Ret_stim_left{sess,run}{trial});
            Results{trial, 5} = char(Parm.Ret_stim_right{sess,run}{trial});
            
            index_left_left = find(strcmp(Parm.Enc_stim_left{sess,run}, Parm.Ret_stim_left{sess,run}{trial}));
            index_left_right = find(strcmp(Parm.Enc_stim_right{sess,run}, Parm.Ret_stim_left{sess,run}{trial}));
            TF_left_left = isempty(index_left_left);
            TF_left_right = isempty(index_left_right);
            index_right_left = find(strcmp(Parm.Enc_stim_left{sess,run}, Parm.Ret_stim_right{sess,run}{trial}));
            index_right_right = find(strcmp(Parm.Enc_stim_right{sess,run}, Parm.Ret_stim_right{sess,run}{trial}));
            TF_right_left = isempty(index_right_left);
            TF_right_right = isempty(index_right_right);
            
            % Spatial Condition
            if TF_left_left == 1
                index_left = index_left_right;
                Results{trial, 6} = index_left;
                Results{trial, 7}= 2;
                Results{trial, 8} = 'left_shuffle';
            elseif TF_left_right == 1
                index_left = index_left_left;
                Results{trial, 6} = index_left;
                Results{trial, 7} = 1;
                Results{trial, 8} = 'left_constant';
            end
            
            Results{trial, 9} = Parm.EncCond_item{sess,run}{index_left};
            
            if TF_right_left == 1
                index_right = index_right_right;
                Results{trial, 10} = index_right;
                Results{trial, 11} = 1;
                Results{trial, 12} = 'right_constant';
            elseif TF_right_right == 1
                index_right = index_right_left;
                Results{trial, 10} = index_right;
                Results{trial, 11} = 2;
                Results{trial, 12} = 'right_shuffle';
            end
            
            Results{trial, 13} = Parm.EncCond_item{sess,run}{index_right};
            
            if index_left == index_right
                Results{trial, 14} = 0;
            elseif index_left > index_right
                Results{trial, 14} = index_left - index_right;
            elseif index_left < index_right
                Results{trial, 14} = index_right - index_left;
            end
            
            if TF_left_right == 1
                if TF_right_left == 1
                    Results{trial, 15} = 1;
                    Results{trial, 16} = 'constant';
                else
                    Results{trial, 15} = 2;
                    Results{trial, 16} = 'shuffle';
                end
            else
                Results{trial, 15} = 2;
                Results{trial, 16} = 'shuffle';
            end
            
            % Question Condition
            if index_left == index_right
                Results{trial, 17} = 2;
                Results{trial, 18} = 'Asso';
                Results{trial, 19} = 'Together';
            elseif index_left < index_right
                Results{trial, 17} = 1;
                Results{trial, 18} = 'TOJ(before)';
                Results{trial, 19} = 'Left';
            elseif index_left > index_right
                Results{trial, 17} = 3;
                Results{trial, 18} = 'TOJ(after)';
                Results{trial, 19} = 'Right';
            end
            
            % Response & RT & Correctness
            
            if response ~= 0
                if response == map.E
                    Results{trial, 20} = response;
                    Results{trial, 21} = 'Left';
                    if index_left < index_right
                        Results{trial, 22} = 1;
                        Results{trial, 23} = 'Correct';
                    else
                        Results{trial, 22} = 0;
                        Results{trial, 23} = 'Incorret';
                    end
                elseif response == map.W
                    Results{trial, 20} = response;
                    Results{trial, 21} = 'Together';
                    if index_left == index_right
                        Results{trial, 22} = 1;
                        Results{trial, 23} = 'Correct';
                    else
                        Results{trial, 22} = 0;
                        Results{trial, 23} = 'Incorret';
                    end
                elseif response == map.N
                    Results{trial, 20} = response;
                    Results{trial, 21} = 'Right';
                    if index_left > index_right
                        Results{trial, 22} = 1;
                        Results{trial, 23} = 'Correct';
                    else
                        Results{trial, 22} = 0;
                        Results{trial, 23} = 'Incorret';
                    end
                end
                Results{trial, 24} = response_time;
            else
                Results{trial, 20} = 0;
                Results{trial, 21} = 'Miss';
                Results{trial, 22} = 0;
                Results{trial, 23} = 'Miss';
                Results{trial, 24} = 0;
            end
            
            % fixation drawing
            clearpict(1);
            settextstyle('HY울릉도B', 30);
            setforecolour(0, 0, 0);
            preparestring('+',1,0,50);
            drawpict(1);
            wait(Parm.ISI_Ret{sess,run}{trial} * 2000);
            clearpict(1);
            
        end
        
        % saving result with mat file
        Results = [Headers_Ret; Results];
        % prevent overwriting
        if ~exist(['BHV_Results_RET\mat\Exp_results_Subj_' num2str(SID) '_Sess_' num2str(sess) '_Run_' num2str(run) '.mat'])
            save(['BHV_Results_RET\mat\Exp_results_Subj_' num2str(SID) '_Sess_' num2str(sess) '_Run_' num2str(run) '.mat'], 'Results');
            xlswrite(['BHV_Results_RET\xls\Exp_results_Subj_' num2str(SID) '_Sess_' num2str(sess) '_Run_' num2str(run) '.xls'], Results);
        else
            save(['BHV_Results_RET\mat\Exp_results_Subj_' num2str(SID) '_Sess_' num2str(sess) '_Run_' num2str(run) '_2.mat'], 'Results');
            xlswrite(['BHV_Results_RET\xls\Exp_results_Subj_' num2str(SID) '_Sess_' num2str(sess) '_Run_' num2str(run) '_2.xls'], Results);
        end
        clear Results
        
        % message after task
        settextstyle('HY울릉도B', 40);
        setforecolour(0, 0, 0);
        preparestring('잠시 후 같은 실험이 반복됩니다', 1, 0, 0);
        drawpict(1);
        wait(3000);
        clearpict(1);
    end
    
    % message after task    (changed 03.02, check if no error)
    if sess == 1
        if run == 4
            settextstyle('HY울릉도B', 40);
            setforecolour(0, 0, 0);
            preparestring('고생하셨습니다', 1, 0, 60);
            preparestring('휴식 후 같은 실험을 반복해주세요', 1, 0, 0);
            drawpict(1);
            waitkeydown(10000);
            clearpict(1);
        else
            settextstyle('HY울릉도B', 40);
            setforecolour(0, 0, 0);
            preparestring('고생하셨습니다', 1, 0, 60);
            preparestring('같은 실험을 반복해주세요', 1, 0, 0);
            drawpict(1);
            waitkeydown(10000);
            clearpict(1);
        end
    elseif sess == 2
        if run == 4
            settextstyle('HY울릉도B', 40);
            setforecolour(0, 0, 0);
            preparestring('고생하셨습니다', 1, 0, 60);
            preparestring('실험이 끝났습니다 휴식을 취해주세요', 1, 0, 0);
            drawpict(1);
            waitkeydown(10000); % scanner sync (s)
            clearpict(1);
        else
            settextstyle('HY울릉도B', 40);
            setforecolour(0, 0, 0);
            preparestring('고생하셨습니다', 1, 0, 60);
            preparestring('같은 실험을 반복해주세요', 1, 0, 0);
            drawpict(1);
            waitkeydown(10000);
            clearpict(1);
        end
    else
        settextstyle('HY울릉도B', 40);
        setforecolour(0, 0, 0);
        preparestring('고생하셨습니다', 1, 0, 60);
        preparestring('휴식 후 같은 실험을 반복해주세요', 1, 0, 0);
        drawpict(1);
        waitkeydown(10000);
        clearpict(1);
    end
            
%     wait(3000);
%     clearpict(1);
    
end

stop_cogent
