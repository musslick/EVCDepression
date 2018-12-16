classdef LearningFnc < EVC.EVCFnc
    
    % this class implements DDM-specific functions
    % most of the functions may reference an EVCDDM instance

    properties
       % output                               % fnc handle{} -> points to dynamic output value 
    end
    
    properties (Constant)
        SIMPLE_SALIENCY_RL = 1;             % learns about saliency based on error feedback
        SIMPLE_STATE_RL = 2;                % learns about any state representation based of experience of that representation
        FORAGING_REWARDS = 3;
        FILL_CLEAR = 4;
        SIMPLE_STATE_RL_UW = 5;
        AUTOMATIC_RL = 6;
        CONTROL_EFFICACY = 8;
        GILZENRAT_2010 = 9;
        APPROACH_AVOIDANCE = 10;
        BRAUN_ARRINGTON_2018 = 11;
        WESTBROOK_BRAVER_2015 = 12;
        
        % holds amount of required parameters for each function 
        paramReqLF = [1, ...  SIMPLE_SALIENCY_RL: 1) learning rate
                      1, ...   SIMPLE_STATE_RL: 1) learning rate
                      1 ...    FORAGING_REWARDS: 1) learning rate ? (maybe not necessary)
                      3, ...    FILL_CLEAR: 1) max reward difference, 2) jump probability
                      1, ...   SIMPLE_STATE_RL_UW: 1) learning rate
                      1, ... AUTOMATIC_RL: 1) learning rate
                      1, ... TRIAL_DELTA: 1) learning rate
                      1, ... CONTROL_EFFICACY: 1) learning rate
                      1, ... GILZENRAT_2010: 1) learning rate
                      1, ... APPROACH_AVOIDANCE: 1) learning rate
                      1, ... BRAUN_ARRINGTON_2018: 1) learning rate
                      1, ... WESTBROOK_BRAVER_2015: 1) learning rate
                     ];
                
        % holds amount of required dynamic input parameters for each function
        inputReqLF = [2, ...  SIMPLE_SALIENCY_RL: 1) actual outcome 2) expected outcome
                      1, ...  SIMPLE_STATE_RL: 1) learned value
                      1 ...  FORAGING_REWARDS: 1) actual outcome
                      1, ...    FILL_CLEAR: 1) previous reward intensity
                      1, ...  SIMPLE_STATE_RL_UW: 1) learned value
                      1, ... AUTOMATIC_RL: 1) actual outcome
                      0, ... TRIAL_DELTA: 1) learning rate
                      2, ... CONTROL_EFFICACY: 1) actual outcome 2) expected outcome
                      1 ...  GILZENRAT_2010: 1) actual outcome
                      1 ...  APPROACH_AVOIDANCE: 1) actual outcome
                      1 ...  BRAUN_ARRINGTON_2018: 1) actual outcome
                      1 ...  WESTBROOK_BRAVER_2015: 1) actual outcome
                     ];
    end
    
    methods
        
        function this = LearningFnc(type, params, varargin)
          
          % call superclass constructor to create an EVC function instance
          this = this@EVC.EVCFnc(type, params, 'paramReq', 'paramReqLF');
         
          [inputVal EVCM] = this.extractInput(varargin);
          this.input = inputVal;
          this.EVCModel = EVCM;

        end
        
    end
    
    methods (Access = public)
        
        % calculates output value dependent on specified function type
        function out = getVal(this, varargin)

            [inputVal, EVCM] = this.extractInput(varargin);
            
            % general learning functions
            
            if(this.type == this.SIMPLE_SALIENCY_RL)
               out = simpleSaliencyRL(this, EVCM); 
            end
            
            if(this.type == this.SIMPLE_STATE_RL)
               out = simpleStateRL(this, EVCM); 
            end
            
            if(this.type == this.SIMPLE_STATE_RL_UW)
               out = simpleStateRL_unweighted(this, EVCM); 
            end
            
            if(this.type == this.FORAGING_REWARDS)
               out = setForagingRewards(this, EVCM); 
            end
            
            if(this.type == this.FILL_CLEAR)
               out = fillClearUpdate(this, EVCM); 
            end
            
            if(this.type == this.AUTOMATIC_RL)
               out = reinforceAutomaticPathway(this, EVCM);
            end
            
            if(this.type == this.CONTROL_EFFICACY)
               out = controlEfficacy(this, EVCM);
            end
            
            if(this.type == this.GILZENRAT_2010)
               out = setGilzenrat2010Trial(this, EVCM); 
            end
            
            if(this.type == this.APPROACH_AVOIDANCE)
               out = setApproachAvoidanceTrial(this, EVCM); 
            end
            
            if(this.type == this.BRAUN_ARRINGTON_2018)
               out = setBraunArrington2018Trial(this, EVCM); 
            end
            
            if(this.type == this.WESTBROOK_BRAVER_2015)
               out = setWestbrookBraver2015Trial(this, EVCM); 
            end
            
            % implementation-specific learning functions may need an
            % additional verification of EVCM class (e.g. EVCDDM)
            
        end
       
    end
    
    methods (Access = private)
        
        %% function definitions
        
        % STIMBIAS: returns a response bias for the specified interval (specified
        % in this.params) given 
        function out = simpleSaliencyRL(this, EVCM)
            
            %% learn about automatic processing bias based on outcome differences between actual and expected state
            
            % EVCFnc parameters
            feedback = this.input{1}() - this.input{2}(); % actual error prob - expected error prob
            learningRate = this.params{1};
            
            % parameters retrieved from EVC model
            stateIdx = EVCM.getExpectedStateIdx();
            stimSalience = EVCM.State.ExpectedSpace(stateIdx).stimSalience();
            currTypeWeight = EVCM.getCurrentTypeWeight();

            % TODO: normalization would interfere with controlSpotlight
            %normSalience = stimSalience / sum(stimSalience);
            normSalience = stimSalience;
            
            stimRespMap = EVCM.State.ExpectedSpace(stateIdx).getStimRespMap();

            % reduce stimulus-response-mapping to one dimension
            stimRespWeights = stimRespMap;
            stimRespWeights(:,2) = stimRespWeights(:,2) * -1;    % val response 2 < 0 < val response 1
            stimRespWeights = transpose(sum(stimRespWeights,2)); % represents how much a given stimulus pushes towards one response (due to stimulus-response mapping) irrespective of it's salience
           
            normSalience = normSalience + stimRespWeights * learningRate * (-1) * feedback * currTypeWeight;
            normSalience = max(0.00001,normSalience);
            % TODO: normalization would interfere with controlSpotlight
            %normSalience = normSalience / sum(normSalience);
            %normSalience = min(0.99, normSalience);
            normSalience = max(0.01, normSalience);
            EVCM.State.ExpectedSpace(stateIdx).setSaliency(normSalience);
            
            %EVC.HelperFnc.disp(source, 'actual prop2', actualProp2, 'expected prop2', expectedProp2, 'old expected saliency', stimSalience, 'new expected saliency', EVCM.State.ExpectedSpace(stateIdx).stimSalience); % DIAGNOSIS
           
            out = normSalience;
        end
        
        function out = simpleStateRL(this, EVCM)
            
            % EVCFnc parameters
            feedback = this.input{1}();
            learningRate = this.params{1};
            
            % parameters retrieved from EVC model
            stateIdx = EVCM.getExpectedStateIdx();
            
            EVCM.State.Expected = EVCM.State.ExpectedSpace(stateIdx);
            %% learn about reward contingencies & stimulus response mappings
            
            EVCM.State.ExpectedSpace(stateIdx).outcomeValues = EVCM.State.Expected.outcomeValues + [1-feedback feedback] .* (EVCM.State.Actual.outcomeValues - EVCM.State.ExpectedSpace(stateIdx).outcomeValues) * learningRate;
            
            currentTaskSet = EVCM.currentTaskSet(1);
            EVCM.State.ExpectedSpace(stateIdx).stimRespMap = EVCM.State.ExpectedSpace(stateIdx).stimRespMap + repmat([1-feedback feedback],size(currentTaskSet,1),1) .* (currentTaskSet - EVCM.State.ExpectedSpace(stateIdx).stimRespMap) * learningRate;
            
            EVCM.State.Expected = EVCM.State.ExpectedSpace(stateIdx);
            
            out(1).learned = EVCM.State.ExpectedSpace(stateIdx).outcomeValues;
            out(2).learned = EVCM.State.ExpectedSpace(stateIdx).stimRespMap;
        end
        
        function out = simpleStateRL_unweighted(this, EVCM)
            
            learningRate = this.params{1};
            
            % parameters retrieved from EVC model
            stateIdx = EVCM.getExpectedStateIdx();

            
            EVCM.State.Expected = EVCM.State.ExpectedSpace(stateIdx);
            %% learn about reward contingencies & stimulus response mappings
            
            EVCM.State.ExpectedSpace(stateIdx).outcomeValues = EVCM.State.Expected.outcomeValues + (EVCM.State.Actual.outcomeValues - EVCM.State.ExpectedSpace(stateIdx).outcomeValues) * learningRate;
            
            currentTaskSet = EVCM.currentTaskSet(1);
            EVCM.State.ExpectedSpace(stateIdx).stimRespMap = EVCM.State.ExpectedSpace(stateIdx).stimRespMap + (currentTaskSet - EVCM.State.ExpectedSpace(stateIdx).stimRespMap) * learningRate;
            
            EVCM.State.Expected = EVCM.State.ExpectedSpace(stateIdx);
            
            out(1).learned = EVCM.State.ExpectedSpace(stateIdx).outcomeValues;
            out(2).learned = EVCM.State.ExpectedSpace(stateIdx).stimRespMap;
        end
          
        function out = fillClearUpdate(this, EVCM)
            
            % parameters retrieved from EVC model
            stateIdx = EVCM.getExpectedStateIdx();
            signalIntensities = [this.input{1}() this.input{2}()]; %[EVCM.System.CtrlSignals(1).Intensity EVCM.System.CtrlSignals(2).Intensity];
            
            maxDiff = this.params{1};
            p = this.params{2};
            rewardIncr = this.params{3};
            jump = 0;
            defaultReward = EVCM.State.TaskEnv.Sequence(1).outcomeValues;
            currentReward = EVCM.State.TaskEnv.Sequence(EVCM.State.TaskEnv.CurrTrialIdx).outcomeValues;
            % switch set 
            if(abs(diff(currentReward)) ~= 0)

                if(abs(diff(currentReward)) >= maxDiff || rand < p)
                    newDiff = round(rand*maxDiff);
                                        
                    if rand < 0.5 %(diff(signalIntensities) < 0)
                        newReward = defaultReward + [0 newDiff];
                    else
                        newReward = defaultReward + [newDiff 0];
                    end
                    jump = 1;
                    if(EVCM.State.TaskEnv.CurrTrialIdx < length(EVCM.State.TaskEnv.Sequence))
                        EVCM.State.TaskEnv.Sequence(EVCM.State.TaskEnv.CurrTrialIdx+1).conditions = 1;
                    end
                end
                
            end
            
            % modify reward depending on choice 
            if(~jump)
                if diff(signalIntensities) < 0 % signal A chosen
                    rewardChange = [rewardIncr 0];
                    newReward = EVCM.State.TaskEnv.Sequence(EVCM.State.TaskEnv.CurrTrialIdx).outcomeValues + rewardChange;
                else    % signal B chosen
                    rewardChange = [0 rewardIncr];
                    newReward = EVCM.State.TaskEnv.Sequence(EVCM.State.TaskEnv.CurrTrialIdx).outcomeValues + rewardChange;
                end
            end
            
            if(EVCM.State.TaskEnv.CurrTrialIdx < length(EVCM.State.TaskEnv.Sequence))
                EVCM.State.TaskEnv.Sequence(EVCM.State.TaskEnv.CurrTrialIdx+1).outcomeValues = newReward;
                EVCM.State.ExpectedSpace(stateIdx).outcomeValues = newReward;
            end
            
            out(1).learned = EVCM.State.ExpectedSpace(stateIdx).outcomeValues;

            
        end
        
        function out = setForagingRewards(this, EVCM)
            
            % this value is 1 if the agent choose to harvest current patch (hit the upper
            % threshold)
            harvestChoice = ~round(this.input{1}());
            
            % these are the previous trial rewards associated with
            % harvesting and switching
            oldHarvestReward = EVCM.State.Actual.outcomeValues(1);
            oldSwitchReward = EVCM.State.Actual.outcomeValues(2);
            
            % calculate new rewards here:
            startVols = (30:15:150)*1;
            %disp(['TRIAL' num2str(EVCM.State.TaskEnv.CurrTrialIdx)]);
            %disp(['harvest choice was ' num2str(harvestChoice) '. Old reward was ' num2str(oldHarvestReward)]);
            if(~harvestChoice)
                newHarvestReward = randsample(startVols, 1);
            else
                p = log10(oldHarvestReward/0.060233);
                newHarvestReward = 10^(p-.06)*0.060233;
                newHarvestReward = max(newHarvestReward, 6);
            end
            newSwitchReward = mean(startVols);
            %disp(['New reward is ' num2str(newHarvestReward)]);
            
            % feed new calculated values to EVC model
            if(EVCM.State.TaskEnv.CurrTrialIdx < length(EVCM.State.TaskEnv.Sequence))
                EVCM.State.TaskEnv.Sequence(EVCM.State.TaskEnv.CurrTrialIdx+1).outcomeValues = [newHarvestReward newSwitchReward];
            end
            stateIdx = EVCM.getExpectedStateIdx();
            EVCM.State.ExpectedSpace(stateIdx).outcomeValues = [newHarvestReward newSwitchReward];
            EVCM.State.Expected = EVCM.State.ExpectedSpace(stateIdx);
            
            out(1).learned = newHarvestReward;
            out(2).learned = newSwitchReward;
        end
        
        function out = reinforceAutomaticPathway(this, EVCM)
            
            import EVC.*;
            
            % EVCFnc parameters
            feedback = this.input{1}();
            learningRate = this.params{1};
            
            % parameters retrieved from EVC model
            stateIdx = EVCM.getExpectedStateIdx();
            
            EVCM.State.Expected = EVCM.State.ExpectedSpace(stateIdx);
            %% learn about reward contingencies & stimulus response mappings
            
            currentTaskSet = EVCM.currentTaskSet(1);
            newStimRespMap = repmat([1-feedback feedback],size(currentTaskSet,1),1) .* repmat(EVCM.State.Actual.outcomeValues,size(currentTaskSet,1),1) .* currentTaskSet * learningRate;
            newSaliency = EVCM.State.Actual.stimSalience + transpose(sum(newStimRespMap,2));
            
            for i = 1:length(EVCM.State.ExpectedSpace)
            EVCM.State.ExpectedSpace(i).stimSalience = newSaliency;
            
            end
            
            for i = EVCM.State.TaskEnv.CurrTrialIdx:length(EVCM.State.TaskEnv.Sequence)
                EVCM.State.TaskEnv.Sequence(i).stimSalience = newSaliency;
            end
            
            EVCM.State.Expected = EVCM.State.ExpectedSpace(stateIdx);

            out(1).learned = newSaliency;
%             switch EVCM.State.Actual.conditions(1)
%                 case 0
%                     disp('congruent');
%                 case 1
%                     disp('incongruent');
%                 case 2
%                     disp('neutral');
%             end
%             disp('old Saliency:');
%             disp(EVCM.State.Actual.stimSalience);
%             disp('reward');
%             disp(EVCM.State.Actual.outcomeValues);
%             disp('feedback:')
%             disp([1-feedback feedback]);
%             disp('new Saliency:');
%             disp(newSaliency);
%             disp('distractor adjustment:');
%             disp(newSaliency(1)-EVCM.State.Actual.stimSalience(1));
%             disp('target adjustment:');
%             disp(newSaliency(2)-EVCM.State.Actual.stimSalience(2));
%             disp('--------------');
        end
        
        % CONTROL EFFICACY: changes efficacy parameter for each control
        % signal
        function out = controlEfficacy(this, EVCM)
            
            %% learn about automatic processing bias based on outcome differences between actual and expected state
            
            % EVCFnc parameters
            feedback = this.input{1}() - this.input{2}(); % actual error prob - expected error prob
            learningRate = this.params{1};
            
            DDMProcesses_expected = EVCM.DDMProcesses_expected;
            learnedControlEfficacy = [];
            
            % loop through each process
            for i = 1:length(DDMProcesses_expected)
                
                % only consider control processes since we are learning the
                % control efficacy
                if(ismember(DDMProc.CONTROL, DDMProcesses_expected(i).type))
                    
                    % implement different learning rules in dependence of
                    % DDM proxy
                    switch DDMProcesses_expected(i).DDMProxy
                        case DDMProc.DRIFT
                             expected_driftEffect = DDMProcesses_expected(i).getVal();
                             % determine direction of adjustment
                             if(expected_driftEffect > 0)
                                 adjustmentDirection = -1;
                             else
                                 adjustmentDirection = 1;
                             end

                        case DDMProc.THRESH
                             % for threshold adjustments there is one
                             % default direction (control increases
                             % threshold)
                             adjustmentDirection = -1;
                              
                        case DDMProc.BIAS
                             % not implemented
                             warning('Control efficacy learning not implementd for bias parameter.');
                        case DDMProc.NOISE
                             % not implemented
                             warning('Control efficacy learning not implementd for noise parameter.');
                        case DDMProc.T0
                             % not implemented
                    end
                    
                    % get applied control intensity
                    controlIntensity =  DDMProcesses_expected(i).input.params{1}.getVal();
                    
                    % get old control efficacy (parameter that
                    % scales the effect of the control signal on
                    % the stimulus)
                    old_efficacy = DDMProcesses_expected(i).input.params{3}.params{2};

                    % learning rule
                    new_efficacy = old_efficacy + adjustmentDirection * feedback * controlIntensity * learningRate;
                    
                    % make sure that efficacy stays > 1
                    new_efficacy = max([new_efficacy, 1e-16]);
                    
                    % apply updated control efficacy
                    learnedControlEfficacy(end+1) = new_efficacy;
                    DDMProcesses_expected(i).input.params{3}.params{2} = new_efficacy;
                end
                
            end
           
            out = learnedControlEfficacy;
        end
        
        function out = setGilzenrat2010Trial(this, EVCM)
            
            % retrieve performance (accuracy)
            error = round(this.input{1}());
            
            % determine which choice was selected
            if(strcmp(EVCM.State.Actual.descr, 'base'))
                explore = 1;
            elseif(strcmp(EVCM.State.Actual.descr, 'exploit'))
                explore = 0;
            else
                error('Description of chosen trials must either be ''base'' or ''exploit''.');
            end
            
            % get trial type indices
            baseTrialIdx = [];
            exploitTrialIdx = [];
            for i = 1:length(EVCM.State.TaskEnv.trialTypes)
                if(strcmp(EVCM.State.TaskEnv.trialTypes{i}.descr, 'base'))
                    baseTrialIdx = i;
                end
                if(strcmp(EVCM.State.TaskEnv.trialTypes{i}.descr, 'exploit'))
                    exploitTrialIdx = i;
                end
            end
            if(isempty(baseTrialIdx) || isempty(exploitTrialIdx))
                error('Description of chosen trials must either be ''base'' or ''exploit''.');
            end
            
            % log points offered on exploit trial
            if(~isfield(EVCM.Log, 'offeredPoints'))
                EVCM.Log.offeredPoints = nan(length(EVCM.State.TaskEnv.Sequence),1);
            end
            EVCM.Log.offeredPoints(EVCM.Log.LogCount, 1) = EVCM.State.TaskEnv.trialTypes{exploitTrialIdx}.outcomeValues(1);
%             disp(EVCM.Log.offeredPoints(EVCM.Log.LogCount, 1));
            
            % if subject selected base trial, need to reset the exploit trial
            if(explore)
                
                EVCM.State.TaskEnv.trialTypes{exploitTrialIdx}.outcomeValues = EVCM.State.TaskEnv.trialTypes{baseTrialIdx}.outcomeValues;
                EVCM.State.TaskEnv.trialTypes{exploitTrialIdx}.stimRespMap(1,1) = EVCM.State.TaskEnv.trialTypes{baseTrialIdx}.stimRespMap(1,1);
                
            end
            
            % if subject made no error, then increase reward (oldReward + 5) and difficulty (oldStimSalience / 2)
            if(error == 0)
                EVCM.State.TaskEnv.trialTypes{exploitTrialIdx}.outcomeValues(1) = EVCM.State.TaskEnv.trialTypes{exploitTrialIdx}.outcomeValues(1) + 5;
                EVCM.State.TaskEnv.trialTypes{exploitTrialIdx}.stimRespMap(1,1) = EVCM.State.TaskEnv.trialTypes{exploitTrialIdx}.stimRespMap(1,1)/2;

                % if minimum difficulty is reached, set it to maximum (stimSalience = 0) after 9 trials
                if(EVCM.State.TaskEnv.trialTypes{exploitTrialIdx}.stimRespMap(1,1) <= (EVCM.State.TaskEnv.trialTypes{baseTrialIdx}.stimRespMap(1,1) / 2^9))
                    EVCM.State.TaskEnv.trialTypes{exploitTrialIdx}.stimRespMap(1,1) = 0;
                end
            else % if subject made an error, then keep difficulty but decrease reward
                EVCM.State.TaskEnv.trialTypes{exploitTrialIdx}.outcomeValues(1) = max(EVCM.State.TaskEnv.trialTypes{exploitTrialIdx}.outcomeValues(1) - 10, 0);
            end
            
%             disp([num2str(error)]);
%             disp('resulting outcomeval:')
%             EVCM.State.TaskEnv.trialTypes{exploitTrialIdx}.outcomeValues
%             EVCM.State.TaskEnv.trialTypes{exploitTrialIdx}.stimRespMap(1,1)
%             disp('---')

%             if(explore)
%                 disp('EXPLORED');
%             else 
%                 disp('exploit');
%             end
            
            out = 0;
            
        end
        
        function out = setApproachAvoidanceTrial(this, EVCM)
            
            minimumSaliency = 0.5;
            stepSize = 0.01;
            minimumTrial = 10;
            
            % retrieve performance (accuracy)
            error = round(this.input{1}());
            
            % determine which choice was selected
            if(strcmp(EVCM.State.Actual.descr, 'taskA'))
                taskA = 1;
            elseif(strcmp(EVCM.State.Actual.descr, 'taskB'))
                taskA = 0;
            else
                error('Description of chosen trials must either be ''taskA'' or ''taskB''.');
            end
            
            % get trial type indices
            taskATrialIdx = [];
            taskBTrialIdx = [];
            for i = 1:length(EVCM.State.TaskEnv.trialTypes)
                if(strcmp(EVCM.State.TaskEnv.trialTypes{i}.descr, 'taskA'))
                    taskATrialIdx = i;
                end
                if(strcmp(EVCM.State.TaskEnv.trialTypes{i}.descr, 'taskB'))
                    taskBTrialIdx = i;
                end
            end
            if(isempty(taskATrialIdx) || isempty(taskBTrialIdx))
                error('Description of chosen trials must either be ''taskA'' or ''taskB''.');
            end
            
            % increase difficulty of task A
            if(isfield(EVCM.State, 'CurrTrlIdx'))
                if(EVCM.State.CurrTrlIdx > minimumTrial)
                    EVCM.State.TaskEnv.trialTypes{taskATrialIdx}.stimSalience = [max([(EVCM.State.TaskEnv.trialTypes{taskATrialIdx}.stimSalience(1)-stepSize) minimumSaliency]) 1];
                end
            end
            
            out = 0;
            
        end
        
        function out = setBraunArrington2018Trial(this, EVCM)
            
            minReward = 1;
            maxReward = 10;
            defaultReward = 5;
            
            % determine which choice was selected
            if(strcmp(EVCM.State.Actual.descr, 'taskA'))
                taskA = 1;
            elseif(strcmp(EVCM.State.Actual.descr, 'taskB'))
                taskA = 0;
            else
                error('Description of chosen trials must either be ''taskA'' or ''taskB''.');
            end
            
            % retrieve performance (accuracy)
            if(taskA)
                MadeAnError = 1-EVCM.State.Actual.performance.probs(1);
            else
                MadeAnError = 1-EVCM.State.Actual.performance.probs(2);
            end
            
            % get trial type indices
            taskATrialIdx = [];
            taskBTrialIdx = [];
            for i = 1:length(EVCM.State.TaskEnv.trialTypes)
                if(strcmp(EVCM.State.TaskEnv.trialTypes{i}.descr, 'taskA'))
                    taskATrialIdx = i;
                end
                if(strcmp(EVCM.State.TaskEnv.trialTypes{i}.descr, 'taskB'))
                    taskBTrialIdx = i;
                end
            end
            if(isempty(taskATrialIdx) || isempty(taskBTrialIdx))
                error('Description of chosen trials must either be ''taskA'' or ''taskB''.');
            end
            
            
            % log block
            if(~isfield(EVCM.Log, 'block'))
                EVCM.Log.block = ones(length(EVCM.State.TaskEnv.Sequence),1);
            end
            if(EVCM.Log.LogCount > 1 & EVCM.Log.block(EVCM.Log.LogCount) <= 1)
                EVCM.Log.block(EVCM.Log.LogCount) = EVCM.Log.block(EVCM.Log.LogCount-1);
            end
            
            % log block trial
            if(~isfield(EVCM.Log, 'blockTrial'))
                EVCM.Log.blockTrial = ones(length(EVCM.State.TaskEnv.Sequence),1);
            end
            if(EVCM.Log.LogCount > 1)
                if(EVCM.Log.block(EVCM.Log.LogCount-1) ~= EVCM.Log.block(EVCM.Log.LogCount))
                    EVCM.Log.blockTrial(EVCM.Log.LogCount) = 1;
                else
                    EVCM.Log.blockTrial(EVCM.Log.LogCount) = EVCM.Log.blockTrial(EVCM.Log.LogCount-1)+1;
                end
            end
            
             % log task choice
            if(~isfield(EVCM.Log, 'taskAchosen'))
                EVCM.Log.taskAchosen = nan(length(EVCM.State.TaskEnv.Sequence),1);
            end
            EVCM.Log.taskAchosen(EVCM.Log.LogCount) = taskA;
            if(~isfield(EVCM.Log, 'taskBchosen'))
                EVCM.Log.taskBchosen = nan(length(EVCM.State.TaskEnv.Sequence),1);
            end
            EVCM.Log.taskBchosen(EVCM.Log.LogCount) = ~taskA;
            
            % log offered points
            if(~isfield(EVCM.Log, 'offeredPoints'))
                EVCM.Log.offeredPoints = nan(length(EVCM.State.TaskEnv.Sequence),2);
            end
            EVCM.Log.offeredPoints(EVCM.Log.LogCount, 1:2) = [EVCM.State.TaskEnv.trialTypes{taskATrialIdx}.outcomeValues(1) EVCM.State.TaskEnv.trialTypes{taskBTrialIdx}.outcomeValues(2)];
            
            % log offered points relative to current choice
            if(~isfield(EVCM.Log, 'currentValue'))
                EVCM.Log.currentValue = nan(length(EVCM.State.TaskEnv.Sequence),1);
            end
            if(taskA)
                EVCM.Log.currentValue(EVCM.Log.LogCount) = EVCM.Log.offeredPoints(EVCM.Log.LogCount, 1);
            else
                EVCM.Log.currentValue(EVCM.Log.LogCount) = EVCM.Log.offeredPoints(EVCM.Log.LogCount, 2);
            end
            
            if(~isfield(EVCM.Log, 'otherValue'))
                EVCM.Log.otherValue = nan(length(EVCM.State.TaskEnv.Sequence),1);
            end
            if(taskA)
                EVCM.Log.otherValue(EVCM.Log.LogCount) = EVCM.Log.offeredPoints(EVCM.Log.LogCount, 2);
            else
                EVCM.Log.otherValue(EVCM.Log.LogCount) = EVCM.Log.offeredPoints(EVCM.Log.LogCount, 1);
            end
            
            % log offered points relative to last choice
            if(~isfield(EVCM.Log, 'currentValueOfLastChoice'))
                EVCM.Log.currentValueOfLastChoice = nan(length(EVCM.State.TaskEnv.Sequence),1);
            end
            
            if(~isfield(EVCM.Log, 'otherValueOfLastChoice'))
                EVCM.Log.otherValueOfLastChoice = nan(length(EVCM.State.TaskEnv.Sequence),1);
            end
            
            if(EVCM.Log.LogCount > 1)
                taskAprev = EVCM.Log.taskAchosen(EVCM.Log.LogCount-1);

                if(taskAprev)
                    EVCM.Log.currentValueOfLastChoice(EVCM.Log.LogCount) = EVCM.Log.offeredPoints(EVCM.Log.LogCount, 1);
                else
                    EVCM.Log.currentValueOfLastChoice(EVCM.Log.LogCount) = EVCM.Log.offeredPoints(EVCM.Log.LogCount, 2);
                end

                if(taskAprev)
                    EVCM.Log.otherValueOfLastChoice(EVCM.Log.LogCount) = EVCM.Log.offeredPoints(EVCM.Log.LogCount, 2);
                else
                    EVCM.Log.otherValueOfLastChoice(EVCM.Log.LogCount) = EVCM.Log.offeredPoints(EVCM.Log.LogCount, 1);
                end
            end
            
            % log accumulated points
            if(~isfield(EVCM.Log, 'accumulatedPoints'))
                EVCM.Log.accumulatedPoints = zeros(length(EVCM.State.TaskEnv.Sequence),1);
            end
            if(EVCM.Log.blockTrial(EVCM.Log.LogCount) > 1)
                EVCM.Log.accumulatedPoints(EVCM.Log.LogCount) = EVCM.Log.accumulatedPoints(EVCM.Log.LogCount-1);
            end
            if(~MadeAnError)
                EVCM.Log.accumulatedPoints(EVCM.Log.LogCount) = EVCM.Log.accumulatedPoints(EVCM.Log.LogCount) + EVCM.Log.currentValue(EVCM.Log.LogCount);
            end
            
            % reset block or adjust points for each task
            resetBlock = 0;
            if(EVCM.Log.accumulatedPoints(EVCM.Log.LogCount) > 500)
                resetBlock = 1;
                if(EVCM.Log.LogCount < length(EVCM.Log.accumulatedPoints))
                    EVCM.Log.accumulatedPoints(EVCM.Log.LogCount+1) = 0;
                end
            end

            
            % if reset, then reset block number to 1
            if(resetBlock)
                if(EVCM.Log.LogCount  < length(EVCM.Log.block))
                    EVCM.Log.block(EVCM.Log.LogCount+1) = EVCM.Log.block(EVCM.Log.LogCount)+1;
                end
             end
            
            % log if value of current task increased
            if(~isfield(EVCM.Log, 'otherValueIncrease'))
                EVCM.Log.otherValueIncrease = nan(length(EVCM.State.TaskEnv.Sequence),1);
            end
            if(EVCM.Log.LogCount > 1)
                taskAprev = EVCM.Log.taskAchosen(EVCM.Log.LogCount-1);
                if(taskAprev == 1)
                    if(EVCM.Log.offeredPoints(EVCM.Log.LogCount-1, 2) < EVCM.Log.offeredPoints(EVCM.Log.LogCount, 2))
                        EVCM.Log.otherValueIncrease(EVCM.Log.LogCount) = 1;
                    else
                        EVCM.Log.otherValueIncrease(EVCM.Log.LogCount) = 0;
                    end
                else
                    if(EVCM.Log.offeredPoints(EVCM.Log.LogCount-1, 1) < EVCM.Log.offeredPoints(EVCM.Log.LogCount, 1))
                        EVCM.Log.otherValueIncrease(EVCM.Log.LogCount) = 1;
                    else
                        EVCM.Log.otherValueIncrease(EVCM.Log.LogCount) = 0;
                    end
                end
            end
            
            % log if value of other task decreased
            if(~isfield(EVCM.Log, 'currentValueDecrease'))
                EVCM.Log.currentValueDecrease = nan(length(EVCM.State.TaskEnv.Sequence),1);
            end
            if(EVCM.Log.LogCount > 1)
                taskAprev = EVCM.Log.taskAchosen(EVCM.Log.LogCount-1);
                if(taskAprev == 1)
                    if(EVCM.Log.offeredPoints(EVCM.Log.LogCount-1, 1) > EVCM.Log.offeredPoints(EVCM.Log.LogCount, 1))
                        EVCM.Log.currentValueDecrease(EVCM.Log.LogCount) = 1;
                    else
                        EVCM.Log.currentValueDecrease(EVCM.Log.LogCount) = 0;
                    end
                else
                    if(EVCM.Log.offeredPoints(EVCM.Log.LogCount-1, 2) > EVCM.Log.offeredPoints(EVCM.Log.LogCount, 2))
                        EVCM.Log.currentValueDecrease(EVCM.Log.LogCount) = 1;
                    else
                        EVCM.Log.currentValueDecrease(EVCM.Log.LogCount) = 0;
                    end
                end
            end
            
            % log EVC of switch option
            if(~isfield(EVCM.Log, 'switchEVC'))
                EVCM.Log.switchEVC = nan(length(EVCM.State.TaskEnv.Sequence),1);
            end
            if(EVCM.Log.LogCount > 1)
                taskAprev = EVCM.Log.taskAchosen(EVCM.Log.LogCount-1);
                if(taskAprev)
                    EVCM.Log.switchEVC(EVCM.Log.LogCount) = EVCM.Log.stateEVCs(EVCM.Log.LogCount,2) - EVCM.Log.stateEVCs(EVCM.Log.LogCount,1);
                else
                    EVCM.Log.switchEVC(EVCM.Log.LogCount) = EVCM.Log.stateEVCs(EVCM.Log.LogCount,1) - EVCM.Log.stateEVCs(EVCM.Log.LogCount,2);
                end
            end
            
            % log if task just switched or repeated
            if(~isfield(EVCM.Log, 'switchTrial'))
                EVCM.Log.switchTrial = nan(length(EVCM.State.TaskEnv.Sequence),1);
            end
            if(EVCM.Log.LogCount > 1)
                if(taskA)
                    if(strcmp(EVCM.Log.ExpectedState(EVCM.Log.LogCount-1).descr, 'taskA'))
                        EVCM.Log.switchTrial(EVCM.Log.LogCount) = 0;
                    else
                        EVCM.Log.switchTrial(EVCM.Log.LogCount) = 1;
                    end
                else
                    if(strcmp(EVCM.Log.ExpectedState(EVCM.Log.LogCount-1).descr, 'taskB'))
                        EVCM.Log.switchTrial(EVCM.Log.LogCount) = 0;
                    else
                        EVCM.Log.switchTrial(EVCM.Log.LogCount) = 1;
                    end
                end
                
            end

            % set points for new trial
            
            if(resetBlock)
                
                EVCM.State.TaskEnv.trialTypes{taskATrialIdx}.outcomeValues = [defaultReward 0];
                EVCM.State.TaskEnv.trialTypes{taskBTrialIdx}.outcomeValues = [0 defaultReward];
                
            else
            
                if(rand < 0.5)
                    taskAChange = 1;
                else
                    taskAChange = 0;
                end

                if(rand < 0.5)
                    taskBChange = 1;
                else
                    taskBChange = 0;
                end

                if(taskA)
                    EVCM.State.TaskEnv.trialTypes{taskATrialIdx}.outcomeValues = [max([EVCM.State.TaskEnv.trialTypes{taskATrialIdx}.outcomeValues(1)-taskAChange  minReward]), 0];
                    EVCM.State.TaskEnv.trialTypes{taskBTrialIdx}.outcomeValues = [0, min([EVCM.State.TaskEnv.trialTypes{taskBTrialIdx}.outcomeValues(2) + taskBChange, maxReward])];
                else
                    EVCM.State.TaskEnv.trialTypes{taskATrialIdx}.outcomeValues = [min([EVCM.State.TaskEnv.trialTypes{taskATrialIdx}.outcomeValues(1)+taskAChange  maxReward]), 0];
                    EVCM.State.TaskEnv.trialTypes{taskBTrialIdx}.outcomeValues = [0, max([EVCM.State.TaskEnv.trialTypes{taskBTrialIdx}.outcomeValues(2) - taskBChange, minReward])];
                end
            
            end
            out = 0;
            
            % if end of 6th block is reached terminate experiment 
            if(resetBlock & EVCM.Log.block(EVCM.Log.LogCount) == 6)
                EVCM.State.terminate = 1;
                disp('TERMINATED EXPERIMENT (all 6 blocks completed)');
            end
            
            disp(['Trial ' num2str(EVCM.Log.blockTrial(EVCM.Log.LogCount)) ' of block ' num2str(EVCM.Log.block(EVCM.Log.LogCount)) ...
                    ', accumulated points: ' num2str(EVCM.Log.accumulatedPoints(EVCM.Log.LogCount))]);
            
        end
        
        function out = setWestbrookBraver2015Trial(this, EVCM)
            
            maximumReward = 200;
            stepSize = 1;
            minimumTrial = 1;
            
            % retrieve performance (accuracy)
            error = round(this.input{1}());
            
            % determine which choice was selected
            if(strcmp(EVCM.State.Actual.descr, 'taskA'))
                taskA = 1;
            elseif(strcmp(EVCM.State.Actual.descr, 'taskB'))
                taskA = 0;
            else
                error('Description of chosen trials must either be ''taskA'' or ''taskB''.');
            end
            
            % get trial type indices
            taskATrialIdx = [];
            taskBTrialIdx = [];
            for i = 1:length(EVCM.State.TaskEnv.trialTypes)
                if(strcmp(EVCM.State.TaskEnv.trialTypes{i}.descr, 'taskA'))
                    taskATrialIdx = i;
                end
                if(strcmp(EVCM.State.TaskEnv.trialTypes{i}.descr, 'taskB'))
                    taskBTrialIdx = i;
                end
            end
            if(isempty(taskATrialIdx) || isempty(taskBTrialIdx))
                error('Description of chosen trials must either be ''taskA'' or ''taskB''.');
            end
            
            % increase reward for task B
            if(isfield(EVCM.State, 'CurrTrlIdx'))
                if(EVCM.State.CurrTrlIdx > minimumTrial)
                    EVCM.State.TaskEnv.trialTypes{taskBTrialIdx}.outcomeValues = [0 min([(EVCM.State.TaskEnv.trialTypes{taskBTrialIdx}.outcomeValues(2)+stepSize) maximumReward])];
                end
            end
            
            out = 0;
            
        end
        
    end
    
end