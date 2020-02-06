% EVC Sim created by A.Shenhav on 3/20/13
% extended by S. Musslick on 12/11/14

classdef EVCModel < handle
    properties (SetAccess = public)
        System                % holds possible control signals, indicies of current control signals
        State                 % represents the actual world & the control system's representation of the actual world (implementation specific)
        learning              % holds learning functions
        rewardFnc             % EVCFnc: represents reward function
        Log                   % log structure
        
        useActualState        % flag if using actual (vs. expected) state for EVC calculations
        RTscale               % degree to which EVC is normalized by RT
        RTconst               % constant added to RTscale*RT -> EVC = EVC / (RTconst + RTscale*RT)
        binaryErrors;         % retrieve binary errors instead of error probabilities?
        trialHistoryFilter      % EVCFnc: function determines weight of a trial based on it's history
        id                          % model ID
    end
    
    properties (Constant)
        defRTscale = 0;
        defRTconst = 1;
        defUseActualState = 1;
        defCapacity = 1e16;
        defID = 0;
        startBinaryErrors = 0;
    end
    
    methods 
       
        % EVCModel constructor
        function this = EVCModel(CtrlSignals, TaskEnv)
            
            % check if dummy instance
            if(isnumeric(CtrlSignals) && isnumeric(TaskEnv))
                if(CtrlSignals == 0 && TaskEnv == 0)
                   return; 
                end
            end
            
            % check input arguments
            if(isa(CtrlSignals, 'EVC.CtrlSignal') == 0)
               error('Please make sure that the control signal array is an instance of the class ''EVC.CtrlSignal''.'); 
            end
            if(isa(TaskEnv, 'EVC.TaskEnv') == 0)
               error('Please make sure that the task sequence is an instance of the class ''EVC.TaskEnv''.'); 
            end
            
            % assign control signal values
            if(length(CtrlSignals) >= 1)
                this.System.CtrlSignals = CtrlSignals;
                this.System.CurrSigIdx = 1:length(this.System.CtrlSignals);
            else
                error('CtrlSignal array has to contain at least one control signal');
            end

            this.State.TaskEnv = TaskEnv;
            this.State.terminate = 0;

            this.System.reconfigurationCost = [];
            this.learning = EVC.LearningFnc.empty(1,0);
            this.System.capacity = this.defCapacity;
            this.RTscale = this.defRTscale;
            this.RTconst = this.defRTconst;
            this.useActualState = this.defUseActualState;
            this.binaryErrors = this.startBinaryErrors;
            this.id = this.defID;
            this.System.ACC_HypotheticalBOLD = 0;
            this.initLog();
            
            % assign default reward function
            this.rewardFnc = EVC.EVCFnc(EVC.EVCFnc.REWRATE, [1 1], this);
            
            % assign default trial history filter
            this.trialHistoryFilter = EVC.EVCFnc(EVC.EVCFnc.VALUE, 1, this);
            
            % index control signals & link this EVCModel instance to ctrl signals
            maxIndex = max([this.System.CtrlSignals(:).id]);
            maxIndex = EVC.HelperFnc.ifelse(maxIndex == EVC.CtrlSignal.defID, 0, maxIndex);
            counter = 1;
            for idx = 1:length(this.System.CtrlSignals)
               
                if(this.System.CtrlSignals(idx).id == EVC.CtrlSignal.defID)
                    this.System.CtrlSignals(idx).id = counter + maxIndex;
                    counter = counter+1;
                end
               
                this.System.CtrlSignals(idx).EVCModel = this;
            end

          % set actual and expected state
          this.setActualState(TaskEnv);
          this.initExpectedStates();   
          this.State.Expected = this.getExpectedState();
        end
        
    end
    
%% Core EVC implementation

    methods
        
        % returns signal combination with maximal EVC 
        % signal* <- max(EVC(signal_i, state))
        function [optIntensities optEVC optStateID signalMap EVCMap optPerformance optEV optCost stateData] = getOptSignals(this)
            
            [ST,I] = dbstack; % track function name for debugging functions
            source = ST.name;
             
            %% first create all possible intensity tubles
            
            % loop through all current control signals
            currentSigIdx = nan;

            for idx = 1:length(this.System.CurrSigIdx)
                currentSigIdx = this.System.CurrSigIdx(idx);
                % build array of all possible intensity combinations
                if(idx == 1)
                    signalCombs = this.System.CtrlSignals(currentSigIdx).IntensityRange;
                else
                    signalCombs = Simulations.combvec(signalCombs, this.System.CtrlSignals(currentSigIdx).IntensityRange);
                end                    
            end

            % find all combinations for which sum of signal intensities
            % exceeds capacity limit...
            deleteIdx = find(sum(signalCombs, 1) > this.System.capacity);
            
            % ... and delete them from the signal search space
            signalCombs(:, deleteIdx) = [];
            
            
            %% calc EVCs
            
            % backup current intensities of current control signals
            backupIntensities = this.getIntensities();

            this.System.previousIntensities = backupIntensities;
                        
            % set intensity search range for each signal
            for idx = 1:length(this.System.CurrSigIdx)
                currentSigIdx = this.System.CurrSigIdx(idx);
                % set new intensity to search range
                this.System.CtrlSignals(currentSigIdx).Intensity = signalCombs(idx, :);
            end
            
            
            %% simulation across expected states
            
            % CHOICE trials: select state + control signal
            if(length(this.State.Actual.cueID) > 1)
                selectState = 1;
            else
            % CUED trials: select control signal
                selectState = 0;
                optStateID = [];
            end
                
            % retrieve Cue: which expected states are relevant?
            expectedStateIdx = this.State.expectedStatesIdx;
            
            if(~selectState)
                % get trial weights based on recent frequency
                typeWeights = this.getTrialTypeWeights(this.State.Actual.cueID);
            end
                
            expectedStateBackup = this.State.Expected;
            useActualStateBackup = this.useActualState;
            
            this.useActualState = 0;
            
            for i = 1:length(expectedStateIdx)
                
                % set current expected state
                this.State.Expected = this.State.ExpectedSpace(expectedStateIdx(i));
            
                % retrieve array of EVCs for all possible signal combinations
                [overallStateEVC,overallStateCost,overallStateEV,overallStatePerformance] = getEVC(this);
                
                if(selectState)
                    % if choice trial, store state EVC's separately
                    stateEVCs(i,:) = overallStateEVC;
                    stateCosts(i,:) = overallStateCost;
                    stateEVs(i,:) = overallStateEV;
                    statePerformance(i) = overallStatePerformance;

                else
                    currentType = this.State.ExpectedSpace(expectedStateIdx(i)).typeID;
                    % else sum EVC across all states
                    if(i == 1) 

                        overallEVC = overallStateEVC;
                        overallEV = overallStateEV;
                        overallCost = overallStateCost;
                        overallPerformance = EVC.HelperFnc.weightStruct(overallStatePerformance, typeWeights(2,typeWeights(1,:)==currentType));
                    else
                        
                        overallEVC = overallEVC + overallStateEVC * typeWeights(2,typeWeights(1,:)==currentType);
                        overallEV = overallEV + overallStateEV * typeWeights(2,typeWeights(1,:)==currentType);
                        overallCost = overallCost + overallStateCost * typeWeights(2,typeWeights(1,:)==currentType);
                        overallPerformance = EVC.HelperFnc.addStruct(overallPerformance, overallStatePerformance, typeWeights(2,typeWeights(1,:)==currentType));

                    end
                end

            end

            %% find optimal signal

            % TODO: this may be the place to add some constraints assuming
            % non-optimality of the system
            % ALSO: think about constraint on number of signals -> search
            % which search space to use

            % get indicies of optimal signal combinations
            if(selectState)
                optStatesIdx = find(max(stateEVCs,[],2) == max(max(stateEVCs,[],2)));
                
                % if multiple optimal states, then select the state that
                % demands the least control adjustment
                for i = 1:length(optStatesIdx)
                    combIdx = find(stateEVCs(optStatesIdx(i),:) == max(stateEVCs(optStatesIdx(i),:)));
                    stateDistances = sum(abs(signalCombs(:,combIdx) - repmat(backupIntensities,1,size(combIdx,2))),1);
                    minDistance(i) = min(stateDistances);
                end
                
                optStatesIdx = optStatesIdx(find(minDistance == min(minDistance)));
                % select random state if still multiple equal options
                optStateIdx = optStatesIdx(randi(numel(optStatesIdx)));
                optTrialTypeID = this.State.ExpectedSpace(expectedStateIdx(optStateIdx)).typeID;
                
                optStateID = [];
                % get optimal state ID (corresponding to the index of this.State.TaskEnv.trialTypes)
                for trialTypeIdx = 1:length(this.State.TaskEnv.trialTypes)
                    if(this.State.TaskEnv.trialTypes{trialTypeIdx}.typeID == optTrialTypeID)
                        optStateID = trialTypeIdx;
                    end
                end
                if(isempty(optStateID))
                    Error('Could not find trial type index for this.State.TaskEnv.trialTypes that matches optimal typeID.');
                end
                
                combIdx = find(stateEVCs(optStateIdx,:) == max(stateEVCs(optStateIdx,:)));
            else
                combIdx = find(overallEVC == max(overallEVC));
            end

            if(isempty(combIdx))       
                warning('Did not find optimal signal (possibly because of NaN values). Will stick with default allocation, selecting first signal.');
                combIdx = 1;
            end
            
            % choose the combination with least distance to current signal combination
            distances = sum(abs(signalCombs(:,combIdx) - repmat(backupIntensities,1,size(combIdx,2))),1);
            combIdxIdx = find(distances == min(distances));

            % select random signal combination if still multiple equal options
            combIdxIdx = combIdxIdx(randi(numel(combIdxIdx)));

            optIntensities = signalCombs(:,combIdx(combIdxIdx));

            % extract optimal outcome values
            if(selectState)
                optEVC = stateEVCs(optStateIdx,combIdx(combIdxIdx));
                optEV = stateEVs(optStateIdx,combIdx(combIdxIdx));
                optCost = stateCosts(optStateIdx,combIdx(combIdxIdx));
                EVCMap = stateEVCs(optStateIdx,:);

                performanceNames = fieldnames(statePerformance(optStateIdx));
                optPerformance = statePerformance(optStateIdx);

                for i = 1:length(performanceNames)
                    eval(strcat('optPerformance.',performanceNames{i},' = statePerformance(optStateIdx).',performanceNames{i},'(:,combIdx(combIdxIdx));'));
                end

                overallEVC = stateEVCs;
            else
                optEVC = overallEVC(combIdx(combIdxIdx));
                optEV = overallEV(combIdx(combIdxIdx));
                optCost = overallCost(combIdx(combIdxIdx));
                EVCMap = overallEVC;

                performanceNames = fieldnames(overallPerformance);
                optPerformance = overallPerformance;

                for i = 1:length(performanceNames)
                    eval(strcat('optPerformance.',performanceNames{i},' = overallPerformance.',performanceNames{i},'(:,combIdx(combIdxIdx));'));
                end
     
            end

            this.State.Expected = expectedStateBackup;
            this.useActualState = useActualStateBackup;
            
            %% calculate ACC BOLD
            
            if(selectState)
                % TODO
            else
                this.System.ACC_HypotheticalBOLD = this.getHypotheticalBOLD(overallEVC, optIntensities, this.getCCost(), this.getReconfCost(), signalCombs, combIdx(combIdxIdx), backupIntensities);
            end
            
            %% store state EVCs
            if(selectState)
                
                stateData.stateEVCs = nan(1, length(expectedStateIdx));
                stateData.stateCosts = nan(1, length(expectedStateIdx));
                stateData.stateEVs = nan(1, length(expectedStateIdx));
                
                for expStateIdx = 1:length(expectedStateIdx)
                    bestEVCIdx = find(stateEVCs(expStateIdx, :) == max(stateEVCs(expStateIdx, :)));
                    stateData.stateEVCs(expStateIdx) = stateEVCs(expStateIdx, bestEVCIdx(1));
                    stateData.stateCosts(expStateIdx) = stateCosts(expStateIdx, bestEVCIdx(1));
                    stateData.stateEVs(expStateIdx) = stateEVs(expStateIdx, bestEVCIdx(1));
                    stateData.statePerformance{expStateIdx} = statePerformance;
                end
            else 
                stateData.stateEVCs = [];
                stateData.stateCosts = [];
                stateData.stateEVs = [];
                stateData.statePerformance = [];
            end
            
            
            %% restore intensities
            this.setIntensities(backupIntensities);
            
            signalMap = signalCombs;
            
            EVC.HelperFnc.disp(source, 'optimal EVC', optEVC, 'optimal intensities', optIntensities, 'overallEVC', overallEVC, 'signalCombs', signalCombs); % DIAGNOSIS
            
            %% warn user if control signal limit is reached
            
            for idx = 1:length(this.System.CurrSigIdx)
                currentSigIdx = this.System.CurrSigIdx(idx);
                % build array of all possible intensity combinations
                if(optIntensities(idx) == this.System.CtrlSignals(currentSigIdx).IntensityRange(end))
                    warning(['Reached upper limit of control signal range when looking for optimal intensities of signal ' num2str(idx)]);
                end
            end
        end
        
        % returns the expected value of control based on the current
        % environmental state and the system's current control signal(s)
        function [EVC_val Cost expectedValue performance] = getEVC(this)
            [ST,I] = dbstack; % track function name for debugging functions
            source = ST.name;
             
            import EVC.*;

            % calculate outcomes from child class (DDM, MSDDM, neural network simulation, etc.)
            [performance] = this.simulateOutcomes();
            
            % calculate expected value
            expectedValue = this.rewardFnc.getVal(performance);

            % substract Control Costs
            Cost = this.getCCost() + this.getReconfCost();
            
            % expected value of control
            EVC_val = expectedValue - Cost; % expected value of control

            HelperFnc.disp(source, 'EVC', EVC_val, 'outcomes', performance, 'Costs'); % DIAGNOSIS
            
        end
       
        % returns control costs for current control identities
       function costs = getCCost(this)
           
           % control signal costs
           costs = zeros(size(this.System.CtrlSignals(1).Intensity));
           
            for idx = 1:length(this.System.CurrSigIdx)
                
                costs = costs + this.System.CtrlSignals(idx).getCCost();
                
            end
           
       end
       
       function reconfigurationCosts = getReconfCost(this)
           
            % reconfiguration costs
            
            if(~isempty(this.System.reconfigurationCost))
                
                if(isfield(this.System, 'previousIntensities'))
                % calculate euklidian distance
                newIntensities = this.getIntensities();
                signalDistance = sqrt(sum((newIntensities-repmat(this.System.previousIntensities,1,size(newIntensities,2))).^2,1));
                
                % get maximal euklidian distance for given signal space
                for idx = 1:length(this.System.CurrSigIdx)
                    currentSigIdx = this.System.CurrSigIdx(idx);
                    rangeSpan(idx) = max(this.System.CtrlSignals(currentSigIdx).IntensityRange)-min(this.System.CtrlSignals(currentSigIdx).IntensityRange);                   
                end
                maxDistance = sqrt(sum(rangeSpan.^2));
                
                % normalize current distance
                signalDistance = signalDistance / maxDistance;
                
                reconfigurationCosts = this.System.reconfigurationCost.getVal(signalDistance);
                
                else
                    reconfigurationCosts = 0; % no costs if no origin
                end

            end
           
       end
       
       % returns control costs for current control identities
       function costs = getCCost_old(this)
            
           costs = [];
           % add control costs for all current control identities
           for CSidx = 1:length(this.System.CurrSigIdx)
                costs(CSidx,:) = this.System.CtrlSignals(this.System.CurrSigIdx(CSidx)).getCCost();
           end
           
           % sum costs across control signals
           costs = sum(costs,1);
            
       end
       
       % sets trial history finction
       function setTrialHistoryFilter(this, historyFilter)
        if(isa(historyFilter, 'EVC.EVCFnc') == 0)
           error('Please make sure that the first argument is an instance of the class ''EVC.EVCFnc''.'); 
        end
           this.trialHistoryFilter = historyFilter;
       end
       
       % calculate ACC BOLD response
       %
       % Note that function is not meant to reflect a specific prediction
       % of the EVC theory as Shenhav et al. (2013) do not stipulate a relationship
       % between EVC variables and dACC activity. However, the user may use this
       % function to identify neural correlates of various EVC variables.
       function hypotheticalBOLD = getHypotheticalBOLD(this, EVCMap, optIntensities, signalCosts, reconfigurationCosts, signalCombs, optIdx, backupIntensities)
           [ST,I] = dbstack; % track function name for debugging functions
            source = ST.name;
           
           % 1) take distance from current EVCMap to recent EVCMap
           if(isfield(this.System, 'EVCMap'))
               if(numel(this.System.EVCMap) == numel(EVCMap))
                    diffEVCMap = abs(EVCMap-this.System.EVCMap);
                    hypotheticalBOLD = sum(diffEVCMap(~isnan(diffEVCMap)));
               else
                 HelperFnc.warning('Cannot calculate BOLD activity: EVCMaps don''t match.');
                 hypotheticalBOLD = 0;
               end
           else
              EVC.HelperFnc.warning('Cannot calculate BOLD activity: No recent EVCMap.'); 
              hypotheticalBOLD = 0;
           end
           
           % 1) take signed differecne between current and recent EVCMap
           if(isfield(this.System, 'EVCMap'))
                hypotheticalBOLD = [hypotheticalBOLD sum(EVCMap-this.System.EVCMap)];
           else
                hypotheticalBOLD = [hypotheticalBOLD 0];
           end
       
           % 2) take EVC of currently implemented signal
           hypotheticalBOLD = [hypotheticalBOLD max(max(EVCMap))];
           
           % 3) take overall control intensity
           hypotheticalBOLD = [hypotheticalBOLD sum(optIntensities)];
           
           % 4) Take current costs
           hypotheticalBOLD = [hypotheticalBOLD signalCosts(optIdx)];
           
           % 5) Take current reconfiguration costs
           hypotheticalBOLD = [hypotheticalBOLD reconfigurationCosts(optIdx)];
           
           % 6) Take absolute difference between current and previous intensities
           hypotheticalBOLD = [hypotheticalBOLD sum(abs(optIntensities-backupIntensities))];
           
           % 7) Take signed difference between current and previous intensities
           hypotheticalBOLD = [hypotheticalBOLD sum(optIntensities-backupIntensities)];
           
           
           this.System.EVCMap = EVCMap;
       end
        
    end
        
 %% Sim controller functions
    
    % these functions aren't implemented in the EVCModel class, they just
    % serve as virtual functions for child classes
    methods
        
        % calculates EVC for future states
        function discountFutureStates(this, Intensities)
            
            futureState = this.getFutureExpectedState();
            
            % need to redo -> think about new averaging method
        end
        
        % sets intensities of current control signals
        function setIntensities(this, Intensities)
            
            if(length(this.System.CurrSigIdx) ~= size(Intensities,1) && length(Intensities) ~= 0)
               error('Number of rows in intensity-vector has to match number of current control signals.');
            end
            
            % loop through current control signals and set intensities
            for idx = 1:length(this.System.CurrSigIdx)
               currentSigIdx = this.System.CurrSigIdx(idx);
                % set new intensity to search range
                this.System.CtrlSignals(currentSigIdx).Intensity = Intensities(idx,:);
            end
        end
        
        % returns intensities of current control signals
        function out = getIntensities(this, varargin)
            if(~isempty(varargin))
                signalIdx = varargin{1};
                
                out = this.System.CtrlSignals(signalIdx).Intensity;
            else
            
            for idx = 1:length(this.System.CurrSigIdx)
                currentSigIdx = this.System.CurrSigIdx(idx);
                out(idx,:) = this.System.CtrlSignals(currentSigIdx).Intensity;
            end
            
            end
        end
        
        % return control-stimulus map
        function out = getCtrlSigStimMap(this, varargin)
             if(~isempty(varargin))
                signalIdx = varargin{1};
                
                out = this.System.CtrlSignals(signalIdx).CtrlSigStimMap;
             else
                for idx = 1:length(this.System.CurrSigIdx)
                    currentSigIdx = this.System.CurrSigIdx(idx);
                    out(idx,:) = this.System.CtrlSignals(currentSigIdx).CtrlSigStimMap;
                end
             end
        end
        
        % sets reconfiguration costs
        function setReconfigurationCost(this, reconfCostFnc)
            
            % check input
            if(isa(reconfCostFnc, 'EVC.EVCFnc') == 0)
               error('Please make sure that the first argument is an instance of the class ''EVC.EVCFnc''.'); 
            end
            
            % set reconfiguration costs
            this.System.reconfigurationCost = reconfCostFnc;
        end
        
        % sets reconfiguration costs
        function setRewardFnc(this, rewardFunction)
            
            % check input
            if(isa(rewardFunction, 'EVC.EVCFnc') == 0)
               error('Please make sure that the first argument is an instance of the class ''EVC.EVCFnc''.'); 
            end
            
            % set reconfiguration costs
            this.rewardFnc = rewardFunction;
            this.rewardFnc.EVCModel = this;
        end
   
        % sets learning functions
        function setLearning(this, learningFnc)
            
            if(isa(learningFnc, 'EVC.EVCFnc') == 0)
               error('Please make sure that the first argument is an instance of the class ''EVC.EVCFnc''.'); 
            end

            for i = 1:length(learningFnc)
                this.learning(i) = learningFnc(i);
                this.learning(i).EVCModel = this;
            end
        end
        
        % update actual state based on task environment
        function setActualState(this, TaskEnv)
            [ST,I] = dbstack; % track function name for debugging functions
            source = ST.name;

            if(exist('TaskEnv','var'))
                currTrial = TaskEnv.currentTrial();
            else
                currTrial = this.State.TaskEnv.currentTrial();
            end

            % update current State
            this.State.Actual = EVC.Trial(currTrial);
                        
        end  
        
        % returns indicies of expected states that match current actual state
        function expectedStateIdx = getExpectedStateIdx(this)
           expectedStateIdx = find([this.State.ExpectedSpace.typeID] == this.State.Actual.typeID);
        end
        
        % returns the current expected state for the current actual state
        function expectedState = getExpectedState(this)
%             expectedState = this.State.ExpectedSpace([this.State.ExpectedSpace.typeID] == this.State.Actual.typeID);
            expectedState = EVC.Trial(this.State.Actual); % set first expected trial to actual expected trial
        end
        
        % cues next trial
        function cueTrial(this)
           % for now: cueing means that the trial identity (cueID) is given 
           
           this.State.expectedStatesIdx = [];
           for i = 1:length(this.State.Actual.cueID)
               
                this.State.expectedStatesIdx = [this.State.expectedStatesIdx find([this.State.ExpectedSpace.cueID] == this.State.Actual.cueID(i))];
                
                
                %for each of the cued expected states, cue specific trial information

                % cue outcome values
                currentStateIdx = this.State.expectedStatesIdx(end);
                if(isfield(this.State.Actual , 'cueOutcomeValues') || isprop(this.State.Actual , 'cueOutcomeValues'))
                    if(this.State.Actual.cueOutcomeValues)
                        this.State.ExpectedSpace(currentStateIdx).outcomeValues = this.State.Actual.outcomeValues;
                    end
                end

                % cue stimulus salience
                if(isfield(this.State.Actual , 'cueStimulusSalience') || isprop(this.State.Actual , 'cueStimulusSalience'))
                    if(this.State.Actual.cueStimulusSalience)
                        this.State.ExpectedSpace(currentStateIdx).stimSalience = this.State.Actual.stimSalience;
                    end
                end

                % cue stimulus response mapping
                if(isfield(this.State.Actual , 'cueStimRespMap') || isprop(this.State.Actual , 'cueStimRespMap'))
                    if(this.State.Actual.cueStimRespMap)
                        this.State.ExpectedSpace(currentStateIdx).stimRespMap = this.State.Actual.stimRespMap;
                    end
                end
                
           end
           
           % set each expected state to actual state if desired
            if(this.State.TaskEnv.cueAvailableTrialTypes)

                for i  = 1:length(this.State.expectedStatesIdx)
                    
                    trialTypeID = this.State.ExpectedSpace(this.State.expectedStatesIdx(i)).typeID;
                    for j = 1: length(this.State.TaskEnv.trialTypes)
                        
                        if(this.State.TaskEnv.trialTypes{j}.typeID == trialTypeID)
                            this.State.ExpectedSpace(this.State.expectedStatesIdx(i)) = EVC.Trial(this.State.TaskEnv.trialTypes{j});
                        end
                        
                    end
                    
                end

            end
           
        end
        
        % perform trial with current signal configuration
        function executeTrial(this)
            
            import EVC.DDM.*;
            [ST,I] = dbstack; % track function name for debugging functions
            source = ST.name;
            EVC.HelperFnc.disp(source, 'Executing trial...'); % DIAGNOSIS

            % cue trial
            this.cueTrial();
            
            % EVC calculation
            this.useActualState = 0;
            [optIntensities,optEVC,optStateID,~,~, expectedPerformance, optEV, optCost stateData] = this.getOptSignals();
            
            EVC.HelperFnc.disp(source, ' chosen intensities', optIntensities); % DIAGNOSIS
            % state selection
            this.State.ActualOrg = EVC.Trial(this.State.Actual);
            if(~isempty(optStateID)) % select state
                this.State.Actual = EVC.Trial(this.State.TaskEnv.trialTypes{optStateID});
                % find corresponding expected state index
                relevantIdx = [];
                for expectedStateIdx = 1:length(this.State.ExpectedSpace)
                    if(this.State.ExpectedSpace(expectedStateIdx).typeID == this.State.Actual.typeID)
                        relevantIdx = expectedStateIdx;
                        break;
                    end
                end
                if(isempty(relevantIdx))
                    error('Could not find index of expected space that matches the typeID of the optimal (selected) state.');
                else
                    this.State.Expected = EVC.Trial(this.State.ExpectedSpace(relevantIdx));
                end
            end
            
            
            % store expected outcomes
            this.State.Expected.performance = expectedPerformance;
            this.State.Expected.performance.EVC = optEVC;
            this.State.Expected.performance.EV = optEV;
            this.State.Expected.performance.Cost = optCost;
            
            % store state data
            this.State.Expected.performance.stateEVCs = stateData.stateEVCs;
            this.State.Expected.performance.stateCosts = stateData.stateCosts;
            this.State.Expected.performance.stateEVs = stateData.stateEVs;
            this.State.Expected.performance.statePerformance = stateData.statePerformance;
            
%             this.System.expectedProbs = expectedProbs;
%             this.System.expectedRTs = expectedRTs;
            
            % task-set inertia
%             if(~isempty(this.System.taskSetInertia))
%                 newIntensities = optIntensities .* (1-this.System.taskSetInertia) + this.System.taskSetInertia * this.getIntensities();
%             else
%                 newIntensities = optIntensities;
%             end
            
            % select signal intensities associated with maximal EVC
            this.setIntensities(optIntensities);
        
            % simulation for actual state
            useActualStateBackup = this.useActualState;
            this.useActualState = 1;
            
            % retrieve actual outcomes
            [EVC_val, Cost, EV, actualPerformance] = this.getEVC();
            
            
            this.State.Actual.performance = actualPerformance;
            this.State.Actual.performance.EVC = EVC_val;
            this.State.Actual.performance.EV = EV;
            this.State.Actual.performance.Cost = Cost;
%             this.State.Actual.probs = actualProbs;
%             this.State.Actual.RTs = actualFullRT;
            
            % log simulation data
            this.log();
            
            
            % update expected state
            this.updateExpectedState();
            
            
            this.useActualState = useActualStateBackup;
        end
        
        % initialize expected state
        function initExpectedStates(this)
            
            % set task sets of expected states
            numTrialTypes = length(this.State.TaskEnv.trialTypes);
            
            % initialize expected state space
           for i = 1:numTrialTypes
                this.State.ExpectedSpace(i) = EVC.Trial(this.State.TaskEnv.trialTypes{i});
           end
           
        end
        
        % update expected state
        function updateExpectedState(this)
            for i = 1:length(this.learning)
               this.learning(i).getVal(); 
            end
        end
        
        % returns future expected state
        function getFutureExpectedState(this)
            % virtual function: needs to be specified by class child
        end
        
        % switches to the next trial of used task sequence
        function nxtTrl = nextTrial(this)
            nxtTrl = this.State.TaskEnv.nextTrial();
            this.State.CurrTrlIdx = this.State.TaskEnv.CurrTrialIdx;
            
            % set actual state
            this.setActualState(this.State.TaskEnv);
            
            % define cued expected states
            this.State.expectedStatesIdx = find([this.State.ExpectedSpace.cueID] == this.State.Actual.cueID);
            
        end
        
        % returns the current trial of used task sequence
        function currTrl = currentTrial(this)
            currTrl = this.State.TaskEnv.currentTrial();
        end
        
        % returns outcome probability of corresponding state
        % TODO: move this to EVCDDM
        function out = getOutcomeProb(this, actualState)
            if(actualState)
                out = 1- this.State.Actual.performance.probs(1); 
            else
                out = 1- this.State.Expected.performance.probs(1); % probability of hitting lower threshold (error)
            end
        end
        
        % returns outcome probability of actual state
        function out = getOutcomeProbActual(this)
           out = this.State.Actual.performance.probs; 
        end
        
        % returns trial type weights based on recent frequency
        % varargin{1}: cue type
        % varargin{2}: specific trial types (if used, cue type param will be ignored)
        function typeWeights = getTrialTypeWeights(this, varargin)
           
           if(isempty(varargin))
               typeWeights(1,:) = unique(this.State.ExpectedSpace.typeID);
           else
               if(length(varargin) > 1)
                   typeWeights(1,:) = varargin{2};
               else
                   typeWeights(1,:) = unique([this.State.ExpectedSpace([this.State.ExpectedSpace.cueID] == varargin{1}).typeID]);
               end
           end
           
           % hold frequencies in 2nd row
           typeWeights(2,:) = 0;
            
           % TEMP: history length = 100
           startTrial = max(1, this.State.TaskEnv.CurrTrialIdx-10);
           
           for i = startTrial:(this.State.TaskEnv.CurrTrialIdx-1)
                currentType = this.Log.Trials(i).typeID;
                typeWeights(2,typeWeights(1,:)==currentType) = typeWeights(2,typeWeights(1,:)==currentType) + 1;
           end

           if any(typeWeights(2,:))
               % normalize
               typeWeights(2,:) = typeWeights(2,:)/sum(typeWeights(2,:));
           else
               % if no information available, assume equal probability for each trial type
               typeWeights(2,:) = 1/size(typeWeights,2);
           end
           
           % REMOVE! (assume all equal)
%            typeWeights(2,:) = 1;
%            typeWeights(2,1) = 1;
       
           
        end
        
        % returns the type weights of the states that match the current actual state
        function currTypeWeight = getCurrentTypeWeight(this)
            typeWeights = this.getTrialTypeWeights(this.State.Actual.cueID, this.State.expectedStatesIdx);
            currTypeWeight = typeWeights(2,typeWeights(1,:) == this.State.Actual.typeID);
        end
        
        
        % returns the task set for the current trial
        function TaskSet = currentTaskSet(this, varargin)
           [ST,I] = dbstack; % track function name for debugging functions
           source = ST.name;
           
           if(length(varargin) > 1)
              boolActualState = varargin; 
           else 
              boolActualState = this.useActualState;
           end
           
           if(boolActualState)
               TaskSet = this.getActualTaskSet();
           else
               TaskSet = this.getExpectedTaskSet();
           end
           
        end
        
        function TaskSet = getActualTaskSet(this)
%             [ST,I] = dbstack; % track function name for debugging functions
%             source = ST.name;
%             import EVC.*;
% 
%             for i = 1:length(this.System.TaskSets)
%                HelperFnc.disp(source, 'TaskSet length', length(this.System.TaskSets), 'System.TaskSets(i).trialType:', this.System.TaskSets(i).trialID, 'currentTrial().descr:', this.currentTrial().descr); % DIAGNOSIS
%                if (all(this.System.TaskSets(i).trialID == this.currentTrial().ID))
%                     TaskSet = this.System.TaskSets(i);
%                   break;
%                end
%             end

        TaskSet = this.State.Actual.stimRespMap;
        
        end
        
        function TaskSet = getExpectedTaskSet(this)
            
            TaskSet = this.State.Expected.stimRespMap;
            
        end
        
        % initializes log data 
        function initLog(this)
            this.Log.LogCount = 0;
            this.Log.id = [];
            this.Log.CtrlSignals = {};
            this.Log.CtrlIntensities = [];
            this.Log.SignalIdxs = [];
            this.Log.EVCs = [];
            this.Log.RTs = [];
            this.Log.ERs = [];
            this.Log.ACC_HypotheticalBOLD = [];
            this.Log.LogIdx = [];
            this.Log.stateEVCs = [];
            this.Log.stateEVs = [];
            this.Log.stateCosts = [];
        end
        
        % logs trial outcomes for current control signals
        % should be called during each model update (updateExpectedState)
        function log(this, varargin)
            
            % log model outcomes for current control identities
            probs = this.State.Expected.performance.probs;
            RTs = this.State.Expected.performance.RTs;
            expectedProbs = probs;
            expectedRT = sum(probs.*RTs);
            expectedEVC = this.State.Expected.performance.EVC;
            expectedEV = this.State.Expected.performance.EV;
            expectedCost = this.State.Expected.performance.Cost;
            
            stateEVCs = this.State.Expected.performance.stateEVCs;
            stateCosts = this.State.Expected.performance.stateCosts;
            stateEVs = this.State.Expected.performance.stateEVs;
            
            probs = this.State.Actual.performance.probs;
            RTs = this.State.Actual.performance.RTs;
            actualProbs = probs;
            actualRT = sum(probs.*RTs);
            actualEVC = this.State.Actual.performance.EVC;
            actualEV = this.State.Actual.performance.EV;
            actualCost = this.State.Actual.performance.Cost;

            EVClog = [expectedEVC actualEVC];
            EVlog = [expectedEV actualEV];
            CostLog = [expectedCost actualCost];
            RTlog = [expectedRT actualRT];
            
            this.Log.LogCount = this.Log.LogCount + 1;
           
            this.Log.id(this.Log.LogCount,1) = this.id;                                                                     % log model id
            for idx = 1:length(this.System.CurrSigIdx)
                SigIdx = this.System.CurrSigIdx(idx);
                this.Log.CtrlSignals{this.Log.LogCount,idx} =  EVC.CtrlSignal(this.System.CtrlSignals(SigIdx));   % log array for control signals used
            end
            this.Log.CtrlIntensities(this.Log.LogCount,:) = transpose(this.getIntensities()); % log array for control intensities used (separate array b/c adapted signal intensities change)
            this.Log.SignalIdxs(this.Log.LogCount,:) = this.System.CurrSigIdx;                                            % log array for indicies of control signals used
            this.Log.EVCs(this.Log.LogCount,:) = EVClog;                                                                  % log array for EVC outcomes produced [expected actual]          
            this.Log.EVs(this.Log.LogCount,:) = EVlog;                                                                  % log array for EVC outcomes produced [expected actual]          
            this.Log.Costs(this.Log.LogCount,:) = CostLog;                                                                  % log array for EVC outcomes produced [expected actual]          
            this.Log.RTs(this.Log.LogCount,:) = RTlog;                                                                    % log array for reaction times [expected actual] 
            this.Log.actualProbs(this.Log.LogCount,:) = transpose(actualProbs);                                                      % log array for actual outcome probabilities
            this.Log.expectedProbs(this.Log.LogCount,:) = transpose(expectedProbs);                                                  % log array for actual outcome probabilities
            this.Log.ACC_HypotheticalBOLD(this.Log.LogCount,:) = this.System.ACC_HypotheticalBOLD;                                                    % log array for ACC activity
            
            this.Log.Trials(this.Log.LogCount,1) = EVC.Trial(this.currentTrial());                                          % log array for trials
            this.Log.TrialsOrg(this.Log.LogCount,1) = EVC.Trial(this.State.ActualOrg);
            this.Log.ExpectedState(this.Log.LogCount,1) = EVC.Trial(this.State.Expected); 
            this.Log.LogIdx(this.Log.LogCount,1) = this.Log.LogCount;                                                       % log steps
            
%             disp(['+++++ Trial ' num2str(this.State.TaskEnv.CurrTrialIdx)]);
%             disp([this.Log.ExpectedState(this.Log.LogCount,1).typeID]);
%             disp(stateEVCs);
%             disp(this.State.Expected.stimSalience);
            if(~isempty(stateEVCs))
                this.Log.stateEVCs(this.Log.LogCount,:)  = stateEVCs';
            end
            if(~isempty(stateEVs))
                this.Log.stateEVs(this.Log.LogCount,:)  = stateEVs';
            end
            if(~isempty(stateCosts))
                this.Log.stateCosts(this.Log.LogCount,:)  = stateCosts';
            end
            
            % log input arguments into Log.params matrix
            if(length(varargin) >= 1)
                for i = 1:length(varargin)
                    this.Log.params(this.State.TaskEnv.CurrTrialIdx, i) = varargin{i}; 
                end
            end
            
        end
        
        % writes log file; varNames argument contains variable names (string) from the
        % Log structure to write. If there are additional input arguments to write, this string structure should
        % also contain their names at the end. All remaining parameters represent
        % the additional arguments to be recorded
        function writeLogFile(this, filename, override, varNames, varargin)
            [ST,I] = dbstack; % track function name for debugging functions
            source = ST.name;
             
            succeeded = true;
            message = '';
            
            % override file?
            if(~exist('override', 'var'))
               override = 0; 
            end
            
            % check input arguments
            
            if(iscell(varNames)==0)
              succeeded = false;
              message = 'VarNames has to contain parameter names (char structure) of log parameters (& input arguments).';
            end
               
            if(~isempty(varargin))
               
               if(length(varNames) < (length(varargin)))
                  succeeded = false;
                  message = 'String structure VarNames has to have at least the same length as the number of additional input parameters.';
               end

              for i = 1:length(varargin);
                 if(length(varargin{i}) ~= this.Log.LogCount)
                    succeeded = false;
                    message = strcat('Length of additional input parameter ',int2str(i-1), ' doesn''t mach the length of log parameters.'); 
                 end
              end

            end
            
            if(succeeded == false)
                error(strcat('Error with regard to input arguments: ', message));
                return;
            end
            
            % check if file already exists
            if fopen(filename, 'rt')~=-1 && override == 0
                
                fclose('all');
                message = 'Result data file already exists! Choose a different file name.';
                succeeded = false;
                disp(message);
                
            else
                % open file
                datafilepointer = fopen(filename,'wt'); % open ASCII file for writing
                
                % write header
                
                for i = 1:length(varNames)
                    try
                        cols = size(eval(strcat('this.Log.',varNames{i})),2);
                    catch
                        message = strcat('log parameter not found: this.Log.', varNames{i});
                        succeeded = false;    
                    end
                    if(cols == 1)
                        fprintf(datafilepointer,'%s ', varNames{i});
                    else
                        for col = 1:cols
                            fprintf(datafilepointer,'%s ', strcat(varNames{i}, '_',int2str(col)));   
                        end
                    end
                    
                end
                
                fprintf(datafilepointer,'\n');
                
                if(succeeded == false)
                    error(strcat('Error writing file (header): ', message));
                    return;
                end
                
                % write data
                for row = 1:this.Log.LogCount
                    
                    vararginCounter = 1;
                    
                    for col = 1:length(varNames)
                        
                        val = NaN;
                        
                        try
                        	val = eval(strcat('this.Log.',varNames{col}));
                            val = val(row,:);
                            EVC.HelperFnc.disp(source, 'row', row, 'val exisiting:', varNames{col}); % DIAGNOSIS
                        catch
                            if(vararginCounter > length(varargin))
                               succeeded = false;
                               message = 'More variable names than parameters. This might happen b/c one of the input variable names is thought to be an additional input parameter rather than being a this.Log parameter (wrong name)';
                            else 
                               val = varargin{vararginCounter}(row);
                               vararginCounter = vararginCounter + 1;
                            end
                        end
                        
                        for valCol = 1:length(val)
                            scalar = val(valCol);
                            if(~ischar(scalar))
                                if(iscell(scalar))
                                    scalar = scalar{1};
                                else
                                    scalar = num2str(scalar);
                                end
                            end
                            fprintf(datafilepointer,'%s ', scalar);  
                        end
                        
                        
                    end
                    fprintf(datafilepointer,'\n');
                end
                

                fclose(datafilepointer);
            end
            
            if(succeeded == false)
               error(strcat('Error writing file (data): ', message));
               return;
            end
    
        end
        
        
    end
    
% for now EVCModel cannot be abstract: a dummy EVCModel instance needs
% to be created in order to make function pointers work in Simulation.m
%     methods (Abstract)
%         [probs RTs] = simulateOutcomes(this)
%     end
    
end
