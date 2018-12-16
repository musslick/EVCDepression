classdef DDM_WestbrookBraver2015 < Simulations.DDMSim
    
    % description of class
    % simulates approach avoidance behavior
    
    % global parameters
    properties
        automaticControlMappingFnc
        automaticControlFnc
        automaticControlProcess
        
        taskAReward = 10;
        taskBReward = 10;
        taskAAutomaticity = 1;
        taskBAutomaticity = 1;
        startSalience = 1.0;
        rewardMinimum = 1.0;
    end
    
    methods
        
        function this = DDM_WestbrookBraver2015()

            import EVC.*;
            import EVC.DDM.*;
            
            % call parent constructor
            this = this@Simulations.DDMSim();
            
            %% general simulation parameters
            
            this.nSubj = 1;
            this.plotSum = true;
            this.binaryErrors = 1;
            
            this.reconfCostFnc.params{1} = 0; 
            this.reconfCostFnc.params{2} = 0;
            this.defaultCostFnc.params{1} = 0.75;                
            this.defaultCostFnc.params{2} = 0;                 
               
            this.defaultDDMParams.c = 1.0;                     
            this.defaultDDMParams.thresh = 2;                
            this.defaultDDMParams.T0 = 0.2;    
            
            temp.rewardFnc.params{1} = 1; % offset for RT disocunt
            temp.rewardFnc.params{2} = 0; % scale on RT discount
            temp.rewardFnc.params{3} = 1; % value scalar
            temp.rewardFnc.params{4} = 0; % value offset
            temp.rewardFnc.type = EVCFnc.REWRATE_EXT;
            this.rewardFnc = EVCFnc(temp.rewardFnc.type, temp.rewardFnc.params);
            
            %% general simulation parameters: learning functions
            
            this.learningFnc = LearningFnc.empty(1,0);
            temp.EVCDummy = EVCModel(0,0);
            
            temp.learningFnc(1).params{1} = 1;                                           % learning rate for other state parameters (may not be necessary for reward updates)
            temp.learningFnc(1).type = LearningFnc.WESTBROOK_BRAVER_2015;       % learning function type
            temp.learningFnc(1).EVCM = temp.EVCDummy;                          % EVCModel dummy 
            this.learningFnc(1) = LearningFnc(temp.learningFnc(1).type, temp.learningFnc(1).params, temp.learningFnc(1).EVCM);
            this.learningFnc(1).input{1} = @() this.learningFnc(1).EVCModel.getOutcomeProb(1);           % dynamic input value for learning function
            
            %% automatic/control process
           
            % generate 2 control signals, one for each task
            
            % signal range
            tmp.IntensityRange = [0:0.05:4];
            
            % Task A
           this.ctrlSignals(1) = CtrlSignal(this.defaultCtrlSignal);
           this.ctrlSignals(1).CtrlSigStimMap  = [1 0];
           this.ctrlSignals(1).IntensityRange = tmp.IntensityRange;
            
           % Task B
           this.ctrlSignals(2) = CtrlSignal(this.defaultCtrlSignal);
           this.ctrlSignals(2).CtrlSigStimMap  = [0 1];
           this.ctrlSignals(2).IntensityRange = tmp.IntensityRange;
           
            %% define DDM processes
            
            % map all control processes to a specific DDM parameter
            
            % target
            temp.taskAControlFnc.type = DDMFnc.INTENSITY2DDM;
            temp.taskAControlFnc.params{1} = this.ctrlSignals(1);
            temp.taskAControlFnc.params{2} = this.defaultControlProxy;
            temp.taskAControlFnc.params{3} = this.defaultControlMappingFnc;
            temp.taskAControlFnc = DDMFnc(temp.taskAControlFnc.type, ...
                                            temp.taskAControlFnc.params);
            
            % flanker 1                            
            temp.taskBControlFnc.type = DDMFnc.INTENSITY2DDM;
            temp.taskBControlFnc.params{1} = this.ctrlSignals(2);
            temp.taskBControlFnc.params{2} = this.defaultControlProxy;
            temp.taskBControlFnc.params{3} = this.defaultControlMappingFnc;
            temp.taskBControlFnc = DDMFnc(temp.taskBControlFnc.type, ...
                                            temp.taskBControlFnc.params);
            
            % define each DDM process
            temp.taskAProcess = DDMProc(DDMProc.CONTROL, ...                  % default controlled DDM process 
                                                    DDMProc.DRIFT, ...
                                                    temp.taskAControlFnc);
           
            temp.taskBProcess = DDMProc(DDMProc.CONTROL, ...                  % default controlled DDM process 
                                        DDMProc.DRIFT, ...
                                        temp.taskBControlFnc);

            % put all DDM processes together
            this.DDMProcesses = DDMProc.empty(4,0);
            this.DDMProcesses(1) =   this.defaultAutomaticProcess;
            this.DDMProcesses(2) =   temp.taskAProcess;
            this.DDMProcesses(3) =   temp.taskBProcess;
            
            %% task environment parameters: task environment
            
            this.nTrials = 92;
            
            temp.stimSalience = [1 1];
            temp.outcomeValues = [5 0];
            w = 1;
            
            % task A trial
            this.trials(1).ID = 1;                                                          % trial identification number (for task set)
            this.trials(1).typeID = 1;                                                      % trial type (defines task context)
            this.trials(1).cueID = 1;                                                       % cued information about trial identity
            this.trials(1).descr = 'taskA';                                                 % trial description
            this.trials(1).conditions    = [1];                                           % set of trial conditions (for logging)
            this.trials(1).outcomeValues = [this.taskAReward 0];                  % reward for correct outcome = 3; no reward/punishment for incorrect outcome
            this.trials(1).stimSalience  = temp.stimSalience;                         % relative salience between stimulus 1 and stimulus 2 defines level of incongruency; here stim 2 is more dominant (like in stroop)
            this.trials(1).stimRespMap   = [this.taskAAutomaticity 0;                                          % responding to stimulus 1 tends to produce second response by 100%
                                            0 0];                                         % responding to stimulus 2 tends to produce second outcome by 100%
            this.trials(1).params = [];                                                     % DDM specific trial parameters
            this.trials(1).cueOutcomeValues = 1;                                    % cue the outcome value for this trial

            % task B trial
            this.trials(2).ID = 2;                                                          % trial identification number (for task set)
            this.trials(2).typeID = 2;                                                      % trial type (defines task context)
            this.trials(2).cueID = 2;                                                       % cued information about trial identity
            this.trials(2).descr = 'taskB';                                                 % trial description
            this.trials(2).conditions    = [2];                                           % set of trial conditions (for logging)
            this.trials(2).outcomeValues = [0 this.taskBReward] ;                  % reward for correct outcome = 3; no reward/punishment for incorrect outcome
            this.trials(2).stimSalience  = temp.stimSalience;                         % relative salience between stimulus 1 and stimulus 2 defines level of incongruency; here stim 2 is more dominant (like in stroop)
            this.trials(2).stimRespMap   = [0 0;                                          % responding to stimulus 1 tends to produce second response by 100%
                                            0 this.taskBAutomaticity];                                         % responding to stimulus 2 tends to produce second outcome by 100%
            this.trials(2).params = [];                                                     % DDM specific trial parameters
            this.trials(1).cueOutcomeValues = 1;                                    % cue the outcome value for this trial

            
            % choice trial
            this.trials(3).ID = [1 2];                                                          % trial identification number (for task set)
            this.trials(3).typeID = 3;                                                      % trial type (defines task context)
            this.trials(3).cueID = [1 2];                                                       % cued information about trial identity
            this.trials(3).descr = 'choice';                                                 % trial description
            this.trials(3).conditions    = [0];                                           % set of trial conditions (for logging)
            this.trials(3).outcomeValues = temp.outcomeValues;                  % reward for correct outcome = 3; no reward/punishment for incorrect outcome
            this.trials(3).stimSalience  = temp.stimSalience;                         % relative salience between stimulus 1 and stimulus 2 defines level of incongruency; here stim 2 is more dominant (like in stroop)
            this.trials(3).stimRespMap   = [0 0;                                          % responding to stimulus 1 tends to produce second response by 100%
                                                            0 0];                                         % responding to stimulus 2 tends to produce second outcome by 100%
            this.trials(3).params = [];                                                     % DDM specific trial parameters

            %% log parameters
            
            this.writeLogFile = 1; 
            this.logFileName = 'DDM_WestbrookBraver2015'; 
            
            this.logAddVars{3} = '[this.EVCM.Log.Trials.conditions]''';
            this.logAddVarNames{3} = 'condition';
        end
        
        function getResults(this)
            
        end
        
        function dispResults(this)
            
              mc = metaclass(this);
              disp(['++++++++++ ' mc.Name ' ++++++++++']);

              
              
        end
        
        function plotSummary(this) 
            
            
        end
        
        
    end
    
    methods (Access = protected)
        
        % generate task environment
        function initTaskEnv(this)
            
            allTrials(1) = EVC.Trial(this.trials(1));
            allTrials(2) = EVC.Trial(this.trials(2));
            allTrials(3) = EVC.Trial(this.trials(3));
            
            % set outcome values
            allTrials(1).outcomeValues = [this.taskAReward 0]; 
            allTrials(2).outcomeValues = [0 this.rewardMinimum]; 
            
            allTrials(1).stimRespMap   = [this.taskAAutomaticity 0;                    
                                                            0 0];                  
            allTrials(2).stimRespMap   = [0 0;                    
                                                            0 this.taskBAutomaticity];      
                                                        

            choiceTrial = EVC.Trial(this.trials(3));
            
            % generate task environment
            this.taskEnv = EVC.TaskEnv(allTrials, this.nTrials);
            
            % fill sequence with choice trials
            for trial = 1:length(this.taskEnv.Sequence)
                this.taskEnv.Sequence(trial) = EVC.Trial(choiceTrial);
                this.taskEnv.cueAvailableTrialTypes = 1;        % set to 1 if the expected state should always be synchronized to the actual state for all choices
            end
            
            
            
        end
        
    end
    
    methods (Access = public)
        
       function initOptimizationTaskEnv(this)
              
       end
                           

       function criterion = getOptimizationCriterion(this)

       end
        
        
    end
    
end

