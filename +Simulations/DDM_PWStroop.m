classdef DDM_PWStroop < Simulations.DDMSim
    
    % description of class
    % simulates sequential control adjustments in response to conflict (see
    % PWStroop 1992)
    
    % global parameters
    properties
        automaticControlMappingFnc
        automaticControlFnc
        automaticControlProcess
        
        fitMainEffects                          % flag used for fitting: determines if the main effects should be fit (vs. interactions)
        reward = 10;                          % amount of reward provided in reward condition
        noReward = 4;                       % amount of (baseline) reward provided in no-reward condition
        frequencyCondition = 1;         % 1 - balanced congruent/incongruent trials, 2 - mostly incongruent trials, 3 - mostly congruent trials
                                                        % 4 - balanced congruent/incongruent trials (for higher trial number), 5 - allmost all incongruent trials, 6 = allmost all congruent trials
        
    end
    
    methods
        
        function this = DDM_PWStroop()

            import EVC.*;
            import EVC.DDM.*;
            
            % call parent constructor
            this = this@Simulations.DDMSim();
            
            %% general simulation parameters
            
            this.nSubj = 1;                                                 % number of subjects
            this.plotSum = true;                                         % plot simulation summary?
            
            this.learningFnc(1).params{1} = 0;%0.68;              % parameter for state learning function 
            this.reconfCostFnc.params{1} = 0;                   % reconfiguration cost parameter 1 (see EVC/EVCFnc.m)
            this.reconfCostFnc.params{2} = 0;                   % reconfiguration cost parameter 2 (see EVC/EVCFnc.m)
            this.defaultCostFnc.params{1} = 0.1561;         % intensity cost parameter 1 (see EVC/EVCFnc.m)
            this.defaultCostFnc.params{2} = 0;                  % intensity cost parameter 2 (see EVC/EVCFnc.m)
            
            this.rewardFnc.params{2} = 0;                          % RT scaling parameter in reward rate function (0 if RT should not be taken into account)
               
            this.defaultDDMParams.c = 0.5;                      % DDM noise parameter              
            this.defaultDDMParams.thresh = 0.3634;        % DDM threshold parameter
            this.defaultDDMParams.T0 = 0.2500;              % DDM non-decision time                 
            
            %% automatic/control process
           
            % generate 3 control signals, one for the target stimulus & two
            % for the flanker stimuli
            
            % signal range
            tmp.IntensityRange = [0:0.2:10];                                        % Range of control signal intensities to be tested
            
            % target signal (PICTURE)
           this.ctrlSignals(1) = CtrlSignal(this.defaultCtrlSignal);        % define control signal for picture stimulus
           this.ctrlSignals(1).CtrlSigStimMap  = [0 1];                          % set control signal map for picture stimulus
           this.ctrlSignals(1).IntensityRange = tmp.IntensityRange;    % assign range of tested control signal intensities for picture stimulus
            
           % distractor signal (WORD)
           this.ctrlSignals(2) = CtrlSignal(this.defaultCtrlSignal);        % define control signal for word stimulus
           this.ctrlSignals(2).CtrlSigStimMap  = [1 0];                         % set control signal map for word stimulus
           this.ctrlSignals(2).IntensityRange = tmp.IntensityRange;    % assign range of tested control signal intensities for word stimulus
           
            %% define DDM processes
            
            % map all control processes to a specific DDM parameter
            
            % target signal (PICTURE)
            temp.targetControlFnc.type = DDMFnc.INTENSITY2DDM;                      % define function type for mapping control onto DDM parameter (see EVC/DDM/DDMFnc.m)
            temp.targetControlFnc.params{1} = this.ctrlSignals(1);                          % assign control signal for picture stimulus
            temp.targetControlFnc.params{2} = this.defaultControlProxy;               % determine control signal proxy (e.g. drift rate, see Simulations/DDMSim.m)
            temp.targetControlFnc.params{3} = this.defaultControlMappingFnc;     % assign function that maps control signal onto drift rate (see Simulations/DDMSim.m)
            temp.targetControlFnc = DDMFnc(temp.targetControlFnc.type, ...
                                            temp.targetControlFnc.params);
            
            % distractor signal (WORD)                            
            temp.flanker1ControlFnc.type = DDMFnc.INTENSITY2DDM;                  % define function type for mapping control onto DDM parameter (see EVC/DDM/DDMFnc.m)
            temp.flanker1ControlFnc.params{1} = this.ctrlSignals(2);                      % assign control signal for word stimulus
            temp.flanker1ControlFnc.params{2} = this.defaultControlProxy;           % determine control signal proxy (e.g. drift rate, see Simulations/DDMSim.m)
            temp.flanker1ControlFnc.params{3} = this.defaultControlMappingFnc; % assign function that maps control signal onto drift rate (see Simulations/DDMSim.m)
            temp.flanker1ControlFnc = DDMFnc(temp.flanker1ControlFnc.type, ...
                                            temp.flanker1ControlFnc.params);
            
            % define each DDM process
            temp.targetProcess = DDMProc(DDMProc.CONTROL, ...                  % default controlled DDM process 
                                                    DDMProc.DRIFT, ...
                                                    temp.targetControlFnc);
           
            temp.flanker1Process = DDMProc(DDMProc.CONTROL, ...               % default controlled DDM process 
                                        DDMProc.DRIFT, ...
                                        temp.flanker1ControlFnc);

            % put all DDM processes together
            this.DDMProcesses = DDMProc.empty(4,0);
            this.DDMProcesses(1) =   this.defaultAutomaticProcess;
            this.DDMProcesses(2) =   temp.targetProcess;
            this.DDMProcesses(3) =   temp.flanker1Process;
            
            %% task environment parameters: task environment
            
            this.nTrials = 216;
            
            temp.stimSalience = [1 1];
            temp.outcomeValues = [100 0];
            w = 0.2754;
            
            % congruent trial,
            % e.g. [ PICTURE = HOUSE, WORD = HOUSE ]
            this.trials(1).ID = 1;                                                          % trial identification number (for task set)
            this.trials(1).typeID = 1;                                                   % trial type (defines task context)
            this.trials(1).cueID = 1;                                                    % cued information about trial identity (this should number be the same across trials for experiments in which there is no cue about the stimulus)
            this.trials(1).descr = 'con';                                                % trial description
            this.trials(1).conditions    = [1 0];                                     % set of trial conditions (for logging)
            this.trials(1).outcomeValues = temp.outcomeValues;       % reward for correct outcome; no reward/punishment for incorrect outcome
            this.trials(1).stimSalience  = temp.stimSalience;               % relative salience between picture stimulus and word stimulus
            this.trials(1).stimRespMap   = [w 0;                                   % responding to picture leads to correct outcome with 100%
                                            w 0];                                   % responding to word leads to correct outcome with 100%
            this.trials(1).params = [];                                                   % potential other prameters
            this.trials(1).cueOutcomeValues = 1;                                 % cue the outcome value for this trial (reward cue)

            
            % incongruent trial,
            % e.g. [ PICTURE = HOUSE, WORD = BUILDING ]
            this.trials(2).ID = 1;                                                             % trial identification number (for task set)
            this.trials(2).typeID = 1;                                                      % trial type (defines task context)
            this.trials(2).cueID = 1;                                                       % cued information about trial identity (this should number be the same across trials for experiments in which there is no cue about the stimulus)
            this.trials(2).descr = 'inc';                                                   % trial description
            this.trials(2).conditions    = [2 0];                                       % set of trial conditions (for logging)
            this.trials(2).outcomeValues = temp.outcomeValues;         % reward for correct outcome; no reward/punishment for incorrect outcome
            this.trials(2).stimSalience  = temp.stimSalience;                % relative salience between picture stimulus and word stimulus
            this.trials(2).stimRespMap   = [0 w;                                    % responding to picture leads to correct outcome with 100%
                                            w 0];                                    % responding to word leads to incorrect outcome with 100%
            this.trials(2).params = [];                                                    % potential other prameters
            this.trials(2).cueOutcomeValues = 1;                                   % cue the outcome value for this trial (reward cue)
            
            % neutral trial,
            % e.g. [ PICTURE = HOUSE, WORD = XXXXX ]
            this.trials(3).ID = 1;                                                          % trial identification number (for task set)
            this.trials(3).typeID = 1;                                                   % trial type (defines task context)
            this.trials(3).cueID = 1;                                                     % cued information about trial identity (this should number be the same across trials for experiments in which there is no cue about the stimulus)
            this.trials(3).descr = 'ntr';                                                 % trial description
            this.trials(3).conditions    = [0 0];                                     % set of trial conditions (for logging)
            this.trials(3).outcomeValues = temp.outcomeValues;       % reward for correct outcome; no reward/punishment for incorrect outcome
            this.trials(3).stimSalience  = temp.stimSalience;              % relative salience between picture stimulus and word stimulus
            this.trials(3).stimRespMap   = [w/2 w/2;                                   % responding to picture leads to correct outcome with 100%
                                            w 0];                                  % responding to word leads neither to correct or incorrect outcome (no learned associated response)
            this.trials(3).params = [];                                                  % potential other prameters
            this.trials(3).cueOutcomeValues = 1;                                 % cue the outcome value for this trial (reward cue)
            
            
            %% log parameters
            
            this.fitMainEffects = 1;
            
            this.writeLogFile = 1; 
            this.logFileName = 'DDM_PWStroop'; 
            
            this.logAddVars{3} = '[this.EVCM.Log.Trials.conditions]''';
            this.logAddVarNames{3} = 'condition';
        end
        
        % computes simulation results -> put analysis code here
        function getResults(this)
            
            % initialize all conditions
            this.results.RT.Ntr = zeros(1, this.nSubj);
            this.results.RT.Con = zeros(1, this.nSubj);
            this.results.RT.Inc = zeros(1, this.nSubj);
            
            this.results.Accuracy.Ntr = zeros(1, this.nSubj);
            this.results.Accuracy.Con = zeros(1, this.nSubj);
            this.results.Accuracy.Inc = zeros(1, this.nSubj);

            this.results.targetIntensity.Ntr = zeros(1, this.nSubj);
            this.results.targetIntensity.Con = zeros(1, this.nSubj);
            this.results.targetIntensity.Inc = zeros(1, this.nSubj);
            
            this.results.flankerIntensity.Ntr = zeros(1, this.nSubj);
            this.results.flankerIntensity.Con = zeros(1, this.nSubj);
            this.results.flankerIntensity.Inc = zeros(1, this.nSubj);
            
            this.results.analysis.RTData = [];
            this.results.analysis.ERData = [];
            this.results.analysis.rewardCondition = [];
            this.results.analysis.congruencyCondition = [];
            this.results.analysis.subject = [];
            
            % loop through all subjects
            for subj = 1:this.nSubj
                
                % get log structure for current subject
                subjLog = this.subjData(subj).Log;
                
                % extract relevant trial conditions fpr current subject
                conditions = reshape([subjLog.TrialsOrg.conditions],[2 this.nTrials])';
                congruency = conditions(:,1);
                rewardCondition = conditions(:,2);
                trial = transpose(1:length(conditions(:,2)));
                
                % extract RTs and Accuracies from log data for current subject
                RT = this.subjData(subj).Log.RTs(:,2)';
                ER = this.subjData(subj).Log.ERs(:,2)';
                CtrlIntensities = this.subjData(subj).Log.CtrlIntensities;
                targetIntensities = CtrlIntensities(:,1);
                flankerIntensities = CtrlIntensities(:,2);
                
                driftControl = sum(this.subjData(subj).Log.ControlParamVal,2)';
                driftAutomatic = this.subjData(subj).Log.ActualStateParam';
                
                this.results.meanIntensity(subj, :) = mean(subjLog.CtrlIntensities(1));
                
                % fetch analysis data
                this.results.analysis.RTData = [this.results.analysis.RTData; RT'];
                this.results.analysis.ERData = [this.results.analysis.ERData; ER'];
                this.results.analysis.rewardCondition = [this.results.analysis.rewardCondition; (rewardCondition+1)];
                this.results.analysis.congruencyCondition = [this.results.analysis.congruencyCondition; (congruency+1)];
                this.results.analysis.subject = [this.results.analysis.subject; repmat(subj, length(RT), 1)];
                
                % get results
%                 this.results.RT.Ntr(subj) = mean(RT(rewardCondition == 0 & congruency == 0 & trial > 1));
                this.results.RT.Con(subj) = mean(RT(rewardCondition == 0 & congruency == 1 & trial > 1));
                this.results.RT.Inc(subj) = mean(RT(rewardCondition == 0 & congruency == 2 & trial > 1));
                
%                 this.results.Accuracy.Ntr(subj) = mean(1-ER(rewardCondition == 0 & congruency == 0 & trial > 1));
                this.results.Accuracy.Con(subj) = mean(1-ER(rewardCondition == 0 & congruency == 1 & trial > 1));
                this.results.Accuracy.Inc(subj) = mean(1-ER(rewardCondition == 0 & congruency == 2 & trial > 1));
                
%                 this.results.targetIntensity.Ntr(subj) = mean(targetIntensities(rewardCondition == 0 & congruency == 0 & trial > 1));
                this.results.targetIntensity.Con(subj) = mean(targetIntensities(rewardCondition == 0 & congruency == 1 & trial > 1));
                this.results.targetIntensity.Inc(subj) = mean(targetIntensities(rewardCondition == 0 & congruency == 2 & trial > 1));
                this.results.targetIntensity.rew(subj) = mean(targetIntensities(rewardCondition == 1 & trial > 1));
                this.results.targetIntensity.norew(subj) = mean(targetIntensities(rewardCondition == 0 & trial > 1));
                
%                 this.results.flankerIntensity.Ntr(subj) = mean(flankerIntensities(rewardCondition == 0 & congruency == 0 & trial > 1));
                this.results.flankerIntensity.Con(subj) = mean(flankerIntensities(rewardCondition == 0 & congruency == 1 & trial > 1));
                this.results.flankerIntensity.Inc(subj) = mean(flankerIntensities(rewardCondition == 0 & congruency == 2 & trial > 1));
                this.results.flankerIntensity.rew(subj) = mean(flankerIntensities(rewardCondition == 1 & trial > 1));
                this.results.flankerIntensity.norew(subj) = mean(flankerIntensities(rewardCondition == 0 & trial > 1));
                
%                 this.results.driftControl.Ntr(subj) = mean(driftControl(rewardCondition == 0 & congruency == 0 & trial > 1));
                this.results.driftControl.Con(subj) = mean(driftControl(rewardCondition == 0 & congruency == 1 & trial > 1));
                this.results.driftControl.Inc(subj) = mean(driftControl(rewardCondition == 0 & congruency == 2 & trial > 1));
                this.results.driftControl.rew(subj) = mean(driftControl(rewardCondition == 1 & trial > 1));
                this.results.driftControl.norew(subj) = mean(driftControl(rewardCondition == 0 & trial > 1));
                
%                 this.results.driftAutomatic.Ntr(subj) = mean(driftAutomatic(rewardCondition == 0 & congruency == 0 & trial > 1));
                this.results.driftAutomatic.Con(subj) = mean(driftAutomatic(rewardCondition == 0 & congruency == 1 & trial > 1));
                this.results.driftAutomatic.Inc(subj) = mean(driftAutomatic(rewardCondition == 0 & congruency == 2 & trial > 1));
                this.results.driftAutomatic.rew(subj) = mean(driftAutomatic(rewardCondition == 1 & trial > 1));
                this.results.driftAutomatic.norew(subj) = mean(driftAutomatic(rewardCondition == 0 & trial > 1));
                
%                 this.results.drift.Ntr(subj) = this.results.driftControl.Ntr(subj) + this.results.driftAutomatic.Ntr(subj);
                this.results.drift.Con(subj) = this.results.driftControl.Con(subj) + this.results.driftAutomatic.Con(subj);
                this.results.drift.Inc(subj) = this.results.driftControl.Inc(subj) + this.results.driftAutomatic.Inc(subj);
                this.results.drift.rew(subj) = this.results.driftControl.rew(subj) + this.results.driftAutomatic.rew(subj);
                this.results.drift.norew(subj) = this.results.driftControl.norew(subj) + this.results.driftAutomatic.norew(subj);
                
                
            end
            
 
        end
        
        % displays simulation results
        function dispResults(this)
            
            mc = metaclass(this);
            disp(['++++++++++ ' mc.Name ' ++++++++++']);
            
            % original data from Padmala (no-reward condition)
            orgData.RT.Ntr = 0.605;
            orgData.RT.Con = 0.577;
            orgData.RT.Inc = 0.645;

            orgData.Accuracy.Ntr = 0.947;
            orgData.Accuracy.Con = 0.96;
            orgData.Accuracy.Inc = 0.875;
            
            % descriptive statistics
            disp(['Ntr: (RT) ' num2str(mean(this.results.RT.Ntr)) ' vs. ' num2str(orgData.RT.Ntr) ]);
            disp(['Con: (RT) ' num2str(mean(this.results.RT.Con)) ' vs. ' num2str(orgData.RT.Con) ]);
            disp(['Inc: (RT) ' num2str(mean(this.results.RT.Inc)) ' vs. ' num2str(orgData.RT.Inc) ]);

            disp('---')
            disp(['Ntr: (Accuracy) ' num2str(mean(this.results.Accuracy.Ntr)) ' vs. ' num2str(orgData.Accuracy.Ntr) ]);
            disp(['Con: (Accuracy) ' num2str(mean(this.results.Accuracy.Con)) ' vs. ' num2str(orgData.Accuracy.Con) ]);
            disp(['Inc: (Accuracy) ' num2str(mean(this.results.Accuracy.Inc)) ' vs. ' num2str(orgData.Accuracy.Inc) ]);

 
        end
        
        % place for summary plots (based on this.results)
        function plotSummary(this) 
             
        end
        
    end
    
    methods (Access = protected)
        
        % generate task environment
        function initTaskEnv(this)
            
            % use the same task environment as used for fitting (parameter optimization), see blow
            this.initOptimizationTaskEnv();
            
        end
        
    end
    
    methods (Access = public)
        
        % generate task environment used for parameter fitting
       function initOptimizationTaskEnv(this)
           
           % create trials
           congruentTrial = EVC.Trial(this.trials(1));
           incongruentTrial = EVC.Trial(this.trials(2));
           neutralTrial = EVC.Trial(this.trials(3));
           
           % CREATE TASK ENVIRONMENT
           
           % 0 - neutral, 1 - congruent, 2 - incongruent
           switch this.frequencyCondition
               case 1
                    congruencyLevels = [1 2]; % balanced conditions
               case 2
                    congruencyLevels = [0 0 0 1 2 2 2 2 2]; % more incongruent than congruent
               case 3
                    congruencyLevels = [0 0 0 1 1 1 1 1 2]; % more congruent than incongruent
               case 4
                    congruencyLevels = [0 0 0 1 1 1 1 1 1 2 2 2 2 2 2]; % balanced conditions
               case 5
                    congruencyLevels = [0 0 0 1 2 2 2 2 2 2 2 2 2 2 2]; % more incongruent than congruent
               case 6
                    congruencyLevels = [0 0 0 1 1 1 1 1 1 1 1 1 1 1 2]; % more congruent than incongruent
           end
           rewardLevels = [0];
           
           crossedLevels = Simulations.combvec(congruencyLevels, rewardLevels);
           trialsPerCell = floor(this.nTrials/size(crossedLevels,2));
           sequence = repmat(crossedLevels, 1, trialsPerCell);

           congruencyConds =  sequence(1,:);
           rewardConds=  sequence(2,:);
           
          if(length(congruencyConds) ~= length(rewardConds))
              warning('initOptimizationTaskEnv: Length of congruency conditions and reward conditions doens''t match.');
          else
              this.nTrials = length(congruencyConds);
          end
          
          % shuffle
          newOrder= randperm(this.nTrials);
          congruencyConds = congruencyConds(newOrder);
          rewardConds = rewardConds(newOrder);
           
           % generate sequence
            this.taskEnv = EVC.TaskEnv(neutralTrial, this.nTrials);
            for trl = 1:(this.nTrials+1)
                
                if(trl == 1)
                    this.taskEnv.Sequence(trl) = EVC.Trial(neutralTrial);
                    this.taskEnv.Sequence(trl).outcomeValues(1) =  this.noReward;
                    this.taskEnv.Sequence(trl).conditions(2) = 0;
                else
                % assign trial
                 switch congruencyConds(trl-1)
                    case 0
                        this.taskEnv.Sequence(trl) = EVC.Trial(neutralTrial);
                    case 1
                        this.taskEnv.Sequence(trl) = EVC.Trial(congruentTrial);
                    case 2
                        this.taskEnv.Sequence(trl) = EVC.Trial(incongruentTrial);
                 end
                 if rewardConds(trl-1) ==1
                    this.taskEnv.Sequence(trl).outcomeValues(1) =  this.reward;
                    this.taskEnv.Sequence(trl).conditions(2) = 1;
                elseif(rewardConds(trl-1) ==0)
                    this.taskEnv.Sequence(trl).outcomeValues(1) =  this.noReward;
                    this.taskEnv.Sequence(trl).conditions(2) = 0;
                 end
                
                end
                
            end

            this.nTrials = length(this.taskEnv.Sequence);
            
                       
       end
                           
       % this function computes the fitting criterion (the lower, the better the simulation results match the simulated results)
       function criterion = getOptimizationCriterion(this)
           
           % PWStroop in Padmala & Pessoa, 2011 (no reward condition)
           
            orgData.RT.Ntr = 0.605;
            orgData.RT.Con = 0.577;
            orgData.RT.Inc = 0.645;

            orgData.Accuracy.Ntr = 0.947;
            orgData.Accuracy.Con = 0.96;
            orgData.Accuracy.Inc = 0.875;
           
           % FIT INTERACTION

           simulation.Intensity = mean(this.results.meanIntensity,2);

           % simulation data
           simulation.RT.Ntr = mean(this.results.RT.Ntr);
           simulation.RT.Con = mean(this.results.RT.Con);
           simulation.RT.Inc = mean(this.results.RT.Inc);

           simulation.Accuracy.Ntr = mean(this.results.Accuracy.Ntr);
           simulation.Accuracy.Con = mean(this.results.Accuracy.Con);
           simulation.Accuracy.Inc = mean(this.results.Accuracy.Inc);
           
           % simulation.Intensity = mean(this.results.meanIntensity,2);

           RT_Criterion =   abs(simulation.RT.Ntr - orgData.RT.Ntr ) + ...
                                   abs(simulation.RT.Con - orgData.RT.Con ) + ...
                                   abs(simulation.RT.Inc - orgData.RT.Inc );

           Accuracy_Criterion =   abs(simulation.Accuracy.Ntr - orgData.Accuracy.Ntr ) + ...
                                   abs(simulation.Accuracy.Con - orgData.Accuracy.Con ) + ...
                                   abs(simulation.Accuracy.Inc - orgData.Accuracy.Inc );

           if(simulation.Intensity == 0)
               intensity_Criterion = 100;
           else
               intensity_Criterion = 0;
           end

           criterion = RT_Criterion + Accuracy_Criterion * 100 + intensity_Criterion;

           disp([RT_Criterion Accuracy_Criterion intensity_Criterion]);                    
           
        end
        
        
    end
    
end

