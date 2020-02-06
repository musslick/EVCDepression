classdef DDM_Gratton < Simulations.DDMSim
    
    % description of class
    % simulates sequential control adjustments in response to conflict (see
    % Gratton 1992)
    
    % global parameters
    properties
        automaticControlMappingFnc
        automaticControlFnc
        automaticControlProcess
        
        fitMainEffects
        taskSequence
        nBalancedTrials
        congruencyBalance
    end
    
    methods
        
        function this = DDM_Gratton()

            import EVC.*;
            import EVC.DDM.*;
            
            % call parent constructor
            this = this@Simulations.DDMSim();
            
            %% general simulation parameters
            
            this.nSubj = 2;
            this.plotSum = true;
            
            this.learningFnc(1).params{1} = 0.5; 
            this.reconfCostFnc.params{1} = 0; % 10
            this.reconfCostFnc.params{2} = 0; % -4
            this.defaultCostFnc.params{1} = 2.0;                % 2.0 % 4
            this.defaultCostFnc.params{2} = 0;                 
            
            this.rewardFnc.params{2} = 0; % reward RT discount
            
            temp.rewardFnc.params{1} = this.rewardFnc.params{1};
            temp.rewardFnc.params{2} = this.rewardFnc.params{2};
            temp.rewardFnc.params{3} = 1;
            temp.rewardFnc.params{4} = 0;
            temp.rewardFnc.type = EVCFnc.REWRATE_EXT;
            this.rewardFnc = EVCFnc(temp.rewardFnc.type, temp.rewardFnc.params);
               
            this.defaultDDMParams.c = 1.0;                     
            this.defaultDDMParams.thresh = 2;                
            this.defaultDDMParams.T0 = 0.2;                    
            
            %% automatic/control process
           
            % generate 3 control signals, one for the target stimulus & two
            % for the flanker stimuli
            
            % signal range
            tmp.IntensityRange = [0:0.2:10]; % [0:0.2:10]
            
            % target signal
           this.ctrlSignals(1) = CtrlSignal(this.defaultCtrlSignal);
           this.ctrlSignals(1).CtrlSigStimMap  = [0 1];
           this.ctrlSignals(1).IntensityRange = tmp.IntensityRange;
            
           % flanker signal 1
           this.ctrlSignals(2) = CtrlSignal(this.defaultCtrlSignal);
           this.ctrlSignals(2).CtrlSigStimMap  = [1 0];
           this.ctrlSignals(2).IntensityRange = tmp.IntensityRange;
               
            %% define DDM processes
            
            % map all control processes to a specific DDM parameter
            
            % target
            temp.targetControlFnc.type = DDMFnc.INTENSITY2DDM_EFFICACY;
            temp.targetControlFnc.params{1} = this.ctrlSignals(1);
            temp.targetControlFnc.params{2} = this.defaultControlProxy;
            temp.targetControlFnc.params{3} = this.defaultControlMappingFnc;
            temp.targetControlFnc.params{4} = 1;
            temp.targetControlFnc = DDMFnc(temp.targetControlFnc.type, ...
                                            temp.targetControlFnc.params);
            
            % flanker 1                            
            temp.flanker1ControlFnc.type = DDMFnc.INTENSITY2DDM_EFFICACY;
            temp.flanker1ControlFnc.params{1} = this.ctrlSignals(2);
            temp.flanker1ControlFnc.params{2} = this.defaultControlProxy;
            temp.flanker1ControlFnc.params{3} = this.defaultControlMappingFnc;
            temp.flanker1ControlFnc.params{4} = 1;
            temp.flanker1ControlFnc = DDMFnc(temp.flanker1ControlFnc.type, ...
                                            temp.flanker1ControlFnc.params);
                                        
            
            % define each DDM process
            temp.targetProcess = DDMProc(DDMProc.CONTROL, ...                  % default controlled DDM process 
                                                    DDMProc.DRIFT, ...
                                                    temp.targetControlFnc);
           
            temp.flanker1Process = DDMProc(DDMProc.CONTROL, ...                  % default controlled DDM process 
                                        DDMProc.DRIFT, ...
                                        temp.flanker1ControlFnc);

            
            % put all DDM processes together
            this.DDMProcesses = DDMProc.empty(4,0);
            this.DDMProcesses(1) =   this.defaultAutomaticProcess;
            this.DDMProcesses(2) =   temp.targetProcess;
            this.DDMProcesses(3) =   temp.flanker1Process;    
            
            %% task environment parameters: task environment
            
            this.congruencyBalance = [0 1];
            this.nBalancedTrials = 100;
%             this.nTrials = 50;    % this variable is set below (after adding filler trials)
            
            temp.stimSalience = [1 1];
            temp.outcomeValues = [100 0];
            w = 1;
            
            % congruent trial,
            % e.g. [ < < < ]
            this.trials(1).ID = 1;                                                          % trial identification number (for task set)
            this.trials(1).typeID = 1;                                                      % trial type (defines task context)
            this.trials(1).cueID = 1;                                                       % cued information about trial identity
            this.trials(1).descr = 'con';                                                 % trial description
            this.trials(1).conditions    = [0 0];                                           % set of trial conditions (for logging)
            this.trials(1).outcomeValues = temp.outcomeValues;                  % reward for correct outcome = 3; no reward/punishment for incorrect outcome
            this.trials(1).stimSalience  = temp.stimSalience;                         % relative salience between stimulus 1 and stimulus 2 defines level of incongruency; here stim 2 is more dominant (like in stroop)
            this.trials(1).stimRespMap   = [w 0;                                          % responding to stimulus 1 tends to produce second response by 100%
                                                        w 0];                                         % responding to stimulus 2 tends to produce second outcome by 100%
            this.trials(1).params = [];                                                     % DDM specific trial parameters

            
            % incongruent trial,
            % e.g. [ < < < ]
            this.trials(2).ID = 1;                                                          % trial identification number (for task set)
            this.trials(2).typeID = 1;                                                      % trial type (defines task context)
            this.trials(2).cueID = 1;                                                       % cued information about trial identity
            this.trials(2).descr = 'inc';                                                 % trial description
            this.trials(2).conditions    = [1 0];                                           % set of trial conditions (for logging)
            this.trials(2).outcomeValues = temp.outcomeValues;                  % reward for correct outcome = 3; no reward/punishment for incorrect outcome
            this.trials(2).stimSalience  = temp.stimSalience;                         % relative salience between stimulus 1 and stimulus 2 defines level of incongruency; here stim 2 is more dominant (like in stroop)
            this.trials(2).stimRespMap   = [0 w;                                          % responding to stimulus 1 tends to produce second response by 100%
                                                        w 0];                                         % responding to stimulus 2 tends to produce second outcome by 100%
            this.trials(2).params = [];                                                     % DDM specific trial parameters
            
            %% log parameters
            
            this.fitMainEffects = 1;
            
            this.writeLogFile = 1; 
            this.logFileName = 'DDM_Gratton'; 
            
            this.logAddVars{3} = '[this.EVCM.Log.Trials.conditions]''';
            this.logAddVarNames{3} = 'condition';
        end
        
        function getResults(this)
            this.results.RT.conINC = zeros(1, this.nSubj);
            this.results.RT.conCON = zeros(1, this.nSubj);
            this.results.RT.incCON = zeros(1, this.nSubj);
            this.results.RT.incINC = zeros(1, this.nSubj);
            this.results.RT.con = zeros(1, this.nSubj);
            this.results.RT.inc = zeros(1, this.nSubj);
            this.results.RT.adaptEffect = zeros(1, this.nSubj);
            this.results.RT.congrEffect = zeros(1, this.nSubj);
            
            this.results.ER.con = zeros(1, this.nSubj);
            this.results.ER.inc = zeros(1, this.nSubj);
            
            this.results.analysis.RTData = [];
            this.results.analysis.ERData = [];
            this.results.analysis.previousCongruency = [];
            this.results.analysis.currentCongruency = [];
            this.results.analysis.subject = [];
            
            % loop through all subjects
            for subj = 1:this.nSubj
                
                 % get log structure for current subject
                 subjLog = this.subjData(subj).Log;
                
                 % get conditions
                 conditions = reshape([subjLog.TrialsOrg.conditions],[2 this.nTrials])';
                 prevInc = conditions(:,2);
                 currInc = conditions(:,1);
                 
                 % extract RT's
                 RT = this.subjData(subj).Log.RTs(:,2)';
                 ER = this.subjData(subj).Log.ERs(:,2)';
                 
                 this.results.meanIntensity(subj, :) = mean(subjLog.CtrlIntensities(:));
                 targetIntensity = subjLog.CtrlIntensities(:,1);
                 flankerIntensity = subjLog.CtrlIntensities(:,2);
                 this.results.targetIntensity.mean(subj, :) = targetIntensity;
                 this.results.flankerIntensity.mean(subj, :) = flankerIntensity;
                 
                 % fetch analysis data
                 analysisTrials = 2:length(RT);
                this.results.analysis.RTData = [this.results.analysis.RTData; RT(analysisTrials)'];
                this.results.analysis.ERData = [this.results.analysis.ERData; ER(analysisTrials)'];
                this.results.analysis.previousCongruency = [this.results.analysis.previousCongruency; (prevInc(analysisTrials)+1)];
                this.results.analysis.currentCongruency = [this.results.analysis.currentCongruency; (currInc(analysisTrials)+1)];
                this.results.analysis.subject = [this.results.analysis.subject; repmat(subj, length((analysisTrials)), 1)];
                 
                 % get results
                 this.results.meanRT(subj) = mean(RT);
                 this.results.RT.conINC(subj) = mean(RT(prevInc == 0 & currInc == 1));
                 this.results.RT.conCON(subj) = mean(RT(prevInc == 0 & currInc == 0));
                 this.results.RT.incCON(subj) = mean(RT(prevInc == 1 & currInc == 0));
                 this.results.RT.incINC(subj) = mean(RT(prevInc == 1 & currInc == 1));
                 this.results.RT.con(subj) = mean(RT(currInc == 0));
                 this.results.RT.inc(subj) = mean(RT(currInc == 1));
                 this.results.RT.congrEffect(subj) = mean(RT(currInc == 1)) - mean(RT(currInc == 0));
                 this.results.RT.adaptEffect(subj) = (this.results.RT.conINC(subj) - this.results.RT.conCON(subj)) - (this.results.RT.incINC(subj) - this.results.RT.incCON(subj));
                 
                 this.results.meanER(subj) = mean(ER);
                 this.results.ER.conINC(subj) = mean(ER(prevInc == 0 & currInc == 1));
                 this.results.ER.conCON(subj) = mean(ER(prevInc == 0 & currInc == 0));
                 this.results.ER.incCON(subj) = mean(ER(prevInc == 1 & currInc == 0));
                 this.results.ER.incINC(subj) = mean(ER(prevInc == 1 & currInc == 1));
                 this.results.ER.con(subj) = mean(ER(currInc == 0));
                 this.results.ER.inc(subj) = mean(ER(currInc == 1));
                 this.results.ER.congrEffect(subj) = mean(ER(currInc == 1)) - mean(ER(currInc == 0));
                 this.results.ER.adaptEffect(subj) = (this.results.ER.conINC(subj) - this.results.ER.conCON(subj)) - (this.results.ER.incINC(subj) - this.results.RT.incCON(subj));
                 
                 this.results.targetIntensity.prevCON(subj) = mean(targetIntensity(prevInc == 0));
                 this.results.targetIntensity.prevINC(subj) = mean(targetIntensity(prevInc == 1));
                 this.results.flankerIntensity.prevCON(subj) = mean(flankerIntensity(prevInc == 0));
                 this.results.flankerIntensity.prevINC(subj) = mean(flankerIntensity(prevInc == 1));
                 
            end
            
             % perform 2x2 repeated measures ANOVA

             this.results.ANOVA_ER.main.stats = rm_anova2(this.results.analysis.ERData, ...
                                                                            this.results.analysis.subject, ...
                                                                            this.results.analysis.previousCongruency, ...
                                                                            this.results.analysis.currentCongruency, ...
                                                                            {'prevCongruent', 'currCongruent'});

            this.results.ANOVA_RT.main.stats = rm_anova2(this.results.analysis.RTData, ...
                                                                                this.results.analysis.subject, ...
                                                                                this.results.analysis.previousCongruency, ...
                                                                                this.results.analysis.currentCongruency, ...
                                                                                {'prevCongruent', 'currCongruent'});

            [h,p,ci,stats] = ttest(this.results.RT.conINC, this.results.RT.incINC, 'Tail', 'right');
            this.results.ttest_RT.incAdaptation.p = p;
            this.results.ttest_RT.incAdaptation.df = stats.df;
            this.results.ttest_RT.incAdaptation.tstat = stats.tstat;

            [h,p,ci,stats] = ttest(this.results.RT.conCON, this.results.RT.incCON, 'Tail', 'left');
            this.results.ttest_RT.conAdaptation.p = p;
            this.results.ttest_RT.conAdaptation.df = stats.df;
            this.results.ttest_RT.conAdaptation.tstat = stats.tstat;

            [h,p,ci,stats] = ttest(this.results.ER.conINC, this.results.ER.incINC, 'Tail', 'right');
            this.results.ttest_ER.incAdaptation.p = p;
            this.results.ttest_ER.incAdaptation.df = stats.df;
            this.results.ttest_ER.incAdaptation.tstat = stats.tstat;

            [h,p,ci,stats] = ttest(this.results.ER.conCON, this.results.ER.incCON, 'Tail', 'left');
            this.results.ttest_ER.conAdaptation.p = p;
            this.results.ttest_ER.conAdaptation.df = stats.df;
            this.results.ttest_ER.conAdaptation.tstat = stats.tstat;
            
        end
        
        function dispResults(this)
            
              mc = metaclass(this);
              disp(['++++++++++ ' mc.Name ' ++++++++++']);
              disp(strcat('congruency Effect: ',num2str(mean(this.results.RT.congrEffect)),'ms'));
              disp(strcat('adaptation Effect: ',num2str(mean(this.results.RT.adaptEffect)),'ms'));
              
               % inference statistics
               disp('---');
               disp(['RT ANOVA (2x2), previous congruency, F(' num2str(this.results.ANOVA_RT.main.stats{2,3}) ... % df for interaction = df(condition 1) * df(condition 2)
                                                                                       ', ' num2str(this.results.ANOVA_RT.main.stats{5,3}) ...    
                                                                                       ') = ' num2str(this.results.ANOVA_RT.main.stats{2,end-1}) ...
                                                                                        ', p = ' num2str(this.results.ANOVA_RT.main.stats{2,end})]);
               disp(['RT ANOVA (2x2), current congruency, F(' num2str(this.results.ANOVA_RT.main.stats{3,3}) ... % df for interaction = df(condition 1) * df(condition 2)
                                                                                       ', ' num2str(this.results.ANOVA_RT.main.stats{6,3}) ...    
                                                                                       ') = ' num2str(this.results.ANOVA_RT.main.stats{3,end-1}) ...
                                                                                        ', p = ' num2str(this.results.ANOVA_RT.main.stats{3,end})]);
               disp(['RT ANOVA (2x2), previous x current congruency, F(' num2str(this.results.ANOVA_RT.main.stats{4,3}) ... % df for interaction = df(condition 1) * df(condition 2)
                                                                                       ', ' num2str(this.results.ANOVA_RT.main.stats{7,3}) ...    
                                                                                       ') = ' num2str(this.results.ANOVA_RT.main.stats{4,end-1}) ...
                                                                                        ', p = ' num2str(this.results.ANOVA_RT.main.stats{4,end})]);
               disp(['ER ANOVA (2x2), previous congruency, F(' num2str(this.results.ANOVA_ER.main.stats{2,3}) ... % df for interaction = df(condition 1) * df(condition 2)
                                                                                       ', ' num2str(this.results.ANOVA_ER.main.stats{5,3}) ...    
                                                                                       ') = ' num2str(this.results.ANOVA_ER.main.stats{2,end-1}) ...
                                                                                        ', p = ' num2str(this.results.ANOVA_ER.main.stats{2,end})]);
               disp(['ER ANOVA (2x2), current congruency, F(' num2str(this.results.ANOVA_ER.main.stats{3,3}) ... % df for interaction = df(condition 1) * df(condition 2)
                                                                                       ', ' num2str(this.results.ANOVA_ER.main.stats{6,3}) ...    
                                                                                       ') = ' num2str(this.results.ANOVA_ER.main.stats{3,end-1}) ...
                                                                                        ', p = ' num2str(this.results.ANOVA_ER.main.stats{3,end})]);
               disp(['ER ANOVA (2x2), previous x current congruency, F(' num2str(this.results.ANOVA_ER.main.stats{4,3}) ... % df for interaction = df(condition 1) * df(condition 2)
                                                                                       ', ' num2str(this.results.ANOVA_ER.main.stats{7,3}) ...    
                                                                                       ') = ' num2str(this.results.ANOVA_ER.main.stats{4,end-1}) ...
                                                                                        ', p = ' num2str(this.results.ANOVA_ER.main.stats{4,end})]);
                                                                                    
               disp('--- single contrasts: ')
               disp(['RT t-test (conINC vs. incINC), t(' num2str(this.results.ttest_RT.incAdaptation.df) ...
                                                                        ') = ' num2str(this.results.ttest_RT.incAdaptation.tstat) ...
                                                                        ', p = ' num2str(this.results.ttest_RT.incAdaptation.p)]);
                disp(['RT t-test (conCON vs. incCON), t(' num2str(this.results.ttest_RT.conAdaptation.df) ...
                                                                        ') = ' num2str(this.results.ttest_RT.conAdaptation.tstat) ...
                                                                        ', p = ' num2str(this.results.ttest_RT.conAdaptation.p)]);
                disp(['ER t-test (conINC vs. incINC), t(' num2str(this.results.ttest_ER.incAdaptation.df) ...
                                                                        ') = ' num2str(this.results.ttest_ER.incAdaptation.tstat) ...
                                                                        ', p = ' num2str(this.results.ttest_ER.incAdaptation.p)]);
                disp(['ER t-test (conCON vs. incCON), t(' num2str(this.results.ttest_ER.conAdaptation.df) ...
                                                                        ') = ' num2str(this.results.ttest_ER.conAdaptation.tstat) ...
                                                                        ', p = ' num2str(this.results.ttest_ER.conAdaptation.p)]);
               

                disp('--- means: ');
                
                orgData.meanRT_incInc = 0.360;
                orgData.meanRT_conInc = 0.368;
                orgData.meanRT_conCon = 0.312;
                orgData.meanRT_incCon = 0.323;

                orgData.meanER_incInc = 0.22;
                orgData.meanER_conInc = 0.32;
                orgData.meanER_conCon = 0.06;
                orgData.meanER_incCon = 0.06;
                
                disp(['congruent RT conINC: ' num2str(mean(this.results.RT.conINC)) ' vs. ' num2str(orgData.meanRT_conInc) ]);
                disp(['congruent RT incINC: ' num2str(mean(this.results.RT.incINC)) ' vs. ' num2str(orgData.meanRT_incInc) ]);
                disp(['congruent RT conCON: ' num2str(mean(this.results.RT.conCON)) ' vs. ' num2str(orgData.meanRT_conCon) ]);
                disp(['congruent RT incCON: ' num2str(mean(this.results.RT.incCON)) ' vs. ' num2str(orgData.meanRT_incCon) ]);

                disp(['congruent ER conINC: ' num2str(mean(this.results.ER.conINC)) ' vs. ' num2str(orgData.meanER_conInc) ]);
                disp(['congruent ER incINC: ' num2str(mean(this.results.ER.incINC)) ' vs. ' num2str(orgData.meanER_incInc) ]);
                disp(['congruent ER conCON: ' num2str(mean(this.results.ER.conCON)) ' vs. ' num2str(orgData.meanER_conCon) ]);
                disp(['congruent ER incCON: ' num2str(mean(this.results.ER.incCON)) ' vs. ' num2str(orgData.meanER_incCon) ]);

        end
        
        function plotSummary(this) 
            
%             exampleSubj = 1;
%             sampleTrials = 1:this.nTrials;
%             f1 = figure(1);
%             set(f1, 'Position', [0 0 600 500])
%             subplot(4,1,1);
%             this.plotER(exampleSubj, 'actual', sampleTrials);
%             subplot(4,1,2);
%             this.plotRT(exampleSubj, 'actual', sampleTrials);
%             subplot(4,1,3);
%             %this.plotEVC(exampleSubj, 'expected', sampleTrials);
%             this.plotDifficulty(exampleSubj, sampleTrials);
%             subplot(4,1,4);
%             this.plotCtrlIntensity(exampleSubj, sampleTrials);
%             ylim([0 0.7]);
%             
%             figure(2);
%             y = [mean(this.results.RT.conCON), mean(this.results.RT.conINC), mean(this.results.RT.incCON), mean(this.results.RT.incINC)];
%             h = bar(y, 'BarWidth', 1,'FaceColor',[0.7 0.7 0.7]);
%             ylim([0, max(y)*1.1]);
%             set(gca,'XTickLabel',{'c-C', 'c-I','i-C','i-I'},'fontsize',this.plotParams.axisFontSize);
%             ylabel('Reaction time (ms)','fontWeight', 'bold', 'fontSize',16);
            
               
        end
        
        
    end
    
    methods (Access = protected)
        
        % generate task environment
        function initTaskEnv(this)
            
            % CREATE TRIALS
           congruentTrial = EVC.Trial(this.trials(1));
           incongruentTrial = EVC.Trial(this.trials(2));
            
           % BALANCE SEQUENCE
           congruency = this.congruencyBalance;
           congruencyTransition = [0 1];

           congruencyIdx = 1;
           congrTransitionIdx = 2;
           
           % counterbalance previous and current congruency
           if(length(congruency) == 2)
               

               trialCombs = Simulations.combvec(congruency, congruencyTransition);
               sequenceRepetitions = ceil(this.nBalancedTrials/size(trialCombs,2));
               trialPool = repmat(trialCombs, 1, sequenceRepetitions);

               trialSequence = [];

               % GENERATE SEQUENCE

               while size(trialSequence, 2) < size(trialPool, 2)

                   selectedTrialPool = trialPool;
                   trialSequence = [];

                   for i = 1:size(trialPool,2)

                       if(i > 1)
                           goodOptions = find( (selectedTrialPool(congruencyIdx,:) == trialSequence(congruencyIdx,i-1) & selectedTrialPool(congrTransitionIdx,:) == 0) ...
                                                        | (selectedTrialPool(congruencyIdx,:) ~= trialSequence(congruencyIdx,i-1) & selectedTrialPool(congrTransitionIdx,:) == 1));
                       else
                           goodOptions = 1:size(selectedTrialPool, 2);

                       end

                       if(isempty(goodOptions))
                           break;
                       end
                       sample = goodOptions(randsample(length(goodOptions), 1));

                       trialSequence = [trialSequence selectedTrialPool(:, sample)];
                       selectedTrialPool(:, sample) = [];

                   end

               end
           
               % add first trial

               if(trialSequence(congrTransitionIdx, 1) == 0)
                   trial = [trialSequence(congruencyIdx, 1)  ; 0];
               else
                   trial = 1- [trialSequence(congruencyIdx, 1)  ; 0];
               end

           else
               trialSequence = repmat(congruency, 1, this.nBalancedTrials/length(congruency));
               trialSequence = trialSequence(randperm(length(trialSequence)));
               trial = randsample(2,1) - 1;
           end

           
           
           trialSequence = [trial trialSequence];

           % CREATE TASK ENVIRONMENT
           this.taskSequence = trialSequence;
           % generate task environment
            this.taskEnv = EVC.TaskEnv(congruentTrial, size(trialSequence,2));
            for trl = 1:size(trialSequence,2)
                % assign trial
                if(trialSequence(congruencyIdx, trl) == 0)
                    this.taskEnv.Sequence(trl) = EVC.Trial(congruentTrial);
                else
                    this.taskEnv.Sequence(trl) = EVC.Trial(incongruentTrial);
                end
                % determine condition for previous trial type
                if(trl > 1)
                        this.taskEnv.Sequence(trl).conditions(2) = this.taskEnv.Sequence(trl-1).conditions(1);
                else
                        this.taskEnv.Sequence(trl).conditions(2) = -1;
                end
            end

            this.nTrials = length(this.taskEnv.Sequence);
            
        end
        
    end
    
    methods (Access = public)
        
       function initOptimizationTaskEnv(this)
           
           % create trials
           congruentTrial = EVC.Trial(this.trials(1));
           incongruentTrial = EVC.Trial(this.trials(2));
           
%            disp(this.trials(1).stimRespMap);
%            disp(congruentTrial.stimRespMap);
                
           if(this.fitMainEffects)

                % generate task environment
                this.taskEnv = EVC.TaskEnv(congruentTrial, 6);

                % set specified sequence
                this.taskEnv.Sequence(1) = congruentTrial; 
    %             this.taskEnv.Sequence(2) = incongruentTrial; 
    %             this.taskEnv.Sequence(3) = congruentTrial; 
    %             this.taskEnv.Sequence(4) = congruentTrial; 
    %             this.taskEnv.Sequence(5) = incongruentTrial; 
    %             this.taskEnv.Sequence(6) = incongruentTrial; 

                this.nTrials = length(this.taskEnv.Sequence);

           else
           
               %sequenceRepetitions = 2;
               
               congruency = [0 1];
               congruencyTransition = [0 1];
               
               congruencyIdx = 1;
               congrTransitionIdx = 2;

               trialCombs = Simulations.combvec(congruency, congruencyTransition);
               sequenceRepetitions = ceil(this.nTrials/size(trialCombs,2));
               trialPool = repmat(trialCombs, 1, sequenceRepetitions);
               
               trialSequence = [];
               
               % GENERATE SEQUENCE
               
               while size(trialSequence, 2) < size(trialPool, 2)
                   
                   selectedTrialPool = trialPool;
                   trialSequence = [];
                   
                   for i = 1:size(trialPool,2)
                       
                       if(i > 1)
                           goodOptions = find( (selectedTrialPool(congruencyIdx,:) == trialSequence(congruencyIdx,i-1) & selectedTrialPool(congrTransitionIdx,:) == 0) ...
                                                        | (selectedTrialPool(congruencyIdx,:) ~= trialSequence(congruencyIdx,i-1) & selectedTrialPool(congrTransitionIdx,:) == 1));
                       else
                           goodOptions = 1:size(selectedTrialPool, 2);
                       
                       end
                       
                       if(isempty(goodOptions))
                           break;
                       end
                       sample = goodOptions(randsample(length(goodOptions), 1));
                       
                       trialSequence = [trialSequence selectedTrialPool(:, sample)];
                       selectedTrialPool(:, sample) = [];
                           
                   end
                  
               end
               
               % add first trial
               
               if(trialSequence(congrTransitionIdx, 1) == 0)
                   trial = [trialSequence(congruencyIdx, 1)  ; 0];
               else
                   trial = 1- [trialSequence(congruencyIdx, 1)  ; 0];
               end
               
               trialSequence = [trial trialSequence];
               
               % CREATE TASK ENVIRONMENT
               this.taskSequence = trialSequence;
               % generate task environment
                this.taskEnv = EVC.TaskEnv(congruentTrial, size(trialSequence,2));
                for trl = 1:size(trialSequence,2)
                    % assign trial
                    if(trialSequence(congruencyIdx, trl) == 0)
                        this.taskEnv.Sequence(trl) = EVC.Trial(congruentTrial);
                    else
                        this.taskEnv.Sequence(trl) = EVC.Trial(incongruentTrial);
                    end
                    % determine condition for previous trial type
                    if(trl > 1)
                            this.taskEnv.Sequence(trl).conditions(2) = this.taskEnv.Sequence(trl-1).conditions(1);
                    else
                            this.taskEnv.Sequence(trl).conditions(2) = -1;
                    end
                end
                
                this.nTrials = length(this.taskEnv.Sequence);
                
           end
          
       end
                           

       function criterion = getOptimizationCriterion(this)
           
           % Laming 1979, Fig. 1, alternate stimulus
           orgData.meanRT_con = 0.320;
           orgData.meanRT_inc = 0.360;

           orgData.meanER_con = 0.06;
           orgData.meanER_inc = 0.26;
           
           orgData.meanRT_incInc = 0.360;
           orgData.meanRT_conInc = 0.368;
           orgData.meanRT_conCon = 0.312;
           orgData.meanRT_incCon = 0.323;
           
           orgData.meanER_incInc = 0.22;
           orgData.meanER_conInc = 0.32;
           orgData.meanER_conCon = 0.06;
           orgData.meanER_incCon = 0.06;
           
           % FIT MAIN EFFECTS
           if(this.fitMainEffects)

               % simulation data
               simulation.RTs = mean(this.results.meanRT);
               simulation.ERs = mean(this.results.meanER);
               simulation.Intensity = mean(this.results.meanIntensity,2);

               % calculate optimization criteria
               if(strcmp(this.trials(1).descr, 'con'))
                   ER_Criterion = abs(sum(simulation.ERs - orgData.meanER_con));
                   RT_Criterion = abs(sum(simulation.RTs - orgData.meanRT_con));
               else
                   ER_Criterion = abs(sum(simulation.ERs - orgData.meanER_inc));
                   RT_Criterion = abs(sum(simulation.RTs - orgData.meanRT_inc));
                    if( simulation.ERs >= 0.5)
                        ER_Criterion = ER_Criterion + (simulation.ERs - 0.5) * 100;
                    end

               end

               if(simulation.Intensity == 0)
                   intensity_Criterion = 100;
               else
                   intensity_Criterion = 0;
               end

               criterion = RT_Criterion + ER_Criterion * 100 + intensity_Criterion;
           
           else
               % FIT INTERACTION
               
               simulation.Intensity = mean(this.results.meanIntensity,2);
               
               % simulation data
               simulation.RTs.conInc = mean(this.results.RT.conINC);
               simulation.RTs.conCon = mean(this.results.RT.conCON);
               simulation.RTs.incInc = mean(this.results.RT.incINC);
               simulation.RTs.incCon = mean(this.results.RT.incCON);
               
               simulation.ERs.conInc = mean(this.results.ER.conINC);
               simulation.ERs.conCon = mean(this.results.ER.conCON);
               simulation.ERs.incInc = mean(this.results.ER.incINC);
               simulation.ERs.incCon = mean(this.results.ER.incCON);
               % simulation.Intensity = mean(this.results.meanIntensity,2);
               
               RT_Criterion =   abs(simulation.RTs.conInc - orgData.meanRT_conInc) + ...
                                       abs(simulation.RTs.incInc - orgData.meanRT_incInc) + ...
                                       abs(simulation.RTs.conCon - orgData.meanRT_conCon) + ...
                                       abs(simulation.RTs.incCon - orgData.meanRT_incCon);
                                   
               ER_Criterion =   abs(simulation.ERs.conInc - orgData.meanER_conInc) + ...
                                       abs(simulation.ERs.incInc - orgData.meanER_incInc) + ...
                                       abs(simulation.ERs.conCon - orgData.meanER_conCon) + ...
                                       abs(simulation.ERs.incCon - orgData.meanER_incCon); 
                                   
               if(simulation.Intensity == 0)
                   intensity_Criterion = 100;
               else
                   intensity_Criterion = 0;
               end
               
               criterion = RT_Criterion + ER_Criterion * 100 + intensity_Criterion;
               
               disp([RT_Criterion ER_Criterion intensity_Criterion]);                    
           end
           
        end
        
        
    end
    
end

