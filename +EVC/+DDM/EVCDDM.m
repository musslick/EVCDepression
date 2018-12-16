classdef EVCDDM < EVC.EVCModel
   
    properties

        DDMProcesses            % DDMproc[]: description of DDM parametrization
        DDMProcesses_expected
        actualStateParam
        expectedStateParam
        
    end
    
    
    methods 
       
       % constructor
       function this = EVCDDM(CtrlSignals, TaskEnv, DDMProcesses)
          
          % call superclass constructor
          % note: this function automatically calls setActualState()

          this = this@EVC.EVCModel(CtrlSignals, TaskEnv);
          
          % reference model to DDM processes
          this.DDMProcesses = DDMProcesses;
          
          this.DDMProcesses_expected = EVC.DDM.DDMProc.empty(length(DDMProcesses),0);
          
          for i = 1:length(this.DDMProcesses)
             this.DDMProcesses(i).input.EVCModel = this;
             
             this.DDMProcesses_expected(i) = EVC.DDM.DDMProc(this.DDMProcesses(i));
          end
        
          
          
       end
       
       
      
 
        
       % calculates probabilities & RTs for each possible outcome depending on the current
        function [performance]= simulateOutcomes(this)
            [ST,I] = dbstack; % track function name for debugging functions
            source = ST.name;
            
            import EVC.DDM.*;
            
            drift = 0;
            bias = 0;
            ddmp.z = 0;
            ddmp.c = 0;
            ddmp.T0 = 0;
            
            if(this.useActualState)
                stateType = DDMProc.ACTUAL_STATE;
                DDMProcesses_current = this.DDMProcesses;
            else
                stateType = DDMProc.EXPECTED_STATE;
                DDMProcesses_current = this.DDMProcesses_expected;
            end
            
            % log state params
            this.actualStateParam = 0;
            this.expectedStateParam = 0;
            
            for i = 1:length(DDMProcesses_current)
               
                % check if ddm parameter is either control, default or
                % relevant state parameter
                if(ismember(DDMProc.DEFAULT,DDMProcesses_current(i).type) ...
                   || ismember(DDMProc.CONTROL,DDMProcesses_current(i).type) ...
                   || ismember(stateType,DDMProcesses_current(i).type))
                        
                    switch this.DDMProcesses(i).DDMProxy
                        case DDMProc.DRIFT
                             drift = drift + DDMProcesses_current(i).getVal();
                        case DDMProc.THRESH
                             ddmp.z = ddmp.z + DDMProcesses_current(i).getVal();
                        case DDMProc.BIAS
                             bias = bias + DDMProcesses_current(i).getVal();
                        case DDMProc.NOISE
                             ddmp.c = ddmp.c + DDMProcesses_current(i).getVal();
                        case DDMProc.T0
                             ddmp.T0 = ddmp.T0 + DDMProcesses_current(i).getVal();
                    end
                    
                    % debugging % this.State.TaskEnv.CurrTrialIdx == 1 &&
%                     if( ~this.useActualState && this.DDMProcesses(i).DDMProxy == DDMProc.DRIFT)
%                         disp(['DDM Drift Process Name: ' num2str(DDMProcesses_current(i).type) ';Value: ' num2str(DDMProcesses_current(i).getVal())])
%                     end
                    
                end
                
                if(ismember(DDMProc.ACTUAL_STATE, DDMProcesses_current(i).type))
                    this.actualStateParam = this.actualStateParam + DDMProcesses_current(i).getVal();
                end
                if(ismember(DDMProc.EXPECTED_STATE, DDMProcesses_current(i).type))
                    this.expectedStateParam = this.expectedStateParam + DDMProcesses_current(i).getVal();
                end
            end
            
            if(this.useActualState)
                EVC.HelperFnc.disp(source, 'drift', drift, 'ddmp.z', ddmp.z, 'bias', bias, 'ddmp.c', ddmp.c, 'ddmp.T0', ddmp.T0); % DIAGNOSIS
                
%                 if(this.State.TaskEnv.Sequence(this.State.TaskEnv.CurrTrialIdx).conditions(3) == 1 & ...
%                    this.State.TaskEnv.Sequence(this.State.TaskEnv.CurrTrialIdx).conditions(2) == 1)
%                     disp(drift)
%                 end
            end
            
%             % debugging
%             if(this.State.TaskEnv.CurrTrialIdx == 1 && ~this.useActualState)
%                 disp('DDM PARAMS:');
%                 disp(ddmp)
%                 bias
%                 drift
%                 disp('-------');
%             end
            
            % call DDM
            %[tmpER,~,~,~,~,~,allFinalRTs_sepCDFs] = AS_ddmSimFRG(drift,bias,0,ddmp,0,0,null(1),[1], 1); % separate RT's
            % [meanERs,~,~,meanRTs] = AS_ddmSimFRG_Mat(drift,bias,ddmp); % use only mean RT for now; tmpER represents probability of hitting bottom threshold
            [meanERs,meanRTs,~,condRTs,condVarRTs, condSkewRTs] = AS_ddmSimFRG_Mat(drift,bias,ddmp,1);
            
%             if(this.useActualState)
%                 disp('RT check');
%                 disp(ddmp.T0);
%                 disp(meanRTs);
%                 disp(condRTs);
%                 disp(drift);
%                 disp('---');
%             end
            
            EVC.HelperFnc.disp(source, 'condRTs', condRTs); % DIAGNOSIS
            % break point condition: this.useActualState == 0 && this.State.TaskEnv.CurrTrialIdx == 6
            
%             if(this.useActualState)
%                 if(this.State.Actual.conditions(1) == 0)
%                 disp('-- congruent');
%                 elseif(this.State.Actual.conditions(1) == 1)
%                 disp('-- incongruent');    
%                 elseif(this.State.Actual.conditions(1) == 2)
%                     disp('-- neutral');
%                 end
%                 disp(this.State.Actual.typeID);
%                 
% %                 disp(this.State.Actual.stimSalience);
% %                 disp(this.State.Actual.stimRespMap);
%                  disp(this.System.CtrlSignals.Intensity);
% %                 disp(this.DDMProcesses(1).getVal());
%                disp(drift);
%                disp(bias);
%                disp(ddmp);
%                disp(meanERs);
%                disp(meanRTs);
%             end
            
            
            % check for weird outcomes
            meanERs(drift==0) = 1-bias;
            meanRTs(drift==0) = 1e+12;
            
            if(~isreal(meanERs))
                disp('drift, bias, thresh');
                disp([drift bias ddmp.z]);
               error('Something imaginary bad happened.');
            end

            % calculate outcomes
            probs(1,:) = (1-meanERs);
            probs(2,:) = meanERs;
            RTs = condRTs; % [meanRTs; meanRTs];
            
            
            % handle binary errors
            if(this.binaryErrors && this.useActualState)
               number = rand;
               if(number < probs(1))
                   probs(1) = 1;
               else
                   probs(1) = 0;
               end
               probs(2) = 1 - probs(1);
            end
            
            % TODO check this... can't update actual state when using binary errors because log()
            % function calls getEVC() again, so the binary probs would
            % change
%             if(this.useActualState)
%                %this.State.Actual.probs = probs;
%                %this.State.Actual.RTs = RTs; 
%             else
%                this.State.Expected.probs = probs;
%                this.State.Expected.RTs = RTs;  
%             end

            performance.probs = probs;
            
            % RT measures
            performance.RTs = RTs;
            performance.meanRTs = meanRTs;
            performance.varRTs = condVarRTs;
            performance.skewRTs = condSkewRTs;
            performance.drift = drift;
            
            if(this.useActualState)
                EVC.HelperFnc.disp(source, 'tmpER', meanERs, 'RTs', RTs, 'meanRTs', meanRTs); % DIAGNOSIS
            end
            
        end
        

        function initLog(this)
            % call superclass initLog function
            initLog@EVC.EVCModel(this);
            
            this.Log.meanRT = [];
            this.Log.expected_condRTs = [];
            this.Log.actual_condRTs = [];
            this.Log.expected_varRTs = [];
            this.Log.actual_varRTs = [];
            this.Log.expected_skewRTs = [];
            this.Log.actual_skewRTs = [];
            
        end
        
        % logs trial outcomes for current control signals
        function log(this, varargin)
            
            % call superclass log function
            if(length(varargin)>0)
                superclass_args = varargin;
                log@EVC.EVCModel(this, superclass_args{:})
            else 
                log@EVC.EVCModel(this);
            end
            
            % log specific parameters from EVCDDM
            
            % - log actual error rate
            expectedProbs = this.State.Expected.performance.probs;
            actualProbs = this.State.Actual.performance.probs;
            
            if(this.State.Actual.outcomeValues(2) > this.State.Actual.outcomeValues(1)) % if second option is defined as the correct outcome, then 1st column must provide error probability
                ERcol = 1;
            else
                ERcol = 2;
            end
            expectedER = expectedProbs(ERcol);
            actualER = actualProbs(ERcol);
            this.Log.ERs(this.Log.LogCount,:) = [expectedER actualER];
 
            % log RT variance & skew
            
            this.Log.meanRT(this.Log.LogCount,:) = [this.State.Expected.performance.meanRTs this.State.Actual.performance.meanRTs];
            this.Log.expected_condRTs(this.Log.LogCount,:) = this.State.Expected.performance.RTs';
            this.Log.actual_condRTs(this.Log.LogCount,:) = this.State.Actual.performance.RTs';
            this.Log.expected_varRTs(this.Log.LogCount,:) = [this.State.Expected.performance.meanRTs this.State.Actual.performance.meanRTs];
            this.Log.expected_varRTs(this.Log.LogCount,:) = this.State.Expected.performance.varRTs';
            this.Log.actual_varRTs(this.Log.LogCount,:) = this.State.Actual.performance.varRTs';
            this.Log.expected_skewRTs(this.Log.LogCount,:) = this.State.Expected.performance.skewRTs';
            this.Log.actual_skewRTs(this.Log.LogCount,:) = this.State.Actual.performance.skewRTs';
            
            % - log control parameter
            controlProc = EVC.DDM.DDMProc.filterProcess(this.DDMProcesses, EVC.DDM.DDMProc.CONTROL);
            proxyArray = NaN(1,length(controlProc));
            valArray = NaN(1,length(controlProc));
            for i = 1:length(controlProc)
               proxyArray(i) = controlProc(i).DDMProxy;
               valArray(i) = controlProc(i).getVal();
            end
            this.Log.ControlParamType(this.Log.LogCount,:) = proxyArray;
            this.Log.ControlParamVal(this.Log.LogCount,:) = valArray;
            this.Log.drift(this.Log.LogCount,:) = this.State.Actual.performance.drift;
            
            % - log expected state parameter
            %expectedProc = EVC.DDM.DDMProc.filterProcess(this.DDMProcesses, EVC.DDM.DDMProc.EXPECTED_STATE);
            this.Log.ExpectedStateParam(this.Log.LogCount,:) = this.expectedStateParam;
            
            
            % - log actual state parameter
            %actualProc = EVC.DDM.DDMProc.filterProcess(this.DDMProcesses, EVC.DDM.DDMProc.ACTUAL_STATE);
            this.Log.ActualStateParam(this.Log.LogCount,:) = this.actualStateParam;
            
            
        end
        
    end

    
end