classdef DDMFnc < EVC.EVCFnc
    
    % this class implements DDM-specific functions
    % most of the functions may reference an EVCDDM instance

    
    properties (Constant)
        STIMBIASRANGE = 1;           % map saliency to range
        INTENSITY2DDM = 2;
        DDMRANDOMUNIFORM = 3;
        INTENSITY2DDM_WEIGHTED = 4;
        CONTROLAUTOMATIC1 = 5;
        SPOTLIGHTCONTROL = 6;
        STIMBIAS = 7;
        INTENSITY2DDM_EXPEFFICACY = 8;
        
        % DDM poxies
        DRIFT = EVC.DDM.DDMProc.DRIFT;
        THRESH = 2;
        BIAS = 3;
        NOISE = 4;
        T0 = 5;
        
        % holds amount of required parameters for each function type
        paramReqDDM = [1, ...   STIMBIASRANGERANGE: 1) lower and upper bound of DDM parameter range, i.e. [min max]
                       3, ... INTENSITY2DDM:
                       1, ... RANDOM UNIFORM: f(x) = rand*(p1(2) - p1(1)) + p1(1)
                       3, ... INTENSITY2DDM_WEIGHTED:
                       2, ... CONTROLAUTOMATIC1:
                       2, ... CONTROLAUTOMATIC1:
                       0, ... STIMBIAS
                       4, ... INTENSITY2DDM_EXPEFFICACY
                       ];         
    end
    
    methods
        
        function this = DDMFnc(type, params, varargin)
            
          % call superclass constructor to create an EVC function instance
          this = this@EVC.EVCFnc(type, params, 'paramReq', 'paramReqDDM');
          
          [inputVal EVCM] = this.extractInput(varargin);
          this.input = inputVal;
          this.EVCModel = EVCM;
        end
        
    end
    
    methods (Access = public)
        
        % calculates output value dependent on specified function type
        function out = getVal(this, varargin)

            [inputVal, EVCM] = this.extractInput(varargin);
            
            if(this.type == this.STIMBIASRANGE)
               out = getStimBiasRange(this, EVCM); 
            end
            
            if(this.type == this.INTENSITY2DDM)
               out =  intensityToDDM(this, EVCM);
            end
                        
            if(this.type == this.DDMRANDOMUNIFORM)
                out = getRandomUniformOut(this, EVCM);
            end
            
            if(this.type == this.INTENSITY2DDM_WEIGHTED)
               out =  intensityToDDM_Weighted(this, EVCM);
            end
            
            if(this.type == this.CONTROLAUTOMATIC1)
               out =  getControlAutomatic(this, EVCM);
            end
                
            if(this.type == this.SPOTLIGHTCONTROL)
               out =  spotlightControl(this, EVCM);
            end
            
            if(this.type == this.STIMBIAS)
               out = getStimBias(this, EVCM); 
            end
            
            if(this.type == this.INTENSITY2DDM_EXPEFFICACY)
               out = intensityToDDM_ExpEfficacy(this, EVCM); 
            end
        end
       
    end
    
    methods (Access = private)
        
        %% function definitions
        
        % DDMRANDOMUNIFORM... unform random value
        function out = getRandomUniformOut(this, EVCM)
            
            interval = this.params{1}; 
            
            if(EVCM.useActualState)
                out = rand*(interval(2) - interval(1)) + interval(1);
            else
                out = (interval(2) - interval(1))/2 + interval(1);
            end
            
        end
   
        % STIMBIASRANGE... returns a response bias for the specified interval
        % param{1}: [b1 b2] (lower and upper bound of output DDM parameter range)
        function out = getStimBiasRange(this, EVCM)
            
            % determine which state to operate on
            if(EVCM.useActualState)
                currTrial = EVCM.State.Actual;
            else
                currTrial = EVCM.State.Expected;
            end

            StimSalRespMap = currTrial.getRespSaliencyMap();
            
            salR1 = sum(StimSalRespMap(:,1));
            % TODO: maybe just use sum(StimSalRespMap(:,2)); and ignore other
            % responses?
            salR2 = sum(sum((StimSalRespMap(:,2:length(StimSalRespMap(1,:)))))); %FIX THIS: all other responses than R1 may belong to the lower threshold
            
            rangeWidth = this.params{1}(2) - this.params{1}(1);

            % how much bias toward upper threshold?
            R1Delta = + rangeWidth/2*salR1;
            
            % how much bias toward lower threshold?
            R2Delta = -rangeWidth/2*salR2;
            
            middleOfRange = this.params{1}(1) + (rangeWidth)/2;

            %EVC.HelperFnc.disp(source, 'middleOfRange', middleOfRange, 'rangeWidth', rangeWidth, 'R1Delta', R1Delta, 'R2Delta', R2Delta); % DIAGNOSIS 
            
            out = middleOfRange + R1Delta + R2Delta;
           
        end
        
        % STIMBIAS... returns an additive response bias
        % no params required
        function out = getStimBias(this, EVCM)
            
            [ST,I] = dbstack; % track function name for debugging functions
            source = ST.name;
            
            % determine which state to operate on
            if(EVCM.useActualState)
                currTrial = EVCM.State.Actual;
            else
                currTrial = EVCM.State.Expected;
            end

            StimSalRespMap = currTrial.getRespSaliencyMap();
            
            salR1 = sum(StimSalRespMap(:,1));
            % TODO: maybe just use sum(StimSalRespMap(:,2)); and ignore other
            % responses?
            salR2 = sum(sum((StimSalRespMap(:,2:length(StimSalRespMap(1,:)))))); %FIX THIS: all other responses than R1 may belong to the lower threshold

            out = salR1 - salR2;
            
                        
            if(EVCM.useActualState)
                EVC.HelperFnc.disp(source, 'salR1', salR1, 'salR2', salR2); % DIAGNOSIS
            end
           
        end
        
        % INTENSITY2DDM... maps control signal intensity to DDM parameter
        % param{1}: EVC.CtrlSignal
        % param{2}: integer (refers to DDM proxy, see constants)
        % param{3}: EVC.EVCFnc (mapping function)
        function out = intensityToDDM(this, EVCM)

            % current control signal
            if(~isa(this.params{1}, 'EVC.CtrlSignal'))
               error('First EVCFnc parameter has to be an instance of EVC.CtrlSignal (control signal instnace with corresponding intensity)');
            else
                currSignal = this.params{1};
            end
            
            % mapping function
            if(~isa(this.params{3}, 'EVC.EVCFnc'))
               error('Second EVCFnc parameter has to be an instance of EVC.EVCFnc (intensity to DDM mapping function)');
            else
                mappingFnc = this.params{3};
            end
            
            % for drift & bias parameters, calculate in which direction control acts
            % (+1 -> upper threshold, -1 -> lower threshold)
            if(this.params{2} == this.DRIFT || this.params{2} == this.BIAS)
            
                % how many stimuli are available?
                ctrlSigStimMap = currSignal.CtrlSigStimMap;
                stimRespMap = EVCM.currentTaskSet();
                [nStimuliEnv,~] = size(stimRespMap);
                [nStimuliCtrl] = length(ctrlSigStimMap);
                
                % check if maps correspond in size
                if(nStimuliCtrl < nStimuliEnv)
                    ctrlSigStimMap = [ctrlSigStimMap repmat(0,1,nStimuliEnv-nStimuliCtrl) ];
                    warning('CtrlSigStimMap and stimRespMap don''t match. Filling missing parts with zero.');
                end
            
                % 
                if(nStimuliCtrl > nStimuliEnv)
                    ctrlSigStimMap = ctrlSigStimMap(1:nStimuliEnv);
                    warning('CtrlSigStimMap and stimRespMap don''t match. Cutting off excess of Stimuli.');
                end
                
                % translate both maps into control-to-response map
                CtrlIToResp = transpose(stimRespMap) * transpose(ctrlSigStimMap);
                
                % make sure that there are only 2 possible responses (DDM default)
                if(length(CtrlIToResp) < 2)
                    CtrlIToResp = [1 0]; % by default the current control intensity acts in favor of the upper response
                end
                
                % translate control-to-response map into scalar value
                % (positive: upper response; negative: bottom response)
                CtrlDirection = CtrlIToResp(1) - sum(CtrlIToResp(2:length(CtrlIToResp)));           
            else
               CtrlDirection = 1;
            end
            
            out = CtrlDirection * mappingFnc.getVal(currSignal.getIntensity());
        end
        
        % INTENSITY2DDM... maps control signal intensity to DDM parameter
        % param{1}: EVC.CtrlSignal
        % param{2}: integer (refers to DDM proxy, see constants)
        % param{3}: EVC.EVCFnc (mapping function)
        function out = intensityToDDM_Weighted(this, EVCM)
            
            % current control signal
            if(~isa(this.params{1}, 'EVC.CtrlSignal'))
               error('First EVCFnc parameter has to be an instance of EVC.CtrlSignal (control signal instnace with corresponding intensity)');
            else
                currSignal = this.params{1};
            end
            
            % mapping function
            if(~isa(this.params{3}, 'EVC.EVCFnc'))
               error('Second EVCFnc parameter has to be an instance of EVC.EVCFnc (intensity to DDM mapping function)');
            else
                mappingFnc = this.params{3};
            end
            
            % for drift & bias parameters, calculate in which direction control acts
            % (+1 -> upper threshold, -1 -> lower threshold)
            if(this.params{2} == this.DRIFT || this.params{2} == this.BIAS)
            
                % how many stimuli are available?
                ctrlSigStimMap = currSignal.CtrlSigStimMap;
                stimRespMap = EVCM.currentTaskSet();
                [nStimuliEnv,~] = size(stimRespMap);
                [nStimuliCtrl] = length(ctrlSigStimMap);
                
                % check if maps correspond in size
                if(nStimuliCtrl < nStimuliEnv)
                    ctrlSigStimMap = [ctrlSigStimMap repmat(0,1,nStimuliEnv-nStimuliCtrl) ];
                    warning('CtrlSigStimMap and stimRespMap don''t match. Filling missing parts with zero.');
                end
            
                % 
                if(nStimuliCtrl > nStimuliEnv)
                    ctrlSigStimMap = ctrlSigStimMap(1:nStimuliEnv);
                    warning('CtrlSigStimMap and stimRespMap don''t match. Cutting off excess of Stimuli.');
                end
                
                % translate both maps into control-to-response map
                CtrlIToResp = transpose(stimRespMap) * transpose(ctrlSigStimMap);
                
                % make sure that there are only 2 possible responses (DDM default)
                if(length(CtrlIToResp) < 2)
                    CtrlIToResp = [1 0]; % by default the current control intensity acts in favor of the upper response
                end
                
                % translate control-to-response map into scalar value
                % (positive: upper response; negative: bottom response)
                CtrlDirection = CtrlIToResp(1) - sum(CtrlIToResp(2:length(CtrlIToResp)));           
            else
               CtrlDirection = 1;
            end
            
            out = CtrlDirection * mappingFnc.getVal(currSignal.getIntensity());
        end
        
        % CONTROLAUTOMATIC... returns a response bias for the specified interval
        % param{1}: integer (refers to DDM proxy, see constants)
        % param{2}: EVC.EVCFnc (mapping function)
        function out = getControlAutomatic(this, EVCM)

            % mapping function
            if(~isa(this.params{2}, 'EVC.EVCFnc'))
               error('Second EVCFnc parameter has to be an instance of EVC.EVCFnc (intensity to DDM mapping function)');
            else
                mappingFnc = this.params{2};
            end
            
            % determine which state to operate on
            if(EVCM.useActualState)
                currTrial = EVCM.State.Actual;
            else
                currTrial = EVCM.State.Expected;
            end

            StimSalRespMap = currTrial.getRespSaliencyMap();
            controlIntensities = transpose(EVCM.getIntensities());
            
            % for drift & bias parameters, calculate in which direction control acts
            % (+1 -> upper threshold, -1 -> lower threshold)
            if(this.params{1} == this.DRIFT || this.params{1} == this.BIAS)

                automaticComponent = StimSalRespMap(:,1) - StimSalRespMap(:,2);
                controlComponent = 1 + mappingFnc.getVal(controlIntensities*EVCM.getCtrlSigStimMap()); % NOTE +1 to account for control signal intensities between 0 and 1
                controlAutomatic = controlComponent*automaticComponent ./ (sum(automaticComponent) + sum(controlComponent,2)); 
                
            else 
               error('No implementation of DDMFnc.getControlAutomatic() for DDM parameters other than drift rate or bias yet.'); 
            end
            
            out = transpose(controlAutomatic);
            
            
%             salR1 = sum(StimSalRespMap(:,1));
%             salR2 = sum(StimSalRespMap(:,2));
%             
%             rangeWidth = this.params{1}(2) - this.params{1}(1);
% 
%             % how much bias toward upper threshold?
%             R1Delta = + rangeWidth./2.*salR1;
%             
%             % how much bias toward lower threshold?
%             R2Delta = -rangeWidth./2.*salR2;
%             
%             middleOfRange = this.params{1}(1) + (rangeWidth)/2;
% 
%             %EVC.HelperFnc.disp(source, 'middleOfRange', middleOfRange, 'rangeWidth', rangeWidth, 'R1Delta', R1Delta, 'R2Delta', R2Delta); % DIAGNOSIS 
%             
%             out = middleOfRange + R1Delta + R2Delta;
           
        end
        
        
        % spotlightControl: control is allocated over both targets and
        % stimuli in a graded fashion
        % param{1}: integer (refers to DDM proxy, see constants)
        % param{2}: EVC.EVCFnc (mapping function)
        function out = spotlightControl(this, EVCM)

            % current control signal
            if(~isa(this.params{1}, 'EVC.CtrlSignal'))
               error('First EVCFnc parameter has to be an instance of EVC.CtrlSignal (control signal instnace with corresponding intensity)');
            else
                currSignal = this.params{1};
            end
            
            % mapping function
            if(~isa(this.params{3}, 'EVC.EVCFnc'))
               error('Second EVCFnc parameter has to be an instance of EVC.EVCFnc (intensity to DDM mapping function)');
            else
                mappingFnc = this.params{3};
            end
            
            % for drift & bias parameters, calculate in which direction control acts
            % (+1 -> upper threshold, -1 -> lower threshold)
            if(this.params{2} == this.DRIFT || this.params{2} == this.BIAS)
            
                % how many stimuli are available?
                ctrlSigStimMap = currSignal.CtrlSigStimMap;
                distrStimMap = 1-ctrlSigStimMap;
                
                % determine which state to operate on
                if(EVCM.useActualState)
                    currTrial = EVCM.State.Actual;
                else
                    currTrial = EVCM.State.Expected;
                end

                StimSalRespMap = currTrial.getRespSaliencyMap();

                [nStimuliEnv,~] = size(StimSalRespMap);
                [nStimuliCtrl] = length(ctrlSigStimMap);
                
                % check if maps correspond in size
                if(nStimuliCtrl < nStimuliEnv)
                    warning('CtrlSigStimMap and stimRespMap don''t match. Filling missing parts with zero.');
                    disp('CtrlSigStimMap:');
                    disp(ctrlSigStimMap);
                    disp('StimSalRespMap:');
                    disp(StimSalRespMap);
                    ctrlSigStimMap = [ctrlSigStimMap repmat(0,1,nStimuliEnv-nStimuliCtrl) ];
                end
            
                % 
                if(nStimuliCtrl > nStimuliEnv)
                    ctrlSigStimMap = ctrlSigStimMap(1:nStimuliEnv);
                    warning('CtrlSigStimMap and stimRespMap don''t match. Cutting off excess of Stimuli.');
                end
                
                % translate both maps into control-to-response map
                CtrlIToResp = transpose(StimSalRespMap) * transpose(ctrlSigStimMap);
                DistrToResp = transpose(StimSalRespMap) * transpose(distrStimMap);
                
                % make sure that there are only 2 possible responses (DDM default)
                if(length(CtrlIToResp) < 2)
                    CtrlIToResp = [1 0]; % by default the current control intensity acts in favor of the upper response boundary
                end
                
                % make sure that there are only 2 possible responses (DDM default)
                if(length(DistrToResp) < 2)
                    DistrToResp = [0 1]; % by default the current distractor acts in favor of the lower response boundary
                end
                
                % translate control-to-response map into scalar value
                % (positive: upper response; negative: bottom response)
                CtrlDirection = CtrlIToResp(1) - sum(CtrlIToResp(2:length(CtrlIToResp)));   
                DistrDirection = DistrToResp(1) - sum(DistrToResp(2:length(DistrToResp)));   
            else
               CtrlDirection = 1;
               DistrDirection = -1;
            end
            
            out = CtrlDirection * mappingFnc.getVal(currSignal.getIntensity()) + (1-mappingFnc.getVal(currSignal.getIntensity())) * DistrDirection;
            
           
        end
        
        % INTENSITY2DDM... maps control signal intensity to DDM parameter
        % param{1}: EVC.CtrlSignal
        % param{2}: integer (refers to DDM proxy, see constants)
        % param{3}: EVC.EVCFnc (mapping function)
        function out = intensityToDDM_ExpEfficacy(this, EVCM)

            % current control signal
            if(~isa(this.params{1}, 'EVC.CtrlSignal'))
               error('First EVCFnc parameter has to be an instance of EVC.CtrlSignal (control signal instnace with corresponding intensity)');
            else
                currSignal = this.params{1};
            end
            
            % mapping function
            if(~isa(this.params{3}, 'EVC.EVCFnc'))
               error('Second EVCFnc parameter has to be an instance of EVC.EVCFnc (intensity to DDM mapping function)');
            else
                mappingFnc = this.params{3};
            end
            
            efficacy = this.params{4};
            
            % for drift & bias parameters, calculate in which direction control acts
            % (+1 -> upper threshold, -1 -> lower threshold)
            if(this.params{2} == this.DRIFT || this.params{2} == this.BIAS)
            
                % how many stimuli are available?
                ctrlSigStimMap = currSignal.CtrlSigStimMap;
                stimRespMap = EVCM.currentTaskSet();
                [nStimuliEnv,~] = size(stimRespMap);
                [nStimuliCtrl] = length(ctrlSigStimMap);
                
                % check if maps correspond in size
                if(nStimuliCtrl < nStimuliEnv)
                    ctrlSigStimMap = [ctrlSigStimMap repmat(0,1,nStimuliEnv-nStimuliCtrl) ];
                    warning('CtrlSigStimMap and stimRespMap don''t match. Filling missing parts with zero.');
                end
            
                % 
                if(nStimuliCtrl > nStimuliEnv)
                    ctrlSigStimMap = ctrlSigStimMap(1:nStimuliEnv);
                    warning('CtrlSigStimMap and stimRespMap don''t match. Cutting off excess of Stimuli.');
                end
                
                % translate both maps into control-to-response map
                CtrlIToResp = transpose(stimRespMap) * transpose(ctrlSigStimMap);
                
                % make sure that there are only 2 possible responses (DDM default)
                if(length(CtrlIToResp) < 2)
                    CtrlIToResp = [1 0]; % by default the current control intensity acts in favor of the upper response
                end
                
                % translate control-to-response map into scalar value
                % (positive: upper response; negative: bottom response)
                CtrlDirection = CtrlIToResp(1) - sum(CtrlIToResp(2:length(CtrlIToResp)));           
            else
               CtrlDirection = 1;
            end
            
            if(EVCM.useActualState) % modulate efficacy of expected state
                out = CtrlDirection * mappingFnc.getVal(currSignal.getIntensity());
            else
                out = CtrlDirection * mappingFnc.getVal(currSignal.getIntensity()) * efficacy;
            end
        end
        
    end
    
end