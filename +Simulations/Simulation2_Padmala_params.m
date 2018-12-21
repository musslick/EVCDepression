% set parameters

% T0                                     = 0.2500;
% c                                       =  0.5000;
% thresh                               = 0.3634;
% wTarget                            = 0.2754;
% TargetDistractorRatio       = 1.000;
% costParam1                       = 0.1561;
% reward                              = 10.000;
% noRewardToRewardRatio  = 0.4;
% learningRate                     = 0.68;

T0                                     = 0.2500;
c                                       =  0.5000;
thresh                               = 0.3634;
wTarget                            = 0.10; %0.2754;
TargetDistractorRatio       = 0.3;
costParam1                       = 0.18;
reward                              = 100.000;
noRewardToRewardRatio  = 0.4;
learningRate                     = 0;% 0.58;

wDistractor = wTarget / TargetDistractorRatio;
noReward = noRewardToRewardRatio * reward;

% set parameters
EVCSim.defaultDDMParams.T0 = T0; 
EVCSim.defaultDDMParams.c = c;
EVCSim.defaultDDMParams.thresh = thresh;

EVCSim.trials(1).stimRespMap   = [wDistractor 0;                                          % responding to stimulus 1 tends to produce second response by 100%
                                  wTarget 0];                                         % responding to stimulus 2 tends to produce second outcome by 100%

EVCSim.trials(2).stimRespMap   = [0 wDistractor;                                          % responding to stimulus 1 tends to produce second response by 100%
                                  wTarget 0];                                         % responding to stimulus 2 tends to produce second outcome by 100%

EVCSim.trials(3).stimRespMap   = [wDistractor/2 wDistractor/2;                                          % responding to stimulus 1 tends to produce second response by 100%
                                  wTarget 0];                                         % responding to stimulus 2 tends to produce second outcome by 100%

EVCSim.defaultCostFnc.params{1} = costParam1; 

EVCSim.reward = reward;
EVCSim.noReward = noReward;

EVCSim.learningFnc(1).params{1} = learningRate;

