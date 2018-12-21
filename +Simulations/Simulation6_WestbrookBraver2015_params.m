T0 = 0.2;
c = 0.8;
thresh = 1.0;
costParam1 = 1.5;
reconfCostParam1 = 0;

% set trial number
EVCSim.nTrials = 100; 

% set parameters
EVCSim.defaultDDMParams.T0 = T0; 
EVCSim.defaultDDMParams.c = c;
EVCSim.defaultDDMParams.thresh = thresh;
EVCSim.reconfCostFnc.params{1} = reconfCostParam1;
EVCSim.defaultCostFnc.params{1} = costParam1; 

EVCSim.taskAReward = 200;
EVCSim.taskAAutomaticity = 0.1;
EVCSim.taskBAutomaticity = 1.2;
EVCSim.rewardMinimum = 1;

