function DepressionSim_COGED_RewardSensitivityExp()

%% RUN SIMULATION

clear all;
close all;
clc;
import Simulations.*;
import EVC.*;

plotOffline = 0;

% simulation settings
nSubj = 1;     % 54
nTrials = 210;
probedRewardSensitivityExp = [0.5 1 1.5];

baselineAutomaticity = 1;
taskAAutomaticityRange = [linspace(0.2, 0.4, 4) baselineAutomaticity ];%[ 0.1 0.11 0.12 0.15 0.9 1]; %linspace(0.01, 1, 6); 0.2:0.2:1
taskAAutomaticityRange(end) = [];
taskAAutomaticityRange = fliplr(taskAAutomaticityRange); % easy to difficult

% plot settings
plotSEM = 0;
showLegend = 0;
fixYLimit = 1;

% limits    
ylimitSV = [0.6 1.6];


probedRewardSensitivityExp = [0.5 1 2];
numCostParams = length(probedRewardSensitivityExp);

% regressor = repmat(probedRewardSensitivityExp, nSubj, 1);
subjectiveValueLog = nan(length(probedRewardSensitivityExp), length(taskAAutomaticityRange));
accuracyLog = nan(length(probedRewardSensitivityExp), length(taskAAutomaticityRange));

for parameterCondition = 1:numCostParams
    
    
    disp(['*************************TESTED REWARD SENSITIVITY EXP: ' num2str(probedRewardSensitivityExp(parameterCondition))]);

    subjectiveValue = nan(1, length(taskAAutomaticityRange));
    controlSignalLog = zeros(1, length(taskAAutomaticityRange));
    accuracy = zeros(1, length(taskAAutomaticityRange));

    for taskAAutomaticityIdx = 1:length(taskAAutomaticityRange);

            EVCSim = DDM_WestbrookBraver2015();
            Simulation6_WestbrookBraver2015_params;
            EVCSim.nSubj = nSubj;
            EVCSim.nTrials = nTrials;
            EVCSim.taskAAutomaticity = taskAAutomaticityRange(taskAAutomaticityIdx);
            EVCSim.taskBAutomaticity = baselineAutomaticity;
            EVCSim.rewardMinimum = 1;
            
            temp.rewardFnc.params{1} = 1;
            temp.rewardFnc.params{2} = 0;
            temp.rewardFnc.params{3} = probedRewardSensitivityExp(parameterCondition);
            temp.rewardFnc.params{4} = 0;
            temp.rewardFnc.type = EVCFnc.REWRATE_EXT;
            EVCSim.rewardFnc = EVCFnc(temp.rewardFnc.type, temp.rewardFnc.params);
            
            EVCSim.run();

            % find equilibrium trial
            t_eq = nan;
            SV = nan;
            for t = 1:length(EVCSim.subjData.Log.ExpectedState)
                if(round(EVCSim.subjData.Log.ExpectedState(t).descr == 'taskB'))
                    SV = EVCSim.subjData.Log.ExpectedState(t).outcomeValues(2);
                    t_eq = t;
                    break;
                end
            end

            subjectiveValue(taskAAutomaticityIdx) = SV;
            controlSignalLog(taskAAutomaticityIdx) = EVCSim.subjData.Log.CtrlIntensities(t_eq, 2);
            accuracy(taskAAutomaticityIdx) = 1-EVCSim.subjData.Log.ERs(t_eq-1, 1);

            disp(['task A automaticity ' num2str(taskAAutomaticityIdx) '/' num2str(length(taskAAutomaticityRange)) '.']);

    end
    
    subjectiveValueLog(parameterCondition, :) = subjectiveValue;
    accuracyLog(parameterCondition, :) = accuracy;
    
end

save(['logfiles/DepressionSim_COGED_nSubj' num2str(nSubj) '_RewardSensitivityExp_' num2str(min(probedRewardSensitivityExp))  ...
                                                                                               '_' num2str(max(probedRewardSensitivityExp)) ...
                                                                                               '_' num2str(probedRewardSensitivityExp(2)-probedRewardSensitivityExp(1)) ...
                                                                                               '.mat']);
                                                                                           
%% perform regression and plot
% load('logfiles/DepressionSim_Padmala_RewardSensitivityExp_0.1_0.2_0.01.mat');

ylimitSV = [0.6 1];

EVCPlotSettings;
DepressionSim_ylimits;

close all;
fig1 = figure(1);
set(fig1, 'Position', [100 100 250 230]);
if(plotOffline)
    set(fig1, 'visible','off');
end

colorGradient = getAlphaGradient([0.7 0.7 0.7], [0 0 0], size(subjectiveValueLog,1));
legendText = {};

for parameterCondition = 1:size(subjectiveValueLog,1)
    
    subjectiveValue = subjectiveValueLog(parameterCondition,:);
    plot([1:length(subjectiveValue)]+1,subjectiveValue/EVCSim.taskAReward, 'LineWidth', lineWidth.line1, 'Color', colorGradient(parameterCondition,:)); hold on;
    legendText{parameterCondition} = ['Reward Sensitivity v = ' num2str(probedRewardSensitivityExp(parameterCondition))];
    
end
hold off;

xlabel('Task Difficulty', 'fontSize', fontSize.xlabel);
ylabel({'Subjective Value'}, 'fontSize', fontSize.ylabel);
leg = legend(legendText, 'Location', 'southwest');
set(leg, 'FontSize', fontSize.xlabel-3);
set(gca, 'fontSize', fontSize.xlabel);
if(fixYLimit)
    ylim(ylimitSV);
end

saveas(fig1,['figures/DepressionSim_COGED_RewardSensitivityExp_Full_nSubj' num2str(nSubj) '_RewardSensitivityExp_' num2str(min(probedRewardSensitivityExp))  ...
                                                                                               '_' num2str(max(probedRewardSensitivityExp)) ...
                                                                                               '_' num2str(probedRewardSensitivityExp(2)-probedRewardSensitivityExp(1)) '.fig'],'fig');
                                                                                           
%% PRINT PARAMETERS

printSimulationParameters(EVCSim);

% openfig('figures/DepressionSim_Padmala_RewardSensitivityExp_Subset_nSubj100_RewardSensitivityExp_0.1_0.2_0.05.fig','new','visible')
end