function DepressionSim_COGED_ControlCostExp()

%% RUN SIMULATION

clear all;
close all;
clc;
import Simulations.*;

plotOffline = 0;

% simulation settings
nSubj = 1;     % 54
nTrials = 210;

baselineAutomaticity = 1;
taskAAutomaticityRange = [linspace(0.2, 0.4, 4) baselineAutomaticity ];%[ 0.1 0.11 0.12 0.15 0.9 1]; %linspace(0.01, 1, 6); 0.2:0.2:1
taskAAutomaticityRange(end) = [];
taskAAutomaticityRange = fliplr(taskAAutomaticityRange); % easy to difficult

% plot settings
plotSEM = 0;
showLegend = 0;
fixYLimit = 1;

% limits    
ylimitSV = [0 1.75];

% RUN REGRESSION

probedControlCostExp = [1 1.5 2];
numCostParams = length(probedControlCostExp);

% regressor = repmat(probedControlCostExp, nSubj, 1);
subjectiveValueLog = nan(length(probedControlCostExp), length(taskAAutomaticityRange));

for parameterCondition = 1:numCostParams
    
    
    disp(['*************************TESTED CONTROL COST EXP: ' num2str(probedControlCostExp(parameterCondition))]);

    subjectiveValue = nan(1, length(taskAAutomaticityRange));
    controlSignalLog = zeros(1, length(taskAAutomaticityRange));

    for taskAAutomaticityIdx = 1:length(taskAAutomaticityRange);

            EVCSim = DDM_WestbrookBraver2015();
            Simulation6_WestbrookBraver2015_params;
            EVCSim.nSubj = nSubj;
            EVCSim.nTrials = nTrials;
            EVCSim.taskAAutomaticity = taskAAutomaticityRange(taskAAutomaticityIdx);
            EVCSim.taskBAutomaticity = baselineAutomaticity;
            EVCSim.rewardMinimum = 1;
            
            EVCSim.defaultCostFnc.params{1} = probedControlCostExp(parameterCondition);
            
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

            disp(['task A automaticity ' num2str(taskAAutomaticityIdx) '/' num2str(length(taskAAutomaticityRange)) '.']);

    end
    
    subjectiveValueLog(parameterCondition, :) = subjectiveValue;
    
end

save(['logfiles/DepressionSim_COGED_nSubj' num2str(nSubj) '_ControlCostExp_' num2str(min(probedControlCostExp))  ...
                                                                                               '_' num2str(max(probedControlCostExp)) ...
                                                                                               '_' num2str(probedControlCostExp(2)-probedControlCostExp(1)) ...
                                                                                               '.mat']);
                                                                                           
%% perform regression and plot
% load('logfiles/DepressionSim_Padmala_ControlCostExp_0.1_0.2_0.01.mat');

EVCPlotSettings;
DepressionSim_ylimits;

ylimitSV = [0.6 1];

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
    legendText{parameterCondition} = ['Control Cost c = ' num2str(probedControlCostExp(parameterCondition))];
    
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

saveas(fig1,['figures/DepressionSim_COGED_ControlCostExp_Full_nSubj' num2str(nSubj) '_ControlCostExp_' num2str(min(probedControlCostExp))  ...
                                                                                               '_' num2str(max(probedControlCostExp)) ...
                                                                                               '_' num2str(probedControlCostExp(2)-probedControlCostExp(1)) '.fig'],'fig')
                                                                                                                                                                            
%% PRINT PARAMETERS

printSimulationParameters(EVCSim);

% openfig('figures/DepressionSim_COGED_ControlCostExp_Full_nSubj1_ControlCostExp_1_2_0.5.fig','new','visible')
end