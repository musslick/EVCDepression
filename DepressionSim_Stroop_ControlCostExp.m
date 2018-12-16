function DepressionSim_Stroop_ControlCostExp()

%% RUN SIMULATION

clear all;
close all;
clc;
import Simulations.*;

plotOffline = 0;

% simulation settings
nSubj = 1;     
nTrials = 200;
Simulation2_Padmala_params;

% plot settings
plotSEM = 0;
showLegend = 0;
fixYLimit = 1;

% RUN REGRESSION

Simulation2_Padmala_params;
probedControlCostExp = [0.14:0.01:0.22]; % [0.5:0.1:1.5]
numCostParams = length(probedControlCostExp);

regressor = repmat(probedControlCostExp, nSubj, 1);
ControlSignalIntensity = nan(nSubj, numCostParams);
RT_Interference = nan(nSubj, numCostParams);
ER_Interference = nan(nSubj, numCostParams);

for parameterCondition = 1:numCostParams
    
    
    disp(['*************************TESTED CONTROL COST EXP: ' num2str(probedControlCostExp(parameterCondition))]);

    EVCSim = DDM_PWStroop();
    EVCSim.printResults = 1;
    EVCSim.nSubj = nSubj;
    EVCSim.nTrials = nTrials;

    % set parameters and learning rate
    Simulation2_Padmala_params;
    EVCSim.defaultCostFnc.params{1} = probedControlCostExp(parameterCondition);

    EVCSim.run();
    
    ControlSignalIntensity(:,parameterCondition) = EVCSim.results.targetIntensity.norew;
    
    RT_Interference(:,parameterCondition) = EVCSim.results.RT.Inc - EVCSim.results.RT.Con;

    EVCSim.results.Accuracy.Interference = EVCSim.results.Accuracy.Inc - EVCSim.results.Accuracy.Con;

    ER_Interference(:,parameterCondition) = - EVCSim.results.Accuracy.Interference * 100;

end

save(['logfiles/DepressionSim_Stroop_nSubj' num2str(nSubj) '_ControlCostExp_' num2str(min(probedControlCostExp))  ...
                                                                                               '_' num2str(max(probedControlCostExp)) ...
                                                                                               '_' num2str(probedControlCostExp(2)-probedControlCostExp(1)) ...
                                                                                               '.mat']);
% perform regression and plot

EVCPlotSettings;
DepressionSim_ylimits;

plotRTER = 2; % set to 1 for plotting RT, set to 2 for plotting error rate

fig1 = figure(1);
set(fig1, 'Position', [100 100 250 300]);
if(plotOffline)
    set(fig1, 'visible','off');
end

% Control Signal Intensity

subplot(2,1,1);
plot(probedControlCostExp, mean(ControlSignalIntensity,1), style.line1, 'LineWidth', lineWidth.line1, 'Color', color.line1);  hold off;
set(gca, 'FontSize', fontSize.TickLabel, 'FontWeight', fontWeight.TickLabel);
ylabel({'Control Signal', 'Intensity'}, 'FontSize', fontSize.ylabel, 'FontWeight', fontWeight.ylabel);
xlabel('Control Cost c', 'FontSize', fontSize.xlabel, 'FontWeight', fontWeight.xlabel);
if(fixYLimit)
     ylim([ylimitIntensity]);
end

% Reaction Times

subplot(2,1,2);
if(plotRTER == 1)
   plot(probedControlCostExp, mean(RT_Interference,1), style.line1, 'LineWidth', lineWidth.line1, 'Color', color.line1);  hold off; 
   ylabel({'Incongruency Cost', 'in Reaction Time (s)'}, 'FontSize', fontSize.ylabel, 'FontWeight', fontWeight.ylabel);
   if(fixYLimit)
       ylim([ylimitRT]);
   end
elseif(plotRTER == 2)
   plot(probedControlCostExp, mean(ER_Interference,1), style.line1, 'LineWidth', lineWidth.line1, 'Color', color.line1);  hold off; 
   ylabel({'Incongruency Cost', 'in Error Rate (%)'}, 'FontSize', fontSize.ylabel, 'FontWeight', fontWeight.ylabel); 
   if(fixYLimit)
       ylim([ylimitER]);
   end
end
set(gca, 'FontSize', fontSize.TickLabel, 'FontWeight', fontWeight.TickLabel);
xlabel('Control Cost c', 'FontSize', fontSize.xlabel, 'FontWeight', fontWeight.xlabel);

saveas(fig1,['figures/DepressionSim_Stroop_ControlCostExp_Full_nSubj' num2str(nSubj) '_ControlCostExp_' num2str(min(probedControlCostExp))  ...
                                                                                               '_' num2str(max(probedControlCostExp)) ...
                                                                                               '_' num2str(probedControlCostExp(2)-probedControlCostExp(1)) '.fig'],'fig')
%% PRINT PARAMETERS

printSimulationParameters(EVCSim);

% openfig('figures/DepressionSim_Padmala_ControlCostExp_Subset_nSubj100_ControlCostExp_0.1_0.2_0.05.fig','new','visible')
end