function DepressionSim_Stroop_ControlEfficacyExp()

%% RUN SIMULATION

clear all;
close all;
clc;
import Simulations.*;
import EVC.DDM.*;

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
probedControlEfficacyExp = [0.5:0.1:1.5];  % [0.5:0.5:5];
numCostParams = length(probedControlEfficacyExp);

regressor = repmat(probedControlEfficacyExp, nSubj, 1);
ControlSignalIntensity = nan(nSubj, numCostParams);
RT_Interference = nan(nSubj, numCostParams);
ER_Interference = nan(nSubj, numCostParams);

for parameterCondition = 1:numCostParams
    
    
    disp(['*************************TESTED CONTROL EFFICACY EXP: ' num2str(probedControlEfficacyExp(parameterCondition))]);

    EVCSim = DDM_PWStroop();
    EVCSim.printResults = 1;
    EVCSim.nSubj = nSubj;
    EVCSim.nTrials = nTrials;

    % set parameters and learning rate
    Simulation2_Padmala_params;
    temp.targetControlFnc.type = DDMFnc.INTENSITY2DDM_EXPEFFICACY;                      % define function type for mapping control onto DDM parameter (see EVC/DDM/DDMFnc.m)
    temp.targetControlFnc.params{1} = EVCSim.ctrlSignals(1);                          % assign control signal for picture stimulus
    temp.targetControlFnc.params{2} = EVCSim.defaultControlProxy;               % determine control signal proxy (e.g. drift rate, see Simulations/DDMSim.m)
    temp.targetControlFnc.params{3} = EVCSim.defaultControlMappingFnc;     % assign function that maps control signal onto drift rate (see Simulations/DDMSim.m)
    temp.targetControlFnc.params{4} = probedControlEfficacyExp(parameterCondition);
    EVCSim.DDMProcesses(2).input = DDMFnc(temp.targetControlFnc.type, ...
                                    temp.targetControlFnc.params);
    
    temp.flanker1ControlFnc.type = DDMFnc.INTENSITY2DDM_EXPEFFICACY;                  % define function type for mapping control onto DDM parameter (see EVC/DDM/DDMFnc.m)
    temp.flanker1ControlFnc.params{1} = EVCSim.ctrlSignals(2);                      % assign control signal for word stimulus
    temp.flanker1ControlFnc.params{2} = EVCSim.defaultControlProxy;           % determine control signal proxy (e.g. drift rate, see Simulations/DDMSim.m)
    temp.flanker1ControlFnc.params{3} = EVCSim.defaultControlMappingFnc; % assign function that maps control signal onto drift rate (see Simulations/DDMSim.m)
    temp.flanker1ControlFnc.params{4} = probedControlEfficacyExp(parameterCondition);
    EVCSim.DDMProcesses(3).input = DDMFnc(temp.flanker1ControlFnc.type, ...
                                    temp.flanker1ControlFnc.params);

    EVCSim.run();
    
    ControlSignalIntensity(:,parameterCondition) = EVCSim.results.targetIntensity.norew;
    
    RT_Interference(:,parameterCondition) = EVCSim.results.RT.Inc - EVCSim.results.RT.Con;

    EVCSim.results.Accuracy.Interference = EVCSim.results.Accuracy.Inc - EVCSim.results.Accuracy.Con;

    ER_Interference(:,parameterCondition) = - EVCSim.results.Accuracy.Interference * 100;

end

save(['logfiles/DepressionSim_Stroop_nSubj' num2str(nSubj) '_ControlEfficacyExp_' num2str(min(probedControlEfficacyExp))  ...
                                                                                               '_' num2str(max(probedControlEfficacyExp)) ...
                                                                                               '_' num2str(probedControlEfficacyExp(2)-probedControlEfficacyExp(1)) ...
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
plot(probedControlEfficacyExp, mean(ControlSignalIntensity,1), style.line1, 'LineWidth', lineWidth.line1, 'Color', color.line1);  hold off;
set(gca, 'FontSize', fontSize.TickLabel, 'FontWeight', fontWeight.TickLabel);
ylabel({'Control Signal', 'Intensity'}, 'FontSize', fontSize.ylabel, 'FontWeight', fontWeight.ylabel);
xlabel('Control Efficacy \epsilon', 'FontSize', fontSize.xlabel, 'FontWeight', fontWeight.xlabel);
if(fixYLimit)
     ylim([ylimitIntensity]);
end

% Reaction Times

subplot(2,1,2);
if(plotRTER == 1)
   plot(probedControlEfficacyExp, mean(RT_Interference,1), style.line1, 'LineWidth', lineWidth.line1, 'Color', color.line1);  hold off; 
   ylabel({'Incongruency Cost', 'in Reaction Time (s)'}, 'FontSize', fontSize.ylabel, 'FontWeight', fontWeight.ylabel);
   if(fixYLimit)
       ylim([ylimitRT]);
   end
elseif(plotRTER == 2)
   plot(probedControlEfficacyExp, mean(ER_Interference,1), style.line1, 'LineWidth', lineWidth.line1, 'Color', color.line1);  hold off; 
   ylabel({'Incongruency Cost', 'in Error Rate (%)'}, 'FontSize', fontSize.ylabel, 'FontWeight', fontWeight.ylabel); 
   if(fixYLimit)
       ylim([ylimitER]);
   end
end
set(gca, 'FontSize', fontSize.TickLabel, 'FontWeight', fontWeight.TickLabel);
xlabel('Control Efficacy \epsilon', 'FontSize', fontSize.xlabel, 'FontWeight', fontWeight.xlabel);

saveas(fig1,['figures/DepressionSim_Stroop_ControlEfficacyExp_Full_nSubj' num2str(nSubj) '_ControlEfficacyExp_' num2str(min(probedControlEfficacyExp))  ...
                                                                                               '_' num2str(max(probedControlEfficacyExp)) ...
                                                                                               '_' num2str(probedControlEfficacyExp(2)-probedControlEfficacyExp(1)) '.fig'],'fig')
%% PRINT PARAMETERS

printSimulationParameters(EVCSim);

% openfig('figures/DepressionSim_Padmala_ControlEfficacyExp_Subset_nSubj100_ControlEfficacyExp_0.1_0.2_0.05.fig','new','visible')
end