function printSimulationParameters(EVCSim)

disp('S-R map of first trial:');
EVCSim.EVCM.State.TaskEnv.trialTypes{1}.stimRespMap

disp(['a_f (Parameter of control implementation cost function) = ' num2str(EVCSim.defaultCostFnc.params{1})]);
disp(['b_f (Parameter of control implementation cost function) = ' num2str(EVCSim.defaultCostFnc.params{2})]);
disp(['a_g (Parameter of control reconfiguration cost function) = ' num2str(EVCSim.reconfCostFnc.params{1})]);
disp(['b_g (Parameter of control reconfiguration cost function) = ' num2str(EVCSim.reconfCostFnc.params{2})]);
disp(['T0 (Non-decision time in seconds) = ' num2str(EVCSim.defaultDDMParams.T0)]);
disp(['c (Standard deviation of DDM noise distribution) = ' num2str(EVCSim.defaultDDMParams.c)]);
disp(['z (DDM threshold) = ' num2str(EVCSim.defaultDDMParams.thresh)]);
disp(['x_0 (Starting point) = ' num2str(EVCSim.defaultDDMParams.bias - 0.5)]);

disp(['R_1 (Reward obtained for responding to the (upper) threshold specified by the target feature) = ' num2str(EVCSim.EVCM.State.TaskEnv.trialTypes{1}.outcomeValues(1))]);
disp(['R_2 (Reward obtained for responding to the (lower) threshold specified by the distractor feature) = ' num2str(EVCSim.EVCM.State.TaskEnv.trialTypes{1}.outcomeValues(2))]);

disp(['\alpha_e (Learning rate based on external feedback) = ' num2str(EVCSim.learningFnc(1).params{1})]);



end