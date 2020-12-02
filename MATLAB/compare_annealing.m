% Import the data
simdata1 = readtable(strcat("/Users/lizzy/BN2MF/Sims/dgp_rep1/sim_dgp_rep1_", num2str(1), ".csv"));
% Convert to output type
simdata1 = table2array(simdata1);

[EWAreg, EHreg, varHreg, alphaHreg, betaHreg] = NPBayesNMF(simdata1);

size(EWAreg)
sum(sum(varHreg))
sum(sum(EHreg))
sum(sum(varHreg ./ EHreg))