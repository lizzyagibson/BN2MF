% Import the data
simdata1 = readtable(strcat("/Users/lizzy/BN2MF/Sims/dgp_rep1/sim_dgp_rep1_", num2str(1), ".csv"));
% Convert to output type
simdata1 = table2array(simdata1);

tic()
[EWAreg, EHreg, varHreg, alphaHreg, betaHreg, varWAreg, scorereg] = NPBayesNMF(simdata1);
toc()

size(EWAreg)
scorereg
mean(mean(varHreg)) % 6.2878e-04
mean(mean(EHreg)) % 1.0550
