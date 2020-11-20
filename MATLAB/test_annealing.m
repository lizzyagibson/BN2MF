%% Setup the Import Options and import the data
% opts = delimitedTextImportOptions("NumVariables", 50);
% % Specify range and delimiter
% opts.DataLines = [2, Inf];
% opts.Delimiter = ",";
% % Specify column names and types
% opts.VariableTypes = ["double", "double", "double", "double", "double", "double", ...
%     "double", "double", "double", "double", "double", "double", "double", "double", "double", ...
%     "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", ...
%     "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", ...
%     "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", ...
%     "double", "double", "double", "double", "double"];
% % Specify file level properties
% opts.ExtraColumnsRule = "ignore";
% opts.EmptyLineRule = "read";
% % Import the data
% simdata1 = readtable(strcat("/Users/lizzy/BN2MF/Sims/dgp_rep1/sim_dgp_rep1_", num2str(1), ".csv"), opts);
% %% Convert to output type
% simdata1 = table2array(simdata1);
% %% Clear temporary variables
% clear opts

%[EWA, EH, varH, alphaH, betaH] = NPBayesNMF(simdata1);
% x = gaminv(p,a,b) returns the icdf of the gamma distribution with 
% shape parameter a and the scale parameter b, evaluated at the values in p.

% shapeH = 1./betaH;
% gaminv(0.025, alphaH, shapeH);
% gaminv(0.975, alphaH, shapeH);
sum(sum( varH ./ EH))
    
[EWAa, EHa, varHa, alphaHa, betaHa] = NPBayesNMF_annealing(simdata1);
sum(sum( varHa ./ EHa))

shapeHa = 1./betaHa;
lower = gaminv(0.025, alphaHa, shapeHa);
upper = gaminv(0.975, alphaHa, shapeHa);

% Hard to compare with truth
% Need to normalize truth and EH and then compare
% ship to R

% save("/Users/lizzy/BN2MF/MATLAB/test_loadings.mat", 'EHa');
% save("/Users/lizzy/BN2MF/MATLAB/test_upper.mat", 'upper');
% save("/Users/lizzy/BN2MF/MATLAB/test_lower.mat", 'lower');


