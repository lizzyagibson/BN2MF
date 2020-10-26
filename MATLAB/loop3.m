function [EWA, EH, varWA, varH, alphaH, betaH] = loop3()
%% Make this whole thing loop

for i = 101:150

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 50);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
simdata1 = readtable(strcat("/Users/lizzy/BN2MF/Sims/dgp_rep1/sim_dgp_rep1_", num2str(i), ".csv"), opts);

%% Convert to output type
simdata1 = table2array(simdata1);

%% Clear temporary variables
clear opts

[EWA, EH, varWA, varH, alphaH, betaH] = NPBayesNMF(simdata1);

save(strcat("/Users/lizzy/BN2MF/MATLAB/dgp_rep1/rep1_ewa_dist_", num2str(i), ".mat"), 'EWA');
save(strcat("/Users/lizzy/BN2MF/MATLAB/dgp_rep1/rep1_eh_dist_",  num2str(i), ".mat"), 'EH');

end

end