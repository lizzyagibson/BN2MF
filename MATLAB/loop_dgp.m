% Get job number
i = getenv('SGE_TASK_ID')

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
simdata1 = readtable(strcat("/ifs/scratch/msph/ehs/eag2186/Data/dgp_rep1/sim_dgp_rep1_", num2str(i), ".csv"), opts);
% /ifs/scratch/msph/ehs/eag2186/Data/

%% Convert to output type
simdata1 = table2array(simdata1);

%% Clear temporary variables
clear opts

[EWA, EH] = NPBayesNMF(simdata1);

save(strcat("/ifs/scratch/msph/ehs/eag2186/Data/npbnmf/dgp_rep1_100/rep1_ewa_dist_", num2str(i), ".mat"), 'EWA');
save(strcat("/ifs/scratch/msph/ehs/eag2186/Data/npbnmf/dgp_rep1_100/rep1_eh_dist_",  num2str(i), ".mat"), 'EH');

