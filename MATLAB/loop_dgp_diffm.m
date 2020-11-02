% Make this whole thing loop

i = getenv('SGE_TASK_ID')

%% Setup the Import Options and import the data
% opts = delimitedTextImportOptions("NumVariables", 50);
% l
% % Specify range and delimiter
% opts.DataLines = [2, Inf];
% opts.Delimiter = ",";
% 
% % Specify column names and types
% opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
% 
% % Specify file level properties
% opts.ExtraColumnsRule = "ignore";
% opts.EmptyLineRule = "read";

% Import the data
simdata1 = readtable(strcat("/ifs/scratch/msph/ehs/eag2186/Data/dgp_local/sim_dgp_local_", num2str(i), ".csv"), opts);

%% Convert to output type
simdata1 = table2array(simdata1);

%% Clear temporary variables
clear opts

[EWA, EH] = NPBayesNMF(simdata1);

save(strcat("/ifs/scratch/msph/ehs/eag2186/Data/npbnmf/dgp_100/ewa_dist_", num2str(i), ".mat"), 'EWA');
save(strcat("/ifs/scratch/msph/ehs/eag2186/Data/npbnmf/dgp_100/eh_dist_",  num2str(i), ".mat"), 'EH');
