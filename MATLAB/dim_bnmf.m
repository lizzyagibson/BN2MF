%% Get job number
i = getenv('SGE_TASK_ID')

%% Import the data
simdata1 = readtable(strcat("/ifs/scratch/msph/ehs/eag2186/Data/sim_dim/sim_dim", num2str(i), ".csv"));

%% Convert to output type
simdata1 = table2array(simdata1);

%% Run BN2MF
[EWA, EH, varH, alphaH, betaH, varWA] = NPBayesNMF(simdata1);

%% Save output
save(strcat("/ifs/scratch/msph/ehs/eag2186/npbnmf/bnmf_dim_out/ewa_dim_", num2str(i), ".mat"), 'EWA');
save(strcat("/ifs/scratch/msph/ehs/eag2186/npbnmf/bnmf_dim_out/eh_dim_",  num2str(i), ".mat"), 'EH');


