%% Get job number
i = getenv('SGE_TASK_ID')

% Import the data
simdata1 = readtable(strcat("/ifs/scratch/msph/ehs/eag2186/Data/dgp_csv/sim_dim", num2str(i), ".csv"));

%% Convert to output type
simdata1 = table2array(simdata1);

%% Clear temporary variables
clear opts

[EWA, EH] = NPBayesNMF(simdata1);

save(strcat("/Users/lizzy/BN2MF/MATLAB/loop_dgp_rep1/sim_dgp_rep1_", num2str(i), ".mat"), 'simdata1');
save(strcat("/Users/lizzy/BN2MF/MATLAB/loop_dgp_rep1/rep1_ewa_dist_", num2str(i), ".mat"), 'EWA');
save(strcat("/Users/lizzy/BN2MF/MATLAB/loop_dgp_rep1/rep1_eh_dist_",  num2str(i), ".mat"), 'EH');

end


