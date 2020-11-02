% Get job number
i = getenv('SGE_TASK_ID')

% Import the data
strcat("Data/sim_dist_", i, "_stand.csv")

simdata_dist = readtable(strcat("Data/sim_dist_", i, "_stand.csv"));
simdata_over = readtable(strcat("Data/sim_over_", i, "_stand.csv"));
simdata_cor = readtable(strcat("Data/sim_cor_", i, "_stand.csv"));

%% Convert to output type
simdata_dist = table2array(simdata_dist);
simdata_over = table2array(simdata_over);
simdata_cor = table2array(simdata_cor);

%% Clear temporary variables
clear opts

%% Run npBNMF
[ewa_dist,eh_dist] = NPBayesNMF(simdata_dist, 50);
[ewa_over,eh_over] = NPBayesNMF(simdata_over, 50);
[ewa_cor,eh_cor] = NPBayesNMF(simdata_cor, 50);

% labels = ["pcb28", "pcb66", "pcb74", "pcb99", "pcb105", "pcb118", "pcb138_158" "pcb146", "pcb153", "pcb156", "pcb167", ...
%    "pcb170", "pcb178", "pcb183", "pcb187", "pcb180", "pcb189", "pcb194", "pcb196_203" "pcb199", "pcb206", "pcb209", "BDE17", ...
%    "BDE28", "BDE47", "BDE66", "BDE85", "BDE99", "BDE100", "BDE153", "BDE154", "BDE183", "BDE209", "MECPP", "MEHHP", "MEOHP", ...
%    "MCPP", "MIBP", "MBP", "MBZP", "MEP", "MEHP", "dcp_24", "dcp_25", "b_pb", "bp_3", "m_pb", "p_pb", "tcs", "bpa"];

%% Save output
save(strcat("ewa_dist_", i, "_stand.mat"), 'ewa_dist');
save(strcat("eh_dist_", i, "_stand.mat"), 'eh_dist');

save(strcat("ewa_over_", i, "_stand.mat"), 'ewa_over');
save(strcat("eh_over_", i, "_stand.mat"), 'eh_over');

save(strcat("ewa_cor_", i, "_stand.mat"), 'ewa_cor');
save(strcat("eh_cor_", i, "_stand.mat"), 'eh_cor');

