%% Get job number
i = getenv('SGE_TASK_ID')

% Import the data
simdata1 = readtable(strcat("/ifs/scratch/msph/ehs/eag2186/Data/dgp_csv/sim_dgp_rep1_", num2str(i), ".csv"));

%% Convert to output type
simdata1 = table2array(simdata1);

[EWA, EH, varH, alphaH, betaH, varWA, finalscore, final_iter] = BN2MF(simdata1);

save(strcat("/ifs/scratch/msph/ehs/eag2186/npbnmf/bnmf_dgp_out/dgp_ewa_",    num2str(i), ".mat"), 'EWA');
save(strcat("/ifs/scratch/msph/ehs/eag2186/npbnmf/bnmf_dgp_out/dgp_eh_",     num2str(i), ".mat"), 'EH');
save(strcat("/ifs/scratch/msph/ehs/eag2186/npbnmf/bnmf_dgp_out/dgp_varwa_",  num2str(i), ".mat"), 'varWA');
save(strcat("/ifs/scratch/msph/ehs/eag2186/npbnmf/bnmf_dgp_out/dgp_alphah_", num2str(i), ".mat"), 'alphaH');
save(strcat("/ifs/scratch/msph/ehs/eag2186/npbnmf/bnmf_dgp_out/dgp_betah_",  num2str(i), ".mat"), 'betaH');



