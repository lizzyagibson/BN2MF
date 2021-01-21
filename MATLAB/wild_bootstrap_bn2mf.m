%% Wild Bootstrap BN2MF
% Run bn2mf, bootstrap residuals
% Then combine results to create empirical distribution for scores
% Want to compare bootstrap CI to VCI

%% Get job number
%j = getenv('SGE_TASK_ID')

%% Choose 1 example of correlated simulations
sim       = table2array(readtable("/Users/lizzy/BN2MF/Results/Main/Corr Ex/corr_sim.csv"));
pre_noise = table2array(readtable("/Users/lizzy/BN2MF/Results/Main/Corr Ex/corr_chem.csv"));
patterns  = table2array(readtable("/Users/lizzy/BN2MF/Results/Main/Corr Ex/corr_patterns.csv"));
scores    = table2array(readtable("/Users/lizzy/BN2MF/Results/Main/Corr Ex/corr_scores.csv"));

% sim       = table2array(readtable("/ifs/scratch/msph/ehs/eag2186/Data/Corr Ex/corr_sim.csv"));
% pre_noise = table2array(readtable("/ifs/scratch/msph/ehs/eag2186/Data/Corr Ex/corr_chem.csv"));
% patterns  = table2array(readtable("/ifs/scratch/msph/ehs/eag2186/Data/Corr Ex/corr_patterns.csv"));
% scores    = table2array(readtable("/ifs/scratch/msph/ehs/eag2186/Data/Corr Ex/corr_scores.csv"));

%% Normalize truth
patterns_denom      = sum(patterns, 2);
patterns_scaled     = patterns ./ patterns_denom;
patterns_denom_diag = diag(patterns_denom);
scores_scaled       = scores * patterns_denom_diag;

%% Run bn2mf
[EWA0, EH0, varH0, alphaH0, betaH0, alphaW0, betaW0, ...
   alphaA0, betaA0, varWA0, finalscore0, final_iter0, init] = BN2MF(sim);

pred       = EWA0 * EH0;
resid      = sim - pred;
true_resid = chem - pred;

%% Take a bootstrapped sample
[n, p] = size(sim);
sam = randsample(1:n, n, true);
sim_sample = sim(sam, :);

%% Normalize EH to L1 norm across chemicals
H_denom   = sum(EH0, 2);
EH_scaled = EH0 ./ H_denom;

%% Scale EWA by corresponding normalization constant
H_denom_diag = diag(H_denom);
EWA_scaled   = EWA0 * H_denom_diag;

%% Rearrange solution to match truth
[~,Pi] = factor_correspondence(patterns_scaled',EH_scaled');
EH_final = (EH_scaled' * Pi)';
EWA_final = EWA_scaled * Pi;

%save(strcat("/ifs/scratch/msph/ehs/eag2186/npbnmf/bootstrap_cor/bs_ewa_", num2str(j), ".mat"), 'EWA_final');
%save(strcat("/ifs/scratch/msph/ehs/eag2186/npbnmf/bootstrap_cor/bs_eh_",  num2str(j), ".mat"), 'EH_final');
