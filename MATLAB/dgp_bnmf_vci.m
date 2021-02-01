%% BN2MF on main sims
% All solutions are rearranged to match patterns in truth
% All H matrices are L1 normed
% All WA matrices are scaled by the corresponding normalization constant

%% Add path with functions
addpath('/ifs/scratch/msph/ehs/eag2186/npbnmf')

%% Run script 1 time for each simulated dataset (1:300)

%% Get job number
j = getenv('SGE_TASK_ID')

% Import the data
simdata1 = readtable(strcat("/ifs/scratch/msph/ehs/eag2186/Data/dgp_csv/sim_dgp_", j, ".csv"));
patterns = readtable(strcat("/ifs/scratch/msph/ehs/eag2186/Data/dgp_csv/patterns_dgp_", j, ".csv"));
scores   = readtable(strcat("/ifs/scratch/msph/ehs/eag2186/Data/dgp_csv/scores_dgp_",   j, ".csv"));

%% Convert to output type
simdata1 = table2array(simdata1);
patterns = table2array(patterns);
scores   = table2array(scores);

%% Run model
[EWA0, EH0, varH0, alphaH0, betaH0, alphaW0, betaW0, ...
    alphaA0, betaA0, varWA0, finalscore0, final_iter0] = BN2MF(simdata1);

% Save pattern number
[kk,~] = size(EH0);

% Normalize  truth
patterns_denom = sum(patterns, 2);
patterns_scaled = patterns ./ patterns_denom;
scores_scaled = scores * diag(patterns_denom);
        
% Rearrange solution matrices to match truth
% If  matrices are the same size, rearrange
% If matrices are not the same size, rearrange by identity
% ie dont rearrange
try
    [~,Pi] = factor_correspondence(patterns',EH0');
catch
    warning('Matrices not the same size. Replacing permutation matrix with identity.');
    Pi = eye(r);
end

EH = (EH0' * Pi)';
EWA = EWA0 * Pi;

alphaW = alphaW0 * Pi;
betaW  = betaW0  * Pi;

alphaH = (alphaH0' * Pi)';
betaH  = (betaH0'  * Pi)';

alphaA = alphaA0 * Pi;
betaA  = betaA0  * Pi;

% Scale is inverse rate
% theta is inverse beta
thetaA = 1 ./ betaA;
thetaW = 1 ./ betaW;
thetaH = 1 ./ betaH;

% Take alphas and thetas
% Take 1000 draws from distributions for W, A, & H
draws = 1000;

% This creates an empirical distribution for each
W_dist  = zeros(1000,    kk, draws);
A_dist  = zeros(1,       kk, draws);
H_dist  = zeros(kk,      50, draws);
WA_dist = zeros(1000,   kk, draws);

for i = 1:draws
    W_dist(:,:,i)  = gamrnd(alphaW, thetaW);
    A_dist(:,:,i)  = gamrnd(alphaA, thetaA);
    H_dist(:,:,i)  = gamrnd(alphaH, thetaH);
    WA_dist(:,:,i) = W_dist(:,:,i) * diag(A_dist(:,:,i));
end

% Normalize all H matrices to L1 norm across chemicals
H_denom  = sum(H_dist, 2);
H_scaled = H_dist ./ H_denom;

% Create CI
upper_ci_H = quantile(H_scaled, 0.975, 3);
lower_ci_H = quantile(H_scaled, 0.025, 3);

% Scale all WA matrices by corresponding normalization constant
% This creates a scaled empirical distribution for scores
WA_scaled = zeros(1000, kk, draws);

for i = 1:draws
    WA_scaled(:,:,i) = WA_dist(:,:,i) * diag(H_denom(:,:,i));
end

% Create CI
upper_ci_WA = quantile(WA_scaled, 0.975, 3);
lower_ci_WA = quantile(WA_scaled, 0.025, 3);

% Normalize  expected values, too
EH_denom = sum(EH, 2);
EH_scaled = EH ./ EH_denom;
EWA_scaled = EWA * diag(EH_denom);

% Save matrices
save(strcat("/ifs/scratch/msph/ehs/eag2186/npbnmf/main/dgp_vci_out/dgp_ewa_NOTscaled", j, ".mat"), 'EWA');
save(strcat("/ifs/scratch/msph/ehs/eag2186/npbnmf/main/dgp_vci_out/dgp_eh_NOTscaled",  j, ".mat"), 'EH');

% Save scaled versions, too
save(strcat("/ifs/scratch/msph/ehs/eag2186/npbnmf/main/dgp_vci_out/dgp_ewa_scaled", j, ".mat"), 'EWA_scaled');
save(strcat("/ifs/scratch/msph/ehs/eag2186/npbnmf/main/dgp_vci_out/dgp_eh_scaled",  j, ".mat"), 'EH_scaled');

% CI are for scaled WA and norm H
save(strcat("/ifs/scratch/msph/ehs/eag2186/npbnmf/main/dgp_vci_out/dgp_upperWA_", j, ".mat"), 'upper_ci_WA');
save(strcat("/ifs/scratch/msph/ehs/eag2186/npbnmf/main/dgp_vci_out/dgp_lowerWA_", j, ".mat"), 'lower_ci_WA');

save(strcat("/ifs/scratch/msph/ehs/eag2186/npbnmf/main/dgp_vci_out/dgp_upperH_", j, ".mat"), 'upper_ci_H');
save(strcat("/ifs/scratch/msph/ehs/eag2186/npbnmf/main/dgp_vci_out/dgp_lowerH_", j, ".mat"), 'lower_ci_H');

% Save variational distribution arrays, too
save(strcat("/ifs/scratch/msph/ehs/eag2186/npbnmf/main/dgp_vci_out/dgp_distWA_", j, ".mat"), 'WA_scaled');
save(strcat("/ifs/scratch/msph/ehs/eag2186/npbnmf/main/dgp_vci_out/dgp_distEH_", j, ".mat"), 'H_scaled');

% Results go into `combine_main.R`, `get_ex_bootstrap.R`, and `bootstrap_combo.R`
