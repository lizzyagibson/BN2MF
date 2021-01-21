%% Compare EWA to bootstrap results

%% Choose 1 example of correlated simulations
sim       = table2array(readtable("/Users/lizzy/BN2MF/Results/Main/Corr Ex/corr_sim.csv"));
pre_noise = table2array(readtable("/Users/lizzy/BN2MF/Results/Main/Corr Ex/corr_chem.csv"));
patterns  = table2array(readtable("/Users/lizzy/BN2MF/Results/Main/Corr Ex/corr_patterns.csv"));
scores    = table2array(readtable("/Users/lizzy/BN2MF/Results/Main/Corr Ex/corr_scores.csv"));

%% Normalize truth
patterns_denom      = sum(patterns, 2);
patterns_scaled     = patterns ./ patterns_denom;
patterns_denom_diag = diag(patterns_denom);
scores_scaled       = scores * patterns_denom_diag;

%% Run bn2mf
[EWA0, EH0, varH0, alphaH0, betaH0, alphaW0, betaW0, ...
   alphaA0, betaA0, varWA0, finalscore0, final_iter0, init] = BN2MF(sim);

pred = EWA0 * EH0;

%% Rearrange solution matrices to match truth
[e,Pi] = factor_correspondence(patterns',EH0');
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

draws = 1000;

%% Create an empirical distribution for each
% Take alphas and thetas
% Take 1000 draws from distributions for W, A, & H
W_dist = zeros(1000,   4,  draws);
A_dist = zeros(1,      4,  draws);
H_dist = zeros(4,      50, draws);
WA_dist = zeros(1000,  4,  draws);
A_dist_diag = zeros(4, 4,  draws);

for i = 1:draws
    W_dist(:,:,i) = gamrnd(alphaW, thetaW);
    A_dist(:,:,i) = gamrnd(alphaA, thetaA);
    H_dist(:,:,i) = gamrnd(alphaH, thetaH);

    A_dist_diag(:,:,i) = diag(A_dist(:,:,i));
    WA_dist(:,:,i)     = W_dist(:,:,i) * A_dist_diag(:,:,i);
end

%% Normalize all H matrices to L1 norm across chemicals
H_denom  = sum(H_dist, 2);
H_scaled = H_dist ./ H_denom;

H_denom_diag = zeros(4, 4, draws);
for i = 1:draws
    H_denom_diag(:, :, i) = diag(H_denom(:,:,i));
end

%% Normalize EH to L1 norm across chemicals
H_denom   = sum(EH, 2);
EH_scaled = EH ./ H_denom;

%% Scale all WA matrices by corresponding normalization constant
% This creates a scaled empirical distribution for scores
WA_scaled = zeros(1000, 4, draws);

for i = 1:draws
    WA_scaled(:,:,i) = WA_dist(:,:,i) * H_denom_diag(:,:,i);
end

%% Scale EWA by corresponding normalization constant
H_denom_diag = diag(H_denom);
EWA_scaled   = EWA * H_denom_diag;

%% Create CI
upper_ci = quantile(WA_scaled, 0.975, 3);
lower_ci = quantile(WA_scaled, 0.025, 3);

%% Variational distribution is WA_scaled
sum(sum(scores_scaled <= upper_ci & scores_scaled >= lower_ci)) / (1000*4)

save(strcat("/Users/lizzy/BN2MF/Bootstrap/compare/cor_EWA.mat"), 'EWA_scaled');
save(strcat("/Users/lizzy/BN2MF/Bootstrap/compare/cor_EH.mat"), 'EH_scaled');
save(strcat("/Users/lizzy/BN2MF/Bootstrap/compare/cor_WA_upper.mat"), 'upper_ci');
save(strcat("/Users/lizzy/BN2MF/Bootstrap/compare/cor_WA_lower.mat"), 'lower_ci');
save(strcat("/Users/lizzy/BN2MF/Bootstrap/compare/cor_var_dist_WA.mat"), 'WA_scaled');
