%% BN2MF on main sims

%% Distinct example
simdata1 = readtable(strcat("/Users/lizzy/BN2MF/Sims/Test/test_sim_", num2str(101), "_c.csv"));
patterns = readtable(strcat("/Users/lizzy/BN2MF/Sims/Test/test_patterns_", num2str(101), "_c.csv"));
scores   = readtable(strcat("/Users/lizzy/BN2MF/Sims/Test/test_scores_", num2str(101), "_c.csv"));

%% Convert to output type
simdata1 = table2array(simdata1);
patterns = table2array(patterns);
scores   = table2array(scores);

%% Run model
[EWA0, EH0, varH0, alphaH0, betaH0, alphaW0, betaW0, ...
    alphaA0, betaA0, varWA0, finalscore0, final_iter0] = BN2MF(simdata1, 10);

% Save pattern number
[k,~] = size(EH0);

% Normalize  truth
patterns_denom = sum(patterns, 2);
patterns_scaled = patterns ./ patterns_denom;
scores_scaled = scores * diag(patterns_denom);

% Rearrange solution matrices to match truth
% If  matrices are the same size, rearrange
% If matrices are not the same size, rearrange by identity
% ie dont rearrange

[~,Pi] = factor_correspondence(patterns',EH0');
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

[n,p] = size(simdata1);

% This creates an empirical distribution for each
W_dist  = zeros(n, k, draws);
A_dist  = zeros(1, k, draws);
H_dist  = zeros(k, p, draws);
WA_dist = zeros(n, k, draws);

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
WA_scaled = zeros(1000, k, draws);

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

prop1 = sum(sum(scores_scaled <= upper_ci_WA & scores_scaled >= lower_ci_WA )) /4000