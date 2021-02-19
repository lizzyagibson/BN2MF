%% BN2MF on main sims

%% Distinct example
simdata101 = readtable(strcat("/Users/lizzy/BN2MF/Sims/Test/test_sim_", num2str(101), ".csv"));
patterns101 = readtable(strcat("/Users/lizzy/BN2MF/Sims/Test/test_patterns_", num2str(101), ".csv"));
scores101   = readtable(strcat("/Users/lizzy/BN2MF/Sims/Test/test_scores_", num2str(101), ".csv"));

%% Convert to output type
simdata101 = table2array(simdata101);
patterns101 = table2array(patterns101);
scores101   = table2array(scores101);

%% Run model
[EWA0101, EH0101, varH0101, alphaH0101, betaH0101, alphaW0101, betaW0101, ...
    alphaA0101, betaA0101, varWA0101, finalscore0101, final_iter0101] = BN2MF(simdata101, 10);

% Save pattern number
[k101,~] = size(EH0101);

% Normalize  truth
patterns_denom101 = sum(patterns101, 2);
patterns_scaled101 = patterns101 ./ patterns_denom101;
scores_scaled101 = scores101 * diag(patterns_denom101);

% Rearrange solution matrices to match truth
% If  matrices are the same size, rearrange
% If matrices are not the same size, rearrange by identity
% ie dont rearrange

[~,Pi101] = factor_correspondence(patterns101',EH0101');
EH101 = (EH0101' * Pi101)';
EWA101 = EWA0101 * Pi101;

alphaW101 = alphaW0101 * Pi101;
betaW101  = betaW0101  * Pi101;

alphaH101 = (alphaH0101' * Pi101)';
betaH101  = (betaH0101'  * Pi101)';

alphaA101 = alphaA0101 * Pi101;
betaA101  = betaA0101  * Pi101;

% Scale is inverse rate
% theta is inverse beta
thetaA101 = 1 ./ betaA101;
thetaW101 = 1 ./ betaW101;
thetaH101 = 1 ./ betaH101;

% Take alphas and thetas
% Take 1000 draws from distributions for W, A, & H
draws = 1000;

[n101,p101] = size(simdata101);

% This creates an empirical distribution for each
W_dist101  = zeros(n101, k101, draws);
A_dist101  = zeros(1, k101, draws);
H_dist101  = zeros(k101, p101, draws);
WA_dist101 = zeros(n101, k101, draws);

for i = 1:draws
    W_dist101(:,:,i)  = gamrnd(alphaW101, thetaW101);
    A_dist101(:,:,i)  = gamrnd(alphaA101, thetaA101);
    H_dist101(:,:,i)  = gamrnd(alphaH101, thetaH101);
    WA_dist101(:,:,i) = W_dist101(:,:,i) * diag(A_dist101(:,:,i));
end

% Normalize all H matrices to L1 norm across chemicals
H_denom101  = sum(H_dist101, 2);
H_scaled101 = H_dist101 ./ H_denom101;

% Create CI
upper_ci_H101 = quantile(H_scaled101, 0.975, 3);
lower_ci_H101 = quantile(H_scaled101, 0.025, 3);

% Scale all WA matrices by corresponding normalization constant
% This creates a scaled empirical distribution for scores
WA_scaled101 = zeros(1000, k101, draws);

for i = 1:draws
    WA_scaled101(:,:,i) = WA_dist101(:,:,i) * diag(H_denom101(:,:,i));
end

% Create CI
upper_ci_WA101 = quantile(WA_scaled101, 0.975, 3);
lower_ci_WA101 = quantile(WA_scaled101, 0.025, 3);

% Normalize  expected values, too
EH_denom101 = sum(EH101, 2);
EH_scaled101 = EH101 ./ EH_denom101;
EWA_scaled101 = EWA101 * diag(EH_denom101);

prop101 = sum(sum(scores_scaled101 <= upper_ci_WA101 & scores_scaled101 >= lower_ci_WA101 )) /4000