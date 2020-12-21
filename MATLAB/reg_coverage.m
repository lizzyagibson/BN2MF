%% Import the data
%simdata1 = table2array(readtable("/Users/lizzy/BN2MF/Results/Main/dgp_csv/sim_dgp_rep1_1.csv"));
%patterns = table2array(readtable("/Users/lizzy/BN2MF/Results/Main/patterns.csv"));
%scores   = table2array(readtable("/Users/lizzy/BN2MF/Results/Main/scores.csv"));

%% Run BN2MF
%[EWAreg, EHreg, varHreg, alphaHreg, betaHreg,  alphaWreg, betaWreg, alphaAreg, betaAreg, varWAreg, finalscorereg] = NPBayesNMF(simdata1);

%% Coverage from reg run


%% Rearrange solution matrices to match truth
[e,Pi] = factor_correspondence(patterns',EHreg');
EH = (EHreg' * Pi)';
EWA = EWAreg * Pi;

alphaW = alphaWreg * Pi;
betaW  = betaWreg  * Pi;

alphaH = (alphaHreg' * Pi)';
betaH  = (betaHreg'  * Pi)';

alphaA = alphaAreg * Pi;
betaA  = betaAreg  * Pi;

% Scale is inverse rate
% theta is inverse beta
thetaA = 1 ./ betaA;
thetaW = 1 ./ betaW;
thetaH = 1 ./ betaH;

draws = 10000;

%% Take 1000 draws from distributions for W, A, & H
% with alphas and thetas
% This creates an empirical distribution for each
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

% Normalize all H matrices to L1 norm across chemicals
H_denom  = sum(H_dist, 2);
H_scaled = H_dist ./ H_denom;

H_denom_diag = zeros(4, 4, draws);
for i = 1:draws
    H_denom_diag(:, :, i) = diag(H_denom(:,:,i));
end

% Scale all WA matrices by corresponding normalization constant
% This creates a scaled empirical distribution for scores
WA_scaled = zeros(1000, 4, draws);

for i = 1:draws
    WA_scaled(:,:,i) = WA_dist(:,:,i) * H_denom_diag(:,:,i);
end

% Create CI
upper_ci = quantile(WA_scaled, 0.975, 3);
lower_ci = quantile(WA_scaled, 0.025, 3);

% Normalize  truth, too
patterns_denom = sum(patterns, 2);
patterns_scaled = patterns ./ patterns_denom;
patterns_denom_diag = diag(patterns_denom);
scores_scaled = scores * patterns_denom_diag;

sum(sum(scores_scaled <= upper_ci & scores_scaled >= lower_ci)) / (1000*4)