simdata1 = readtable("/Users/lizzy/BN2MF/Sims/Main/sim_dgp_1.csv");
patterns = readtable("/Users/lizzy/BN2MF/Sims/Main/patterns_dgp_1.csv");
scores   = readtable("/Users/lizzy/BN2MF/Sims/Main/scores_dgp_1.csv");
true   = readtable("/Users/lizzy/BN2MF/Sims/Main/chem_dgp_1.csv");

%% Convert to output type
simdata1 = table2array(simdata1);
patterns = table2array(patterns);
scores   = table2array(scores);
true    = table2array(true);

% tic()
% [EWA1, EH1, varH1, alphaH1, betaH1, alphaW1, betaW1, ...
%     alphaA1, betaA1, varWA1, finalscore1, final_iter1] = BN2MF(simdata1);
% toc()

norm(simdata1 - (EWA1 * EH1), 'fro')/norm(simdata1, 'fro')
norm(true - (EWA1 * EH1), 'fro')/norm(simdata1, 'fro')

%% Run model
tic()
[EWA0, EH0, varH0, alphaH0, betaH0, alphaW0, betaW0, ...
    alphaA0, betaA0, varWA0, finalscore0, final_iter0] = BN2MF_sparse_h(simdata1);
toc()

norm(simdata1 - (EWA0 * EH0), 'fro')/norm(simdata1, 'fro')
norm(true - (EWA0 * EH0), 'fro')/norm(simdata1, 'fro')

%% Rearrange
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

size(patterns)
size(EH0)

try
    [~,Pi] = factor_correspondence(patterns',EH0');
catch
    warning('Matrices not the same size. Replacing permutation matrix with identity.');
    [r,~] = size(EH0)
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

%% Rearrange
% Save pattern number
[kk,~] = size(EH1);

% Rearrange solution matrices to match truth
% If  matrices are the same size, rearrange
% If matrices are not the same size, rearrange by identity
% ie dont rearrange

size(patterns)
size(EH1)

try
    [~,Pi1] = factor_correspondence(patterns',EH1');
catch
    warning('Matrices not the same size. Replacing permutation matrix with identity.');
    [r,~] = size(EH0)
    Pi1 = eye(r);
end

EH11 = (EH1' * Pi1)';
EWA11 = EWA1 * Pi1;

alphaW11 = alphaW1 * Pi1;
betaW11  = betaW1  * Pi1;

alphaH11 = (alphaH1' * Pi1)';
betaH11  = (betaH1'  * Pi1)';

alphaA11 = alphaA1 * Pi1;
betaA11  = betaA1  * Pi1;

% Scale is inverse rate
% theta is inverse beta
thetaA11 = 1 ./ betaA11;
thetaW11 = 1 ./ betaW11;
thetaH11 = 1 ./ betaH11;

% Take alphas and thetas
% Take 1000 draws from distributions for W, A, & H

% This creates an empirical distribution for each
W_dist11  = zeros(1000,    kk, draws);
A_dist11  = zeros(1,       kk, draws);
H_dist11  = zeros(kk,      50, draws);
WA_dist11 = zeros(1000,   kk, draws);

for i = 1:draws
    W_dist11(:,:,i)  = gamrnd(alphaW11, thetaW11);
    A_dist11(:,:,i)  = gamrnd(alphaA11, thetaA11);
    H_dist11(:,:,i)  = gamrnd(alphaH11, thetaH11);
    WA_dist11(:,:,i) = W_dist11(:,:,i) * diag(A_dist11(:,:,i));
end

% Normalize all H matrices to L1 norm across chemicals
H_denom11  = sum(H_dist11, 2);
H_scaled11 = H_dist11 ./ H_denom11;

% Create CI
upper_ci_H11 = quantile(H_scaled11, 0.975, 3);
lower_ci_H11 = quantile(H_scaled11, 0.025, 3);

% Scale all WA matrices by corresponding normalization constant
% This creates a scaled empirical distribution for scores
WA_scaled11 = zeros(1000, kk, draws);

for i = 1:draws
    WA_scaled11(:,:,i) = WA_dist11(:,:,i) * diag(H_denom11(:,:,i));
end

% Create CI
upper_ci_WA11 = quantile(WA_scaled11, 0.975, 3);
lower_ci_WA11 = quantile(WA_scaled11, 0.025, 3);

% Normalize  expected values, too
EH_denom11 = sum(EH, 2);
EH_scaled11 = EH11 ./ EH_denom11;
EWA_scaled11 = EWA11 * diag(EH_denom11);

norm(patterns - EH, 'fro')/norm(patterns, 'fro')
norm(patterns - EH11, 'fro')/norm(patterns, 'fro')

norm(scores - EWA, 'fro')/norm(scores, 'fro')
norm(scores - EWA11, 'fro')/norm(scores, 'fro')

norm(patterns_scaled - EH_scaled, 'fro')/norm(patterns_scaled, 'fro')
norm(patterns_scaled - EH_scaled11, 'fro')/norm(patterns_scaled, 'fro')

norm(scores_scaled - EWA_scaled, 'fro')/norm(scores_scaled, 'fro')
norm(scores_scaled - EWA_scaled11, 'fro')/norm(scores_scaled, 'fro')

sum(sum(patterns_scaled >= lower_ci_H & patterns_scaled <= upper_ci_H))/200
sum(sum(patterns_scaled >= lower_ci_H11 & patterns_scaled <= upper_ci_H11))/200

sum(sum(scores_scaled >= lower_ci_WA & scores_scaled <= upper_ci_WA))/4000
sum(sum(scores_scaled >= lower_ci_WA11 & scores_scaled <= upper_ci_WA11))/4000