simdata1 = readtable("/Users/lizzy/BN2MF/Sims/Main/sim_dgp_201.csv");
patterns = readtable("/Users/lizzy/BN2MF/Sims/Main/patterns_dgp_201.csv");
scores   = readtable("/Users/lizzy/BN2MF/Sims/Main/scores_dgp_201.csv");
true   = readtable("/Users/lizzy/BN2MF/Sims/Main/chem_dgp_201.csv");

%% Convert to output type
simdata1 = table2array(simdata1);
patterns = table2array(patterns);
scores   = table2array(scores);
true    = table2array(true);

% tic()
% [EWA1, EH1, varH1, alphaH1, betaH1, alphaW1, betaW1, ...
%     alphaA1, betaA1, varWA1, finalscore1, final_iter1] = BN2MF(simdata1);
% toc()
% 
% norm(simdata1 - (EWA1 * EH1), 'fro')/norm(simdata1, 'fro')
% %    0.1558
% norm(true - (EWA1 * EH1), 'fro')/norm(simdata1, 'fro')
% %    0.0582

%% Run model
tic()
[EWA0, EH0, varH0, alphaH0, betaH0, alphaW0, betaW0, ...
    alphaA0, betaA0, varWA0, finalscore0, final_iter0] = BN2MF_sparse_h(simdata1);
toc()

norm(simdata1 - (EWA0 * EH0), 'fro')/norm(simdata1, 'fro')
%    0.1626 @ 10^-3
%    0.1614 @ 10^-4
%    0.1605 @ 10^-5
%    0.1618 @ 10^-2
norm(true - (EWA0 * EH0), 'fro')/norm(simdata1, 'fro')
%    0.0743
%    0.0714
%    0.0692
%    0.0726

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
W_dist  = zeros(1000, kk, draws);
A_dist  = zeros(1,    kk, draws);
H_dist  = zeros(kk,   50, draws);
WA_dist = zeros(1000, kk, draws);

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

%% Same steps for reg BN2MF
% Rearrange
% Save pattern number
[kk,~] = size(EH1);

% Rearrange solution matrices to match truth
% If  matrices are the same size, rearrange
% If matrices are not the same size, rearrange by identity
% ie dont rearrange

try
    [~,Pi1] = factor_correspondence(patterns',EH1');
catch
    warning('Matrices not the same size. Replacing permutation matrix with identity.');
    [r,~] = size(EH1)
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
EH_denom11 = sum(EH11, 2);
EH_scaled11 = EH11 ./ EH_denom11;
EWA_scaled11 = EWA11 * diag(EH_denom11);

%% Measure
norm(patterns - EH, 'fro')/norm(patterns, 'fro')
%    0.6972
%    1.5614
%    0.6975
%    0.6968
norm(patterns - EH11, 'fro')/norm(patterns, 'fro')
%    1.3452

norm(scores - EWA, 'fro')/norm(scores, 'fro')
%    2.1840
%    0.5932
%    2.1965
%    2.1845
norm(scores - EWA11, 'fro')/norm(scores, 'fro')
%    0.5743
%norm(scores - EWAp, 'fro')/norm(scores, 'fro')
%    0.0713
    
norm(patterns_scaled - EH_scaled, 'fro')/norm(patterns_scaled, 'fro')
%    0.6002
%    0.5971
%    0.5686
%    0.5898
norm(patterns_scaled - EH_scaled11, 'fro')/norm(patterns_scaled, 'fro')
%    0.4435

norm(scores_scaled - EWA_scaled, 'fro')/norm(scores_scaled, 'fro')
%    0.3695
%    0.3715
%    0.3561
%    0.3656
norm(scores_scaled - EWA_scaled11, 'fro')/norm(scores_scaled, 'fro')
%    0.3160
norm(scores_scaled - EWA_scaledp, 'fro')/norm(scores_scaled, 'fro')
%    0.0714

sum(sum(patterns_scaled >= lower_ci_H & patterns_scaled <= upper_ci_H))/200
%    0.4950
%    0.5 @ 10^-2
sum(sum(patterns_scaled >= lower_ci_H11 & patterns_scaled <= upper_ci_H11))/200
%    0

sum(sum(scores_scaled >= lower_ci_WA & scores_scaled <= upper_ci_WA))/4000
%    0.5600
%    0.5723
%    0.5665 @ 10^-2
sum(sum(scores_scaled >= lower_ci_WA11 & scores_scaled <= upper_ci_WA11))/4000
%    0.6012

lower_ci_H(:,1:5)
patterns_scaled(:,1:5)
EH_scaled(:,1:5)
upper_ci_H(:,1:5)

%% when patterns are known
tic()
[EWAp, EHp, varHp, alphaHp, betaHp, alphaWp, betaWp, ...
    alphaAp, betaAp, varWAp, finalscorep, final_iterp] = BN2MF_patterns(simdata1, patterns);
toc()

% Normalize  expected values
EH_denomp = sum(EHp, 2);
EH_scaledp = EHp ./ EH_denomp;
EWA_scaledp = EWAp * diag(EH_denomp);
