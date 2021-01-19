sim = table2array(readtable("/Users/lizzy/BN2MF/Results/Main/sim1.csv"));
pre_noise = table2array(readtable("/Users/lizzy/BN2MF/Results/Main/chem1.csv"));
patterns = table2array(readtable("/Users/lizzy/BN2MF/Results/Main/patterns1.csv"));
scores   = table2array(readtable("/Users/lizzy/BN2MF/Results/Main/scores1.csv"));

% Normalize truth
patterns_denom = sum(patterns, 2);
patterns_scaled = patterns ./ patterns_denom;
patterns_denom_diag = diag(patterns_denom);
scores_scaled = scores * patterns_denom_diag;

%% Run with known patterns
[EWA0, EH0, varH0, alphaH0, betaH0, alphaW0, betaW0, ...
   alphaA0, betaA0, varWA0, finalscore0, final_iter0, init] = BN2MF_patterns(sim, patterns);

pred = EWA0 * EH0;

% Scale is inverse rate
% theta is inverse beta
thetaA = 1 ./ betaA0;
thetaW = 1 ./ betaW0;

draws = 1000;

% Take alphas and thetas
% Take 1000 draws from distributions for W, A, & H
% This creates an empirical distribution for each
W_dist = zeros(1000,   4,  draws);
A_dist = zeros(1,      4,  draws);

WA_dist = zeros(1000,  4,  draws);
A_dist_diag = zeros(4, 4,  draws);

for i = 1:draws
    W_dist(:,:,i) = gamrnd(alphaW0, thetaW);
    A_dist(:,:,i) = gamrnd(alphaA0, thetaA);

    A_dist_diag(:,:,i) = diag(A_dist(:,:,i));
    WA_dist(:,:,i)     = W_dist(:,:,i) * A_dist_diag(:,:,i);
end

% Normalize all H matrices to L1 norm across chemicals
H_denom  = sum(EH0, 2);
H_scaled = EH0 ./ H_denom;

H_denom_diag = diag(H_denom);
EWA_scaled = EWA0 * H_denom_diag;

% Scale all WA matrices by corresponding normalization constant
% This creates a scaled empirical distribution for scores
WA_scaled = zeros(1000, 4, draws);

for i = 1:draws
    WA_scaled(:,:,i) = WA_dist(:,:,i) * H_denom_diag;
end

% Create CI
upper_ci = quantile(WA_scaled, 0.975, 3);
lower_ci = quantile(WA_scaled, 0.025, 3);

sum(sum(scores_scaled <= upper_ci & scores_scaled >= lower_ci)) / (1000*4)

%% Same steps without known patterns
[EWA0_reg, EH0_reg, varH0_reg, alphaH0_reg, betaH0_reg, alphaW0_reg, betaW0_reg, ...
      alphaA0_reg, betaA0_reg, varWA0_reg, finalscore0_reg, final_iter0_reg, INIT_reg] = BN2MF(sim);

pred_reg = EWA0_reg * EH0_reg;

[~,Pi_reg] = factor_correspondence(patterns',EH0_reg');
EH_reg = (EH0_reg' * Pi_reg)';
EWA_reg = EWA0_reg * Pi_reg;

alphaW_reg = alphaW0_reg * Pi_reg;
betaW_reg  = betaW0_reg  * Pi_reg;

alphaH_reg = (alphaH0_reg' * Pi_reg)';
betaH_reg  = (betaH0_reg'  * Pi_reg)';

alphaA_reg = alphaA0_reg * Pi_reg;
betaA_reg  = betaA0_reg  * Pi_reg;

% Scale is inverse rate
% theta is inverse beta
thetaA_reg = 1 ./ betaA_reg;
thetaW_reg = 1 ./ betaW_reg;
thetaH_reg = 1 ./ betaH_reg;

draws = 1000;

% Take alphas and thetas
% Take 1000 draws from distributions for W, A, & H
% This creates an empirical distribution for each
W_dist_reg = zeros(1000,   4,  draws);
A_dist_reg = zeros(1,      4,  draws);
H_dist_reg = zeros(4,      50, draws);
WA_dist_reg = zeros(1000,  4,  draws);
A_dist_diag_reg = zeros(4, 4,  draws);

for i = 1:draws
    W_dist_reg(:,:,i) = gamrnd(alphaW_reg, thetaW_reg);
    A_dist_reg(:,:,i) = gamrnd(alphaA_reg, thetaA_reg);
    H_dist_reg(:,:,i) = gamrnd(alphaH_reg, thetaH_reg);

    A_dist_diag_reg(:,:,i) = diag(A_dist_reg(:,:,i));
    WA_dist_reg(:,:,i)     = W_dist_reg(:,:,i) * A_dist_diag_reg(:,:,i);
end

% Normalize all H matrices to L1 norm across chemicals
H_denom_reg  = sum(H_dist_reg, 2);
H_scaled_reg = H_dist_reg ./ H_denom_reg;

H_denom_diag_reg = zeros(4, 4, draws);
for i = 1:draws
    H_denom_diag_reg(:, :, i) = diag(H_denom_reg(:,:,i));
end

EH_denom_reg  = sum(EH_reg, 2);
EH_scaled_reg = EH_reg ./ EH_denom_reg;
EWA_scaled_reg = EWA_reg * diag(EH_denom_reg);

% Scale all WA matrices by corresponding normalization constant
% This creates a scaled empirical distribution for scores
WA_scaled_reg = zeros(1000, 4, draws);

for i = 1:draws
    WA_scaled_reg(:,:,i) = WA_dist_reg(:,:,i) * H_denom_diag_reg(:,:,i);
end

% Create CI
upper_ci_reg = quantile(WA_scaled_reg, 0.975, 3);
lower_ci_reg = quantile(WA_scaled_reg, 0.025, 3);

sum(sum(scores_scaled <= upper_ci_reg & scores_scaled >= lower_ci_reg)) / (1000*4)

%% Compare results

% Results without known patterns
disp("Results without known patterns")
sum(sum(scores_scaled <= upper_ci_reg & scores_scaled >= lower_ci_reg)) / (1000*4)
% 0.8183
norm(sim - pred_reg, "fro") / norm(sim, "fro")
% 0.0637
norm(pre_noise - pred_reg, "fro") / norm(pre_noise, "fro")
% 0.0276
norm(patterns_scaled - EH_scaled_reg, "fro") / norm(patterns_scaled, "fro")
% 0.0735
norm(scores_scaled - EWA_scaled_reg, "fro") / norm(scores_scaled, "fro")
% 0.0759

% Results with known patterns
disp("Results with known patterns")
sum(sum(scores_scaled <= upper_ci & scores_scaled >= lower_ci)) / (1000*4)
% 0.9575
norm(sim - pred, "fro") / norm(sim, "fro")
% 0.0645
norm(pre_noise - pred, "fro") / norm(pre_noise, "fro")
% 0.0204
norm(patterns_scaled - H_scaled, "fro") / norm(patterns_scaled, "fro")
% 0
norm(scores_scaled - EWA_scaled, "fro") / norm(scores_scaled, "fro")
% 0.0314
norm(patterns - EH0, "fro") / norm(patterns, "fro")
% 0
norm(scores - EWA0, "fro") / norm(scores, "fro")
% 0.0315

figure(1);
subplot(2,2,1);
stem(EH0(1,:));
subplot(2,2,2);
stem(EH0(2,:));
subplot(2,2,3);
stem(EH0(3,:));
subplot(2,2,4);
stem(EH0(4,:));

figure(2);
subplot(2,2,1);
stem(EH_reg(1,:));
subplot(2,2,2);
stem(EH_reg(2,:));
subplot(2,2,3);
stem(EH_reg(3,:));
subplot(2,2,4);
stem(EH_reg(4,:));


