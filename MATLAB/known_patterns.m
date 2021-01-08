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

% Normalize H matrix to L1 norm across chemicals
H_denom  = sum(EH0, 2);
H_scaled = EH0 ./ H_denom;
H_denom_diag = diag(H_denom);

% Scale WA matrix by corresponding normalization constant
WA_scaled = EWA0 * H_denom_diag;

% Scale is inverse rate
% theta is inverse beta
thetaA0 = 1 ./ betaA0;
thetaW0 = 1 ./ betaW0;
thetaH0 = 1 ./ betaH0;

% Create CI with inverse cumulative distribution
% x = gaminv(p,a,b) returns the icdf of the gamma distribution 
% with shape parameter a and the scale parameter b, evaluated at the values in p.
upper_w = gaminv(0.975, alphaW0, thetaW0);
lower_w = gaminv(0.025, alphaW0, thetaW0);

upper_a = gaminv(0.975, alphaA0, thetaA0);
lower_a = gaminv(0.025, alphaA0, thetaA0);

upper_ci = upper_w * diag(upper_a) * H_denom_diag;
lower_ci = lower_w * diag(lower_a) * H_denom_diag;

% Rearrange solution matrices to match truth
[e,Pi] = factor_correspondence(patterns_scaled',H_scaled');
EH = (H_scaled' * Pi)';
EWA = WA_scaled * Pi;

upper  = upper_ci * Pi;
lower  = lower_ci * Pi;

%% Same steps without known patterns
[EWA0_reg, EH0_reg, varH0_reg, alphaH0_reg, betaH0_reg, alphaW0_reg, betaW0_reg, ...
      alphaA0_reg, betaA0_reg, varWA0_reg, finalscore0_reg, final_iter0_reg, INIT_reg] = BN2MF_ones(sim);

pred_reg = EWA0_reg * EH0_reg;

% Normalize H matrix to L1 norm across chemicals
H_denom_reg  = sum(EH0_reg, 2);
H_scaled_reg = EH0_reg ./ H_denom_reg;
H_denom_diag_reg = diag(H_denom_reg);

% Scale WA matrix by corresponding normalization constant
WA_scaled_reg = EWA0_reg * H_denom_diag_reg;

% Scale is inverse rate
% theta is inverse beta
thetaA0_reg = 1 ./ betaA0_reg;
thetaW0_reg = 1 ./ betaW0_reg;
thetaH0_reg = 1 ./ betaH0_reg;

% Create CI with inverse cumulative distribution
% x = gaminv(p,a,b) returns the icdf of the gamma distribution 
% with shape parameter a and the scale parameter b, evaluated at the values in p.
upper_w_reg = gaminv(0.975, alphaW0_reg, thetaW0_reg);
lower_w_reg = gaminv(0.025, alphaW0_reg, thetaW0_reg);

upper_a_reg = gaminv(0.975, alphaA0_reg, thetaA0_reg);
lower_a_reg = gaminv(0.025, alphaA0_reg, thetaA0_reg);

upper_ci_reg = upper_w_reg * diag(upper_a_reg) * H_denom_diag_reg;
lower_ci_reg = lower_w_reg * diag(lower_a_reg) * H_denom_diag_reg;

% Rearrange solution matrices to match truth
[e_reg,Pi_reg] = factor_correspondence(patterns_scaled',H_scaled_reg');
EH_reg  = (H_scaled_reg' * Pi_reg)';
EWA_reg = WA_scaled_reg  * Pi_reg;

upper_reg  = upper_ci_reg * Pi_reg;
lower_reg  = lower_ci_reg * Pi_reg;

%% Compare results

% Results without known patterns
disp("Results without known patterns")
sum(sum(scores_scaled <= upper_reg & scores_scaled >= lower_reg)) / (1000*4)
% 0.8257
norm(sim - pred_reg, "fro") / norm(sim, "fro")
% 0.0637
norm(pre_noise - pred_reg, "fro") / norm(pre_noise, "fro")
% 0.0276
norm(patterns_scaled - EH_reg, "fro") / norm(patterns_scaled, "fro")
% 0.0735
norm(scores_scaled - EWA_reg, "fro") / norm(scores_scaled, "fro")
% 0.0759

% Results with known patterns
disp("Results with known patterns")
sum(sum(scores_scaled <= upper & scores_scaled >= lower)) / (1000*4)
% 0.9605
norm(sim - pred, "fro") / norm(sim, "fro")
% 0.0645
norm(pre_noise - pred, "fro") / norm(pre_noise, "fro")
% 0.0204
norm(patterns_scaled - EH, "fro") / norm(patterns_scaled, "fro")
% 7.3008e-09
norm(scores_scaled - EWA, "fro") / norm(scores_scaled, "fro")
% 0.0314

figure(1);
subplot(2,2,1);
stem(patterns_scaled(1,:));
subplot(2,2,2);
stem(patterns_scaled(2,:));
subplot(2,2,3);
stem(patterns_scaled(3,:));
subplot(2,2,4);
stem(patterns_scaled(4,:));

figure(2);
subplot(2,2,1);
stem(EH(1,:));
subplot(2,2,2);
stem(EH(2,:));
subplot(2,2,3);
stem(EH(3,:));
subplot(2,2,4);
stem(EH(4,:));

figure(3);
subplot(2,2,1);
stem(EH_reg(1,:));
subplot(2,2,2);
stem(EH_reg(2,:));
subplot(2,2,3);
stem(EH_reg(3,:));
subplot(2,2,4);
stem(EH_reg(4,:));


