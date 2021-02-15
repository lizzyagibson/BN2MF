simdata1 = readtable("/Users/lizzy/BN2MF/Sims/Main/sim_dgp_201.csv");
patterns = readtable("/Users/lizzy/BN2MF/Sims/Main/patterns_dgp_201.csv");
scores   = readtable("/Users/lizzy/BN2MF/Sims/Main/scores_dgp_201.csv");

%% Convert to output type
simdata1 = table2array(simdata1);
patterns = table2array(patterns);
scores   = table2array(scores);

%% Normalize truth
patterns_denom      = sum(patterns, 2);
patterns_scaled     = patterns ./ patterns_denom;
scores_scaled       = scores * diag(patterns_denom);

bootstraps = 1:150;

[ewa_scaled, ewa, upper_wa, lower_wa] = deal(zeros(1000, 5, 150));
[eh_scaled, eh, upper_h, lower_h] = deal(zeros(4, 50, 150));

for j = bootstraps
    [n, p] = size(simdata1);
    sam = randsample(1:n, n, true);
    sim_sample = simdata1(sam, :);

    %% Run bn2mf
    [EWA0, EH0, varH0, alphaH0, betaH0, alphaW0, betaW0, ...
       alphaA0, betaA0, varWA0, finalscore0, final_iter0, init] = BN2MF(sim_sample, 10);

    % Save pattern number
    [kk,~] = size(EH0);

    % Rearrange solution matrices to match truth
    % If  matrices are the same size, rearrange
    % If matrices are not the same size, rearrange by identity
    % ie dont rearrange
    try
        [~,Pi] = factor_correspondence(patterns',EH0');
    catch
        warning('Matrices not the same size. Replacing permutation matrix with identity.');
        [r,~] = size(EH0);
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
    
    EWA_scaled_sam = [sam' EWA_scaled];
    EWA_sam = [sam' EWA];
    lower_ci_WA_sam = [sam' lower_ci_WA];
    upper_ci_WA_sam = [sam' upper_ci_WA];

    % Save matrices
    ewa_scaled(:,:,j) = EWA_scaled_sam;
    ewa(:,:,j) = EWA_sam;
    upper_wa(:,:,j) = upper_ci_WA_sam;
    lower_wa(:,:,j) = lower_ci_WA_sam;
    
    eh_scaled(:,:,j) = EH_scaled;
    eh(:,:,j) = EH;
    upper_h(:,:,j) = upper_ci_H; 
    lower_h(:,:,j) = lower_ci_H;
    
    if j == 1
        WA_dist_1 = WA_scaled;
        H_dist_1  = H_scaled;
    end
    if j == 2
        WA_dist_2 = WA_scaled;
        H_dist_2  = H_scaled;
    end
    if j == 3
        WA_dist_3 = WA_scaled;
        H_dist_3  = H_scaled;
    end
    
    disp(["Bootstrap number:", num2str(j)])
end

%% Save
save("/Users/lizzy/BN2MF/Bootstrap/vci_on_bs/ewa.mat", 'ewa');
save("/Users/lizzy/BN2MF/Bootstrap/vci_on_bs/ewa_scaled.mat", 'ewa_scaled');

save("/Users/lizzy/BN2MF/Bootstrap/vci_on_bs/eh.mat", 'eh');
save("/Users/lizzy/BN2MF/Bootstrap/vci_on_bs/eh_scaled.mat", 'eh_scaled');

% CI are for scaled WA and norm H
save("/Users/lizzy/BN2MF/Bootstrap/vci_on_bs/upper_wa.mat", 'upper_wa');
save("/Users/lizzy/BN2MF/Bootstrap/vci_on_bs/ower_wa.mat", 'lower_wa');

save("/Users/lizzy/BN2MF/Bootstrap/vci_on_bs/upper_h.mat", 'upper_h');
save("/Users/lizzy/BN2MF/Bootstrap/vci_on_bs/lower_h.mat", 'lower_h');

% Save variational distribution arrays, too
save("/Users/lizzy/BN2MF/Bootstrap/vci_on_bs/WA_dist_1.mat", 'WA_dist_1');
save("/Users/lizzy/BN2MF/Bootstrap/vci_on_bs/H_dist_1.mat", 'H_dist_1');

save("/Users/lizzy/BN2MF/Bootstrap/vci_on_bs/WA_dist_2.mat", 'WA_dist_2');
save("/Users/lizzy/BN2MF/Bootstrap/vci_on_bs/H_dist_2.mat", 'H_dist_2');

save("/Users/lizzy/BN2MF/Bootstrap/vci_on_bs/WA_dist_3.mat", 'WA_dist_3');
save("/Users/lizzy/BN2MF/Bootstrap/vci_on_bs/H_dist_3.mat", 'H_dist_3');

