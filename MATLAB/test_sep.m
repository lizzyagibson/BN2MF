%% Check noise effect

% Same dataset with overlapping patterns
% increasing noise level 0--1
% Known patterns plugged into algorithm
% Does inference get worse as noise increases?

%% BN2MF on main sims
% All solutions are rearranged to match patterns in truth
% All H matrices are L1 normed
% All WA matrices are scaled by the corresponding normalization constant

out = zeros(11,7);

for j = 1:11

    simdata1 = readtable(strcat("/Users/lizzy/BN2MF/Sims/sep_csv/sim_sep_", num2str(j), ".csv"));
    patterns = readtable(strcat("/Users/lizzy/BN2MF/Sims/sep_csv/patterns_sep_", num2str(j), ".csv"));
    scores   = readtable(strcat("/Users/lizzy/BN2MF/Sims/sep_csv/scores_sep_", num2str(j), ".csv"));
    
    %% Convert to output type
    simdata1 = table2array(simdata1);
    patterns = table2array(patterns);
    scores   = table2array(scores);
    chem = scores * patterns;

    %% Run model
    [EWA, EH, varH, alphaH, betaH, alphaW, betaW, ...
        alphaA, betaA, varWA, finalscore, final_iter] = BN2MF_patterns(simdata1, patterns);

    pred = EWA * EH;
    
    % Save pattern number
    [k,~] = size(EH);

    % Normalize  truth
    patterns_denom = sum(patterns, 2);
    patterns_scaled = patterns ./ patterns_denom;
    scores_scaled = scores * diag(patterns_denom);

    % Rearrange solution matrices to match truth
    % If  matrices are the same size, rearrange
    % If matrices are not the same size, rearrange by identity
    % ie dont rearrange

    thetaA = 1 ./ betaA;
    thetaH = 1 ./ betaH;
    thetaW = 1 ./ betaW;
    
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

    % Compare error and coverage across noise levels
    prop = sum(sum(scores_scaled <= upper_ci_WA & scores_scaled >= lower_ci_WA))/4000
    pred_err = norm(chem - pred, 'fro')/norm(chem, 'fro')
    load_err = norm(patterns - EH, 'fro')/norm(patterns, 'fro')
    score_err = norm(scores - EWA, 'fro')/norm(scores, 'fro')
    load_err2 = norm(patterns_scaled - EH_scaled, 'fro')/norm(patterns_scaled, 'fro')
    score_err2 = norm(scores_scaled - EWA_scaled, 'fro')/norm(scores_scaled, 'fro')
    
    out(j,:) = [j prop pred_err load_err score_err load_err2 score_err2] ;
end    