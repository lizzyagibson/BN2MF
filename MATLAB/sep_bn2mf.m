%% BN2MF on main sims
% All solutions are rearranged to match patterns in truth
% All H matrices are L1 normed
% All WA matrices are scaled by the corresponding normalization constant
j = getenv('SGE_TASK_ID')

    simdata1 = readtable(strcat("/ifs/scratch/msph/ehs/eag2186/Data/sep_csv/sim_q_", num2str(j), ".csv"));
    patterns = readtable(strcat("/ifs/scratch/msph/ehs/eag2186/Data/sep_csv/patterns_q_", num2str(j), ".csv"));
    scores   = readtable(strcat("/ifs/scratch/msph/ehs/eag2186/Data/sep_csv/scores_q_", num2str(j), ".csv"));

    %% Convert to output type
    simdata1 = table2array(simdata1);
    patterns = table2array(patterns);
    scores   = table2array(scores);

    %% Run model
    [EWA0, EH0, varH0, alphaH0, betaH0, alphaW0, betaW0, ...
        alphaA0, betaA0, varWA0, finalscore0, final_iter0] = BN2MF(simdata1);

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

    % Save matrices
    %save(strcat("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_vci_out/q_ewa_NOTscaled", num2str(j), ".mat"), 'EWA');
    %save(strcat("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_vci_out/q_eh_NOTscaled",  num2str(j), ".mat"), 'EH');

    % Save scaled versions, too
    save(strcat("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_vci_out/q_ewa_scaled", num2str(j), ".mat"), 'EWA_scaled');
    save(strcat("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_vci_out/q_eh_scaled",  num2str(j), ".mat"), 'EH_scaled');

    % CI are for scaled WA and norm H
    save(strcat("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_vci_out/q_upperWA_", num2str(j), ".mat"), 'upper_ci_WA');
    save(strcat("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_vci_out/q_lowerWA_", num2str(j), ".mat"), 'lower_ci_WA');

    save(strcat("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_vci_out/q_upperH_", num2str(j), ".mat"), 'upper_ci_H');
    save(strcat("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_vci_out/q_lowerH_", num2str(j), ".mat"), 'lower_ci_H');

    % Save variational distribution arrays, too
    save(strcat("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_vci_out/q_distWA_", num2str(j), ".mat"), 'WA_scaled');
    save(strcat("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_vci_out/q_distEH_", num2str(j), ".mat"), 'H_scaled');

    % Save pred
    save(strcat("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_vci_out/q_pred_", num2str(j), ".mat"), 'pred');

    disp(j)
