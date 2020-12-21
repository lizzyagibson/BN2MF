% Import the data
simdata1 = table2array(readtable("/Users/lizzy/BN2MF/Results/Sup Noise/noise_sim1.csv"));
patterns = table2array(readtable("/Users/lizzy/BN2MF/Results/Sup Noise/noise_patterns1.csv"));
scores   = table2array(readtable("/Users/lizzy/BN2MF/Results/Sup Noise/noise_scores1.csv"));

% Normalize  truth, too
patterns_denom = sum(patterns, 2);
patterns_scaled = patterns ./ patterns_denom;
patterns_denom_diag = diag(patterns_denom);
scores_scaled = scores * patterns_denom_diag;
        
[EWAreg, EHreg, varHreg, alphaHreg, betaHreg,  alphaWreg, ...
    betaWreg, alphaAreg, betaAreg, varWAreg, finalscorereg] = NPBayesNMF(simdata1);

values = 0:0.001:0.02;
[~, sz] = size(values);
out_m = zeros(sz,7);

for j = 1:21
        % Run BN2MF
        [EWA0, EH0, varH0, alphaH0, betaH0, alphaW0, betaW0, ...
           alphaA0, betaA0, varWA0, final_score, final_iter] = BN2MF(simdata1, values(j));
        
        % Rearrange solution matrices to match truth
        [e,Pi] = factor_correspondence(patterns',EH0');
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

        draws = 10000;

        % Take alphas and thetas
        % Take 1000 draws from distributions for W, A, & H
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

        [N,rank] = size(EWA);

        out_m(j, 1) = rank;
        out_m(j, 2) = final_score;
        out_m(j, 3) = final_iter;
        out_m(j, 4) = mean(mean(varH0)) - mean(mean(varHreg));
        out_m(j, 5) = mean(mean(varH0));
        out_m(j, 6) = sum(sum(varH0));
        out_m(j, 7) = sum(sum(scores_scaled <= upper_ci & scores_scaled >= lower_ci)) / (1000*4);
        disp(j)
end

mean(mean(varHreg))
sum(sum(varHreg))

out_m = [values' out_m];
save("/Users/lizzy/BN2MF/MATLAB/test_t0/t0_vci_l1_n.mat", 'out_m')
