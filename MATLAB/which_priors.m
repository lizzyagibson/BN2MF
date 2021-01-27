j = '130';

%% Import the data
simdat   = table2array(readtable(strcat("/Users/lizzy/BN2MF/Sims/Main/sim_dgp_", j, ".csv")));
patterns = table2array(readtable(strcat("/Users/lizzy/BN2MF/Sims/Main/patterns_dgp_", j, ".csv")));
scores   = table2array(readtable(strcat("/Users/lizzy/BN2MF/Sims/Main/scores_dgp_",   j, ".csv")));

% Running these with diff prior values
tic()
[EWA0, EH0, varH0, alphaH0, betaH0, alphaW0, betaW0, ...
    alphaA0, betaA0, varWA0, finalscore0, final_iter0] = BN2MF(simdat);
toc()

tic()
[EWA1, EH1, varH1, alphaH1, betaH1, alphaW1, betaW1, ...
    alphaA1, betaA1, varWA1, finalscore1, final_iter1] = BN2MF_ones(simdat);
toc()

%% Normalize  truth
patterns_denom = sum(patterns, 2);
patterns_scaled = patterns ./ patterns_denom;
scores_scaled = scores * diag(patterns_denom);
        
%% Rearrange solution matrices to match truth
[~,Pi0] = factor_correspondence(patterns',EH0');
EH00 = (EH0' * Pi0)';
EWA00 = EWA0 * Pi0;

[~,Pi1] = factor_correspondence(patterns',EH1');
EH11 = (EH1' * Pi1)';
EWA11 = EWA1 * Pi1;

%% Normalize H to L1 norm across chemicals
EH0_denom  = sum(EH00, 2);
EH0_scaled = EH00 ./ EH0_denom;
EWA0_scaled = EWA00 * diag(EH0_denom);

EH1_denom  = sum(EH11, 2);
EH1_scaled = EH11 ./ EH1_denom;
EWA1_scaled = EWA11 * diag(EH1_denom);

norm(patterns_scaled - EH0_scaled, 'fro')/norm(patterns_scaled, 'fro')
norm(patterns_scaled - EH1_scaled, 'fro')/norm(patterns_scaled, 'fro')

norm(scores_scaled - EH0_scaled, 'fro')/norm(scores_scaled, 'fro')
norm(scores_scaled - EH1_scaled, 'fro')/norm(scores_scaled, 'fro')
