sim = table2array(readtable("/Users/lizzy/BN2MF/Results/Main/sim1.csv"));
pre_noise = table2array(readtable("/Users/lizzy/BN2MF/Results/Main/chem1.csv"));
patterns = table2array(readtable("/Users/lizzy/BN2MF/Results/Main/patterns1.csv"));
scores   = table2array(readtable("/Users/lizzy/BN2MF/Results/Main/scores1.csv"));
        
[EWA0, EH0, varH0, alphaH0, betaH0, alphaW0, betaW0, ...
   alphaA0, betaA0, varWA0, finalscore0, final_iter0, init] = BN2MF(sim);

[EWA1, EH1, varH1, alphaH1, betaH1, alphaW1, betaW1, ...
   alphaA1, betaA1, varWA1, finalscore1, final_iter1, init1] = BN2MF_ones(sim);

% Same error
norm(sim - EWA0*EH0, 'fro')/norm(sim, 'fro')
norm(sim - EWA1*EH1, 'fro')/norm(sim, 'fro')

[e0, pi0] = factor_correspondence(patterns', EH0');
[e1, pi1] = factor_correspondence(patterns', EH1');

EH0_re = (EH0' * pi0)';
EH1_re = (EH1' * pi1)';

EWA0_re = (EWA0 * pi0);
EWA1_re = (EWA1 * pi1);

% not the same
norm(EWA0_re - EWA1_re, 'fro')
norm(EH0_re  - EH1_re,  'fro')

% scores -- ones better
norm(scores - EWA0_re, 'fro')/norm(scores, 'fro')
norm(scores - EWA1_re, 'fro')/norm(scores, 'fro')

% patterns -- ones better
norm(patterns - EH0_re, 'fro')/norm(patterns, 'fro')
norm(patterns - EH1_re, 'fro')/norm(patterns, 'fro')
