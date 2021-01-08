%done = readtable("/ifs/scratch/msph/ehs/eag2186/npbnmf/noise_done.csv");
%done = table2array(done(:,2));

% if any(done == j)
%     exit()
% end

%% Get job number
%j = getenv('SGE_TASK_ID')

% Try script locally
j = 1;
simdata1 = table2array(readtable("/Users/lizzy/BN2MF/Results/Sup Noise/noise_sim1.csv"));
patterns = table2array(readtable("/Users/lizzy/BN2MF/Results/Sup Noise/noise_patterns1.csv"));
scores   = table2array(readtable("/Users/lizzy/BN2MF/Results/Sup Noise/noise_scores1.csv"));

% Import the data
% simdata1 = readtable(strcat("/ifs/scratch/msph/ehs/eag2186/Data/dgp_csv/sim_dgp_rep1_", num2str(j), ".csv"));
% patterns = readtable(strcat("/ifs/scratch/msph/ehs/eag2186/Data/dgp_csv/dgp_patterns_", num2str(j), ".csv"));
% scores   = readtable(strcat("/ifs/scratch/msph/ehs/eag2186/Data/dgp_csv/dgp_scores_",   num2str(j), ".csv"));

%% Convert to output type
% simdata1 = table2array(simdata1);
% patterns = table2array(patterns);
% scores   = table2array(scores);

% Normalize truth
patterns_denom = sum(patterns, 2);
patterns_scaled = patterns ./ patterns_denom;
patterns_denom_diag = diag(patterns_denom);
scores_scaled = scores * patterns_denom_diag;
        
[EWA0, EH0, varH0, alphaH0, betaH0, alphaW0, betaW0, ...
    alphaA0, betaA0, varWA0, finalscore0, final_iter0, init_seed_struct] = BN2MF(simdata1);

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

% Create CI with inverse cumulative distribution
% x = gaminv(p,a,b) returns the icdf of the gamma distribution 
% with shape parameter a and the scale parameter b, evaluated at the values in p.
% upper_w = gaminv(0.975, alphaW0, thetaW0);
% lower_w = gaminv(0.025, alphaW0, thetaW0);
% 
% upper_a = gaminv(0.975, alphaA0, thetaA0);
% lower_a = gaminv(0.025, alphaA0, thetaA0);
% 
% upper_ci = upper_w * diag(upper_a) * H_denom_diag;
% lower_ci = lower_w * diag(lower_a) * H_denom_diag;

% alpha = mean^2/var
% theta = var/mean
% beta = mean/var
upper_ci = gaminv(0.975, EWA0.^2 ./ varWA0, varWA0 ./ EWA0); 
lower_ci = gaminv(0.025, EWA0.^2 ./varWA0, varWA0 ./ EWA0);

% Scale 
upper_ci_scaled = upper_ci * H_denom_diag;
lower_ci_scaled = lower_ci * H_denom_diag;

% Rearrange solution matrices to match truth
[e,Pi] = factor_correspondence(patterns_scaled',H_scaled');
EH = (H_scaled' * Pi)';
EWA = WA_scaled * Pi;

% upper  = upper_ci * Pi;
% lower  = lower_ci  * Pi;

upper  = upper_ci_scaled  * Pi;
lower  = lower_ci_scaled  * Pi;

prop = sum(sum(scores_scaled <= upper & scores_scaled >= lower)) / (1000*4)
save_prop = [str2num(j), prop]

% Save prop
% save(strcat("/ifs/scratch/msph/ehs/eag2186/npbnmf/prop_out_cuml/cuml_prop",  num2str(j), ".mat"), 'save_prop');


