%% Read Data
% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 17);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["mecpp", "mehhp", "meohp", "mcpp", "mibp", "mbp", "mbzp", "mep", "mehp", "dcp_24", "dcp_25", "b_pb", "bp_3", "m_pb", "p_pb", "tcs", "bpa"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
edc = readtable("/Users/lizzy/BN2MF/Data/isee_2020_dat.csv", opts);

%% Clear temporary variables
clear opts

%% Clean up Data
% Convert to output type
edc = table2array(edc);

% Calculate std for each column
sd = std(edc); 

% Subtract mean and divide by std
dataNorm = (edc) ./ sd; 

%% Run Model
tic()
[EWA, EH, varH, alphaH, betaH, alphaW, betaW, ...
    alphaA, betaA, varWA, finalscore, final_iter, init_seed_struct] = BN2MF(dataNorm, 50);
toc()

%% Plot Loadings
labels = ["mecpp", "mehhp", "meohp", "mcpp", "mibp", "mbp", "mbzp", "mep", "mehp", ...
    "dcp_24", "dcp_25", "b_pb", "bp_3", "m_pb", "p_pb", "tcs", "bpa"];
 
denom = sum(EH, 2);
EHp = EH ./ denom;

%PLOT
figure(3);
subplot(3,1,1);
stem(EHp(1,:));
set(gca,'XTick',1:size(EHp,2));
set(gca,'XTickLabels',labels);
subplot(3,1,2);
stem(EHp(2,:));
set(gca,'XTick',1:size(EHp,2));
set(gca,'XTickLabels',labels);
subplot(3,1,3);
stem(EHp(3,:));
set(gca,'XTick',1:size(EHp,2));
set(gca,'XTickLabels',labels);
% 2: Final Iter: 352. Final Iter Number: 2533.9419

%% Create VCI

% theta is inverse beta
thetaA = 1 ./ betaA;
thetaW = 1 ./ betaW;
thetaH = 1 ./ betaH;

% Take alphas and thetas
% Take 1000 draws from distributions for W, A, & H
draws = 1000;
[k,~] = size(EH);
[n,p] = size(edc);

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
WA_scaled = zeros(n, k, draws);

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

%% SAVE OUTPUT
save("/Users/lizzy/BN2MF/Apply/mn_bn2mf.mat");
% this saves whole workspace

save("/Users/lizzy/BN2MF/Apply/mn_ewa.mat", 'EWA');
save("/Users/lizzy/BN2MF/Apply/mn_eh.mat", 'EH');

% Save scaled versions, too
save("/Users/lizzy/BN2MF/Apply/mn_ewa_scaled.mat", 'EWA_scaled');
save("/Users/lizzy/BN2MF/Apply/mn_eh_scaled.mat", 'EH_scaled');

% CI are for scaled WA and norm H
save("/Users/lizzy/BN2MF/Apply/mn_upperWA_.mat", 'upper_ci_WA');
save("/Users/lizzy/BN2MF/Apply/mn_lowerWA_.mat", 'lower_ci_WA');

save("/Users/lizzy/BN2MF/Apply/mn_upperH_.mat", 'upper_ci_H');
save("/Users/lizzy/BN2MF/Apply/mn_lowerH_.mat", 'lower_ci_H');

% Save variational distribution arrays, too
save("/Users/lizzy/BN2MF/Apply/mn_distWA_.mat", 'WA_scaled');
save("/Users/lizzy/BN2MF/Apply/mn_distEH_.mat", 'H_scaled');


