%% Bootstrap BN2MF
% Each script should take a random subsample with replacement, runn bn2mf
% Then combine results to create empirical distribution for scores
% Want to compare bootstrap CI to VCI

%% Add path with functions
addpath('/ifs/scratch/msph/ehs/eag2186/npbnmf')

% Run script 1:1773

%% Get job number
j = getenv('SGE_TASK_ID')

grid = (readtable("/ifs/scratch/msph/ehs/eag2186/Data/bs_not_done.csv"));

grid = table2array(grid);
    
boot  = grid(str2num(j),2)
runn  = grid(str2num(j),1)

path = "/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/";

if isfile(strcat(path, "bootstrap_ewa/ewa_bs_sim_", num2str(runn), "_bs_", num2str(boot),  ".mat"))
    quit
end

%% Choose 1 simulation
sim_data = table2array(readtable(strcat("/ifs/scratch/msph/ehs/eag2186/Data/sep_csv_sim/sim_sep_", num2str(runn), ".csv")));
patterns = table2array(readtable(strcat("/ifs/scratch/msph/ehs/eag2186/Data/sep_csv_patterns/patterns_sep_", num2str(runn), ".csv")));
scores   = table2array(readtable(strcat("/ifs/scratch/msph/ehs/eag2186/Data/sep_csv_scores/scores_sep_",   num2str(runn), ".csv")));

%% Normalize truth
patterns_denom      = sum(patterns, 2);
patterns_scaled     = patterns ./ patterns_denom;
scores_scaled       = scores * diag(patterns_denom);

%% Random vector for RNG
rand_vec = [4282307143 1806115923 1299418625 1014640254 3486320955 3475949215 2438382271 2749185709 3478516535 1495918571 ...
4214076998 2967289476  936457185 3576372199  314821303 2595507463 1063843859 1543561191  824609675 1712896119 ...
2275510221  906902719 1102768280 1565382719 2458103832 1934259547 4110530065 3172261655 1756069398 4117109288 ...
 751378139 3358094859 2951244163 3749366427 1204141114 1014653924 2916334112 3916729194 2084175308 3106886505 ...
1762218979  901393025  590875514  458515468 3831944298 1885301291 1988275923 3401519643  546660693 1579230951 ...
2818904091 3320964346 1742591982  823234263   21775846  914306680 2058656664 1517194933  985642925 2328230523 ...
2057029767  104176713 3902126681 1122338440 1694204895 1066444215 4146081778 2894090010 3757562589 3300392865 ...
 867410024  180046225 3685856552  653264723 1408870447 2766712039   97505463   68259393 2604995099 3350156526 ...
 681648631 3504080718 2497860746  685664700 2317772749 3596614000 3309305597 1571092971 3395861944 3950927040 ...
4249714110 4143734009 3200663752 3511488724 1674961676 1935407036 2367228338  370450842 1281451333 2528876247 ...
3134382177 2239935549 3796876962 4281018602 1266389374 2588202377   11665238  330352859 2287076630 2191092998 ...
2817540797  333458979 1552900832  985229275 4008917355 1508945993 2389526045 1795703665 3810290076  465843161 ...
2094271746  954032925 4075942401 1290954149 2612522934 2899208188 4055866823 3823914385 2959294343 2807454478 ...
3679669798 2537314458 3518769387 4247397712 4244591706 2196559874 3809840576 1648399904  522653552 1621644349 ...
 367017685 1867240871 2470122043 4248258747  409378864 1268791057 3541264323 2066601472 1704569774 3426526461];

% Set up rng
rng(rand_vec(boot))
init_seed_struct = rng; % bc default is seed = 0
randn(1000); % Warming up the mersenne twister rng

%% Take a bootstrapped sample
[n, p] = size(sim_data)
sam = randsample(1:n, n, true);
sim_sample = sim_data(sam, :);

%% Run bn2mf
[EWA0, EH0, varH0, alphaH0, betaH0, alphaW0, betaW0, ...
   alphaA0, betaA0, varWA0, finalscore0, final_iter0, init_seed] = BN2MF(sim_sample);

size(EH0)

pred = EWA0 * EH0;

%% Normalize EH to L1 norm across chemicals
H_denom   = sum(EH0, 2);
EH_scaled = EH0 ./ H_denom;

%% Scale EWA by corresponding normalization constant
EWA_scaled   = EWA0 * diag(H_denom);

%% Rearrange solution to match truth
% [~,Pi]    = factor_correspondence(patterns_scaled',EH_scaled');
% EH_final  = (EH_scaled' * Pi)';
% EWA_perm = EWA_scaled * Pi;

% EWA_sam = [sam' EWA_perm];

pred_sam = [sam' pred];

path = "/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/";
save(strcat(path, "boot_temp/sam_bs_sim_",   num2str(runn), "_bs_", num2str(boot),  ".mat"), 'sam');
save(strcat(path, "boot_temp/ewa_bs_sim_",   num2str(runn), "_bs_", num2str(boot),  ".mat"), 'EWA_scaled');
save(strcat(path, "boot_temp/eh_bs_sim_",   num2str(runn), "_bs_", num2str(boot),  ".mat"), 'EH_scaled');

save(strcat(path, "bootstrap_pred/pred_bs_sim_",   num2str(runn), "_bs_", num2str(boot),  ".mat"), 'pred_sam');
% Results go into `get_ex_bootstrap.R`, and `bootstrap_combo.R`

disp("Finished!")
