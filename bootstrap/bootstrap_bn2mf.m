%% Bootstrap BN2MF
% Each script should take a random subsample with replacement, runn bn2mf
% Then combine results to create empirical distribution for scores
% Want to compare bootstrap CI to VCI

%% Add path with functions
addpath('/ifs/scratch/msph/ehs/eag2186/npbnmf')

% Run script 1:45000 (150 bootstraps * 300 datasets)

%% Get job number
j = getenv('SGE_TASK_ID')

% bs_ids = table2array(readtable("/ifs/scratch/msph/ehs/eag2186/Data/bs_ids.csv"));
grid = (readtable("/ifs/scratch/msph/ehs/eag2186/Data/bs_not_done.csv"));

% size(grid)
% grid(1:5,:)

grid = table2array(grid);

% bootstrap = 1:150;
% grid = combvec(bs_ids', bootstrap)';
    
boot  = grid(str2num(j),2)
runn  = grid(str2num(j),1)

path = "/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/";

if isfile(strcat(path, "bootstrap_ewa/ewa_bs_sim_", num2str(runn), "_bs_", num2str(boot),  ".mat"))
    quit
end

cvx_setup

%% Choose 1 simulation
sim_data = table2array(readtable(strcat("/ifs/scratch/msph/ehs/eag2186/Data/sep_csv_sim/sim_sep_", num2str(runn), ".csv")));
patterns = table2array(readtable(strcat("/ifs/scratch/msph/ehs/eag2186/Data/sep_csv_patterns/patterns_sep_", num2str(runn), ".csv")));
scores   = table2array(readtable(strcat("/ifs/scratch/msph/ehs/eag2186/Data/sep_csv_scores/scores_sep_",   num2str(runn), ".csv")));

%% Normalize truth
patterns_denom      = sum(patterns, 2);
patterns_scaled     = patterns ./ patterns_denom;
scores_scaled       = scores * diag(patterns_denom);

%% Random vector for RNG
rand_vec = [3501228050 1488771025 2731990861   60509156 1037662741 2352646350 3193004673 2346028629 3986961559 1369621208 ...
             2750016082  128153259  211018766 3697941586 1507705599 2442503730 1364482864 2933618539 3948601767 3589406724 ...
             1203967945  803240999 3996458864 2216572004 2555752827 1748335415 3400037340 1195681205  475904722 3510860186 ...
             2901950830 3748236915 2277064864  461389852  603627806 1020428948 2991905064 1831909180  307451952 2927595011 ...
             2980770365 3904224325 2434286117 1002912196 2628122789 2771125337  972673140 3190203559  405649737  760353336 ...
             2066140896 3514513613 3767636615 2407708960 2691887280 3574658937 1152178504 2748580682 2485107360 2285305131 ...
             2427791059 2651603182 2985994218 2197762559 1200226507 1522538244 1318151763  544334623  728527698 3916398212 ...
             1463116818 1354135440 1308028101 3094616549 3878679330 2325577263 1619494124 3013871835  461820161  579198174 ...
             3144508501 1016267344  784931344 4201097273 2073721923 3934667271 3483665572 3573299237 2611322491 1330461909 ...
             3689820424 2977324801 2768129786 3547904863 2316243868  103965476 1864545527 1632041539  156838467 4010711129 ...
             3346744680 2492466462 1175565523 1976840071 1194581023 1340228632 1238147042 3791466126  110431693 4004828744 ...
             3228401615 1818677192    8042860 4243606376  913927848 1301029335  255307966 1320900958 3569406960 2915419372 ...
              956360536  363571424 3265706894 1030998255  746935384 3937621499 1388583926   56203674 1831571796 3779753323 ...
             2709149406 2780791150  248027367 3539999537 1296905559 2828531492 3753131488  897692140 2222191584  711315105 ...
             3227002434 3755969300 2233364172 1966259028  244069829   84326699 2607789365 1894590942 3784745101 3796372477];

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
[~,Pi]    = factor_correspondence(patterns_scaled',EH_scaled');
EH_final  = (EH_scaled' * Pi)';
EWA_perm = EWA_scaled * Pi;

EWA_sam = [sam' EWA_perm];
pred_sam = [sam' pred];

path = "/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/";
save(strcat(path, "bootstrap_ewa/ewa_bs_sim_",   num2str(runn), "_bs_", num2str(boot),  ".mat"), 'EWA_sam');
save(strcat(path, "bootstrap_eh/eh_bs_sim_",   num2str(runn), "_bs_", num2str(boot),  ".mat"), 'EH_scaled');
save(strcat(path, "bootstrap_pred/pred_bs_sim_",   num2str(runn), "_bs_", num2str(boot),  ".mat"), 'pred_sam');
% Results go into `get_ex_bootstrap.R`, and `bootstrap_combo.R`

disp("Finished!")
