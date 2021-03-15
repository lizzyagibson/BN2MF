%% Bootstrap BN2MF
% Each script should take a random subsample with replacement, runn bn2mf
% Then combine results to create empirical distribution for scores
% Want to compare bootstrap CI to VCI

%% Add path with functions
addpath('/ifs/scratch/msph/ehs/eag2186/npbnmf')

% Run script 1:45000 (150 bootstraps * 300 datasets)

%% Get job number
j = getenv('SGE_TASK_ID')

bs_ids = table2array(readtable("/ifs/scratch/msph/ehs/eag2186/Data/bs_ids.csv"));

bootstrap = 1:150;
grid = combvec(bs_ids', bootstrap)';
    
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
rand_vec = [1943390030  746196695   17852396  940589715   42646086 1692962385 1885869580  174393542 1996151257  500383961 1542277254 ...
 2047864733 1801761505  472352621 1901164064  276781489 2144237700 1287261066 1177070671  508867770   54906098 1710116507 ...
   72460519  706288188 2144008053 1522012788  136404331  947831887   63053959  236697809  480539041  900919009 2050901654 ...
 1746339750 1685586903  133255189 1392031761  705707892  545613736 1279444778 1625786061  207200580 1175861432  707538740 ...
  188018141  278412317 1977970720 1341316833 1585172448 1923173550  964272245 1676828299  194795111    1425488  226202817 ...
 1425405539  227303675  400216140  489674218 1883984130 1646108785 1431359662  521739754  759657866  287866418 1330176548 ...
  897105491  580231842  349185096  929438376 2122283525 1509328056  955285452 1654803512 1098486776 1644870105  239994384 ...
 1339596098 1676930598 1216258819  919448435 1462586332 1020334188 1696099450  794891513  640851235 1641996051  314307348 ...
 1950877610  566719717  774764512 1656803790  111177208 1455327775  684682520 1436512824  864953845  308190087 1517416267 ...
 1705328924  215290646 1637507135 1604953147  648376159  610397044  104610036 1062123019 1757689370  616135341  610141128 ...
585529044 1981931409   12184767  611918834  416469288  609065173  777877805  721691688 1106721534 1232428093  139588859 ...
 68379987 1960042482 1531004038 1596145758  381345935  978806661 1026651933  375121438  484083057 1099967675  748655156 ...
 78323434  131801479   66893448  575151865 1653904622 1282574519   12772149  282507755  632920032 1125519556 1209234991 ...
200430331  403780097 1956014448 1281647281 1465010206  724085778  501929009];

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
save(strcat(path, "bootstrap_ewa/eh_bs_sim_",   num2str(runn), "_bs_", num2str(boot),  ".mat"), 'EH_scaled');
save(strcat(path, "bootstrap_ewa/pred_bs_sim_",   num2str(runn), "_bs_", num2str(boot),  ".mat"), 'pred_sam');
% Results go into `get_ex_bootstrap.R`, and `bootstrap_combo.R`
