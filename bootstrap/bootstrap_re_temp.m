%% Bootstrap BN2MF
% Rearrange output from temp
% and output EH

%% Add path with functions
addpath('/ifs/scratch/msph/ehs/eag2186/npbnmf')

% job number = dataset 1:600
%% Get job number
j = getenv('SGE_TASK_ID')

bs_ids = table2array(readtable("/ifs/scratch/msph/ehs/eag2186/Data/bs_ids.csv"));
    
runn  = bs_ids(str2num(j))

%% Choose 1 simulation
patterns = table2array(readtable(strcat("/ifs/scratch/msph/ehs/eag2186/Data/sep_csv_patterns/patterns_sep_", num2str(runn), ".csv")));

%% Normalize truth
patterns_denom      = sum(patterns, 2);
patterns_scaled     = patterns ./ patterns_denom;

path = "/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/";

for boot = 1:150

    disp(strcat("Boot number:", num2str(boot)))

    if isfile(strcat(path, "bootstrap_eh/eh_bs_sim_", num2str(runn), "_bs_", num2str(boot),  ".mat"))
        load(strcat(path, "bootstrap_eh/eh_bs_sim_",   num2str(runn), "_bs_", num2str(boot),  ".mat"));
        % EH_scaled

        %% Rearrange solution to match truth
        [~,Pi]    = factor_correspondence(patterns_scaled',EH_scaled');
        EH_final  = (EH_scaled' * Pi)';

        save(strcat(path, "bootstrap_eh/eh_bs_sim_",   num2str(runn), "_bs_", num2str(boot),  ".mat"), 'EH_final');
    end

    if isfile(strcat(path, "boot_temp/sam_bs_sim_", num2str(runn), "_bs_", num2str(boot),  ".mat"))
        load(strcat(path, "boot_temp/sam_bs_sim_",   num2str(runn), "_bs_", num2str(boot),  ".mat"));
        load(strcat(path, "boot_temp/ewa_bs_sim_",   num2str(runn), "_bs_", num2str(boot),  ".mat"));
        load(strcat(path, "boot_temp/eh_bs_sim_",   num2str(runn), "_bs_", num2str(boot),  ".mat"));
        % sam
        % EWA_scaled
        % EH_scaled

        [k,p] = size(EH_scaled)

        if k == 4
            %% Rearrange solution to match truth
            [~,Pi]    = factor_correspondence(patterns_scaled',EH_scaled');
            EH_final  = (EH_scaled' * Pi)';
            EWA_perm = EWA_scaled * Pi;

            EWA_sam = [sam' EWA_perm];

            save(strcat(path, "bootstrap_ewa/ewa_bs_sim_",   num2str(runn), "_bs_", num2str(boot),  ".mat"), 'EWA_sam');
            save(strcat(path, "bootstrap_eh/eh_bs_sim_",   num2str(runn), "_bs_", num2str(boot),  ".mat"), 'EH_final');
        end
    end
end

disp("Finished!")

