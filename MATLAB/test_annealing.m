% Import the data
simdata1 = readtable(strcat("/Users/lizzy/BN2MF/Sims/dgp_rep1/sim_dgp_rep1_", num2str(1), ".csv"));
% Convert to output type
simdata1 = table2array(simdata1);

tnot = [0:.0005:0.01 0.01:.001:0.05];
[~, ni] = size(tnot);
out = zeros(ni,6);

[EWAreg, EHreg, varHreg, alphaHreg, betaHreg, varWAreg, finalscorereg] = NPBayesNMF(simdata1);
save("/Users/lizzy/BN2MF/MATLAB/test_t0/ewa_reg.mat", 'EWAreg')
save("/Users/lizzy/BN2MF/MATLAB/test_t0/varwa_reg.mat", 'varWAreg')
save("/Users/lizzy/BN2MF/MATLAB/test_t0/eh_reg.mat", 'EHreg')
save("/Users/lizzy/BN2MF/MATLAB/test_t0/varh_reg.mat", 'varHreg')

for i = 1:ni
                                       
    [EWA, EH, varH, alphaH, betaH, varWA, final_score, final_iter] = ...
            NPBayesNMF_annealing(simdata1, tnot(i));
    save(strcat("/Users/lizzy/BN2MF/MATLAB/test_t0/ewa_t0", num2str(tnot(i)), ".mat"), 'EWA')
    save(strcat("/Users/lizzy/BN2MF/MATLAB/test_t0/varwa_t0", num2str(tnot(i)), ".mat"), 'varWA')
    save(strcat("/Users/lizzy/BN2MF/MATLAB/test_t0/eh_t0", num2str(tnot(i)), ".mat"), 'EH')
    save(strcat("/Users/lizzy/BN2MF/MATLAB/test_t0/varh_t0", num2str(tnot(i)), ".mat"), 'varH')

    [N,rank] = size(EWA);

    out(i, 1) = rank;
    out(i, 2) = final_score;
    out(i, 3) = final_iter;
    out(i, 5) = mean(mean(varH));
    out(i, 6) = mean(mean(varHreg));
    
    try
        varHdiff = varH - varHreg; % difference matrix, positive means annealed is wider
        out(i, 4) = mean(mean(varHdiff));
    catch ME
    end
    
    disp(i)
end

tnot_long = reshape(tnot, 201, 1);
grid_out = [tnot_long,out];
save("/Users/lizzy/BN2MF/MATLAB/test_t0/t0_grid_out.mat", 'grid_out')
