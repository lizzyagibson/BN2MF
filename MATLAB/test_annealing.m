% Import the data
simdata1 = readtable(strcat("/Users/lizzy/BN2MF/Sims/dgp_rep1/sim_dgp_rep1_", num2str(1), ".csv"));
% Convert to output type
simdata1 = table2array(simdata1);

tnot = 0:.0001:0.02;
[~, ni] = size(tnot);
out = zeros(ni,4);

[EWAreg, EHreg, varHreg, alphaHreg, betaHreg] = NPBayesNMF(simdata1);

for i = 1:ni
                                       
    [EWA, EH, varH, alphaH, betaH, final_score, final_iter] = ...
            NPBayesNMF_annealing(simdata1, tnot(i));
    [N,rank] = size(EWA);

    out(i, 1) = rank;
    out(i, 2) = final_score;
    out(i, 3) = final_iter;
    
    try
        varHdiff = varH - varHreg; % difference matrix, positive means annealed is wider
        out(i, 4) = mean(mean(varHdiff));
    catch ME
    end
    
    disp(i)
end

tnot_long = reshape(tnot, 201, 1);
grid_out = [tnot_long,out];
save("/Users/lizzy/BN2MF/MATLAB/grid_out4.m", 'grid_out')
