%% Make this whole thing loop

for i = 1:100

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 50);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["pcb28", "pcb66", "pcb74", "pcb99", "pcb105", "pcb118", "pcb138_158", "pcb146", "pcb153", "pcb156", "pcb167", "pcb170", "pcb178", "pcb183", "pcb187", "pcb180", "pcb189", "pcb194", "pcb196_203", "pcb199", "pcb206", "pcb209", "BDE17", "BDE28", "BDE47", "BDE66", "BDE85", "BDE99", "BDE100", "BDE153", "BDE154", "BDE183", "BDE209", "MECPP", "MEHHP", "MEOHP", "MCPP", "MIBP", "MBP", "MBZP", "MEP", "MEHP", "dcp_24", "dcp_25", "b_pb", "bp_3", "m_pb", "p_pb", "tcs", "bpa"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
simdata2 = readtable(strcat("/Users/lizzy/OneDrive/Columbia/Spring 2020/nmf/Sims/Iterate/sim_over_", num2str(i), ".csv"), opts);

%% Convert to output type
simdata2 = table2array(simdata2);

%% Clear temporary variables
clear opts

[ewa2,eh2, a2] = NPBayesNMF(simdata2, 50, 10000);

%labels = {"pcb28", "pcb66", "pcb74", "pcb99", "pcb105", "pcb118", "pcb138_158" "pcb146", "pcb153", "pcb156", "pcb167", "pcb170", "pcb178", "pcb183", "pcb187", "pcb180", "pcb189", "pcb194", "pcb196_203" "pcb199", "pcb206", "pcb209", "BDE17", "BDE28", "BDE47", "BDE66", "BDE85", "BDE99", "BDE100", "BDE153", "BDE154", "BDE183", "BDE209", "MECPP", "MEHHP", "MEOHP", "MCPP", "MIBP", "MBP", "MBZP", "MEP", "MEHP", "dcp_24", "dcp_25", "b_pb", "bp_3", "m_pb", "p_pb", "tcs", "bpa"};

%%PLOT
% figure;
% subplot(2,2,1);
% stem(eh2(1,:));
% set(gca,'XTick',1:size(eh2,2));
% set(gca,'XTickLabels',labels);
% subplot(2,2,2);
% stem(eh2(2,:));
% set(gca,'XTick',1:size(eh2,2));
% set(gca,'XTickLabels',labels);
% subplot(2,2,3);
% stem(eh2(3,:));
% set(gca,'XTick',1:size(eh2,2));
% set(gca,'XTickLabels',labels);
% subplot(2,2,4);
% stem(eh2(4,:));
% set(gca,'XTick',1:size(eh2,2));
% set(gca,'XTickLabels',labels);

save(strcat("/Users/lizzy/OneDrive/Columbia/Spring 2020/nmf/Sims/Iterate_Out/ewa_over_", num2str(i), ".mat"), 'ewa2');
save(strcat("/Users/lizzy/OneDrive/Columbia/Spring 2020/nmf/Sims/Iterate_Out/eh_over_", num2str(i), ".mat"), 'eh2');

end
