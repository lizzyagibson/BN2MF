%% Import data from text file
% Script for importing data from the following text file:
%
%    filename: /Users/lizzy/OneDrive/Columbia/Spring 2020/nmf/Data/edc_pred.csv
%
% Auto-generated by MATLAB on 06-Mar-2020 15:11:47

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 23);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["P1x", "P2x", "P3x", "sid", "P1y", "P2y", "P3y", "eth", "edu", "age", "care", "lotion", "perfume", "liquid_soap", "hair_gel", "hair_spray", "nail_polish", "makeup", "lipstick", "eye_makeup", "sunscreen", "shampoo", "creme_rinse"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["eth", "edu", "lotion", "perfume", "liquid_soap", "hair_gel", "hair_spray", "nail_polish", "makeup", "lipstick", "eye_makeup", "sunscreen", "shampoo", "creme_rinse"], "TrimNonNumeric", true);
opts = setvaropts(opts, ["eth", "edu", "lotion", "perfume", "liquid_soap", "hair_gel", "hair_spray", "nail_polish", "makeup", "lipstick", "eye_makeup", "sunscreen", "shampoo", "creme_rinse"], "ThousandsSeparator", ",");

% Import the data
edcpred = readtable("/Users/lizzy/OneDrive/Columbia/Spring 2020/nmf/Data/edc_pred.csv", opts);

%% Convert to output type
edcpred = table2array(edcpred);

%% Clear temporary variables
clear opts

[ewa_pred,eh_pred] = NPBayesNMF(edc_pred, 51, 10000);

labels = ["pcb28", "pcb66", "pcb74", "pcb99", "pcb105", "pcb114", "pcb118", "pcb138_158" "pcb146", "pcb153", ...
    "pcb156", "pcb157", "pcb167", "pcb170", "pcb178", "pcb180", "pcb183", "pcb187", "pcb189", "pcb194", ...
    "pcb196_203" "pcb199", "pcb206", "pcb209", "BDE28", "BDE47", "BDE66", "BDE85", "BDE99", "BDE100", "BDE153", ...
    "BDE154", "BDE183", "BDE209", "MECPP", "MEHHP", "MEOHP", "MCPP", "MIBP", "MBP", "MBZP", "MEP", "MEHP", ...
    "dcp_24", "dcp_25", "b_pb", "bp_3", "m_pb", "p_pb", "tcs", "bpa", "care", "lotion", "perfume", ...
    "liquid_soap", "hair_gel", "hair_spray", "nail_polish", "makeup", "lipstick", "eye_makeup", ...
    "sunscreen", "shampoo", "creme_rinse"];
       
%PLOT
figure;
subplot(3,1,1);
stem(eh_pred(1,:));
set(gca,'XTick',1:size(eh_pred,2));
set(gca,'XTickLabels',labels);
subplot(3,1,2);
stem(eh_pred(2,:));
set(gca,'XTick',1:size(eh_pred,2));
set(gca,'XTickLabels',labels);
subplot(3,1,3);
stem(eh_pred(3,:));
set(gca,'XTick',1:size(eh_pred,2));
set(gca,'XTickLabels',labels);

save("/Users/lizzy/OneDrive/Columbia/Spring 2020/nmf/Sims/Iterate_Out/ewa_edc_pred.mat", 'ewa_pred');
save("/Users/lizzy/OneDrive/Columbia/Spring 2020/nmf/Sims/Iterate_Out/eh_edc_pred.mat", 'eh_pred');