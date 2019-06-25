dat0 = readtable('sim_data.csv');
dat = table2array(dat0);

[W, H, L] = NMF(dat, 5, 'divergence', 100);