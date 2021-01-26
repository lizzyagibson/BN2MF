a = rand(10,2) * [rand(1,2) 0 0; 0 0 rand(1,2)]
b = [a ; a(1:2,:)]

[EWA, EH, varH, alphaH, betaH, alphaW, betaW, ...
    alphaA, betaA, varWA, finalscore, final_iter, init_seed_struct] = BN2MF(b);

EWA