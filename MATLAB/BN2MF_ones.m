function [EWA, EH, varH, alphaH, betaH, alphaW, betaW, ...
    alphaA, betaA, varWA, finalscore, final_iter, init_seed_struct] = BN2MF_ones(X)
% X is the data matrix of non-negative values
% Kinit is the maximum allowable factorization (initial). The algorithm tries to reduce this number.
%       the size of EWA and EH indicate the learned factorization.
% EWA and EH are the left and right matrices of the factorization. Technically, they're the expected
%           values of these matrices according to their approximate posterior variational distributions.
% num_iter is the number of iterations to run. The code terminates based on convergence.

rng('shuffle')
init_seed_struct = rng; % bc default is seed = 0
randn(1000); % Warming up the mersenne twister rng

[dim,N] = size(X);
Kinit = N;

reps = 10;
end_score = zeros(reps, 1);

for i = 1:reps % Choose best of 10 runs
    
h01 = 1; %/Kinit;
h02 = 1;

w01 = 1; %/dim;
w02 = 1;
W1 = gamrnd(ones(dim,Kinit),w02); % changed
W2 = ones(dim,Kinit);             % changed

a01 = 1/Kinit;
a02 = 1;
A1 = ones(1,Kinit)/Kinit; % changed
A2 = ones(1,Kinit);       % changed

H1 = ones(Kinit,N);
H2 = ones(Kinit,N);

K = Kinit;
num_iter = 100000;
score = zeros(num_iter, 1);

for iter = 1:num_iter

    T = 1 + 0.75^(iter-1); % deterministic annealing temperature

    EW = W1./W2;
    X_reshape = repmat(reshape(X',[1 N dim]),[K 1 1]);
    ElnWA = psi(W1) - log(W2) + repmat(psi(A1)-log(A2),dim,1);
    ElnWA_reshape = repmat(reshape(ElnWA',[K 1 dim]),[1 N 1]);
    t1 = max(ElnWA_reshape,[],1);
    ElnWA_reshape = ElnWA_reshape - repmat(t1,[K 1 1]);
    
    ElnH = psi(H1) - log(H2);
    P = bsxfun(@plus,ElnWA_reshape/T,ElnH/T);
    P = exp(P);
    P = bsxfun(@rdivide,P,sum(P,1));
    H1 = 1 + (h01 + sum(P.*X_reshape,3) - 1)/T;
    H2 = (h02 + repmat(sum(EW.*repmat(A1./A2,dim,1),1),N,1)')/T;
    
    W1 = 1 + (w01 + reshape(sum(X_reshape.*P,2),[K dim])' - 1)/T;
    W2 = (w02 + repmat(sum((H1./H2).*repmat((A1./A2)',1,N),2)',dim,1))/T;
    
    A1 = 1 + (a01 + sum(sum(X_reshape.*P,3),2)' - 1)/T;
    A2 = (a02 + sum(W1./W2,1).*sum(H1./H2,2)')/T; 
    
    idx_prune = find(A1./A2 < 10^-3);
    if ~isempty(idx_prune)
      W1(:,idx_prune) = [];
      W2(:,idx_prune) = [];
      A1(idx_prune) = [];
      A2(idx_prune) = [];
      H1(idx_prune,:) = [];
      H2(idx_prune,:) = [];
    end
    K = length(A1);
    
%    stem(A1./A2); colorbar; pause(.5);
%    disp([num2str(iter) ' : ' num2str(sum(sum(abs(X-(W1./W2)*diag(A1./A2)*(H1./H2)))))]);
     score(iter) = sum(sum(abs(X-(W1./W2)*diag(A1./A2)*(H1./H2))));
     
     if iter > 1 && abs(score(iter-1)-score(iter)) < 1e-5  
     break
     end
 
end

end_score(i) = score(find(score,1,'last'));  
disp(['Run Number: ' num2str(i) '. Iter Number: ' num2str(iter) '. Iter Score: ' num2str(end_score(i))]); 

% Among the results, use the fitted variational parameters that achieve the HIGHEST ELBO
if i == 1 || (i > 1 && (end_score(i) >= max(end_score)))
    EWA = (W1./W2)*diag(A1./A2);
    EH = H1./H2;
    varWA = ((W1 .* A1) .* (W1 + A1 + 1)) ./ (W2.^2 .* A2.^2);
    varH = H1 ./ H2.^2;
    alphaH = H1;
    betaH = H2;
    alphaW = W1;
    betaW = W2;
    alphaA = A1;
    betaA = A2;
    finalscore = end_score(i);
    final_iter = iter;
end

end
