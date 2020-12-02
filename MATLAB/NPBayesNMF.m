function [EWA, EH, varH, alphaH, betaH, varWA] = NPBayesNMF(X)
    %EWA_low, EH_low, varWA_low, varH_low, alphaH_low, betaH_low] = NPBayesNMF(X)
% X is the data matrix of non-negative values
% Kinit is the maximum allowable factorization (initial). The algorithm tries to reduce this number.
%       the size of EWA and EH indicate the learned factorization.
% EWA and EH are the left and right matrices of the factorization. Technically, they're the expected
%           values of these matrices according to their approximate posterior variational distributions.
% num_iter is the number of iterations to run. The code terminates based on convergence currently.

bnp_switch = 1;  % this turns on/off the Bayesian nonparametric part. On now.

[dim,N] = size(X);

end_score = zeros(10, 1);

%Not sure if i need this
EWA = [];
EH = [];
varWA = [];
varH = [];
alphaH = [];
betaH = [];

Kinit = N;

for i = 1:10 % Switch back to 1:10

K = Kinit;
 
%h01 = 1/Kinit;
h01 = 1; % Try h with non-sparse prior, yes!
h02 = 1;

w01 = 1;
%w01 = 1/dim;
w02 = 1;
W1 = gamrnd(dim*ones(dim,Kinit),1/dim);
% gamrnd(A,B) generates random numbers from the gamma distribution 
% with shape parameters in A and scale parameters in B.
W2 = dim*ones(dim,Kinit);

a01 = bnp_switch*1/Kinit + (1-bnp_switch);
a02 = 1;
A1 = a01 + bnp_switch*1000*ones(1,Kinit)/Kinit;
A2 = a02 + bnp_switch*1000*ones(1,Kinit);

H1 = ones(Kinit,N);
H2 = ones(Kinit,N);
% initializing variables

num_iter = 100000;
score = zeros(num_iter, 1);

for iter = 1:num_iter
  
    EW = W1./W2;
    % element-wise multiplication
    % this is the expectation of W
    X_reshape = repmat(reshape(X',[1 N dim]),[K 1 1]);
    % B = repmat(A,n) returns an array containing n copies of A in the row and 
    % column dimensions. The size of B is size(A)*n when A is a matrix.
    % Transpose X and repeat it k times as tensor.

    ElnWA = psi(W1) - log(W2) + repmat(psi(A1)-log(A2),dim,1);
    ElnWA_reshape = repmat(reshape(ElnWA',[K 1 dim]),[1 N 1]);
    t1 = max(ElnWA_reshape,[],1);
    ElnWA_reshape = ElnWA_reshape - repmat(t1,[K 1 1]);
    % Expected value of log W * expected value of log A
    % why subtract t1?

    ElnH = psi(H1) - log(H2);
    % expected value of log H

    P = bsxfun(@plus,ElnWA_reshape,ElnH);
    P = exp(P);
    P = bsxfun(@rdivide,P,sum(P,1));
    % P is a probability to put a lower bound on the ELBO
    % expected value of log WAH is concave
    % include P to make expectation tractable

    % These are update steps from optimizing the ELBO 
    % take the gradient with respect to parameter, set to zero, solve
    H1 = h01 + sum(P.*X_reshape,3);
    H2 = h02 + repmat(sum(EW.*repmat(A1./A2,dim,1),1),N,1)';
    
    W1 = w01 + reshape(sum(X_reshape.*P,2),[K dim])';
    W2 = w02 + repmat(sum((H1./H2).*repmat((A1./A2)',1,N),2)',dim,1);
    
    A1 = a01 + bnp_switch*sum(sum(X_reshape.*P,3),2)';
    A2 = a02 + bnp_switch*sum(W1./W2,1).*sum(H1./H2,2)'; 
    
    % This is the sparse prior on A, pushing A to zero
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

end

end
