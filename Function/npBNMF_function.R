# npBNMF
# % X is the data matrix of non-negative values
# % Kinit is the maximum allowable factorization (initial). The algorithm tries to reduce this number.
# %       the size of EWA and EH indicate the learned factorization.
# % EWA and EH are the left and right matrices of the factorization. Technically, they're the expected
# %           values of these matrices according to their approximate posterior variational distributions.
# % num_iter is the number of iterations to run. The code terminates based on convergence currently.

X = matrix(runif(50), 10, 5)

library(matlab) # This is for repmat

NPBayesNMF <- function(X) {

  bnp_switch = 1  # % this turns on/off the Bayesian nonparametric part. On now.
  dim = nrow(X)
  N = ncol(X)
  Kinit = ncol(X)

  nruns = 100
  end_score = matrix(rep(0, times = nruns))
  
  EA = matrix()
  EWA = matrix()
  EH = matrix()
  EW = matrix()

for (i in 1:nruns) {

    K = Kinit
  
    h01 = 1 # % Try h with non-sparse prior, yes!
  # h01 = 1/Kinit
    h02 = 1
  
    w01 = 1
    w02 = 1
  
  # W1 = gamrnd(dim*ones(dim,Kinit), 1/dim);
    W1 = apply(dim*matrix(1, nrow = dim, ncol = Kinit), 2, function(x) rgamma(x, 1/dim))
  
  # % gamrnd(A,B) generates random numbers from the gamma distribution 
  # % with shape parameters in A and scale parameters in B.
  # W2 = dim*ones(dim,Kinit);
    W2 = dim*matrix(1, nrow = dim, ncol = Kinit)
    
    a01 = bnp_switch*1/Kinit + (1-bnp_switch)
    a02 = 1
  
  # A1 = a01 + bnp_switch*1000*ones(1,Kinit)/Kinit;
  # A2 = a02 + bnp_switch*1000*ones(1,Kinit);
    A1 = a01 + bnp_switch*1000*matrix(1, 1, Kinit)/Kinit
    A2 = a02 + bnp_switch*1000*matrix(1, 1, Kinit)
  
  # H1 = ones(Kinit,N);
  # H2 = ones(Kinit,N);
  # % initializing variables
    H1 = matrix(1, Kinit, N)
    H2 = matrix(1, Kinit, N)
  
    num_iter = 100000 # INCREASE
    
  # score = zeros(num_iter, 1);
    score = vector("numeric", length = num_iter)
  
  for (iter in 1:num_iter) {
  
    # EW = W1./W2;
      EW = W1 / W2
    # % element-wise multiplication
    # % this is the expectation of W
    # X_reshape = repmat(reshape(X', [1 N dim]), [K 1 1]);
    # The equivalent of Matlab's repmat(a,2,3) in base R is kronecker(matrix(1,2,3),a)
      X_reshape = kronecker(array(1, dim = c(K,1,1)), array(t(X), c(1, N, dim)))
      
    # % B = repmat(A,n) returns an array containing n copies of A in the row and 
    # % column dimensions. The size of B is size(A)*n when A is a matrix.
    # % Transpose X and repeat it k times as tensor.
    
    # ElnWA = psi(W1) - log(W2) + repmat(psi(A1)-log(A2),dim,1);
      ElnWA = digamma(W1) - log(W2) + repmat(matrix(digamma(A1)-log(A2), ncol = 1), dim, 1)

    # ElnWA_reshape = repmat(reshape(ElnWA',[K 1 dim]),[1 N 1]);
      ElnWA_reshape = kronecker(array(1, dim = c(1,N,1)), array(t(ElnWA), c(K, 1, dim)))
      
    # t1 = max(ElnWA_reshape,[],1);
    # max returns the maximum element along dimension dim. For example, if A is a matrix, 
    # then max(A,[],2) is a column vector containing the maximum value of each row.
      t1 = array(apply(ElnWA_reshape,c(2,3), max), c(1, 5, 10))
    
    # ElnWA_reshape = ElnWA_reshape - repmat(t1,[K 1 1]);
    # % Expected value of log W * expected value of log A
      ElnWA_reshape = ElnWA_reshape - kronecker(array(1, dim = c(K,1,1)), t1)
    
      # ElnH = psi(H1) - log(H2);
    # % expected value of log H
      ElnH = digamma(H1) - log(H2)                               
                                   
    # P = bsxfun(@plus,ElnWA_reshape,ElnH);
      P = ElnWA_reshape + kronecker(array(1, dim=c(1,1,dim)), ElnH)
      P = exp(P)
    # P = bsxfun(@rdivide,P,sum(P,1));
      P = P / kronecker(array(1, dim = c(K,1,1)), array(apply(P, 3, colSums), dim = c(1, K, dim)))
      
      # % P is a probability to put a lower bound on the ELBO
      # % expected value of log WAH is concave
      # % include P to make expectation tractable
      # % These are update steps from optimizing the ELBO 
      # % take the gradient with respect to parameter, set to zero, solve
        H1 = h01 + apply(P * X_reshape, c(1,2), sum)

        H2 = h02 + t(kronecker(matrix(1, nrow = N, ncol = 1), 
                     matrix(colSums(EW * kronecker(matrix(1, nrow = dim, ncol = 1), A1/A2)), nrow=1)))

        W1 = w01 + t(array(apply(X_reshape*P, c(1,3), sum), dim = c(K, dim)))

        W2 = w02 + kronecker(matrix(1, dim, 1), t(rowSums((H1/H2) * kronecker(matrix(1, 1, N), t(A1/A2)))))

        A1 = a01 + bnp_switch * t(rowSums(apply(X_reshape*P, c(1,2), sum)))
        
        A2 = a02 + bnp_switch * (colSums(W1/W2) * t(rowSums(H1/H2)))
    
    # % This is the sparse prior on A, pushing A to zero
    idx_prune = which(A1/A2 < 10^-3)
    
    if (length(idx_prune) >= 1) {
            W1 = W1[,-idx_prune]
            W2 = We[,-idx_prune]
            A1 = A1[,-idx_prune]
            A2 = A2[,-idx_prune]
            H1 = H1[,-idx_prune]
            H2 = H2[,-idx_prune] 
            }
    
    K = length(A1)
  
   score[iter] = base::sum(abs(X- (W1/W2) %*% diag(as.vector(A1/A2)) %*% (H1/H2)))
    
   if (iter %% 100 == 0) {print(paste0("Run Number: ", i, "; Iter Number: ", iter, "; Iter Score: ", round(score[iter], 5)))}
   
    if (iter > 1 && abs(score[iter-1] - score[iter]) < 1e-5) {
      break} # % Convergence criteria!
}

end_score[i] = score[tail(which(score != 0),1)]
    
print(paste0("Run Number: ", i, "; Final Score: ", round(end_score[i], 5)))

# % Among the results, use the fitted variational parameters that achieve the highest ELBO
  if (i == 1 | (i > 1 && (end_score[i] >= max(end_score)))) {
      EA = A1/A2
      EWA = (W1/W2)*diag(A1/A2)
      EH = H1/H2
      EW = (W1/W2)
      varA = A1/(A2^2)
      varW = W1/(W2^2)
      varH = H1/(H2^2)
      }
}

  list(EWA, EH)
  }

NPBayesNMF(X)
