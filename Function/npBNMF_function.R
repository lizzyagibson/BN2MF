# npBNMF
# % X is the data matrix of non-negative values
# % Kinit is the maximum allowable factorization (initial). The algorithm tries to reduce this number.
# %       the size of EWA and EH indicate the learned factorization.
# % EWA and EH are the left and right matrices of the factorization. Technically, they're the expected
# %           values of these matrices according to their approximate posterior variational distributions.
# % num_iter is the number of iterations to run. The code terminates based on convergence currently.

X = matrix(runif(50), 10, 5)

NPBayesNMF <- function(X) {

  bnp_switch = 1
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
    
      w01 = 1
      w02 = 1
  
      W1 = apply(dim*matrix(1, nrow = dim, ncol = Kinit), 2, function(x) rgamma(x, 1/dim))
      W2 = dim*matrix(1, nrow = dim, ncol = Kinit)
      
      a01 = bnp_switch*1/Kinit + (1-bnp_switch)
      a02 = 1
  
      A1 = a01 + bnp_switch*1000*matrix(1, 1, Kinit)/Kinit
      A2 = a02 + bnp_switch*1000*matrix(1, 1, Kinit)

      h01 = 1 # Non-sparse prior
      # h01 = 1/Kinit # Sparse prior
      h02 = 1
      
      H1 = matrix(1, Kinit, N)
      H2 = matrix(1, Kinit, N)
    
      num_iter = 100000
      
      score = vector("numeric", length = num_iter)
    
      for (iter in 1:num_iter) {
      
          EW = W1 / W2
          
          # The equivalent of Matlab's repmat(a,2,3) in base R is kronecker(matrix(1,2,3),a)
          X_reshape = kronecker(array(1, dim = c(K,1,1)), array(t(X), c(1, N, dim)))
          ElnWA = digamma(W1) - log(W2) + kronecker(array(1, dim), digamma(A1)-log(A2))
          ElnWA_reshape = kronecker(array(1, dim = c(1,N,1)), array(t(ElnWA), c(K, 1, dim)))
            
          t1 = array(apply(ElnWA_reshape,c(2,3), max), c(1, 5, 10))
          
          # % Expected value of log W * expected value of log A
          ElnWA_reshape = ElnWA_reshape - kronecker(array(1, dim = c(K,1,1)), t1)
          
          # % expected value of log H
          ElnH = digamma(H1) - log(H2)                               
                                         
          P = ElnWA_reshape + kronecker(array(1, dim=c(1,1,dim)), ElnH)
          P = exp(P)
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
          if (iter > 1 && abs(score[iter-1] - score[iter]) < 1e-5) {break} # Convergence criteria!
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
          alphaH = H1
          betaH = H2
      }
      H_CI_low  <- qgamma(0.025, shape = alphaH, rate = betaH)
      H_CI_high <- qgamma(0.975, shape = alphaH, rate = betaH)
      }
  list(EWA = EWA, EH = EH, H_CI_low = H_CI_low, H_CI_high = H_CI_high)
  }

NPBayesNMF(X)
