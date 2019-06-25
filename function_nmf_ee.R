##################
## NMF Function ##
## 6/25/2019 #####
##################

# MATLAB

# function [W,H,L] = NMF(X,K,penalty,num_iter)
# % NMF implements non-negative matrix factorization of the matrix X
# % according to a predefined penalty.
# %
# % X : D x N matrix of nonnegative numbers
# % K : rank of factorization, takes value in {1,2,3,...,min(D,N)}
# % penalty : string taking one of two values, 'L2' or 'divergence'
# % num_iter : number of iterations
# 
# [D,N] = size(X);
# 
# W = 1 + rand(D,K);
# H = 1 + rand(K,N);
# 
# temp = sum(X(:))/sum(sum(W*H));
# W = W*sqrt(temp);
# H = H*sqrt(temp);
# 
# L = [];
# for iter = 1:num_iter
# if strcmp(penalty,'L2')
# H = H.*(W'*X)./((W'*W)*H + eps);
# W = W.*(X*H')./(W*(H*H') + eps);
# L(iter) = sum(sum((X - W*H).^2));
# elseif strcmp(penalty,'divergence')
# H = H.*((W./repmat(sum(W,1),D,1))'*(X./(W*H + eps)));
#         W = W.*((X./(W*H + eps))*(H./repmat(sum(H,2),1,N))');
# L(iter) = -sum(sum(X.*log(W*H))) + sum(sum(W*H));
# else
#   disp('Invalid penalty input');
# break;
# end
# end

# % NMF implements non-negative matrix factorization of the matrix X
# % according to a predefined penalty.
# %
# % X : D x N matrix of nonnegative numbers
# % K : rank of factorization, takes value in {1,2,3,...,min(D,N)}
# % penalty : string taking one of two values, 'L2' or 'divergence'
# % num_iter : number of iterations

X <- chem_n
K <- 5
num_iter <- 10

nmf_ee <- function(X, K, penalty, num_iter) {

D <- nrow(X)
N <- ncol(X)

# Initiate W and H matrices with random non-zero values
W <- matrix(1, nrow = D, ncol = K) + matrix(runif(D*K), nrow = D, ncol = K)
H <- matrix(1, nrow = K, ncol = N) + matrix(runif(K*N), nrow = K, ncol = N)

temp <- sum(X)/sum(W %*% H)
W <- W*sqrt(temp)
H <- H*sqrt(temp)

L <- matrix()

# Define 'eps'
# eps returns the distance from 1.0 to the next larger double-precision number, that is, 2-52.
eps <- 2.2204e-16

# Iterate to minimize loss function
for (i in 1:num_iter) {

  if (identical(penalty,'L2')) {
        H <- H * (t(W) %*% X) / ((t(W) %*% W) %*% H + eps)
        W <- W * (X %*% t(H)) / (W %*% (H %*% t(H)) + eps)
        
  # Loss function -- Gaussian MLE
  L(i) <- sum((X - W %*% H)^2)
  }

  else if (identical(penalty,'divergence')) {
        H <- H * (t(W/kronecker(matrix(1,D,1), matrix(colSums(W), nrow=1))) %*% (X / (W %*% H + eps)))
        W <- W * ((X / (W %*% H + eps)) %*% t(H / kronecker(matrix(1,1,N), matrix(rowSums(H), ncol = 1))))

  # Loss function -- Poisson MLE
  L(i) <- -sum(X * log(W %*% H)) + sum(W %*% H)
  }
  
  else {print('Invalid penalty input')}

  }
}