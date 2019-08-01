function [W,H,L] = NMF(X,K,penalty,num_iter)
% NMF implements non-negative matrix factorization of the matrix X
% according to a predefined penalty.
%
% X : D x N matrix of nonnegative numbers
% K : rank of factorization, takes value in {1,2,3,...,min(D,N)}
% penalty : string taking one of two values, 'L2' or 'divergence'
% num_iter : number of iterations

[D,N] = size(X);

W = 1 + rand(D,K);
H = 1 + rand(K,N);

temp = sum(X(:))/sum(sum(W*H));
W = W*sqrt(temp);
H = H*sqrt(temp);

L = [];
for iter = 1:num_iter
    if strcmp(penalty,'L2')
        H = H.*(W'*X)./((W'*W)*H + eps);
        W = W.*(X*H')./(W*(H*H') + eps);
        L(iter) = sum(sum((X - W*H).^2));
    elseif strcmp(penalty,'divergence')
        H = H.*((W./repmat(sum(W,1),D,1))'*(X./(W*H + eps)));
        W = W.*((X./(W*H + eps))*(H./repmat(sum(H,2),1,N))');
        L(iter) = -sum(sum(X.*log(W*H))) + sum(sum(W*H));
    else
        disp('Invalid penalty input');
        break;
    end
end
