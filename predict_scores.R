# y = A*x

# y = p chemicals x n individuals
# A = p x k
# x = k x n

# matrix multiply
A = matrix(exp(rnorm(50)), ncol = 5)
dim(A)

x = matrix(exp(rnorm(100)), nrow = 5)
dim(x)

y = A %*% x
dim(y)

# MSE
mean((x - (Pinv(A) %*% y))^2)
mean((x - (t(A) %*% y))^2)

####
####
####

# NMF pred does not work on new data

library(NMF)
y <- rmatrix(20, 10)
# fit an NMF model
n_out <- nmf(y, 5)

dim(basis(n_out))
# basis are scores
dim(coef(n_out))
#coef are loadings

pred_basis_t <- y %*% t(coef(n_out))
pred_basis_i <- y %*% Pinv(coef(n_out))

mean((basis(n_out) - pred_basis_t)^2)
mean((basis(n_out) - pred_basis_i)^2)

####
####
####

y <- rmatrix(20, 10)
y_scaled <- scale(y)
y_scaled_F <- scale(y, center = TRUE, scale = FALSE)
dim(y)

pca_y <- prcomp(y)

# predict scores
pred_sF <- y_scaled_F %*% pca_y$rotation

# y centered
mean((pca_y$x - pred_sF)^2)

####
####
####

# FA predict works

x <- sim.item(12,50)
dim(x)

f <- fa(x, 2, scores="regression")

# #find the predicted scores
# # Standardize by old data
# p <- predict(f, old.data = x,  data = xx) 
# dim(p)

#test how well these predicted scores match the factor scores from the second set
# f2 <- fa(xx,2,scores="regression")
# round(cor(f2$scores,p),2)

scores_x <- x %*% matrix(f$loadings, ncol = 2)
round(cor(f$scores, scores_x),2)

mean((p - f2$scores)^2)
mean((scores_x - f2$scores)^2)

