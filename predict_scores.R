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

pred_x = Pinv(A) %*% y

# MSE
mean((x - pred_x)^2)

####
####
####

# NMF pred does not work on new data

library(NMF)

# fit an NMF model
nn <- nmf(y, 5)

dim(basis(nn))
# basis are scores
dim(coef(nn))
#coef are loadings

# predicted column and row clusters
predict(nn)
# this gives the column clusters
predict(x, 'rows')
# this gives the row clusters

pred_basis <- y %*% t(coef(nn))
pred_basis <- y %*% Pinv(coef(nn))

mean((basis(nn) - pred_basis)^2)
cor(basis(nn), pred_basis)

####
####
####



####
####
####

# FA predict works

x <- sim.item(12,50)
xx <- sim.item(12,50)
dim(xx)

f <- fa(x,2,scores="regression")  # a two factor solution
#find the predicted scores (The B set)
# Standardize by old data
p <- predict(f, old.data = x,  data = xx) 
dim(p)

#test how well these predicted scores match the factor scores from the second set
f2 <- fa(xx,2,scores="regression")
round(cor(f2$scores,p),2)

scores_xx <- xx %*% matrix(f2$loadings, ncol = 2)
round(cor(f2$scores, scores_xx),2)

mean((p - f2$scores)^2)
mean((scores_xx - f2$scores)^2)

