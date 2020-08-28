library(R.matlab)

##### 
##### CI WA

varwa <- readMat(here::here("./MATLAB/Output/sim_over_100_varwa.mat"))[[1]]
ewa <- readMat(here::here("./MATLAB/Output/sim_over_100_ewa.mat"))[[1]]
head(ewa)
head(varwa)

wa_low <- ewa - 1.96*varwa
wa_up <- ewa + 1.96*varwa
head(wa_low)
head(wa_up)

# Truth
load("./R/Sims/Iterate/sim_over_old.RDA")
true_score <- sim_over[100,2][[1]]
head(true_score)

sum(true_score > wa_low & true_score < wa_up)/nrow(true_score)*ncol(true_score)

#####
##### CI H

alphaH <- readMat(here::here("./MATLAB/Output/sim_over_100_alphah.mat"))[[1]]
betaH <- readMat(here::here("./MATLAB/Output/sim_over_100_betah.mat"))[[1]]
eh <- readMat(here::here("./MATLAB/Output/sim_over_100_eh.mat"))[[1]]

# normalize loadings
row_sum_h <- apply(eh, 1, sum)
eh_norm <- apply(eh, 2, function(x) x / row_sum_h)
head(eh_norm)
rowSums(eh_norm)

# create CI
h_low <- qgamma(0.025, shape = alphaH, rate =betaH)
h_up <- qgamma(0.975, shape = alphaH, rate =betaH)

h_low_norm <- apply(h_low, 2, function(x) x / row_sum_h)
h_up_norm <- apply(h_up, 2, function(x) x / row_sum_h)

# Truth
true_load <- sim_over[100,3][[1]]
head(true_load)[, 1:10]
# normalize loadings
row_sum <- apply(true_load, 1, sum)
true_load_norm <- apply(true_load, 2, function(x) x / row_sum)
head(true_load_norm)
rowSums(true_load_norm)

sum(true_load_norm > low_norm & true_load_norm < up_norm)
sum(true_load > h_low & true_load < h_up)

sum(eh > h_low & eh < h_up)

head(h_low)[,1:5]
head(true_load)[,1:5]
head(eh)[,1:5]
head(h_up)[,1:5]

head(h_low_norm)[,1:5]
head(true_load_norm)[,1:5]
head(eh_norm)[,1:5]
head(h_up_norm)[,1:5]

rowSums(true_load_norm)
rowSums(eh_norm)
rowSums(h_low_norm)
rowSums(h_up_norm)
##########
########## Gamma Dist


#####
alphaw = 2
betaw = 3
gamw <- rgamma(100000, shape = alphaw, rate = betaw)
mean(gamw)
# 2/3
var(gamw)
# 2/3^2

# var term 1 = a(a+1)/b^2
# var term 2 = a^2/b^2
# gamma var = term 1 - term 2

# e(w^2) * e(a^2) - e(w)^2 * e(a)^2

var_1 = function (alpha, beta) {(alpha*(alpha+1))/(beta^2)}
var_2 = function (alpha, beta) {(alpha^2)/(beta^2)}
  
var_1(alphaw, betaw) - var_2(alphaw, betaw)

alphaa = 3
betaa = 4

gama <- rgamma(100000, shape = alphaa, rate = betaa)
mean(gama)
# 3/4
var(gama)
# 3/4^2

var_1(alphaa, betaa) - var_2(alphaa, betaa)

product <- gamw*gama
mean(product)
var(product)

var_1(alphaw, betaw)*var_1(alphaa, betaa) - var_2(alphaw, betaw)*var_2(alphaa, betaa)

newbeta = mean(product)/var(product)
newtheta = var(product)/mean(product)
newalpha = mean(product)*newbeta
newbeta
newalpha

gam3 <- rgamma(100000, shape = newalpha, rate = newbeta)
gam4 <- rgamma(100000, shape = newalpha, scale = newtheta)

mean(gam3)
mean(gam4)
# 1/2
var(gam3)
var(gam4)
# 1/2^2
var_1(newalpha, newbeta) - var_2(newalpha, newbeta)

# varWA = (W1 .* A1) .* (W1 + A1 + 1) ./ (W2.^2 .* A2.^2);

