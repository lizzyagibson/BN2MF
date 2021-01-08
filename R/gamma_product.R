library(Bessel)
library(tidyverse)
library(gsl)

#the pdf of the product z = x*y
log_gam_gam_pdf <- function(alpha1, theta1, alpha2, theta2, z) {

  # (2 * (((z/(theta1 * theta2))^((alpha1 + alpha2)/2)) * 
  #  besselK(2*sqrt(z/(theta1*theta2)), alpha2 - alpha1))) /
  #  (z * gamma(alpha1) * gamma(alpha2))
  
  ((log(2) + ((log(z) - (log(theta1) + log(theta2))) * ((alpha1+ alpha2)/2))) + 
    log(besselK(2*sqrt(z/(theta1*theta2)), alpha2-alpha1))) -
    (log(z) + log(gamma(alpha1)) + log(gamma(alpha2)))
}

# gamma product pdf
gam_gam_pdf <- function(alpha1, theta1, alpha2, theta2, z) {
  exp(log_gam_gam_pdf(alpha1,theta1,alpha2,theta2,z))
}

# gamma pdf
gam_pdf <- function(k, theta, x) {(x^(k-1)*exp(-x/theta))/(gamma(k)*theta^k)}

# gamma product cdf
gam_gam_cdf <- function(y, alpha1, theta1, alpha2, theta2) {
  integrate(function(x) gam_gam_pdf(alpha1, theta1, alpha2, theta2, x), 
            lower=0, upper=y, subdivisions = 100000)[[1]]
}

# against theoretical solution derived above when 𝑎=2, 𝑏=3, 𝛼=4, and 𝛽=1.1
x = rgamma(100000, shape = 2, scale = 3)
y = rgamma(100000, shape = 4, scale = 1.1)
z = x * y

# e[x] is x * density
mean = integrate(function(x) x * gam_gam_pdf(2, 3, 4, 1.1, x), lower=0, upper=1000, subdivisions = 1000)
mean[[1]]
2*3*4*1.1
var = integrate(function(x){ (x- mean[[1]])^2 * gam_gam_pdf(2, 3, 4, 1.1, x)}
                , lower=0, upper=1000, subdivisions = 1000)
var[[1]]
# varWA = ((W1 .* A1) .* (W1 + A1 + 1)) ./ (W2.^2 .* A2.^2);
((2*4) * (2+4+1))/((1/3)^2 * (1/1.1)^2)

aa = 100
ta = 1/100
ab = 100
tb = 1/100
a = rgamma(100000, shape = aa, scale = ta)
b = rgamma(100000, shape = ab, scale = tb)
c = a * b

# Plug this into functions
theCDF <- function(x) {gam_gam_cdf(x, aa, ta, ab, tb)}

gam_gam_cdf(10, aa, ta, ab, tb)

find_upper <- function(thecdf){
  old_x = 1
  new_x = 2 
  old_prob = thecdf(old_x)
  new_prob = thecdf(new_x)
  
  while (new_prob >= old_prob & new_prob < 1 ) {
    old_prob = new_prob
    old_x= new_x 
    new_x = 2 * new_x 
    new_prob = thecdf(new_x)
  }
  return(old_x)
}

find_upper(function(x) gam_gam_cdf(x, 2, 3, 4, 1.1))
find_upper(function(x) gam_gam_cdf(x, aa, ta, ab, tb))

estTile <- function(prob,
                    theCDF, 
                    upper=find_upper(theCDF), 
                    lower=0.0001, tol=0.000001, 
                    limit= 10000){

  mid = (upper - lower)/2
  midProb <- theCDF(mid)
  step <- 1
  
  while(abs(midProb-prob)> tol && step < limit) {
    if (midProb > prob){
      upper = mid
    } else {
      lower = mid 
    }
  
    step = step + 1 
    mid <- (upper + lower) /2 
    midProb <- theCDF(mid)
    
   print(paste0("step ", step))
   print(midProb)
   print(mid)
  }
  return(mid) 
  }

gam_gam_cdf(1, aa, ta, ab, tb)
quantile(c, c(0.025, 0.975))
estTile(0.025,
        function(x) gam_gam_cdf(x, aa, ta, ab, tb))

mycdfGAM = function(x) gam_gam_cdf(x, aa, ta, ab, tb) 
find_upper(mycdfGAM)

#upper = 20* max(qgamma(0.99,aa,ta),qgamma(0.999,ab,tb)))

estTile(0.025, function(x) gam_gam_cdf(x, 2, 3, 4, 1.1))
estTile(0.975, function(x) gam_gam_cdf(x, 2, 3, 4, 1.1))

quantile(z, c(0.025, 0.975)) # Pretty close!
quantile(x, c(0.025, 0.975)) * quantile(y, c(0.025, 0.975)) # Nope!

getVCI <- function(alpha1, theta1, alpha2, theta2) {
  lower = estTile(0.025, function(x) gam_gam_cdf(x, alpha1, theta1, alpha2, theta2))
  upper = estTile(0.975, function(x) gam_gam_cdf(x, alpha1, theta1, alpha2, theta2))
  return(list(lower = lower, uppper = upper))
}

# Visualize
par(mfrow=c(2,4))
plot(x, gam_pdf(2, 3, x))
plot(y, gam_pdf(4, 1.1, y))
plot(z, gam_gam_pdf(2, 3, 4, 1.1, z))
grid =  seq(0.001, 100, by = .1)
plot(grid,gam_gam_pdf(2, 3.01, 4, 1.1, grid))
hist(x)
hist(y)
hist(z)