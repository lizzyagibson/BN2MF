library(parallel)
n <- 300 # observations
p <- 20000 # covariates
# > ## Different sized matrices as filter inputs
#   > ## Matrix A and B form smaller work loads
#   > ## while matrix C forms a bigger workload (2*p)
A <- matrix(replicate( p, rnorm(n, sd = runif(1, 0.1, 10))), n, p)
B <- matrix(replicate( p, rnorm(n, sd = runif(1, 0.1, 10))), n, p)
C <- matrix(replicate(2*p, rnorm(n, sd = runif(1, 0.1, 10))), n, 2*p)

varFilter <- function (X, nSim = 20) {
   for (i in 1:nSim) {
    train <- sample(nrow(X), 2 / 3 * nrow(X))
    colVars <- apply(X[train, ], 2, var)
    keep <- names(head(sort(colVars, decreasing = TRUE), 100))
    # myAlgorithm(X[, keep])
       }
  }

## Runtime comparison -----------------------------------
  # > ## mclapply with affinity.list
  # > ## CPU mapping: A and B run on CPU 1 while C runs on CPU 2: 
affinity <- c(1,1,2)

system.time(
    mclapply(X = list(A,B,C), FUN = varFilter,
               mc.preschedule = FALSE, affinity.list = affinity))
## elapsed
## 32.695 
  
## mclapply without affinity.list
system.time(
      mclapply(X = list(A,B,C), FUN = varFilter, mc.cores = 2,
                 mc.preschedule = FALSE) )
## elapsed
## 42.442 

## mclapply with prescheduling
system.time(
      mclapply(X = list(A,B,C), FUN = varFilter, mc.cores = 2,
                 mc.preschedule = TRUE) )
## elapsed
## 37.970 

system.time(
  lapply(X = list(A,B,C), FUN = varFilter))
## elapsed
## 48.023

#########
library(future)

plan(multiprocess)
pid <- Sys.getpid()
pid

a %<-% {
       pid <- Sys.getpid()
       cat("Future 'a' ...\n")
       3.14
   }

b %<-% {
       rm(pid)
       cat("Future 'b' ...\n")
       Sys.getpid()
   }

c  %<-% {
       cat("Future 'c' ...\n")
       2 * a
   }
# Running c prints a's cat bc it includes a

b

c

a

pid


plan(multiprocess)
system.time(demo("mandelbrot", package = "future", ask = FALSE))

library("graphics")

plot_what_is_done <- function(counts) {
      for (kk in seq_along(counts)) {
        f <- counts[[kk]]

         ## Already plotted?
         if (!inherits(f, "Future")) next

         ## Not resolved?
         if (!resolved(f)) next

         message(sprintf("Plotting tile #%d of %d ...", kk, n))
       counts[[kk]] <- value(f)
        screen(kk)
        plot(counts[[kk]])
      }

      counts
   }

# > ## Options
region <- getOption("future.demo.mandelbrot.region", 1L)
 
if (!is.list(region)) {
     if (region == 1L) {
         region <- list(xmid = -0.75, ymid = 0.0, side = 3.0)
      } else if (region == 2L) {
           region <- list(xmid = 0.283, ymid = -0.0095, side = 0.00026)
         } else if (region == 3L) {
             region <- list(xmid = 0.282989, ymid = -0.01, side = 3e-8)
           }
  }

nrow <- getOption("future.demo.mandelbrot.nrow", 3L)
resolution <- getOption("future.demo.mandelbrot.resolution", 400L)
delay <- getOption("future.demo.mandelbrot.delay", interactive())

if (isTRUE(delay)) {
     delay <- function(counts) Sys.sleep(rexp(1, rate = 2))
   } else if (!is.function(delay)) {
       delay <- function(counts) {}
     }
 
# > ## Generate Mandelbrot tiles to be computed
Cs <- mandelbrot_tiles(xmid = region$xmid, ymid = region$ymid,
                      side = region$side, nrow = nrow,
                      resolution = resolution)

if (interactive()) {
     dev.new()
     plot.new()
     split.screen(dim(Cs))
     for (ii in seq_along(Cs)) {
         screen(ii)
         par(mar = c(0, 0, 0, 0))
         text(x = 1 / 2, y = 1 / 2, sprintf("Future #%d\nunresolved", ii), cex = 2)
       }
   } else {
      split.screen(dim(Cs))
    }

# > ## Create all Mandelbrot tiles via lazy futures
n <- length(Cs)
 
message(sprintf("Creating %d Mandelbrot tiles:", n), appendLF = FALSE)

# Creating 9 Mandelbrot tiles:
counts <- lapply(seq_along(Cs), FUN=function(ii) {
       message(" ", ii, appendLF = FALSE)
       C <- Cs[[ii]]
       future({
           message(sprintf("Calculating tile #%d of %d ...", ii, n), appendLF = FALSE)
           fit <- mandelbrot(C)
       
            ## Emulate slowness
            delay(fit)
      
          message(" done")
           fit
       }, lazy = TRUE)
     })

message(".")
# > ## Calculate and plot tiles
repeat {
        counts <- plot_what_is_done(counts)
        if (!any(sapply(counts, FUN = inherits, "Future"))) break
}

###############

# 1 .... 100
# [1 ... 25 ], [26 : 50]

nruns <- list(1:25, 26:50, 51:75, 76:100)


f <- function(listofruns){

    for (i in listofruns){ # list of runs is 1:25
      
    
  }
  
}
##

for (i in 1:4) {
  res[i]= future(f(i))
}

agg <- c(1:10)
for (i in 1:10) {
  agg[i]<- value(res[i])
}

normal stuff

