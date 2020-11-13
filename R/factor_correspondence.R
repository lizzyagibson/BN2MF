# %       A, B matrices of the same size, with unit L2 norm columns
# %       nn   boolean, optional, default true. If nn, we search for a
# %       *permutation* matrix Pi. If nn is false, we allow *signed
# %       permutations* which can multiply the factors by -1. 
# %
# %   Outputs
# %       e -- the sum of squared errors under the best calibration \Pi
# %       Pi -- the permutation or signed permutation matrix
# %
# %    June 2020 John Wright, jw2966@columbia.edu

library(CVXR)
# library(tidyverse)

# A <- matrix(rnorm(50), nrow = 10)
# Alength <- apply(A, 2, function(y) sqrt(sum(y^2)))
# Al2 <- matrix(NA, nrow = nrow(A), ncol = ncol(A))
# 
# for (i in 1:ncol(A)) {
#   Al2[,i] <- A[,i]/Alength[i]
# }
# 
# Bl2 <- matrix(NA, nrow = nrow(A), ncol = ncol(A))
# Bl2[,1:5] <- Al2[,sample(1:5)]

###

factor_correspondence <- function (A, B, nn = TRUE) {

    G <- t(B) %*% A
    n <- nrow(G)

    # Step 1. Define the variable to be estimated
    # Pi -- the permutation or signed permutation matrix
    Pi <- Variable(n,n)
    
    # Step 2. Define the objective to be optimized
    objective <- Maximize(base::sum(Pi * G))
    
    if (nn) {
          # Step 2.5. Subject to these constraints
          constX = list()
          for (i in 1:nrow(G)) {
            constX <-  c(constX, list(base::sum(Pi[,i]) == 1))
            constX <-  c(constX, list(base::sum(Pi[i,]) == 1))
           }
          constX <- c(constX, Pi >= 0)
    } else {
          # % allow sign flips 
          # Step 3. vector l1 norms along rows and columns
          constX = list()
          for (i in 1:nrow(G)) {
            constX <-  c(constX, list(base::sum(abs(Pi[,i])) <= 1))
            constX <-  c(constX, list(base::sum(abs(Pi[i,])) <= 1))
          }
    }  
        # Step 3. Create a problem to solve
        problem <- Problem(objective, constraints = constX)
        
        # Step 4. Solve it!
        result <- solve(problem)
        
        # Step 5. Extract solution and objective value
        perm <- round(result$getValue(Pi), 0)
        
        e <- norm(B,'f')^2 + norm(A,'f')^2 - 2 * base::sum(perm * G)
        # e -- the sum of squared errors under the best calibration \Pi
        
        # New matrix with best order
        newB <- B %*% perm
        
        return(list(rearranged = newB, permutation_matrix = perm))
    }

# Bout <- factor_correspondence(Al2, Bl2)$rearranged
# factor_correspondence(Al2, Bl2)
# 
# ## Viz 
# Al2 %>% as_tibble() %>% 
#   pivot_longer(V1:V5) %>% mutate(matrix = "Al2", ind = 1:nrow(.)) %>% 
#   rbind(.,
#         Bl2 %>% as_tibble() %>% 
#           pivot_longer(V1:V5) %>% mutate(matrix = "Bl2", ind = 1:nrow(.))) %>% 
#   rbind(.,
#         Bout %>% as_tibble() %>% 
#           pivot_longer(V1:V5) %>% mutate(matrix = "Bl2*perm", ind = 1:nrow(.))) %>% 
#   ggplot(aes(x = value)) +
#   geom_histogram() + 
#   facet_grid(name~matrix)
