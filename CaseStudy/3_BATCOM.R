
dat <- readRDS("RCTD_HMC_dist3_estrho0.5_dat.rds")

### overall

# temp <- matrix(unlist(dat$Y), ncol=length(dat$Y))
# Y <- rowSums(temp)

###

X <- dat$X
spot <- dat$spot
celltype <- dat$celltype
celltype_interaction <- dat$celltype_interaction

Cspots <- dat$Y
ligands <- dat$ligands
receptors <- dat$receptors

rm(dat)

library(Rcpp)
library(RcppArmadillo)
library(RcppGSL)
library(Matrix)

sourceCpp("MainCode.cpp")

nspot <- length(unique(spot[, 1]))
ncelltype <- length(celltype)
burnin <- 1:10000

n <- nrow(X)
q <- ncol(X) + nspot*2
beta.cov <- c(10000, rep(0.01, ncol(X)-1), rep(1, nspot*2))

phi <- 1
p <- 1.5
logphi <- log(phi)
theta <- log((p-1)/(2-p))

# Construct the 0-1 matrix for indicating random effects

newX1 <- matrix(0, ncol=nspot*2, nrow=n)
for (i in 1:nspot) {
    newX1[spot[, 1] == i, i] <- 1
    newX1[spot[, 2] == i, i+nspot] <- 1
}

newX <- as(cbind(X, newX1), "sparseMatrix")

reslist <- vector(mode="list", length=length(Cspots))

for (k in 1:length(Cspots)) {
    
    print(paste0("LRPair: ", k))
    
    Y <- Cspots[[k]]
    
    # Estimate initial values of betas, v^L, and v^R with fixed phi=1 and p=1.5
    begin <- Sys.time()
    beta <- as.vector(newton_raphson_beta_re_sp(newX, Y, beta.cov, logphi, p, nspot, 
                                             spot[, 1]-1, spot[, 2]-1, 0.000001, 100))
    times <- difftime(Sys.time(), begin, units="secs")
    params <- c(logphi, theta, beta)
    
    flag <- 1
    
    # Estimate initial values of betas, v^L, v^R, T, and logphi, and theta(=log((p-1)/(2-p)))
    
    begin <- Sys.time()
    nr_res <- tryCatch(newton_raphson_re_sp(params, newX, Y, beta.cov, nspot, 
                                         spot[, 1]-1, spot[, 2]-1, 0.000001, 100, 1),
                       error=function(err) list(params_new=params, TT=ifelse(Y > 0, 1, 0)))
    params1 <- as.vector(nr_res$params_new)
    # In some extreme cases, this method won't converge. Then we can change the learning rate from 1 to 0.1 and try again
    if(all(params1[-c(1:3)] == 0) | all(params1 == params) |
       abs(params1[1]) > 10 | abs(params1[2]) > 10) {
        nr_res <- tryCatch(newton_raphson_re_sp(params, newX, Y, beta.cov, nspot,
                                             spot[, 1]-1, spot[, 2]-1, 0.000001, 100, 0.1),
                           error=function(err) list(params_new=params, TT=ifelse(Y > 0, 1, 0)))
        params1 <- as.vector(nr_res$params_new)
        flag <- 2
    }
    # If it still does not converge, then we just use the estimated values of betas and arbitrary values of T, phi, and p
    if(all(params1[-c(1:3)] == 0) | all(params1 == params) |
       abs(params1[1]) > 10 | abs(params1[2]) > 10) {
        flag <- 3
        TT <- ifelse(Y > 0, 1, 0)
    }
    if(flag != 3) {
        params <- params1
        TT <- as.vector(nr_res$TT)
    }
    times <- c(times, difftime(Sys.time(), begin, units="secs"))
    
    re <- params[-c(1:(ncol(X)+2))]
    reL <- re[1:(length(re)/2)][spot[, 1]]
    reR <- re[(length(re)/2+1):length(re)][spot[, 2]]
    params <- params[1:(ncol(X)+2)]
    
    # main MCMC algorithm
    
    iter <- 1
    hmcflag <- TRUE
    begin <- Sys.time()
    # The MCMC algorithm sometimes may not provide any acceptable samples, but it just depends on random seeds.
    # So we can try it again and again until we get the results.
    while( hmcflag & (iter <= 10) ) {
        hmcres <- hmc(Y, X, params, TT, reL, reR, 20000)
        paramstemp <- hmcres$params_store[-burnin, ]
        TTtemp <- hmcres$TT
        
        # If the acceptance rate is very low, we will run hmc function again.
        if( sum(paramstemp[, 1] == 0) < 1000 ) {
            hmcflag <- FALSE
        } else {
            iter <- iter + 1
        }
    }
    times <- c(times, difftime(Sys.time(), begin, units="secs"))
    print(difftime(Sys.time(), begin, units="hours"))
    
    reslist[[k]] <- list(params.store=paramstemp, ligand=ligands[[k]], receptor=receptors[[k]], times=times)
    
}

##################################################################

# save results
saveRDS(reslist, "RCTD_HMC_dist3_estrho0.5_res.rds")
