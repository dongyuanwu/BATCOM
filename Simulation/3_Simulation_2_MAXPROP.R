
simdat <- readRDS("simdat1.rds")   # or "simdat2.rds": data from pseudo hurdle gamma model

library(Rcpp)
library(RcppArmadillo)
library(RcppGSL)
library(HDInterval)
library(coda)
library(Matrix)

sourceCpp("MainCode.cpp")
sourceCpp("get_Ms.cpp")

get_measure <- function(records) {
    
    mainres <- c(summary(records), sd(records), 
                 quantile(records, 0.025), quantile(records, 0.975), 
                 quantile(records, 0.005), quantile(records, 0.995), 
                 quantile(records, 0.01), quantile(records, 0.99),
                 hdi(records, 0.95), hdi(records, 0.98), hdi(records, 0.99))
    
    zval <- geweke.diag(mcmc(records), frac1=0.49, frac2=0.49)$z
    pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
    
    return(c(mainres, zval, pval))
    
}

get_DIC <- function(paramstemp, TT, lp) {
    
    beta_cov <- c(10000, rep(1, ncol(X)-1))
    lp_est <- sum(loglik_ind(colMeans(paramstemp), X, Y, TT, beta_cov, reL, reR))
    
    DIC <- 2*lp_est - 4*mean(lp)
    
    return(DIC)
    
}

finalsim <- vector(mode="list", length=length(simdat))
nsim <- 100
nspot <- 100
burnin <- 1:10000

generate_panel_posi <- function(row=10, col=10) {
    
    columnid <- rep(1:col, times=row)
    rowid <- rep(1:row, each=col)
    return(data.frame(spotid=1:(row*col), spotX=columnid, spotY=rowid))
    
}
spotlist <- generate_panel_posi(row=10, col=10)

spot_orig <- cbind(rep(1:nspot, each=nspot), rep(1:nspot, times=nspot))
dist_orig <- sqrt((spotlist$spotX[spot_orig[, 1]] - spotlist$spotX[spot_orig[, 2]])^2 + 
                      (spotlist$spotY[spot_orig[, 1]] - spotlist$spotY[spot_orig[, 2]])^2)

n <- nspot*nspot
phi <- 1
p <- 1.5
logphi <- log(phi)
theta <- log((p-1)/(2-p))

adj_M_row <- function(M_row) {
    
    M_row <- ifelse(M_row == max(M_row), 1, 0)
    return (M_row)
    
}

for (nset in 1:length(simdat)) {
    
    simres <- vector("list", length=nsim)
    simset <- simdat[[nset]]$simset
    
    
    ncelltype <- simset$ncelltype
    q <- ncelltype*ncelltype+nspot*2+1
    
    for(k in 1:nsim) {
        
        print(paste0("rho=", simset$rho, "; phi=", simset$phi, "; p=", simset$p, "; sprate=", simset$sprate, "; iter=", k))
        
        M <- simdat[[nset]]$simdat[[k]]$M
        # Treat the cell type with highest proportion as 1, others as 0.
        M <- t(apply(M, 1, adj_M_row))
        rmposi <- simdat[[nset]]$simdat[[k]]$rmposi
        
        if(length(rmposi) > 0) {
            spot <- spot_orig[-rmposi, ]
            dist <- dist_orig[-rmposi]
        } else {
            spot <- spot_orig
            dist <- dist_orig
        }
        
        Ms <- get_Ms(M, spot)
        
        nonzero_posi <- which(colSums(Ms) != 0)
        Ms <- Ms[, nonzero_posi]
        
        # can adjust 0.5 to other values for the tuning parameter rho-hat
        X <- scale(Ms * exp(-0.5*dist), center=TRUE, scale=TRUE)
        X <- cbind(rep(1, nrow(X)), X)
        
        Y <- simdat[[nset]]$simdat[[k]]$Y
        
        newX1 <- matrix(0, ncol=nspot*2, nrow=nrow(spot))
        for (i in 1:nspot) {
            newX1[spot[, 1] == i, i] <- 1
            newX1[spot[, 2] == i, i+nspot] <- 1
        }
        newX <- as(cbind(X, newX1), "sparseMatrix")
        
        beta.cov <- c(10000, rep(0.01, ncol(X)-1), rep(1, nspot*2))
        
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
        
        # Parameter transformation

        phistore <- exp(paramstemp[, 1])
        paramstemp[, 2] <- ifelse(paramstemp[, 2] > 10, 10, paramstemp[, 2])
        paramstemp[, 2] <- ifelse(paramstemp[, 2] < -10, -10, paramstemp[, 2])
        pstore <- (2*exp(paramstemp[, 2])+1)/(exp(paramstemp[, 2])+1)
        
        # Descriptive statistics for our inference samples of parameters.

        phires <- get_measure(phistore)
        pres <- get_measure(pstore)

        betares <- t(apply(paramstemp[, -c(1:2)], 2, function(x) get_measure(x)))

        simres[[k]] <- list(phires=phires, pres=pres, betares=betares, TT=TTtemp,
                            init=params, initflag=flag, hmciter=iter, times=times,
                            arate=length(unique(paramstemp[, 1]))/10000,
                            arate1=length(unique(paramstemp[1:3333, 1]))/3333,
                            arate2=length(unique(paramstemp[3334:6666, 1]))/3333,
                            arate3=length(unique(paramstemp[6667:10000, 1]))/3334,
                            DIC=get_DIC(paramstemp, TTtemp, hmcres$lp),
                            WAIC=as.vector(hmcres$WAIC),
                            nonzero_posi=nonzero_posi)
        
    }
    
    finalsim[[nset]] <- list(simset=simset, simres=simres)
    
}

##################################################################

# save results
saveRDS(finalsim, "simres1_maxprop.rds")
