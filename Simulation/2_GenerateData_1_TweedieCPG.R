
library(Rcpp)
library(tweedie)

sourceCpp("get_Ms.cpp")
betaset <- readRDS("simdat_beta.rds")

nspot <- 100

# Generate a spot panel with 10*10 grid

generate_panel_posi <- function(row=10, col=10) {
    
    columnid <- rep(1:col, times=row)
    rowid <- rep(1:row, each=col)
    return(data.frame(spotid=1:(row*col), spotX=columnid, spotY=rowid))
    
}
spotlist <- generate_panel_posi(row=10, col=10)

spot <- cbind(rep(1:nspot, each=nspot), rep(1:nspot, times=nspot))
dist <- sqrt((spotlist$spotX[spot[, 1]] - spotlist$spotX[spot[, 2]])^2 + 
                 (spotlist$spotY[spot[, 1]] - spotlist$spotY[spot[, 2]])^2)

nsim <- 100

# Simulation parameters

ncelltypes <- c(5, 10, 15)
rhos <- c(0.2, 0.4, 0.6, 0.8)
phis <- c(0.6, 0.8, 1, 3, 5, 10)
ps <- c(1.1, 1.3, 1.5, 1.7, 1.9)
sprates <- c(0.2, 0.4, 0.6, 0.8)

finalsim <- vector(mode="list", length=length(ncelltypes)*length(rhos)*length(phis)*length(ps)*length(sprates))
ncomb <- 1


for (ncelltype in ncelltypes) {
    
    for (rho in rhos) {
        
        for (phi in phis) {
            
            for (p in ps) {
                
                if ((p == 1.9) & (phi <= 3)) next
                
                for (sprate in sprates) {
                    
                    setposi <- which(sapply(betaset, function(x) (x$sprate == sprate) & (x$ncelltype == ncelltype)))
                    betaset1 <- betaset[[setposi]]
                    truebeta <- betaset1$truebeta
                    truebeta[1] <- truebeta[1] - 2
                    reL <- betaset1$reL
                    reR <- betaset1$reR
                    
                    simset <- list(rho=rho, phi=phi, p=p, sprate=sprate, ncelltype=ncelltype, 
                                   truebeta=truebeta, reL=reL, reR=reR)
                    simdat <- vector("list", length=nsim)
                    
                    for(k in 1:nsim) {
                        
                        print(paste0("rho=", rho, "; phi=", phi, "; p=", p, "; sprate=", sprate, "; iter=", k))
                        
                        # Generate cell type proportion matrix
                        
                        M <- NULL
                        for (i in 1:nspot) {
                            M <- rbind(M, prop.table(runif(ncelltype)))
                        }
                        
                        # Calculate M*M
                        
                        Ms <- get_Ms(M, spot)
                        
                        X <- scale(Ms * exp(-rho*dist), center=TRUE, scale=TRUE)
                        X <- cbind(rep(1, nrow(X)), X)
                        
                        # eta = X * beta + v_i^L + v_j^R
                        
                        eta <- c(X %*% truebeta + reL[spot[, 1]] + reR[spot[, 2]])
                        
                        # Remove some extreme values of eta
                        
                        threshold <- max(7, mean(eta) + 3 * sd(eta))
                        rmposi <- which(eta > threshold)
                        
                        if (length(rmposi) > 0) {
                            eta <- eta[-rmposi]
                        }
                        
                        eta[eta > 10] <- 10
                        mu <- c(exp(eta))
                        
                        #########
                        
                        # Generate Y using Tweedie distribution
                        
                        n <- length(mu)
                        Y <- array(dim = n, NA)
                        lambda <- mu^(2 - p)/(phi * (2 - p))
                        alpha <- (2 - p)/(1 - p)
                        gam <- phi * (p - 1) * mu^(p - 1)
                        N <- rpois(n, lambda = lambda)
                        for (i in (1:n)) {
                            Y[i] <- rgamma(1, shape = -N[i] * alpha, scale = gam[i])
                        }
                        
                        #########
                        
                        simdat[[k]] <- list(Y=as.vector(Y), M=M, N=N, rmposi=rmposi)
                        
                    }
                    
                    finalsim[[ncomb]] <- list(simset=simset, simdat=simdat)
                    ncomb <- ncomb + 1
                    
                }
                
            }
            
        }
        
    }
    
}

finalsim[which(sapply(finalsim, is.null))] <- NULL

saveRDS(finalsim, "simdat1.rds")
