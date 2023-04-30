
library(Rcpp)
library(tweedie)

sourceCpp("get_Ms.cpp")
simdat <- readRDS("simdat1.rds")

nspot <- 100

generate_panel_posi <- function(row=10, col=10) {
    
    columnid <- rep(1:col, times=row)
    rowid <- rep(1:row, each=col)
    return(data.frame(spotid=1:(row*col), spotX=columnid, spotY=rowid))
    
}
spotlist <- generate_panel_posi(row=10, col=10)

spot_orig <- cbind(rep(1:nspot, each=nspot), rep(1:nspot, times=nspot))
dist_orig <- sqrt((spotlist$spotX[spot_orig[, 1]] - spotlist$spotX[spot_orig[, 2]])^2 + 
                      (spotlist$spotY[spot_orig[, 1]] - spotlist$spotY[spot_orig[, 2]])^2)

nsim <- 100

newsim <- vector(mode="list", length=length(simdat))

for (nset in 1:length(simdat)) {
    
    simnewdat <- vector("list", length=nsim)
    simset <- simdat[[nset]]$simset
    
    for (k in 1:nsim) {
        
        print(paste0("rho=", simset$rho, "; phi=", simset$phi, "; p=", simset$p, "; sprate=", simset$sprate, "; iter=", k))
        
        Y1 <- simdat[[nset]]$simdat[[k]]$Y
        M <- simdat[[nset]]$simdat[[k]]$M
        rmposi <- simdat[[nset]]$simdat[[k]]$rmposi
        
        if(length(rmposi) > 0) {
            spot <- spot_orig[-rmposi, ]
            dist <- dist_orig[-rmposi]
        } else {
            spot <- spot_orig
            dist <- dist_orig
        }
        
        Ms <- get_Ms(M, spot)
        
        X <- scale(Ms * exp(-simset$rho*dist), center=TRUE, scale=TRUE)
        X <- cbind(rep(1, nrow(X)), X)
        
        eta <- X %*% simset$truebeta + simset$reL[spot[, 1]] + simset$reR[spot[, 2]]
        
        eta[eta > 10] <- 10
        mu <- c(exp(eta))
        
        p1 <- mu / (1 + mu)
        alpha <- runif(length(Y1), 0.5, 5)
        gam <- mu / p1 / alpha
        
        Y <- rbinom(length(gam), prob=p1, size=1) * rgamma(length(gam), shape=alpha, scale=gam)
        Y[is.na(Y)] <- 0
        
        simnewdat[[k]] <- list(Y=as.vector(Y), M=M, rmposi=rmposi)
        
    }
    
    newsim[[nset]] <- list(simset=simset, simdat=simnewdat)
    
}


saveRDS(newsim, "simdat2.rds")