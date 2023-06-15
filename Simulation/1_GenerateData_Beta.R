
nspot <- 100
i <- 1
simset <- vector(mode="list", length=12)


# Generate the random effects v_i^L and v_j^R
set.seed(2023)
reL <- rnorm(nspot, 0, 0.4)
reR <- rnorm(nspot, 0, 0.4)

for(ncelltype in c(5, 10, 15)){
    
    set.seed(2023+20-ncelltype)
    
    # Generate the original positive betas
    truebeta <- runif(ncelltype*ncelltype+1, 0.1, 0.5)
    
    # Randomly treat half of them be negative
    truebeta[sample(length(truebeta), length(truebeta)/2)] <- truebeta[sample(length(truebeta), length(truebeta)/2)]*(-1)
    
    # Generate zeros (for sparsity of the betas)
    orig <- 1:(ncelltype^2)
    zero_0.2 <- sample(orig, round(ncelltype^2)*0.2)
    zero_0.4 <- c(zero_0.2, sample(orig[!(orig %in% zero_0.2)], round(ncelltype^2)*0.2))
    zero_0.6 <- c(zero_0.4, sample(orig[!(orig %in% zero_0.4)], round(ncelltype^2)*0.2))
    zero_0.8 <- c(zero_0.6, sample(orig[!(orig %in% zero_0.6)], round(ncelltype^2)*0.2))
    
    for(sprate in c(0.2, 0.4, 0.6, 0.8)) {
        
        if(sprate == 0.2) {
            truebeta[zero_0.2+1] <- 0
        }else if(sprate == 0.4) {
            truebeta[zero_0.4+1] <- 0
        }else if(sprate == 0.6) {
            truebeta[zero_0.6+1] <- 0
        }else if(sprate == 0.8) {
            truebeta[zero_0.8+1] <- 0
        }
        
        simset[[i]] <- list(sprate=sprate, ncelltype=ncelltype, 
                            truebeta=truebeta, reL=reL, reR=reR)
        i <- i + 1
        
    }
    
}


saveRDS(simset, "simdat_beta.rds")
