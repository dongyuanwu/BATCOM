
dat <- readRDS("simres1_batcom.rds")


##################################################################

get_performance_beta_qval_auc <- function(summarybeta, truebeta) {
    
    coefest <- summarybeta[, 4]   # summarybeta[, 4] is the mean of inference samples
    Z <- summarybeta[, 4] / summarybeta[, 7]   # summarybeta[, 7] is the SD of inference samples
    pval <- 2 * pnorm(abs(Z), lower.tail=FALSE)
    qval <- p.adjust(pval, method="fdr")
    
    sens <- spec <- NULL
    for(threshold in seq(0, 1, 0.001)) {
        
        pred_pos <- which(qval < threshold & coefest > 0)
        pred_neg <- which(qval < threshold & coefest < 0)
        pred_zero <- which(qval > threshold)
        
        # sensivitity = TP / (TP+FN):
        # specificity = TN / (TN+FP):
        
        #                   Pred
        #                sig   insig
        # Actual sig     TP      FN
        #      insig     FP      TN
        
        TP <- sum(truebeta[pred_pos] > 0) + sum(truebeta[pred_neg] < 0)
        TN <- sum(truebeta[pred_zero] == 0)
        FN <- sum(truebeta[pred_zero] != 0)
        FP <- sum(truebeta[pred_pos] <= 0) + sum(truebeta[pred_neg] >= 0)
        
        senstemp <- TP / (TP + FN)
        senstemp <- ifelse(is.na(senstemp), 0, senstemp)
        spectemp <- TN / (TN + FP)
        spectemp <- ifelse(is.na(spectemp), 0, spectemp)
        
        sens <- c(sens, senstemp)
        spec <- c(spec, spectemp)
        
    }
    
    return(pracma::trapz(1-spec, sens))
    
}



##

get_performance_beta_qval <- function(summarybeta, truebeta) {
    
    coefest <- summarybeta[, 4]   # summarybeta[, 4] is the mean of inference samples
    Z <- summarybeta[, 4] / summarybeta[, 7]   # summarybeta[, 7] is the SD of inference samples
    pval <- 2 * pnorm(abs(Z), lower.tail=FALSE)
    qval <- p.adjust(pval, method="fdr")
    
    pred_pos <- which(qval < 0.05 & coefest > 0)
    pred_neg <- which(qval < 0.05 & coefest < 0)
    pred_zero <- which(qval > 0.05)
    
    # sensivitity = TP / (TP+FN):
    # specificity = TN / (TN+FP):
    
    #                   Pred
    #                sig   insig
    # Actual sig     TP      FN
    #      insig     FP      TN
    
    TP <- sum(truebeta[pred_pos] > 0) + sum(truebeta[pred_neg] < 0)
    TN <- sum(truebeta[pred_zero] == 0)
    FN <- sum(truebeta[pred_zero] != 0)
    FP <- sum(truebeta[pred_pos] <= 0) + sum(truebeta[pred_neg] >= 0)
    
    PP <- length(pred_neg) + length(pred_pos)
    
    sens <- TP / (TP + FN)
    sens <- ifelse(is.na(sens), 0, sens)
    spec <- TN / (TN + FP)
    spec <- ifelse(is.na(spec), 0, spec)
    fdr <- FP / (FP + TP)
    fdr <- ifelse(is.na(fdr), 0, fdr)
    FPR <- FP / (FP + TN)
    FPR <- ifelse(is.na(FPR), 0, FPR)
    
    true_zero <- which(truebeta == 0)
    true_nonzero <- which(truebeta != 0)
    
    mse_nonzero <- mean((truebeta[true_nonzero] - coefest[true_nonzero])^2)
    mse_zero <- mean((truebeta[true_zero] - coefest[true_zero])^2)
    mse_total <- mean((truebeta - coefest)^2)
    
    return(c(sens, spec, fdr, FPR, PP, mse_nonzero, mse_zero, mse_total))
    
}

###################### maxprop version ##########################


get_performance_beta_qval_auc <- function(summarybeta, truebeta, nonzero_posi) {
    
    coefest <- rep(0, length(truebeta))
    sdest <- rep(1, length(truebeta))
    coefest[c(1, nonzero_posi+1)] <- summarybeta[, 4]
    sdest[c(1, nonzero_posi+1)] <- summarybeta[, 7]
    Z <- coefest / sdest
    pval <- 2 * pnorm(abs(Z), lower.tail=FALSE)
    qval <- p.adjust(pval, method="fdr")
    
    sens <- spec <- NULL
    for(threshold in seq(0, 1, 0.001)) {
        
        pred_pos <- which(qval < threshold & coefest > 0)
        pred_neg <- which(qval < threshold & coefest < 0)
        pred_zero <- which(qval > threshold)
        
        # sensivitity = TP / (TP+FN):
        # specificity = TN / (TN+FP):
        
        #                   Pred
        #                sig   insig
        # Actual sig     TP      FN
        #      insig     FP      TN
        
        TP <- sum(truebeta[pred_pos] > 0) + sum(truebeta[pred_neg] < 0)
        TN <- sum(truebeta[pred_zero] == 0)
        FN <- sum(truebeta[pred_zero] != 0)
        FP <- sum(truebeta[pred_pos] <= 0) + sum(truebeta[pred_neg] >= 0)
        
        senstemp <- TP / (TP + FN)
        senstemp <- ifelse(is.na(senstemp), 0, senstemp)
        spectemp <- TN / (TN + FP)
        spectemp <- ifelse(is.na(spectemp), 0, spectemp)
        
        sens <- c(sens, senstemp)
        spec <- c(spec, spectemp)
        
    }
    
    return(pracma::trapz(1-spec, sens))
    
}

##

get_performance_beta_qval <- function(summarybeta, truebeta, nonzero_posi) {
    
    coefest <- rep(0, length(truebeta))
    sdest <- rep(1, length(truebeta))
    coefest[c(1, nonzero_posi+1)] <- summarybeta[, 4]
    sdest[c(1, nonzero_posi+1)] <- summarybeta[, 7]
    Z <- coefest / sdest
    pval <- 2 * pnorm(abs(Z), lower.tail=FALSE)
    qval <- p.adjust(pval, method="fdr")
    
    pred_pos <- which(qval < 0.05 & coefest > 0)
    pred_neg <- which(qval < 0.05 & coefest < 0)
    pred_zero <- which(qval > 0.05)
    
    # sensivitity = TP / (TP+FN):
    # specificity = TN / (TN+FP):
    
    #                   Pred
    #                sig   insig
    # Actual sig     TP      FN
    #      insig     FP      TN
    
    TP <- sum(truebeta[pred_pos] > 0) + sum(truebeta[pred_neg] < 0)
    TN <- sum(truebeta[pred_zero] == 0)
    FN <- sum(truebeta[pred_zero] != 0)
    FP <- sum(truebeta[pred_pos] <= 0) + sum(truebeta[pred_neg] >= 0)
    
    PP <- length(pred_neg) + length(pred_pos)
    
    sens <- TP / (TP + FN)
    sens <- ifelse(is.na(sens), 0, sens)
    spec <- TN / (TN + FP)
    spec <- ifelse(is.na(spec), 0, spec)
    fdr <- FP / (FP + TP)
    fdr <- ifelse(is.na(fdr), 0, fdr)
    FPR <- FP / (FP + TN)
    FPR <- ifelse(is.na(FPR), 0, FPR)
    
    true_zero <- which(truebeta == 0)
    true_nonzero <- which(truebeta != 0)
    
    mse_nonzero <- mean((truebeta[true_nonzero] - coefest[true_nonzero])^2)
    mse_zero <- mean((truebeta[true_zero] - coefest[true_zero])^2)
    mse_total <- mean((truebeta - coefest)^2)
    
    return(c(sens, spec, fdr, FPR, PP, mse_nonzero, mse_zero, mse_total))
    
}

#############################################################

nsim <- 100
restab <- NULL
for (nset in 1:length(dat)) {
    
    print(nset)
    
    setting <- cbind(rep(dat[[nset]]$simset$rho, nsim), 
                     rep(dat[[nset]]$simset$phi, nsim), 
                     rep(dat[[nset]]$simset$p, nsim), 
                     rep(dat[[nset]]$simset$sprate, nsim), 
                     rep(dat[[nset]]$simset$ncelltype, nsim),
                     1:nsim)
    truebeta <- dat[[nset]]$simset$truebeta
    
    betapf <- t(sapply(dat[[nset]]$simres, function(x) get_performance_beta_qval(x$betares, truebeta)))
    # betapf <- t(sapply(dat[[nset]]$simres, function(x) get_performance_beta_qval(x$betares, truebeta, x$nonzero_posi)))
    
    betaauc <- sapply(dat[[nset]]$simres, function(x) get_performance_beta_qval_auc(x$betares, truebeta))
    # betaauc <- sapply(dat[[nset]]$simres, function(x) get_performance_beta_qval_auc(x$betares, truebeta, x$nonzero_posi))
    
    dic <- sapply(dat[[nset]]$simres, function(x) x$DIC)
    
    waic <- t(sapply(dat[[nset]]$simres, function(x) x$WAIC))
    
    restab <- rbind(restab, cbind(setting, betapf, betaauc, dic, waic))
    
}
restab <- as.data.frame(restab)
colnames(restab) <- c("truerho", "phi", "p", "sprate", "ncelltype", "nsim",
                      "sens", "spec", "fdr", "FPR", "PP", "mse_nonzero", "mse_zero", 
                      "mse_total", "AUC", "DIC", "WAIC1", "WAIC2")

saveRDS(restab, "restab1_batcom.rds")
