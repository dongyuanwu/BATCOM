
library(data.table)
library(Rcpp)
library(Seurat)

# load the SRT data

stdat <- t(read.table("GSM4284317_P2_ST_rep2_stdata.tsv"))
stmeta <- as.data.frame(colnames(stdat))
stmeta$x <- as.numeric(sapply(strsplit(stmeta[, 1], "x"), function(x) x[1]))
stmeta$y <- as.numeric(sapply(strsplit(stmeta[, 1], "x"), function(x) x[2]))
colnames(stmeta)[1] <- "spot"
stmeta$spot <- paste0("spot", 1:nrow(stmeta))
colnames(stdat) <- stmeta$spot
rownames(stmeta) <- stmeta$spot

# load the cell type proportion matrix

celltypeinfo <- readRDS("./RCTD/RCTDProportion_level2.rds")
colnames(celltypeinfo) <- gsub("\\.", "_", colnames(celltypeinfo))
celltype <- gsub(" ", "_", colnames(celltypeinfo))
celltype_interaction <- paste0(rep(celltype, each=length(celltype)), "_to_", rep(celltype, times=length(celltype)))

naposi <- which(!(rownames(stmeta) %in% rownames(celltypeinfo)))
stdat <- stdat[, -naposi]
stmeta <- stmeta[-naposi, ]

# filtering

keepspot <- which(apply(stdat, 2, function(x) sum(x != 0)) > 100)
stdat <- stdat[, keepspot]
stmeta <- stmeta[keepspot, ]
celltypeinfo <- celltypeinfo[keepspot, ]
rmgene <- which(apply(stdat, 1, function(x) sum(x!=0)) < (ncol(stdat) * 0.025))
stdat <- stdat[-rmgene, ]

# normalization using Seurat

dat <- CreateSeuratObject(counts=stdat, meta.data=stmeta)
dat <- NormalizeData(dat, assay = "RNA", scale.factor=10000)

stdat <- dat@assays$RNA@data

# We can use any of the LR pairs database, below is an example for CellChatDB

###################### cellchatDB ####################################

# library(CellChat)
# lrdb <- CellChatDB.human
# ligands <- receptors <- vector(mode="list", length=nrow(lrdb$interaction))
# 
# for(i in 1:length(ligands)) {
#     
#     if(!(lrdb$interaction$ligand[i] %in% rownames(lrdb$complex))) {
#         ligands[[i]] <- lrdb$interaction$ligand[i]
#     } else {
#         ltemp <- c(lrdb$complex[rownames(lrdb$complex) == lrdb$interaction$ligand[i], ])
#         ltemp <- ltemp[nzchar(ltemp)]
#         ligands[[i]] <- ltemp
#     }
#     
#     if(!(lrdb$interaction$receptor[i] %in% rownames(lrdb$complex))) {
#         receptors[[i]] <- lrdb$interaction$receptor[i]
#     } else {
#         ltemp <- c(lrdb$complex[rownames(lrdb$complex) == lrdb$interaction$receptor[i], ])
#         ltemp <- ltemp[nzchar(ltemp)]
#         receptors[[i]] <- ltemp
#     }
#     
# }
# ligands <- lapply(ligands, unlist)
# receptors <- lapply(receptors, unlist)

# We can also use our own LR pairs

######################### Self LRPairs #############################

lrpairs <- readRDS("lrpairs.rds")
lrpairs <- subset(lrpairs, species == "Human")

ligands <- lrpairs$ligand
receptors <- lrpairs$receptor

######################################################################

# check existing pairs
checkl <- sapply(ligands, function(x) any(x %in% rownames(stdat)))
checkr <- sapply(receptors, function(x) any(x %in% rownames(stdat)))
ligands <- ligands[checkl & checkr]
receptors <- receptors[checkl & checkr]

sourceCpp("get_Ck.cpp")
sourceCpp("get_Ck_v.cpp")
sourceCpp("get_Ms.cpp")
sourceCpp("if_neighbors_samesample.cpp")
sourceCpp("get_dist.cpp")

start <- Sys.time()

spot <- data.frame(spoti=rep(1:ncol(stdat), each=ncol(stdat)), 
                   spotj=rep(1:ncol(stdat), times=ncol(stdat)))
# set a threshold to control the number of spot pairs will be included in the model fitting
check_neighbor <- if_neighbors_samesample(stmeta$x, stmeta$y, dist=3)
spot <- as.matrix(spot[check_neighbor, ])

Sys.time() - start


############################## each lr pairs ######################################

# Calculate communication scores C_ij for each LR pairs

expr1 <- as.matrix(stdat)
genenames <- rownames(stdat)

Cspots <- vector(mode="list", length(ligands))

start <- Sys.time()
for(k in 1:length(ligands)) {
    
    expr <- expr1[genenames %in% c(ligands[[k]], receptors[[k]]), ]
    
    if(class(expr)[1] != "matrix") {
        
        Cspots[[k]] <- get_Ck_v(expr, spot)
        
    } else {
        
        ligand <- which(rownames(expr) %in% ligands[[k]])
        receptor <- which(rownames(expr) %in% receptors[[k]])
        
        Cspots[[k]] <- get_Ck(expr, spot, ligand, receptor)
        
    }
    
}
Sys.time() - start

# Calculate spot distances

start <- Sys.time()
D <- get_dist(stmeta$x, stmeta$y, spot)
posi <- sapply(Cspots, function(x) sum(x) != 0)
ligands <- ligands[posi]
receptors <- receptors[posi]
Cspots <- Cspots[posi]
Sys.time() - start

# Filtering out some LR pairs that do not have enough information

nonzero_prop <- sapply(Cspots, function(x) sum(x != 0) / nrow(spot))
summary(nonzero_prop)
which.max(nonzero_prop)
nrow(spot)*0.01
sum(nonzero_prop >= 0.02)
ligands <- ligands[nonzero_prop >= 0.02]
receptors <- receptors[nonzero_prop >= 0.02]
Cspots <- Cspots[nonzero_prop >= 0.02]

###############################################################

start <- Sys.time()
M <- as.matrix(celltypeinfo)

adj_M_row <- function(M_row) {
    
    M_row <- ifelse(M_row == max(M_row), 1, 0)
    return (M_row)
    
}

# M <- t(apply(M, 1, adj_M_row))   # If using MAXPROP version
Ms <- get_Ms(M, spot)   # Calculate M*M
Sys.time() - start

X <- Ms * exp(-0.5*D)  # can change 0.5 to other values for the tuning parameter rho
celltype_interaction <- celltype_interaction[-which(apply(X, 2, max) < 0.01)]
X <- X[, -which(apply(X, 2, max) < 0.01)]
X <- scale(X, center=TRUE, scale=TRUE)
X <- cbind(rep(1, nrow(X)), X)

dat <- list(X=X, Y=Cspots, ligands=ligands, receptors=receptors, 
            celltype=celltype, celltype_interaction=celltype_interaction, spot=spot)
saveRDS(dat, "RCTD_HMC_dist3_estrho0.5_dat.rds")
