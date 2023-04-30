
library(data.table)

scdat <- as.data.frame(fread("merge10pts_counts.txt"))

scdat <- scdat[, c(1, which(grepl("P2_", colnames(scdat))))]
rownames(scdat) <- scdat[, 1]
scdat <- scdat[-c(1:2), -1]

scmeta <- read.table("patient_metadata_new.txt")
scmeta <- scmeta[scmeta$patient == "P2", ]

stdat <- t(read.table("GSM4284317_P2_ST_rep2_stdata.tsv"))
stmeta <- as.data.frame(colnames(stdat))
stmeta$x <- as.numeric(sapply(strsplit(stmeta[, 1], "x"), function(x) x[1]))
stmeta$y <- as.numeric(sapply(strsplit(stmeta[, 1], "x"), function(x) x[2]))
colnames(stmeta)[1] <- "spot"
stmeta$spot <- paste0("spot", 1:nrow(stmeta))
colnames(stdat) <- stmeta$spot
rownames(stmeta) <- stmeta$spot

rmgene <- which(apply(scdat, 1, function(x) sum(x == 0)) > 0.975*ncol(scdat))
scdat <- scdat[-rmgene, ]

commongene <- rownames(stdat)[rownames(stdat) %in% rownames(scdat)]
scdat <- scdat[rownames(scdat) %in% commongene, ]
stdat <- stdat[rownames(stdat) %in% commongene, ]
scdat <- scdat[order(match(rownames(scdat), commongene)), ]


## RCTD

library(spacexr)
celltype <- as.factor(scmeta$level2_celltype)
names(celltype) <- rownames(scmeta)
reference <- Reference(scdat, celltype)
puck <- SpatialRNA(stmeta[, -1], stdat)
myRCTD <- create.RCTD(puck, reference, CELL_MIN_INSTANCE = 5, max_cores=6)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
results <- myRCTD@results
# normalize the cell type proportions to sum to 1.
norm_weights <- as.data.frame(as.matrix(normalize_weights(results$weights)))
rowSums(norm_weights)

saveRDS(norm_weights, "RCTDProportion_level2.rds")
