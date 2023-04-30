
simdat <- readRDS("RCTD_HMC_dist3_estrho0.5_dat.rds")
celltype <- simdat$celltype
celltype_interaction <- simdat$celltype_interaction
celltype_interaction_split <- strsplit(celltype_interaction, "_to_")

simres <- readRDS("RCTD_HMC_dist3_estrho0.5_res.rds")

params.store <- simres[[1]]$paramstemp   # need to specify the detailed LR pair k, here k=1
phi_median <- format(round(median(exp(params.store[, 1])), 2), nsmall=2)
phi_ci <- format(round(c(quantile(exp(params.store[, 1]), 0.025), quantile(exp(params.store[, 1]), 0.975)), 2), nsmall=2)
p_store <- (2*exp(params.store[,2])+1)/(exp(params.store[,2])+1)
p_median <- format(round(median(p_store), 2), nsmall=2)
p_ci <- format(round(c(quantile(p_store, 0.025), quantile(p_store, 0.975)), 2), nsmall=2)

phiexpr <- paste0(phi_median, " (", phi_ci[1], ", ", phi_ci[2], ")")
pexpr <- paste0(p_median, " (", p_ci[1], ", ", p_ci[2], ")")

beta_mean <- colMeans(params.store)[-c(1:3)]
beta_sd <- apply(params.store, 2, sd)[-c(1:3)]
Z <- beta_mean / beta_sd
pval <- 2 * pnorm(abs(Z), lower.tail=FALSE)
qval <- p.adjust(pval, method="fdr")


### Heatmap

heatmat <- matrix(NA, nrow=length(celltype), ncol=length(celltype))
rownames(heatmat) <- celltype
colnames(heatmat) <- celltype

celltype_interaction_split <- strsplit(celltype_interaction, "_to_")
for (i in 1:length(celltype_interaction)) {
    if (qval[i] < 0.2) {
        heatmat[rownames(heatmat) == celltype_interaction_split[[i]][1],
                colnames(heatmat) == celltype_interaction_split[[i]][2]] <- beta_mean[i]
    } else {
        heatmat[rownames(heatmat) == celltype_interaction_split[[i]][1],
                colnames(heatmat) == celltype_interaction_split[[i]][2]] <- 0
    }
}

library(ComplexHeatmap)

par(mar=c(1, 1, 1, 3), mfrow=c(1, 2))
selfcolor <- circlize::colorRamp2(c(min(beta_mean), 0, max(beta_mean)), 
                                  c("#86aaeb", "#FFFFFF", "#fc8383"))
htmp <- Heatmap(heatmat, cluster_rows=FALSE, cluster_columns=FALSE, col=selfcolor,
                row_title="Sender Cell Type", column_title="Receiver Cell Type",
                heatmap_legend_param=list(title="Communication",
                                          title_position = "leftcenter-rot",
                                          labels_gp = gpar(fontsize = 5),
                                          title_gp = gpar(fontsize = 8, fontface = "bold"),
                                          legend_heigth = unit(1, "cm")),
                rect_gp = gpar(col = "gray30", lwd = 1),
                column_names_gp = grid::gpar(fontsize = 8),
                row_names_gp = grid::gpar(fontsize = 8),
                column_title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
                row_title_gp = grid::gpar(fontsize = 10, fontface = "bold"))
# draw(htmp, heatmap_legend_side="bottom", annotation_legend_side="right",
#      legend_grouping = "original")






### Network

restab <- data.frame("Celltype_Sender"=sapply(celltype_interaction_split, function(x) x[1]),
                     "Celltype_Receiver"=sapply(celltype_interaction_split, function(x) x[2]),
                     "Estimate"=beta_mean,
                     "Pvalue"=pval,
                     "Qvalue"=qval)

library(igraph)

restab <- restab[restab$Qvalue < 0.05, ]

netdf <- restab[, -c(4:5)]

net1 <- graph_from_data_frame(netdf, directed=TRUE)
E(net1)$width <- abs(netdf$Estimate) * 10

E(net1)$color[netdf$Estimate < 0] <- '#a0bbeb'
E(net1)$color[netdf$Estimate > 0] <- '#FF9797'

# V(net1)$size

set.seed(20230421)
test.layout1 <- layout_nicely(net1)



### Combine two plots


library(ggplotify)

a <- as.grob(htmp)
b <- base2grob(~plot(net1, vertex.color="#C5E0B4", edge.color = E(net1)$color,
                     vertex.label.color="black", vertex.size=10,
                    vertex.label.font=1.8, vertex.label.family="sans",
                    vertex.frame.color="gray50", edge.color="gray", 
                    vertex.label.cex=0.7, 
                    edge.arrow.size=1, edge.curved=0.1,
                    layout=test.layout1))

library(cowplot)
pdf("RealDataPlot.pdf", width=10, height=4.5)
plot_grid(a, b, ncol=2, labels=LETTERS[1:2])
dev.off()

