library(SingleCellExperiment)
library(cellassign)
library(readr)

pbmc <- readRDS("data/pbmc.rds")
marker_mat <- read_csv("data/FL_celltype.csv")
marker_mat <- marker_mat[-c(24),]
col <- marker_mat[,1]
rownames <- col[['Gene']]
counts(pbmc) <- assay(pbmc, "X")
mm <- data.matrix(marker_mat)
mm <- mm[,-1]
rownames(mm) <- rownames
times <- c()
n_obs <- c()
n <- 68579
i <- 1
s <- colSums(counts(pbmc))
while (n > 1000){
start <- Sys.time()
fit <- cellassign(exprs_obj = pbmc[rownames,1:n], 
                  marker_gene_info = mm,
                  s=s[1:n],
                  learning_rate = 1e-2, 
                  shrinkage = TRUE,
                  verbose = TRUE)
times[i] <- difftime(Sys.time(), start, units='secs')
n_obs[i] <- n
i <- i + 1
n <- n %/% 2
}
df <- data.frame("n_obs" = n_obs, "time" = times)
write.csv(df, "data/cell_assign_r_runtime.csv")