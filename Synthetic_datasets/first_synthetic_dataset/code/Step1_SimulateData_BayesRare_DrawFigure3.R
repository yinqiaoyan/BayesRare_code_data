##############################################################
####################### Simulation 1 #########################
##############################################################

########################## Note ##############################
## Choose the working directory as the directory where 
## this code file locates.
##############################################################

########################## Note ##############################
## Please first install required R packages.
## R: System_preparation.R
##############################################################


rm(list = ls())
library(mvtnorm)
library(ggplot2)
library(Seurat)
library(BayesRare)


## Set the working directory to the source file location
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_dir)
print(getwd())


## Data generation
set.seed(123)
n_mat = matrix(NA, nrow=8, ncol=8)
# controls
n_mat[1, ]  = c(80, 80, 80,  0, 15, 1000, 1500, 2000)
n_mat[2, ]  = c(80,  0, 79,  0,  0, 1000, 1500, 2000)
n_mat[3, ]  = c(80, 81, 79,  0, 15, 1000, 1500, 2000)
n_mat[4, ]  = c( 0, 79, 80,  0,  0, 1000, 1500, 2000)
# patients
n_mat[5, ]  = c(81, 79,  0, 15,  0, 1000, 1500, 2000)
n_mat[6, ]  = c(81, 79,  0, 15, 15, 1000, 1500, 2000)
n_mat[7, ]  = c(81,  0, 80, 15,  0, 1000, 1500, 2000)
n_mat[8, ]  = c( 0, 79, 80, 15, 15, 1000, 1500, 2000)

#                  1    2    3    4    5    6    7    8
true_mean1   = c( 10,   2,   1,  -3, -10,  10,   0,  12)
true_mean2   = c( 10,   7,   3,  -4,  -3,   0,   9,   9)
true_sigmasq = c(1.2, 1.2, 1.2, 0.1, 0.1, 1.5, 1.5, 1.5)

D = nrow(n_mat)
K = length(true_mean1)
p_dim = 2

true_labels_list = vector("list", D)
for (d in 1:D) {
  true_labels_list[[d]] = rep(1:length(n_mat[d,]), n_mat[d,])
}

K_each_mat = n_mat
K_each_mat[n_mat > 0] = 1


X_list <- X_all_list <- vector("list", D)
expr_pc_data_nonrare = NULL
for (d in 1:D) {
  for (ii in which(K_each_mat[d,] > 0)) {
    tmpmat = rmvnorm(n_mat[d, ][ii], mean = c(true_mean1[ii], true_mean2[ii]), sigma = diag(rep(true_sigmasq[ii], p_dim)))
    if (ii <= 5) {
      X_list[[d]] = rbind(X_list[[d]], tmpmat)
    }
    if (ii > 5) {
      expr_pc_data_nonrare = rbind(expr_pc_data_nonrare, tmpmat)
    }
    X_all_list[[d]] = rbind(X_all_list[[d]], tmpmat)
  }
}

expr_pc_data = do.call(rbind, X_list)


## Use Seurat
pca_mat <- expr_pc_data_nonrare
pca_mat <- as.matrix(pca_mat)

ncells <- nrow(pca_mat)
npc    <- ncol(pca_mat)

rownames(pca_mat) <- paste0("Cell_", seq_len(ncells))
colnames(pca_mat) <- paste0("PC_", seq_len(npc))

feature_name <- "fakegene"  
placeholder_counts <- matrix(0, nrow=npc, ncol=ncells)

# Create Seurat object
obj <- CreateSeuratObject(counts = placeholder_counts, project = "FromPCA")
DefaultAssay(obj) <- "RNA"

pca_dr <- CreateDimReducObject(
  embeddings = pca_mat,
  key = "PC_",
  assay = DefaultAssay(obj)
)
obj[["pca"]] <- pca_dr

use_dims <- 1:npc
set.seed(42)
obj <- FindNeighbors(obj, reduction = "pca", dims = use_dims)
obj <- FindClusters(obj, resolution = 0.1) 

mu_gk_fixed = NULL
sigmasq_gk_fixed = NULL
n_k_abundant = NULL

res_seurat = as.numeric(obj@meta.data[["seurat_clusters"]])
uniq_seurat = sort(unique(res_seurat))

for (k_abund in seq_len(length(uniq_seurat))) {
  tmp_ids = which(res_seurat == uniq_seurat[k_abund])
  tmp_expr = expr_pc_data_nonrare[tmp_ids, ]
  mu_gk_fixed = cbind(mu_gk_fixed, colMeans(tmp_expr))
  sigmasq_gk_fixed = cbind(sigmasq_gk_fixed, apply(tmp_expr, 2, var))
  n_k_abundant = c(n_k_abundant, length(tmp_ids))
}

## abundant clusters
centers_abundant = t(mu_gk_fixed)  # Kp x G

## Run BayesRare
res_BayesRare = BayesRare_train(X_list=X_list, 
                                K_init=5, 
                                mu_gk_fixed=mu_gk_fixed, 
                                sigmasq_gk_fixed=sigmasq_gk_fixed,
                                patient_group_ids=5:8, control_group_ids=1:4,
                                random_seed=19,
                                do_test=TRUE,
                                verbose_time=TRUE, verbose_em=TRUE)




##############################################################
#####         Figure 3(a) Cell type annotations          #####
##############################################################
cell_type_labels = NULL
for (d in 1:8) {
  cell_type_labels = c(cell_type_labels, rep(paste0("A", 1:8), n_mat[d, ]))
}
df_all = data.frame(rbind(do.call(rbind, X_all_list)))
colnames(df_all) <- c("X1", "X2")
df_all$true_type = cell_type_labels

plot_color <- c(
  "A1" = "#c63a32", 
  "A2" = "#3977af", 
  "A3" = "#f08536", 
  "A4" = "#4f9b6c", 
  "A5" = "#84584e", 
  "A6" = "#f7d8ed", 
  "A7" = "#E9EBD0", 
  "A8" = "#EEE2E0" 
)

p <- ggplot(df_all, aes(x = X1, y = X2, color = factor(true_type))) +
  geom_point(
    data = subset(df_all, true_type %in% c("A6","A7","A8")),
    aes(color = factor(true_type)),
    size = 0.5, alpha = 0.8
  ) +
  geom_point(
    data = subset(df_all, true_type %in% c("A1","A2","A3")),
    aes(color = factor(true_type)),
    size = 1.5, alpha = 0.8
  ) +
  geom_point(
    data = subset(df_all, true_type %in% c("A4","A5")),
    aes(color = factor(true_type)),
    size = 1.5, alpha = 0.8
  ) +
  theme_classic(base_size = 37) +
  labs(
    x = "Dimension 1",
    y = "Dimension 2"
  ) +
  scale_color_manual(
    name   = NULL,
    breaks = paste0("A", 1:8),  
    values = plot_color
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    plot.background  = element_blank(),
    legend.text      = element_text(size = 25),
    legend.key.size  = unit(2.5, "line")
  ) +
  guides(color = guide_legend(title = NULL, override.aes = list(size = 7))) +
  coord_cartesian(xlim = c(-11, 15), ylim = c(-5, 15))

ggsave("../figures/figure3a.png", p, width = 11, height = 9, dpi = 100)




##############################################################
#####    Figure 3(b) Initial rare cell identification    #####
##############################################################
df_abundant = as.data.frame(centers_abundant)
df_all_abund = data.frame(rbind(do.call(rbind, X_list), df_abundant))
colnames(df_all_abund) <- c("X1", "X2")

n_mat_vec = as.numeric(t(n_mat[, 1:5]))
lab <- rep(1:5, length.out = length(n_mat_vec))
res <- unlist(Map(rep, lab, n_mat_vec), use.names = FALSE)
res <- c(res, rep("Abundant", nrow(centers_abundant)))
df_all_abund$type = res
df_all_abund$type[df_all_abund$type %in% c("1","2","3")] = "FP"
df_all_abund$type[df_all_abund$type %in% c("4","5")] = "True rare"

df_all_abund$type <- factor(df_all_abund$type,
                            levels = c("True rare", "FP", "Abundant"))

cols_map <- c(
  "True rare" = "#c63a32",   # red
  "FP"        = "#c2b2d2",   # light purple
  "Abundant"  = "#000000"    # black
)

p <- ggplot(df_all_abund, aes(x = X1, y = X2, color = type)) +
  geom_point(
    data = subset(df_all_abund, type %in% c("True rare","FP")),
    aes(color = factor(type)),
    size = 1.5, alpha = 0.8
  ) +
  geom_point(
    data = subset(df_all_abund, type == "Abundant"),
    aes(color = factor(type)),
    size = 4, alpha = 0.8
  ) +
  theme_classic(base_size = 37) +
  labs(x = "Dimension 1", y = "Dimension 2") +
  scale_color_manual(
    name   = NULL,
    breaks = c("True rare", "FP", "Abundant"), 
    labels = c("True rare", "False positive\ncells", "Abundant\ncenters"),
    values = cols_map
  ) +
  theme(
    legend.position = "right",
    plot.title      = element_text(hjust = 0.5),
    panel.background= element_blank(),
    plot.background = element_blank(),
    legend.text     = element_text(size = 25),
    legend.key.size = unit(4, "line")
  ) +
  guides(color = guide_legend(override.aes = list(size = 7))) + 
  coord_cartesian(xlim = c(-11, 15), ylim = c(-5, 15))

ggsave("../figures/figure3b.png", p, width = 12, height = 9, dpi = 100)




##############################################################
#####            Figure 3(c) BayesRare results           #####
##############################################################
umap_df_plot = as.data.frame(expr_pc_data)
tmpLabels = unlist(res_BayesRare$z_res)

tmpc2 = tmpLabels
tmpc2[tmpLabels == 1] = 5
tmpc2[tmpLabels == 2] = 4
tmpc2[tmpLabels == 3] = 3
tmpc2[tmpLabels == 4] = 1
tmpc2[tmpLabels == 5] = 2
tmpLabels = tmpc2

umap_df_plot$Group = paste0("C", tmpLabels)

plot_color <- c(
  "#c63a32", "#3977af", "#f08536", "#4f9b6c", "#84584e", 
  "#d57fbe", "#b6bc6d", "#be9e96"
)

p <- ggplot(umap_df_plot, aes(x = V1, y = V2, color = factor(Group))) +
  geom_point(size = 1.5, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = plot_color) +
  theme_classic(base_size = 37) +
  labs(
    x = "Dimension 1",
    y = "Dimension 2"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    plot.background  = element_blank(),
    legend.text      = element_text(size = 25),
    legend.key.size  = unit(2.5, "line")
  ) +
  guides(color = guide_legend(title = NULL, override.aes = list(size = 7))) +
  coord_cartesian(xlim = c(-11, 15), ylim = c(-5, 15))

ggsave("../figures/figure3c.png", p, width = 11, height = 9, dpi = 100)




