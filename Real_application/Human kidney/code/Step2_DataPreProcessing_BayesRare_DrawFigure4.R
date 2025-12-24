##############################################################
####################     Human kidney    #####################
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
library(BayesRare)
library(ggplot2)

## Set the working directory to the source file location
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_dir)
print(getwd())


## Read in data
expr_pc = read.csv("../input_data/X_pca_py.csv")  # python preprocessing code
meta_data_all = read.csv("../input_data/raw_metadata.tsv", sep = "\t")
clust_res_all = read.csv("../input_data/scCAD_detection_results.tsv", sep = "\t")

colnames(meta_data_all)[colnames(meta_data_all) == "donor_id"]  <- "donor_id"
colnames(meta_data_all)[colnames(meta_data_all) == "cell_type"] <- "cell_type"

rare_candidates_id = which(clust_res_all$rare_subcluster_id != -1)

out_proc = data_preprocess(expr_pc, meta_data_all, clust_res_all)
res_BayesRare = BayesRare_train(X_list = out_proc$X_list, 
                                K_init = out_proc$K_init, 
                                mu_gk_fixed = out_proc$mu_gk_fixed, 
                                sigmasq_gk_fixed = out_proc$sigmasq_gk_fixed,
                                random_seed = 175,
                                verbose_time = TRUE, verbose_em = TRUE)


true_rare_types = c("B cell", "conventional dendritic cell", "mast cell",
                    "mature NK T cell", "mononuclear phagocyte", "neural cell",
                    "neutrophil", "non-classical monocyte", "papillary tips cell",
                    "parietal epithelial cell", "plasma cell", "plasmacytoid dendritic cell, human")
total_labels = meta_data_all$cell_type[res_BayesRare$rare_ids]
cat("Number of detected rare cells:", length(res_BayesRare$rare_ids), " \n")
cat("Number of true rare cells:", sum(total_labels %in% true_rare_types), " \n")



##############################################################
#####         Figure 4(a) Cell type annotations          #####
##############################################################
plot_color <- c(
  "#c2b2d2", "#c63a32", "#3977af", "#f08536", "#f6bd82", "#4f9b6c", "#84584e", 
  "#d57fbe", "#b6bc6d", "#be9e96", "#54b0c2", "#9d51f3", "#a8dc93",
  "#ff7f50", "#4682b4", "#dda0dd", "#ffb6c1", "#20b2aa", "#ffd700", 
  "#87cefa", "#ff6347", "#90ee90", "#ff69b4", "#708090", "#cd853f", "#00ced1", "#8b0000"
)

umap_df = read.csv("../input_data/X_umap_py.csv")
umap_df_plot = umap_df[, -1]

ann_df = data.frame(types=meta_data_all$cell_type, ann_labels="0")
uniq_cell_type = sort(unique(meta_data_all$cell_type))
for (ll in seq_along(uniq_cell_type)) {
  ann_df$ann_labels[ann_df$types == uniq_cell_type[ll]] = paste0("A", ll)
}

ord <- paste0("A", 1:26)
umap_df_plot$Group = ann_df$ann_labels
umap_df_plot$Group = factor(umap_df_plot$Group, levels = ord)

p <- ggplot(umap_df_plot, aes(x = UMAP1, y = UMAP2, color = factor(Group))) +
  geom_point(size = 0.5, alpha = 0.8) +
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
    legend.text      = element_text(size = 20),
    legend.key.size  = unit(2.3, "line")
  ) +
  guides(color = guide_legend(title = NULL, override.aes = list(size = 6)))

ggsave("../figures/figure4a.png", p, width = 12, height = 9, dpi = 100)




##############################################################
#####                   Figure 4(b)-(f)                  #####
##############################################################
giniclust_res = read.csv("../input_data/res_GiniClust.csv")
scissors_res  = read.csv("../input_data/res_SCISSORS.csv")
tmpLetter = letters[2:6]

for (vv in tmpLetter){
  if (vv == "b") {
    ## true rare cells (b)
    umap_df_plot$Group = "abundant cells"
    umap_df_plot$Group[meta_data_all$cell_type %in% true_rare_types] = "rare cells"
  } else if (vv == "c") {
    ## giniclust (c)
    umap_df_plot$Group = "abundant cells"
    umap_df_plot$Group[umap_df$index %in% giniclust_res$cell_name] = "rare cells" 
  } else if (vv == "d") {
    ## SCISSORS (d)
    umap_df_plot$Group = "abundant cells"
    umap_df_plot$Group[umap_df$index %in% scissors_res$Cell] = "rare cells" 
  } else if (vv == "e") {
    ## scCAD (e)
    umap_df_plot$Group = "abundant cells"
    umap_df_plot$Group[rare_candidates_id] = "rare cells"
  } else if (vv == "f") {
    ## Our (f)
    umap_df_plot$Group = "abundant cells"
    umap_df_plot$Group[res_BayesRare$rare_ids] = "rare cells"
  }
  
  ## abundant vs rare
  p <- ggplot(umap_df_plot, aes(x = UMAP1, y = UMAP2)) +
    # abundant cells
    geom_point(
      data = subset(umap_df_plot, Group == "abundant cells"),
      aes(color = factor(Group)),
      size = 0.1, alpha = 0.8
    ) +
    # rare cells
    geom_point(
      data = subset(umap_df_plot, Group == "rare cells"),
      aes(color = factor(Group)),
      size = 1, alpha = 0.8
    ) +
    theme_classic(base_size = 37) +
    labs(
      x = "Dimension 1",
      y = "Dimension 2"
    ) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5),
      panel.background = element_blank(), 
      plot.background  = element_blank()  
    ) +
    scale_color_manual(
      name   = NULL,
      breaks = c("abundant cells", "rare cells"),
      labels = c("Ab", "Ra"),
      values = plot_color  
    ) +
    guides(color = guide_legend(override.aes = list(size = 10))) +
    theme(
      legend.position   = "right",
      legend.text       = element_text(size = 34), # 图例字体
      legend.key.size   = unit(3, "line")
    )
  
  ggsave(paste0("../figures/figure4", vv, ".png"), p, width = 12, height = 9, dpi = 100)
}




