import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import umap
import matplotlib.pyplot as plt
import scanpy as sc

## Read in raw data
expr_df = pd.read_csv("../raw_data/raw_gene_expression_matrix.tsv", sep="\t", index_col=0)

expr_logged = np.log1p(expr_df)
top_genes = expr_logged.var(axis=0).sort_values(ascending=False).head(1000).index
expr_top = expr_logged[top_genes]

X_pca = PCA(n_components=50, random_state=42).fit_transform(expr_top)
embedding = umap.UMAP(random_state=42).fit_transform(X_pca)

umap_df = pd.DataFrame(
    embedding,
    columns=["UMAP1", "UMAP2"],
    index=expr_top.index
)

## Save results
cell_ids = expr_top.index 
df_pca = pd.DataFrame(X_pca, 
                      index=cell_ids, 
                      columns=[f"PC{i+1}" for i in range(X_pca.shape[1])])

df_pca.to_csv("../input_data/X_pca_py.csv")
umap_df.to_csv("../input_data/X_umap_py.csv")