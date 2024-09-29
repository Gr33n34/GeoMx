import scanpy as sc

# Log-normalize the data
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

### Find relevant genes for each group
# Assuming that 'state' in adata.obs contains the cancer states (e.g., state1, state2, state3, state4)
sc.tl.rank_genes_groups(adata, 'Gleason_grade', method='wilcoxon')

# View the results
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)

### Perform Pseudotime Analysis
import scanpy.external as sce

# Pseudotime analysis using Palantir
sce.tl.palantir(adata, n_comps=30)

# Visualize pseudotime trajectory
sc.pl.palantir_results(adata)
# Alternative Monocle3

### Mixed effect model
import statsmodels.formula.api as smf
# For a single gene example, you can loop through multiple genes as needed
genes = adata.var.index  # List of all gene names
gene_expression = adata.X.toarray()
df_expression = pd.DataFrame(gene_expression, index=adata.obs.index, columns=adata.var.index)
for gene in genes:
    metadata['expression'] = df_expression[gene]  # Update the expression for the current gene
    result = smf.mixedlm("gene_expression ~ state", data=adata.obs, groups=adata.obs['patient_id']).fit()
    #result = smf.mixedlm("gene_expression ~ state", data=metadata, groups=adata.obs['patient_id']).fit()
    print(f"Results for {gene}:")
    print(result.summary())