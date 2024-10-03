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

### Best practice gseapy
# Extract ranked genes for a particular state (state1 vs others)
# Assuming you have performed differential expression analysis
ranked_genes = adata.uns['rank_genes_groups']['names']['state1']  # All genes ranked for state1

# Perform GSEA with the ranked gene list
gsea_results = gp.prerank(rnk=ranked_genes, gene_sets='KEGG_2021', organism='Human')

### Assuming transitions are of most importance
# Assuming you've performed pairwise comparisons
sc.tl.rank_genes_groups(adata, 'state', groups=['state1'], reference='state2', method='wilcoxon')

# Extract DEGs for the pairwise comparison
degs_state1_vs_state2 = adata.uns['rank_genes_groups']['names']['state1']

# Perform GSEA with pairwise DEGs
gsea_results = gp.enrichr(gene_list=degs_state1_vs_state2, gene_sets='KEGG_2021', organism='Human')

### Formatting data for immune cell deconvolution
# Select all columns except the first one for processing
df = df_values[df_values.columns]  # Use all columns

# Sort the DataFrame by the 'TargetName' column
df = df.sort_values(by='TargetName', ascending=True)

# Automatically rename columns (except the first one)
column_mapping = {df.columns[0]: 'Symbol'}  # First column remains 'Symbol'
column_mapping.update({df.columns[i]: f'Sample{i}' for i in range(1, len(df.columns))})

# Rename columns using the generated column_mapping
df = df.rename(columns=column_mapping)

# Perform log normalization on all columns except the first one (starting from index 1)
df.iloc[:, 1:] = np.log(df.iloc[:, 1:] + 1)

# Save the resulting DataFrame as a tab-separated file
df.to_csv('TestFile.txt', sep='\t', index=False)
