import pandas as pd
import numpy as np
import scipy as sp
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import os
import math
import anndata as ad

def GeoMx_Reader(file_path, metadata_columns):
    """
    Reads a GeoMx-formatted Excel file and converts it to an AnnData object.

    Parameters:
    file_path (str): Path to the Excel file.

    Returns:
    AnnData: Annotated data object containing gene counts and metadata.
    """
    # Read the TargetCountMatrix sheet
    target_count_df = pd.read_excel(file_path, sheet_name="TargetCountMatrix")

    # Extract gene names and count data
    gene_names = target_count_df["TargetName"].values
    count_data = target_count_df.drop(columns=["TargetName"]).set_index(gene_names).T

    # Read the SegmentProperties sheet
    segment_properties_df = pd.read_excel(file_path, sheet_name="SegmentProperties")

    # Filter the metadata to include relevant columns
    metadata = segment_properties_df[metadata_columns]

    # Ensure sample names in count data match the SegmentDisplayName in metadata
    metadata = metadata.set_index("SegmentDisplayName")
    count_data = count_data.loc[metadata.index]

    # Create the AnnData object
    adata = ad.AnnData(X=count_data.values, var=pd.DataFrame(index=gene_names), obs=metadata)

    return adata

def filter_by_median(adata, required_columns):
    # Ensure the columns exist in adata.obs
    for col in required_columns:
        if col not in adata.obs:
            raise ValueError(f"Column '{col}' is missing in adata.obs.")

    # Initialize a list to collect indices of rows to keep
    indices_to_keep = []

    # Group by "Patient" and "AOI_Diagnose" and process each group
    grouped = adata.obs.groupby(required_columns)
    for (patient, diagnose), group in grouped:
        if len(group) == 0:
            continue
        
        # Calculate the median of the group
        median_row = group.median(numeric_only=True)
        
        # Find the row in the original group that is closest to the median values
        numeric_columns = group.select_dtypes(include=[np.number]).columns
        if numeric_columns.empty:
            continue
        
        group_numeric = group[numeric_columns]
        median_values = median_row[numeric_columns]
        distances = (group_numeric - median_values).abs().sum(axis=1)
        
        closest_idx = distances.idxmin()
        indices_to_keep.append(closest_idx)
    
    # Create a new AnnData object with the rows to keep
    new_adata = adata[indices_to_keep]
    return new_adata

def vulcano_plot(data, plotsize_x=10, plotsize_y=6, log2fc_threshold=2.0, padj_threshold=0.05,
                 title='Vulcano plot', legend_loc='upper right', grid=True, save_as_svg=False, 
                 svg_filename='vulcano_plot.svg', label=True, downsample=False, pv=None):
    """
    Reads DataFrame with log10 p-values and plots a volcano plot.

    Parameters:
    data (pd.DataFrame): Contains columns 'Target', 'log2FoldChange', 'padj'
    save_as_svg (bool): If True, saves the plot as an SVG file.
    svg_filename (str): Filename for saving the plot if save_as_svg is True.
    """
    # Enables manually choosing another column for pvalues 
    if pv is None:
        pv = 'padj'

    # Convert padj to -log10(padj), setting NaN values to 1.0
    data['negLog10padj'] = -np.log10(data[pv].fillna(1.0))
    
    # Compute significance based on padj and log2FC thresholds
    data['Significant'] = (data[pv] < padj_threshold) & ((data['log2FoldChange'] > log2fc_threshold) | (data['log2FoldChange'] < -log2fc_threshold))
    
    # Figure dimensions
    plt.figure(figsize=(plotsize_x, plotsize_y))
    
    # Helper function to decrease file size
    if downsample == True:
        significant = data[data['Significant'] == True]
        non_significant = data[data['Significant'] == False]
        # Randomly select 10% of non-significant points
        non_significant_sample = non_significant.sample(frac=0.1, random_state=42)
        # Concatenate the sampled non-significant points with the significant ones
        data = pd.concat([significant, non_significant_sample])
        
    # Scatter plot with coloring based on significance
    plt.scatter(
        data['log2FoldChange'], 
        data['negLog10padj'], 
        c=data['Significant'].map({True: 'red', False: 'gray'}),
        alpha=0.7,
        label='Genes'
    )

    # Add threshold lines
    plt.axhline(-np.log10(padj_threshold), color='blue', linestyle='--', label=f'padj = {padj_threshold}')
    plt.axvline(log2fc_threshold, color='green', linestyle='--', label=f'log2FC = {log2fc_threshold}')
    plt.axvline(-log2fc_threshold, color='green', linestyle='--')

    # Add gene labels only for significant points
    if label == True:
        for _, row in data[data['Significant']].iterrows():
            plt.text(
                row['log2FoldChange'], 
                row['negLog10padj'], 
                row['Target'], 
                fontsize=8,
                color='red'
            )

    # Add labels and title
    plt.title(title)
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-log10(FDR-corrected p-value)')
    plt.legend(loc=legend_loc)
    plt.grid(grid)
    
    # Save plot as SVG if enabled
    if save_as_svg:
        plt.savefig(svg_filename, format='svg')
    
    plt.show()
