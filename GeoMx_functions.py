import pandas as pd
import numpy as np
import scipy as sp
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import os
import math
import anndata as ad

def GeoMx_Reader(file_path):
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
    metadata_columns = ["SegmentDisplayName", "Age", "Gender", "Phenotype", "Treatment", "AreaClass","Entity", "Patient"]
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
