#=========================================================================
# Script Name: Hierarchical clustering analysis of averaged single-cell FTIR spectra.py

# Description:This script performs hierarchical clustering analysis (HCA) on averaged
# single-cell FTIR spectral data. It reads a combined spectral dataset, computes the mean
# spectrum for each cell label (group), and applies Ward’s method with Euclidean distance
# to construct a hierarchical cluster tree. The resulting dendrogram visually represents
# spectral similarity among different cell populations and is exported as a high-resolution image.

# Input Data Format:
# File: cellspec_for_HCA.csv
# - Columns 1–2: cell identifiers and group labels (id, label).
# - Columns 3–end: spectral intensities across wavenumbers.
# - Row 1: Corresponding wavenumber of cell spectra.
# The script groups spectra by label and computes class-averaged spectra for clustering.

# Corresponding paper:
#  "PRIMED: High-throughput single-cell profiling of proteome turnover dynamics"
#   Authors: Yuchen Sun+, Minqian Wei+, Lixue Shi*
#
# Dependencies:
#   - Python 3.12 (Python Software Foundation)
#   - Required packages: pandas, matplotlib, scipy

# Citation:
#   If you use this code, please cite the above article.
#
# Copyright (c) 2025 Lixue Shi and collaborators
# Licensed under the Apache License, Version 2.0 (http://www.apache.org/licenses/LICENSE-2.0)
# =========================================================================
import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib.pyplot as plt

# ==== Parameter Settings ====
fig_size = (8, 8)                 # Figure size (width, height)
fig_dpi = 300                     # Figure resolution (DPI)
output_format = 'pdf'             # Output format: 'png' or 'pdf'

# ==== Load Data ====
filename = "cellspec_for_HCA.csv"
data_df = pd.read_csv(filename)

# ==== Compute the mean spectrum for each label ====
# Group the spectral data (from the 3rd column onward) by the 'label' column and take the mean
mean_data_df = data_df.iloc[:, 2:].groupby(data_df['label']).mean()
mean_labels = mean_data_df.index.values
mean_spectral_data = mean_data_df.values

# ==== Perform Hierarchical Clustering ====
# Using Ward’s method with Euclidean distance
linkage_matrix = linkage(mean_spectral_data, method='ward', metric='euclidean')

# ==== Create Figure ====
fig, ax = plt.subplots(figsize=fig_size)

# ==== Plot Dendrogram ====
dendrogram(linkage_matrix,
           labels=mean_labels,         # Use class labels as leaf names
           leaf_rotation=45,           # Rotate leaf labels for readability
           leaf_font_size=10,          # Font size of leaf labels
           color_threshold=0.7,        # Threshold for color segmentation
           above_threshold_color='gray',
           ax=ax)

# ==== Figure Formatting ====
ax.set_xlabel('Categories', fontsize=12, weight='bold')
ax.set_ylabel('Distance', fontsize=12, weight='bold')

# ==== Layout Optimization ====
plt.tight_layout()

# ==== Save Figure ====
output_filename = f"dendrogram.{output_format}"
plt.savefig(output_filename, dpi=fig_dpi, bbox_inches='tight')
print(f"Dendrogram saved as {output_filename}")

# ==== Close Plot ====
plt.close()

