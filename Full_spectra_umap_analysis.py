#=========================================================================
# Copyright (c) 2025 Lixue Shi and collaborators
# Licensed under the Apache License, Version 2.0 (http://www.apache.org/licenses/LICENSE-2.0)

# Script Name: Dimensionality reduction and visualization of FTIR single-cell spectra.py

# Description:This script performs dimensionality reduction and visualization
# of single-cell FTIR spectral data. It aligns spectra from multiple .csv files to a common
# wavelength reference, applies normalization and standardization, and projects them into
# low-dimensional space for comparative analysis. The script outputs 2D/3D scatter plots and
# embedding results for downstream interpretation.

# Input Data Format:
# Directory: path/to/your/data containing multiple spectral .csv files, each representing
# single-cell spectra from one experimental group.
# Wavenumber reference: wave_total_WT.csv, where each column corresponds to one input .csv file
# and defines its wavelength axis.
# Each .csv file should contain:
#  - Columns → individual cell spectra
#  - Rows → spectral intensities corresponding to specific wavenumbers
# The script automatically aligns all spectra to the reference wavelengths before performing normalization,
# standardization, and dimensionality reduction.

# Corresponding paper:
#  "PRIMED: High-throughput single-cell profiling of proteome turnover dynamics"
#   Authors: Yuchen Sun+, Minqian Wei+, Lixue Shi*
#
# Dependencies:
#   - Python 3.12 (Python Software Foundation)
#   - Required packages: numpy, pandas, matplotlib, scikit-learn, scipy, umap-learn

# Citation:
#   If you use this code, please cite the above article.
#
# Copyright (c) 2025 Lixue Shi and collaborators
# Licensed under the Apache License, Version 2.0 (http://www.apache.org/licenses/LICENSE-2.0)
# =========================================================================

#   parameter for analysis
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
import os
import pandas as pd
import numpy as np
from scipy.interpolate import CubicSpline
from sklearn.preprocessing import StandardScaler
import umap as aligned_umap
from matplotlib.widgets import Button

base_path = 'path/to/your/data'


range_selected = [[1000,2235]]

#   parameter for data preprocessing
normalization = True
standardization = True
with_labels = True

#   parameter for color_list
color_list = [
    [135 / 256, 206 / 256, 235 / 256],
    [255 / 256, 105 / 256, 180 / 256],
    [255 / 256, 215 / 256, 0 / 256],
    [154 / 256, 205 / 256, 50 / 256]
]

my_cmap = LinearSegmentedColormap.from_list('rain', color_list)
plt.colormaps.register(cmap=my_cmap)
cmap = plt.get_cmap("rain")

#   parameter for figure
figure_width = 6
figure_height = 6
axis_linewidth = 1
axis_linelength = 8
scatter_size = 10
scatter_alpha = 0.6
scatter_linewidth = 0
scatter_linecolor = "black"
font_size = 1
jitter = 0
ids_show = False
colorbar = False

#   parameter for ellipse
ellipse_95 = False
fill = False
transparency = 1
ellipse_linewidth = 2
confidence = 0.95  # confidence of the ellipse: 0.90; 0.95; 0.99.


# Step 1: Load the data files and extract the wavelengths
data_path = os.path.join(base_path, 'data_filename')
all_files = os.listdir(data_path)
file_names = [file for file in all_files if file.endswith('.csv')]
print("The list of data files is: ", file_names)
wavenum_df = pd.read_csv('wavenumber_file.csv')

# Step 2: Map wavelengths to files using column names
wavelengths = {file_name: wavenum_df[file_name.replace('.csv', '')].values for file_name in file_names}
for file_name, w_array in wavelengths.items():
    cleaned_array = w_array[~np.isnan(w_array)]
    wavelengths[file_name] = cleaned_array

# Step 3: Find the reference file with the most wavelength points
reference_file = None
max_points = 0
for file_name, w_array in wavelengths.items():
    num_points = len(w_array)
    if num_points > max_points:
        max_points = num_points
        reference_file = file_name
print(f"The reference file is {reference_file} with {max_points} wavelength points.")
reference_wavelengths = wavelengths[reference_file]

data = []
label = []

def align_spectrum_to_reference(current_wavelengths, spectrum, reference_wavelengths):
    # print(len(current_wavelengths), len(spectrum))
    interpolator = CubicSpline(current_wavelengths, spectrum)
    aligned_spectrum = interpolator(reference_wavelengths)
    return aligned_spectrum

for i, file_name in enumerate(file_names):
    df = pd.read_csv(os.path.join(data_path, file_name), header=None)
    print(f"Loading file {file_name} with shape of {df.shape} and label of {i}.")
    count = 0
    # Iterate through each row in the DataFrame
    for index, row in df.iterrows():
        spectrum = row.values
        # Check if the spectrum needs to be aligned
        if file_name != reference_file:
            # Assuming the wavelengths for the current file are stored in a variable named `current_wavelengths`
            current_wavelengths = wavelengths[file_name]
            aligned_spectrum = align_spectrum_to_reference(current_wavelengths, spectrum, reference_wavelengths)
            data.append(aligned_spectrum)
        else:
            # If it's the reference file, append the spectrum directly
            data.append(spectrum)
        label.append(file_name.replace('.csv', ''))
        count += 1

data_df = pd.DataFrame(data)
if normalization:
    data_df = data_df.apply(lambda x: x/x.max(), axis=1)
labels = np.array(label)

print("Data preprocessing finished.")
print(f"The shape of data: {data_df.shape}")

data_df = np.array(data_df.values)
if standardization:
    data_df = StandardScaler().fit_transform(data_df)

# Filter data_df to keep only the columns corresponding to the specified wavelength ranges
indices = set()
for wavelength_range in range_selected:
    min_wavelength, max_wavelength = wavelength_range
    range_indices = [i for i, wavelength in enumerate(reference_wavelengths) if min_wavelength <= wavelength <= max_wavelength]
    indices.update(range_indices)

indices = sorted(indices)
data_df = data_df[:, indices]
print(data_df.shape)

#   extract the subscript for each class and set the color for each class
label_class = list(set(labels))
indexes_class = [np.where(np.array(labels) == buffer) for buffer in label_class]
class_number = len(label_class)
total_colors = len(color_list)  # total number of colors in cmap
colors_buffer = [cmap(i / (total_colors - 1)) for i in range(min(class_number, total_colors))]
print("The mapping between label and color is: ", dict(zip(file_names, colors_buffer)))


reducer = aligned_umap.UMAP(random_state=1, n_components=3)
analysis_result = reducer.fit_transform(data_df)

fig, ax = plt.subplots(figsize=(figure_width, figure_height))
ax.spines['top'].set_linewidth(axis_linewidth)
ax.spines['right'].set_linewidth(axis_linewidth)
ax.spines['bottom'].set_linewidth(axis_linewidth)
ax.spines['left'].set_linewidth(axis_linewidth)
plt.tick_params(labelsize=15, width=axis_linewidth, length=axis_linelength)

for i in range(class_number):
    data = analysis_result[indexes_class[i][0], :].T
    x, y ,z = data
    x = x + jitter * np.random.rand(*x.shape)
    y = y + jitter * np.random.rand(*y.shape)
    z = z + jitter * np.random.rand(*z.shape)
    plt.scatter(x, y, s=scatter_size, color=colors_buffer[i], alpha=scatter_alpha, linewidths=scatter_linewidth, edgecolors=scatter_linecolor)
    if ellipse_95:
        make_ellipses(data, ax, confidence=confidence, color=colors_buffer[i], alpha=transparency, eigv=False,
                      fill=fill, linewidth=ellipse_linewidth)


# Create color mapping for different classes
    unique_labels = np.unique(labels)
    colors = [cmap(i / (len(unique_labels) - 1)) for i in range(len(unique_labels))]
    color_dict = dict(zip(unique_labels, colors))

    # Print color mapping information
    print("Color and Label Correspondence:")
    for label, color in color_dict.items():
        print(f"Label: {label}, Color: {color}")

def create_plot():
    """
    Create 2D scatter plot of the analysis results.

    This function creates a 2D visualization of the dimensionality reduction results,
    with optional confidence ellipses and labels for different classes.

    Returns
    -------
    tuple
        Figure and axes objects for the created plot
    """
    fig, ax = plt.subplots(figsize=(figure_width, figure_height))

    # Configure axis appearance
    ax.spines['top'].set_linewidth(axis_linewidth)
    ax.spines['right'].set_linewidth(axis_linewidth)
    ax.spines['bottom'].set_linewidth(axis_linewidth)
    ax.spines['left'].set_linewidth(axis_linewidth)
    plt.tick_params(labelsize=5, width=axis_linewidth, length=axis_linelength)

    # Plot data points for each class
    for label in unique_labels:
        mask = labels == label
        x = analysis_result[mask, 0] + jitter * np.random.randn(np.sum(mask))
        y = analysis_result[mask, 1] + jitter * np.random.randn(np.sum(mask))

        plt.scatter(x, y,
                    s=scatter_size,
                    color=color_dict[label],
                    alpha=scatter_alpha,
                    linewidths=scatter_linewidth,
                    edgecolors=scatter_linecolor,
                    label=label if with_labels else None)

        # Add confidence ellipses if requested
        if ellipse_95:
            make_ellipses(analysis_result[mask].T,
                          ax,
                          confidence=confidence,
                          color=color_dict[label],
                          alpha=transparency,
                          eigv=False,
                          fill=fill,
                          linewidth=ellipse_linewidth)

    # Add legend if labels are enabled
    if with_labels:
        plt.legend()

    # Set plot title and labels
    plt.title(f"{analysis_type} Analysis")
    plt.xlabel(f"{analysis_type} 1")
    plt.ylabel(f"{analysis_type} 2")

    return fig, ax


def create_3d_plot():
    """
    Create 3D scatter plot of the analysis results.

    This function creates a 3D visualization of the dimensionality reduction results,
    with optional labels for different classes.

    Returns
    -------
    tuple
        Figure and axes objects for the created plot
    """
    fig = plt.figure(figsize=(figure_width, figure_height))
    ax = fig.add_subplot(111, projection='3d')
    plt.tick_params(labelsize=5, width=axis_linewidth, length=axis_linelength)

    # Plot data points for each class
    for label in unique_labels:
        mask = labels == label
        x = analysis_result[mask, 0] + jitter * np.random.randn(np.sum(mask))
        y = analysis_result[mask, 1] + jitter * np.random.randn(np.sum(mask))
        z = analysis_result[mask, 2] + jitter * np.random.randn(np.sum(mask))

        ax.scatter(x, y, z,
                   s=scatter_size,
                   color=color_dict[label],
                   alpha=scatter_alpha,
                   linewidths=scatter_linewidth,
                   edgecolors=scatter_linecolor,
                   label=label if with_labels else None)

    # Add legend if labels are enabled
    if with_labels:
        ax.legend()

    # Set plot title and labels
    ax.set_title(f"{analysis_type} Analysis")
    ax.set_xlabel(f"{analysis_type} 1")
    ax.set_ylabel(f"{analysis_type} 2")
    ax.set_zlabel(f"{analysis_type} 3")

    return fig, ax


# Create and display plots
fig, ax = create_plot()
plt.show()

fig, ax = create_3d_plot()
plt.show()


# Save analysis results to CSV
results_df = pd.DataFrame({
    id_col: id,
    label_col: labels,
    f"{analysis_type}-1": analysis_result[:, 0],
    f"{analysis_type}-2": analysis_result[:, 1],
    f"{analysis_type}-3": analysis_result[:, 2]
})
results_df.to_csv(output_results.format(analysis_type), index=True)

print(f"Analysis complete. Results saved to {output_results.format(analysis_type)}")
def resize_and_save(event):
    ax_button.set_visible(False)
    fig.set_size_inches(figure_width, figure_height)
    plt.savefig(analysis_type + '.jpg', dpi=1200)
    ax_button.set_visible(True)
    print('Image resized and saved as ' + analysis_type + '.jpg')





