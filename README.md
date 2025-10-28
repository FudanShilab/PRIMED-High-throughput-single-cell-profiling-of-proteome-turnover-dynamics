# PRIMED-High-throughput-single-cell-profiling-of-proteome-turnover-dynamics
*Authors: Yuchen Sun+, Minqian Wei+, Lixue Shi**
-This package provides MATLAB and Python scripts for baseline correction, data format conversion, protein turnover rate quantification, outlier removal, spectral normalization, single-cell spectrum extraction, and full-spectrum dimensionality reduction analysis (UMAP and HCA) described in the paper.

---

## MATLAB Scripts
- **`bc_rubber.m`** – Baseline correction (rubber-band) for spectral data.  
- **`corratio.m`** - Outliers removal by three-sigma criterion.
- **`data_convert.m`** - Converts .dat files into .tif images for visualization and segmentation.
- **`dataprocessing_baselinecorrection.m`** - Baseline correction on preprocessed FTIR spectral cubes using the bc_rubber algorithm.
- **`IR_ratio_anal.m`** - Single-cell IR peak ratios calculation to quantify protein synthesis and degradation activities (¹²C-Amide/(¹³C-Amide + ¹²C-Amide) & CD/((¹³C-Amide + ¹²C-Amide)).
- **`singlecell_spec.m`** - Single-cell spectra extraction from FTIR data using spec_export algorithm.
- **`spcnormalize.m`** – Normalization of single-cell spectra.  
- **`spec_export.m`** - Exports infrared spectra for each segmented cell.

## Python Scripts
- **`Full_spectra_umap_analysis.py`** – Performs dimensionality reduction and visualization of single-cell FTIR spectral data.  
- **`HCA_Dendrogram.py`** — Performs hierarchical clustering analysis (HCA) on averaged FTIR spectra grouped by cell label.  
  Generates and saves a dendrogram (PDF/PNG) visualizing spectral similarity across experimental groups.


---

## Data Format
- **IR Image (MATLAB):**  
  `Infrared data matrices (.dat) obtained from Cytospec at specific wavenumbers.`
- **IR Spectrum (MATLAB):**
  `4-D .mat file exported from Cytospec, where the first three dimensions represent the x-coordinate, y-coordinate, and wavenumber, and the fourth dimension indicates whether the data have been preprocessed.`
- **UMAP Analysis (Python):**  
  `Multiple spectral .csv files from path/to/your/data, each representing single-cell spectra from one experimental group.` 
- **HCA Dendrogram Analysis (Python):**  
  `Columns 1–2: cell identifiers and group labels (id, label). Columns 3–end: spectral intensities across wavenumbers. Row 1: Corresponding wavenumber of cell spectra.` 

---

## Dependencies
- MATLAB R2023b or later.  
- Python 3.12, with `numpy`, `pandas`, `matplotlib`, `scikit-learn`, `scipy`, `umap-learn`.  

---

## Citation
If you use these scripts, please cite:  
**"PRIMED: High-throughput single-cell profiling of proteome turnover dynamics"**
