%% =========================================================================
% Copyright (c) 2025 Lixue Shi and collaborators
% Licensed under the Apache License, Version 2.0 (http://www.apache.org/licenses/LICENSE-2.0)
%
% Script Name: data_convert.m
% Description:
%   Convert a single FTIR imaging matrix stored in ASCII .dat format into a
%   16-bit grayscale TIFF image for visualization and downstream analysis.
%   The script reads a numeric matrix (e.g., from Cytospec export at a
%   specific wavenumber such as 1614 cm⁻¹ or 2212 cm⁻¹), optionally rotates
%   it to the correct orientation, replaces missing values (NaN) with zeros,
%   applies an intensity scaling strategy, and saves the result as a .tif file.
%
% Input:
%   - A single infrared data matrix in ASCII .dat format (e.g., Cytospec export)
%     for one selected wavenumber.
%
% Output:
%   - A 16-bit grayscale TIFF image (.tif) corresponding to the input matrix.
%
% Notes:
%   - NaN values are replaced with 0 before scaling.
%   - Intensity can be scaled by an empirical multiplier (mf) or mapped to the
%     full uint16 dynamic range if auto-stretch is enabled (see script body).
%   - Orientation can be adjusted via rot90 with a user-defined rotation count.
%   - For binary .dat (non-ASCII) exports, replace the reader with fread and
%     specify the data shape and precision accordingly.

% Corresponding paper:
%   "PRIMED: High-throughput single-cell profiling of proteome turnover dynamics"
%   Authors: Yuchen Sun+, Minqian Wei+, Lixue Shi*
%
% Dependencies:
%   - MATLAB R2023b (MathWorks, Natick, MA)
%
% Citation:
%   If you use this code, please cite the above article.
%
% =========================================================================
%%
clc;clear;close all;

%% Output image filename
write_file=strcat('Output_image','.tif');

%% Measure image size
data = load('Input_file.dat');
[n,m] = size(data);
clear data

%% Read .dat files & Basic parse
fileID = fopen('Input_file.dat','r');
f = fscanf(fileID,'%f');
fclose(fileID);
s = (size(f));
im = reshape(f,m,n);
im(isnan(im))=0;
raw_dat = rot90(im);

%% Format conversion and output
mf=1; 
% mf is a scaling factor that depends on the selected wavenumber, 
% intended to enhance the contrast of the generated TIFF files for 
% improved cell segmentation performance.
raw_image = uint16(raw_dat*mf);
imwrite(raw_image,write_file);



