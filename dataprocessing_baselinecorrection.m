%% =========================================================================
% Copyright (c) 2025 Lixue Shi and collaborators
% Licensed under the Apache License, Version 2.0 (http://www.apache.org/licenses/LICENSE-2.0)
%
% Script Name: dataprocessing_baselinecorrection.m
% Description:
% This script performs baseline correction on preprocessed FTIR hyperspectral data.
% It loads the processed spectral cube and metadata, reshapes the data into a 
% two-dimensional matrix, and applies the bc_rubber algorithm to remove baseline 
% distortions. The corrected spectra are then reconstructed into their original 3D 
% format and saved as a new .mat file for subsequent quantitative and comparative 
% spectral analysis.
%
% Input: 
% A 4-D .mat file exported from Cytospec, where the first three dimensions 
% represent the x-coordinate, y-coordinate, and wavenumber, and the fourth 
% dimension indicates whether the data have been preprocessed.
%
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
%% Improt data
load('unprocessed_data.mat');
data =C(:,:,:,2); 
lwn = str2double(Minfo.File.LWN);
uwn = str2double(Minfo.File.UWN);
nod = str2double(Minfo.File.NOD);
wvs = str2double(Minfo.File.WVS);
clear C

%% Reshape the spec
[a,b,c] = size(data);
pspec=reshape(data,a*b,c);
wave = linspace(lwn,uwn,nod);

%% rubber-band baseline correction
pspec_bc=bc_rubber(pspec);

%% Save data
C=reshape(pspec_bc,a,b,c);
save('processed_data.mat','C','Minfo');




