%% =========================================================================
% Copyright (c) 2025 Lixue Shi and collaborators
% Licensed under the Apache License, Version 2.0 (http://www.apache.org/licenses/LICENSE-2.0)

% Script Name: IR_ratio_anal.m
% Description:
% This script processes segmented FTIR images to quantify single-cell protein 
% synthesis and degradation activities. It reads IR intensity maps (e.g., 
% 1614 cm⁻¹, 1645 cm⁻¹, 2212 cm⁻¹) and corresponding segmentation masks, 
% calculates average absorption intensities for each cell region, performs spectral 
% unmixing to separate ¹²C- and ¹³C-amide components, and computes derived 
% ratios such as  ¹²C-to-total amide ratios. Outlier cells are then removed using 
% filters (corratio) for reliable downstream quantitative analysis.
% The absorbance of cells at 1614, 1645, and 2212 cm⁻¹ can be obtained from 
% the 194# and 204# columns of the single-cell spectral data file spec_fp, and 
% the 52# column of spec_cd.
%
% Input: 
% segment: Cell segmentation mask obtained from CellProfiler.
% C13_raw, C12_raw, CD_raw: Infrared intensity images at selected wavenumbers.
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
%% Import data
segment=imread('Cell_segment.tiff');
C13_raw=imread('Cell_C13amide.tif');
C12_raw=imread('Cell_C12amide.tif');
CD_raw=imread('Cell_CD.tif');

%% Single-cell data extraction
N=max(max(segment));
for i = 1:N
    mask = segment;
    mask(mask~=i)=0;
    mask(mask>0)=1;
    C13_sum(i)=sum(sum(C13_raw.*mask));
    C12_sum(i)=sum(sum(C12_raw.*mask));
    CD_sum(i)=sum(sum(CD_raw.*mask));
    Area(i)=sum(sum(mask)); 
end

C13_avg=C13_sum./Area;
C12_avg=C12_sum./Area;
CD_avg=CD_sum./Area;

%% Unmixing of ¹²C- and ¹³C-amide peaks
% Unmixing for FTIR datasets
C12_unmix = C12_avg*1.281+C13_avg*(-0.548);
C13_unmix = C12_avg*(-0.657)+C13_avg*1.281;

% Unmixing for DFIR datasets 
C12_unmix=C12_avg*(1.094)+C13_avg*(-0.413);
C13_unmix=C13_avg*1.094+C12_avg*(-0.250);

%% Ratios calculation and filtering
sum=C13_unmix+C12_unmix;
ratio_CAA=C12_unmix./sum;
ratio_DAA=100*CD_avg./sum;
[ratio_CAA_f,ratio_DAA_f]=corratio(ratio_CAA,ratio_DAA);

%% Data visualization (optional)
figure
scatter(ratio_CAA_f,ratio_DAA_f);
