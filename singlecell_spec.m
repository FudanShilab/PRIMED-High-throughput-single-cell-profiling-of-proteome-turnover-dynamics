% Script Name: singlecell_spec.m
% Description: This script extracts single-cell infrared spectra from 
% segmented FTIR data in the fingerprint (1000–1800 cm⁻¹) and C–D 
% stretching (2050–2235 cm⁻¹) regions. It loads spectral cubes and 
% corresponding segmentation masks, reconstructs wavenumber axes, 
% and outputs for each cell the summed spectrum, obtained by 
% integrating the spectra of all pixels within the cell region for 
% subsequent metabolic analysis.

%
% Input: 
% segment: Cell segmentation mask obtained from CellProfiler.
% data_fp, data_cd: Infrared absorption spectral data in the ranges 
%   of 1000–1800 cm⁻¹ and 2050–2235 cm⁻¹
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
% License:For academic use only.
% =========================================================================
%%
clc;close all;clear
%% Import data & Basic parse
load('Celldata_Fingerprint_region.mat');
data_fp=C;
lwn_fp = str2double(Minfo.File.LWN);
uwn_fp= str2double(Minfo.File.UWN);
nod_fp = str2double(Minfo.File.NOD);
wvs_fp = str2double(Minfo.File.WVS);
clear C Minfo
wave_fp=linspace(lwn_fp,uwn_fp,nod_fp);

load('Celldata_CD_region.mat');
data_cd=C;
lwn_cd=str2double(Minfo.File.LWN);
uwn_cd=str2double(Minfo.File.UWN);
nod_cd=str2double(Minfo.File.NOD);
clear C Minfo
wave_cd=linspace(lwn_cd,uwn_cd,nod_cd);

%% Single-cell spectra extraction
l=imread('Cell_segment.tiff');
l=rot90(l,3);
spec_fp=spec_export(data_fp,l,wave_fp);
spec_cd=spec_export(data_cd,l,wave_cd);

%% Data visualization (optional)
figure
plot(wave_fp,spec_fp);

figure
plot(wave_cd,spec_cd);

