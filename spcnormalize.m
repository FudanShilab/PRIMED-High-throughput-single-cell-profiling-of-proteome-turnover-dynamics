%% =========================================================================
% Script Name: spcnormalize.m
% Description:
% This function performs row-wise min–max normalization on input spectral 
% data to scale each spectrum to its maximum intensity. By dividing each row 
% (spectrum) by its own maximum value, the function removes absolute intensity 
% differences while preserving relative spectral shapes.
%
% Input: 
% spectrainput: a matrix where each row represents one spectrum and each 
%    column corresponds to a wavenumber.
% spectra: The normalized spectral matrix, where each spectrum is scaled 
%    such that its maximum intensity equals 1.
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

function [ spectra ] = spcnormalize(spectrainput)
% min-max normalize
spectra = spectrainput;
spectra = spectra ./ max(spectra,[],2);
end


