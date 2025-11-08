%% =========================================================================
% Copyright (c) 2025 Lixue Shi and collaborators
% Licensed under the Apache License, Version 2.0 (http://www.apache.org/licenses/LICENSE-2.0)
%
% Script Name: corratio.m
% Description:
% This function removes outlier points from two corresponding datasets A and B 
% based on a three-sigma rule. For each dataset, elements that deviate more 
% than ±3 × standard deviation from the mean are identified and removed from 
% both vectors to maintain pairwise correspondence. The cleaned datasets are 
% then returned as CORA and CORB, representing the outlier-filtered versions 
% of A and B respectively.
%
% Input: 
% -A,B: Ratios of IR absorbance intensities of cells at selected wavenumbers.
%
% Output
% -COR_A, COR_B: Corresponding corrected datasets after outlier removal.
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
function [COR_A,COR_B]=corratio(A,B)
avg1=mean(A);avg2=mean(B);
err1=std(A);err2=std(B);
[row1,col1]=find(A<=(avg1-3*err1));
[row2,col2]=find(A>=(avg1+3*err1));
[row3,col3]=find(B<=(avg2-3*err2));
[row4,col4]=find(B>=(avg2+3*err2));
col12=union(col1,col2);
col123=union(col12,col3);
col=union(col123,col4);
[~,n]=size(col);
k=n;
for i = 1:1:n
    d=col(k);
    A(d)=[];
    B(d)=[];
    k=k-1;
end

COR_A=A;

COR_B=B;

