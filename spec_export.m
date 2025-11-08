%% =========================================================================
% Script Name: spec_export.m
% Description:
% This function extracts and exports summed single-cell infrared spectra 
% from a 3D FTIR data cube. For each segmented cell in the input mask, it 
% identifies all pixel coordinates belonging to that cell, retrieves their 
% corresponding spectra from the data matrix, and sums them to generate 
% one integrated spectrum per cell. The function then plots both all individual 
% cell spectra and their mean spectrum for visualization.
%
% Input: 
% A: 3D FTIR data matrix (x, y, wavenumber).
% I: Cell segmentation mask obtained from CellProfiler.
% W: Vector of wavenumbers corresponding to spectral channels.
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
function CELL_SPEC = spec_export(A,I,W)
m_seg= max(I,[],'all');
CELL_SPEC=[];
for x=1:m_seg
    [row,col]=find(I==x);
    M=[row,col];
    [m,~]= size(M);
    [~,c2]=size(W);
    H=zeros(m,c2);
    for b1 = 1:m
        if A(M(b1,1))==NAN
            continue
        else
            H(b1,:)= A(M(b1,1),M(b1,2),:);
            CELL=sum(H,1);
        end
        CELL_SPEC=[CELL_SPEC; CELL]; 
    end
end
Meanspec = nanmean(CELL_SPEC);
figure
plot(W,CELL_SPEC,'linewidth',1)
xlabel('wavenumbers (cm-1)')
ylabel('Absorbance')
title('all spec')
figure
plot(W,Meanspec,'linewidth',1)
xlabel('wavenumbers (cm-1)')
ylabel('Absorbance')
title('mean spec')


    
