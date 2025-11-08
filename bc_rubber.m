%% =========================================================================
% Copyright (c) 2025 Lixue Shi and collaborators
% Licensed under the Apache License, Version 2.0 (http://www.apache.org/licenses/LICENSE-2.0)
%
% Script Name: bc_rubber.m
% Description:
% This function performs baseline correction on input spectral data using 
% a convex rubberband algorithm, inspired by the Rubberband Baseline 
%
% Input: 
% -X: a matrix of size no× nf, where each row represents a single spectrum
%  to be baseline-corrected.
%
% Output
% -varargout: [Y] or [Y,L]. Y represents the baseline-corrected spectrum,
%  while L are the baselines.
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
function varargout = bc_rubber(X)

msgstring = nargoutchk(1, 2, nargout);
if ~isempty(msgstring)
    error(msgstring);
end;


[no, nf] = size(X);

Y = zeros(no, nf);
L = zeros(no, nf);

for i = 1:no
    if nf > 0
        l = [];
        x = X(i, :);
        if length(x) > 1
            l2 = rubber(x);
        else
            l2 = [];
        end;
        l = [x(1) l2];
        
        Y(i, :) = x-l;
        L(i, :) = l;
    end;
end;



if nargout == 1
    varargout = {Y};
elseif nargout == 2
    varargout = {Y, L};
end;
    
%> @cond
%---------------------------------------------------------------------
% returns a "rubber" vector with one element less than the length of x
function y = rubber(x)

nf = length(x); % number of points

l = linspace(x(1), x(end), nf);

xflat = x-l;
[val, idx] = min(xflat);
if ~isempty(val) && val < 0
    y = [rubber(x(1:idx)), rubber(x(idx:end))];
else
    y = l(2:end);
end;
%> @endcond

