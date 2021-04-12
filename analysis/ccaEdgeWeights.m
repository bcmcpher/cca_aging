function [ bmat ] = ccaEdgeWeights(grotAAd, ccf)
%[ bmat ] = ccaEdgeWeights(grotAAd, ccf);
%   Convert a set of edge-based loadings back into an adjacency matrix of
%   individual edge contributions.
%   
%   INPUTS:
%       grotAAd - dat.dat?.loading; deconvolved loadings from an edge-based CCA inputs
%       ccf     - the canonical correlation to use for reconstructing the loadings
%   OUTPUTS:
%       bmat - the edge-wise contribution matrix that was requested
%
% Copyright (c) Brent McPherson (Indiana University), 2019. All rights reserved.
%

if size(grotAAd, 1) > 377
    
    % preallocate output
    bmat = zeros(376);
    
    % create indices for unique values
    [ x, y ] = find(triu(ones(376), 1));
    
    % create matrix filled with contribution weights of edges
    for ii = 1:size(x, 1)
        bmat(x(ii), y(ii)) = grotAAd(ii, ccf);
        bmat(y(ii), x(ii)) = grotAAd(ii, ccf);
    end
    
    % plot the values to see
    %plotAdjacencyMatrix(bmat);

end

