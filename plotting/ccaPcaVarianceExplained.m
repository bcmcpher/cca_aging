function [ fh ] = ccaPcaVarianceExplained(inp, Ncca, label)
%[ fh ] = ccaPcaVarianceExplained(inp, Ncca, label);
%   Using an input form the CCA preprocessing, this will show how much
%   variability is explained by the PCA components taken during
%   preprocessing.
%
% The inputs will change when cca structure is defined
%
% THIS IS NOT USED / REPORTED IN THE PAPER. IT MAY BE USEFUL.
%
%   INPUTS:
%       inp   - dat.dat?.ss?; which PCA to inspect
%       Ncca  - the number of significant CCA components; draws a line to
%               identify the contribution of significant components
%       label - text label for plot title
%   OUTPUTS:
%       fh - figure handle holding the plot
%
% Copyright (c) Brent McPherson (Indiana University), 2019. All rights reserved.
%

% input either: dd, ss1

% behavior / network variance explained
dat = normalize(diag(inp));

fh = figure; hold on; 
plot(1:20, dat(1:20));
plot(1:20, dat(1:20), 'o', 'MarkerEdgeColor', [ 0 0 0 ], 'MarkerFaceColor', [ 0 0 1 ]);
line([ (Ncca + .5) (Ncca + .5) ], [ 0 1 ], 'color', 'red');
set(gca, 'YLim', [ 0 1 ], 'XLim', [ 0 20 ]);
title([ label ' - Variance Explained' ]);
hold off;

end

