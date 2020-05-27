function [ n, fh ] = plotVarHist(vals, nbin, name)
%[ n, fh ] = plotVarHist(val, nbin, name);
%   Plot a histogram as a density line for a provided variable. Also returns the
%   count for the requested number of bins.
%
%   INPUTS:
%       vals - subj x 1 vector of data to compute the histogram on
%       nbin - the number of bins to split the data into
%       name - a text field of the name of the variable for the plot title
%   OUTPUTS:
%       n  - a nbin x 1 vector of the count of vals
%       fh - the figure handle of the generated plot
%
% Copyright (c) Brent McPherson (Indiana University), 2019. All rights reserved.
%

% pull data

% plot histogram
fh = figure; hold on;

% use two-output version of hist to get values and normalize to unit area
[ n, x ] = histcounts(vals, nbin); 
%n_normalized = n / numel(val) / (x(2) - x(1));

% plot the line
%plot(1:nbin, n_normalized, 'r');
plot(1:nbin, n, 'blue');
title([ 'Histogram of ' name ], 'Interpreter', 'None');
set(gca, 'XTick', 1:nbin, 'XTickLabels', round(x, 2));
xtickangle(90);
xlabel('Value');
ylabel('Count');

end

