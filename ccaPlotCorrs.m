function [ fh ] = ccaPlotCorrs(cca, age)
%[ fh ] = ccaPlotCorrs(cca);
%
%   Plot the cross-validated correlation and error bars for every estimated
%   canonical factor.
%
%   INPUTS:
%       cca - a fit cca structure
%       age - the vector of ages for every participant
%   OUTPUTS:
%       fh - figure handle of the plot produced
%
% Copyright (c) Brent McPherson (Indiana University), 2020. All rights reserved.
%

%% pull the data

% pull the hold-out data
hvl = cca.cca.hocorrs;
hse = cca.cca.hocorrs_se;

% get the number of correlations
nval = size(hvl, 2);

%% make the plot

figure; 

% top panel
subplot(2, 1, 1); hold on;

% for every correlation
for ii = 1:nval

    % pull the points values
    hpt = hvl(ii);
    heb = 3*hse(ii);
    
    % plot the errorbars with the point
    plot([ ii, ii ], [ hpt+heb hpt-heb ], 'color', 'black');
    plot(ii, hpt, 'o', 'MarkerFaceColor', [ 0.1 0.2 0.7 ], ...
         'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'MarkerSize', 6);

end

set(gca, 'XTick', []);
title('Correlation Between Brain and Behavior');
xlabel('Canonical Correlations');
ylabel('Correlation');
hold off

% bottom panel
subplot(2, 1, 2); hold on;

% for every correlation
for ii = 1:nval

    % pull the points values
    [ hpt, heb ] = ccaLinRegCorr(cca, ii, age, 1000);
    
    % 3 standard errors
    heb = 3*heb;
    
    % plot the errorbars with the point
    plot([ ii, ii ], [ hpt+heb hpt-heb ], 'color', 'black');
    plot(ii, hpt, 'o', 'MarkerFaceColor', [ 0.7 0.2 0.1 ], ...
         'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'MarkerSize', 6);

end

set(gca, 'XTick', []);
title('Correlation Between Age and Canonical Axis');
xlabel('Canonical Correlations');
ylabel('Correlation');
hold off

end
