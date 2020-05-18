function [ fh ] = ccaPlotCorrs(cca, age, nse, nperm)
%[ fh ] = ccaPlotCorrs(cca);
%
%   Plot the cross-validated correlation and error bars for every estimated
%   canonical factor.
%
%   INPUTS:
%       cca - a fit cca structure
%       age - the vector of ages for every participant
%       nse - how many unit of standard error to plot around points
%       nperm - the number of permutations to estimate correlation with age
%   OUTPUTS:
%       fh - figure handle of the plot produced
%
% Copyright (c) Brent McPherson (Indiana University), 2020. All rights reserved.
%

%% pull the data

if(~exist('nse', 'var') || isempty(nse))
    nse = 3;
end

if(~exist('nperm', 'var') || isempty(nperm))
    nperm = 1000;
end

% pull the hold-out data
hvl = cca.cca.hocorrs;
hse = cca.cca.hocorrs_se;

% get the number of correlations
nval = size(hvl, 2);

%% pull the null distribution if it exists

if isfield(cca.cca, 'hocorrs_null')
    
    % preallocate null shape
    lbub = nan(nval, 2);
    
    % pull the 5th / 95th percentiles
    for ii = 1:nval
        lbub(ii, :) = prctile(squeeze(cca.cca.hocorrs_null(1, ii, :)), [ 5 95 ]);
    end
    
end

%% make the plot

figure; 

% top panel
subplot(2, 1, 1); hold on;

% if the null exists, shade in the background
if isfield(cca.cca, 'hocorrs_null')
    nll = fill([ 1:nval, fliplr(1:nval) ], [ lbub(:, 1)', lbub(:, 2)' ], [ .75 .75 .75 ]);
    nll.FaceAlpha = 0.5;
end

% for every correlation
for ii = 1:nval

    % pull the points values / units of standard error
    hpt = hvl(ii);
    heb = nse*hse(ii);
    
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

disp('Creating estimate of canonical axes with age...');

% pull the null distribution for every age
for ii = 1:nval
    [ hpt(ii), heb(ii), ~, hnl ] = ccaLinRegCorr(cca, ii, age, nperm);
    anl(ii, :) = prctile(squeeze(hnl(:, 1)), [ 5 95 ]);
end

% fill in null distribution shape from permutation test
anll = fill([ 1:nval, fliplr(1:nval) ], [ anl(:, 1)', anl(:, 2)' ], [ .75 .75 .75 ]);
anll.FaceAlpha = 0.5;

% for every correlation
for ii = 1:nval

    % pull the points values
    %[ hpt, heb, hpv, hnl ] = ccaLinRegCorr(cca, ii, age, 1000);
    
    % units of standard error
    hebp = nse*heb(ii);
    
    % plot the errorbars with the point
    plot([ ii, ii ], [ hpt(ii)+hebp hpt(ii)-hebp ], 'color', 'black');
    plot(ii, hpt(ii), 'o', 'MarkerFaceColor', [ 0.7 0.2 0.1 ], ...
         'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'MarkerSize', 6);

end

set(gca, 'XTick', []);
title('Correlation Between Age and Canonical Axis');
xlabel('Canonical Correlations');
ylabel('Correlation');
hold off

end
