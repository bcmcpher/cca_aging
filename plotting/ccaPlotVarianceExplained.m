function [ fh ] = ccaPlotVarianceExplained(cca, train)
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

% assign default plot option
if(~exist('train', 'var') || isempty(train))
    train = false;
end

%% pull the data

% pull the hold-out data
trv = cca.cca.trRv;
tpc = cca.cca.trPc;

hrv = cca.cca.hoRv;
hpc = cca.cca.hoPc;

% get the number of correlations
nval = size(hrv, 2);

%% make the plot

fh = figure; 

% top panel
subplot(2, 1, 1); hold on;

% for every correlation
for ii = 1:nval
    
    % plot the hold out point
    plot(ii, hrv(ii), 'o', 'MarkerFaceColor', [ 0.1 0.2 0.7 ], ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'MarkerSize', 6);
    
    if train
        plot(ii, trv(ii), 'o', 'MarkerFaceColor', [ 0.2 0.4 1 ], ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'MarkerSize', 6);
    end
    
end

set(gca, 'XTick', []);
title([ 'Variability Explained by CCA Axes: ' num2str(sum(hrv)) ]);
xlabel('Canonical Correlations');
ylabel('Variability Explained');
hold off

% bottom panel
subplot(2, 1, 2); hold on;

% for every correlation
for ii = 1:nval
    
    % plot the hold out percent
    plot(ii, hpc(ii), 'o', 'MarkerFaceColor', [ 0.7 0.2 0.1 ], ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'MarkerSize', 6);
    
    if train
        plot(ii, tpc(ii), 'o', 'MarkerFaceColor', [ 1 0.4 0.2 ], ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'MarkerSize', 6);
    end
    
end

set(gca, 'XTick', []);
title('Percent Explained by CCA Axes');
xlabel('Canonical Correlations');
ylabel('Percent');
hold off

end
