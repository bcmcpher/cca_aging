function [ fh ] = ccaPlotRankedTrendPoints(dat, cca, mod, ccf)
%[ wght, fh ] = ccaPlotRankedTrends(dat, cca, age, mod, type, ccf, nshow);
%   Estimate and plot the ranked trends of the requested data sets as a
%   line plot with the labels, indicating the highest / lowest loading
%   values.
%
%   INPUTS:
%       dat   - the preprocessed data from the CCA
%       cca   - the findings from the CCA
%       mod   - 'brain' or 'beahavior'; selects the data to plot
%       ccf   - the cc to plot data from
%       nshow - the number of loadings to show on both top and bottom.
%   OUTPUTS:
%       fh - the figure handle of the created images
%
% Copyright (c) Brent McPherson (Indiana University), 2020. All rights reserved.
%

if strcmp(mod, 'brain')
    
    % load brain data
    varLabs = dat.dat1.names;
    vals = cca.dat1.loading(:, ccf);
    valv = cca.dat1.loading_se(:, ccf);
    %valp = cca.full.fgrotAAp(:, ccf);
    
end

if strcmp(mod, 'behavior')
    
    % load behavior data
    varLabs = dat.dat2.names;
    vals = cca.dat2.loading(:, ccf);
    valv = cca.dat2.loading_se(:, ccf);
    %valp = cca.full.fgrotBBp(:, ccf);
    
end

%% extract and sort weights

% sort the values high to low
[ sw, si ] = sort(vals, 'descend');
sv = valv(si);

% set variability to 2 standard errors
sv = 2*sv;

% ADD BACK NSHOW OPTION?

% start plot
fh = figure; hold on;

% for every observation
for ii = 1:size(sw, 1)
    
    % plot 2 standard error bar line w/ point
    plot([ ii ii ], [ (sw(ii) - sv(ii)) (sw(ii) + sv(ii)) ], 'color', 'black');
    plot(ii, sw(ii), 'o', 'MarkerFaceColor', 'blue', ...
         'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'MarkerSize', 6);
    
end

% plot zero line in red
plot([ -1 size(sw, 1)+1 ], [ 0 0 ], 'red'); 

% set axis limits and assign labels
set(gca, 'Xlim', [ -1 size(sw, 1)+1 ], 'XTick', 1:size(sw, 1), ...
    'XTickLabel', varLabs(si), 'XTickLabelRotation', 35, 'TickLabelInterpreter', 'none');

hold off;

end