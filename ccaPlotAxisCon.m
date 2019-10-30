function [ fh ] = ccaPlotAxisCon(cca, ccf, nage, cmap)
%[ fh ] = ccaPlotAxisCon(grotU, grotV, ccf, nage, cmap);
%   Plot the CCA of a specific factor, using a color mapped age value for
%   each point.
%
%   The inputs of this data fxn will change when CCA object is defined.
%
%   INPUTS:
%       cca - CCA analysis structure containing loadings
%       %   grotU - first cca scaled dataset
%       %   grotV - second cca scaled dataset
%       ccf   - the cc to plot
%       nage  - scaled age vector to correctly assign colors
%       cmap  - color map scaled to align w/ nage for plotting
%   OUTPUTS:
%       fh - figure handle of requested plot
%
% Copyright (c) Brent McPherson (Indiana University), 2019. All rights reserved.
%

grotU = cca.dat1.factor;
grotV = cca.dat2.factor;

% plot data
fh = figure('Position', [ 450 500 750 675 ]); hold on;

% for every observation
for ii = 1:size(grotU, 1)
    
    % plot the CCA feature, coloring the point be age
    plot(grotU(ii, ccf), grotV(ii, ccf), 'o', 'MarkerFaceColor', cmap(nage(ii), :), ...
         'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'MarkerSize', 6);
     
end

% final plot clean up
ccv = corr(grotU(:, ccf), grotV(:, ccf));
%ccv = grotRv(ccf);
title(['Canonical Correlation of Principal Factor ' num2str(ccf) ': ' num2str(ccv)]);
xlabel('Network Values');
ylabel('Behavioral Scores');
axis equal; axis square;
cb = colorbar; caxis([ 18 88 ]);
ylabel(cb, 'Age');
hold off;

end

