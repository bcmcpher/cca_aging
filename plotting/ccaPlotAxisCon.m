function [ fh ] = ccaPlotAxisCon(cca, ccf, nage, cmap, ebars, full)
%[ fh ] = ccaPlotAxisCon(grotU, grotV, ccf, nage, cmap);
%   Plot the CCA of a specific factor, using a color mapped age value for
%   each point.
%
%   INPUTS:
%       cca - CCA analysis structure containing loadings
%       ccf   - the cc to plot
%       nage  - scaled age vector to correctly assign colors
%       cmap  - color map scaled to align w/ nage for plotting
%       full  - use the full set of observations w/o cross-validation. This
%               will override a request for errorbars b/c there is not
%               variability estimate in a full model.
%   OUTPUTS:
%       fh - figure handle of requested plot
%
% Copyright (c) Brent McPherson (Indiana University), 2019. All rights reserved.
%

if(~exist('ebars', 'var') || isempty(ebars))
    ebars = false;
end

if(~exist('full', 'var') || isempty(full))
    full = false;
end

disp('Determining correlation of main axis with age...');

if full
    ebars = false;
    grotU = cca.full.grotU;
    grotV = cca.full.grotV;
    ccv = corr(grotU(:, ccf), grotV(:, ccf));
    [ ageR, ageS ] = ccaLinRegCorr(cca, ccf, nage, 1000, true);
else
    grotU = cca.dat1.factor;
    grotV = cca.dat2.factor;
    ccv = cca.cca.hocorrs(ccf);
    ccs = cca.cca.hocorrs_se(ccf);
    [ ageR, ageS ] = ccaLinRegCorr(cca, ccf, nage, 1000);
end

if isfield(cca.dat1, 'factor_se')
    grotUs = cca.dat1.factor_se;
    grotVs = cca.dat2.factor_se;
else
    grotUs = zeros(size(cca.dat1.factor));
    grotVs = zeros(size(cca.dat2.factor));
end

% plot data
fh = figure('Position', [ 450 500 750 675 ]); hold on;

% for every observation, put all the error bars in the background
for ii = 1:size(grotU, 1)
    
    % grab the centers and variability of the points
    xc = grotU(ii, ccf);
    yc = grotV(ii, ccf);
    xv = 2 * grotUs(ii, ccf);
    yv = 2 * grotVs(ii, ccf);
    
    % plot the error bars around the points
    if ebars
        plot([ (xc-xv) (xc+xv) ], [ yc yc ], 'color', 'black');
        plot([ xc xc ], [ (yc-yv) (yc+yv) ], 'color', 'black');
    end
    
end

% for every observation, put all the points on top
for ii = 1:size(grotU, 1)
    
    % grab the centers and variability of the points
    xc = grotU(ii, ccf);
    yc = grotV(ii, ccf);
    
    % plot the CCA point, coloring the point be age
    plot(xc, yc, 'o', 'MarkerFaceColor', cmap(nage(ii), :), ...
         'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'MarkerSize', 6);

end

% if the full model is passed, there is no variability on the estimate
if full
    ccatxt = ['Canonical Correlation of Factor ' num2str(ccf) ': ' num2str(ccv) ];
else
    ccatxt = ['Canonical Correlation of Factor ' num2str(ccf) ': ' num2str(ccv) ' +/- ' num2str(ccs) ];
end

% catch the correlation with age and format it into the title
agetxt = [ 'Correlation of Factor ' num2str(ccf) ' with Age: ' num2str(ageR) ' +/- ' num2str(ageS) ];

% final plot clean up
title({ccatxt, agetxt});
xlabel('Network Values');
ylabel('Behavioral Scores');
axis equal; axis square;
cb = colorbar; caxis([ 18 88 ]);
ylabel(cb, 'Age');
hold off;

end

