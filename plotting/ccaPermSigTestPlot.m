function [ fh ] = ccaPermSigTestPlot(cca, ccf)
%[ fh ] = ccaPermSigTestPlot(grotU, grotV, grotRp, ccf);
%   Plot all canonical correlations to visualize the relative strengths
%
%   INPUTS:
%       cca - the values from a CCA analysis
%       ccf - the canonical correlation to emphasize and print
%   OUTPUTS:
%       fh - figure handle of the requested plot
%
% Copyright (c) Brent McPherson (Indiana University), 2019. All rights reserved.
%

grotU = cca.dat1.factor;
grotV = cca.dat2.factor;
grotRp = squeeze(cca.cca.hocorrs_null(1, ccf, :));

grotR = corr(grotU(:, ccf), grotV(:, ccf));

fh = figure; hold on;
hist(grotRp, 64);
N = hist(grotRp, 64);
line([ grotR grotR ], [ 0 1.5*max(N) ], 'color', 'red', 'LineWidth', 3);
title([ 'Null Distribution of Permuted CCA Correlations for CC: ' num2str(ccf) ]);
hold off;

end

