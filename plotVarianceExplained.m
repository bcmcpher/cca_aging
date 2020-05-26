function [ fh ] = plotVarianceExplained(cca, dat)
%[ fh ] = plotVarianceExplained(cca, dat);
%   Plot the variance explained by each factor in the CCA for the requested
%   data set.
%
%   INTPUTS:
%       cca - the cca analysis object
%       dat - 'dat1' or 'dat2'; which data field to look at, brain or behavior
%   OUTPUTS:
%       fh - figure handle of the plot
%
% Copyright (c) Brent McPherson (Indiana University), 2019. All rights reserved.
%

%
% TRY AND GET AXIS SCALED TO PERCENT CORRECTLY
%

% variability = mean(fgrotAAd(:, ii).^2, 'omitnan');
obs = mean(cca.(dat).variability;
perm = cca.(dat).var;

% check that the data is the same size
if size(obs, 1) ~= size(perm, 2)
    error('The input data sets have mismatched x-axes.');
end

% define the size of the x-axis
Nkeep = size(obs, 1);

% boundaries of background data 
lb = min(perm(2:end, :));
ub = max(perm(2:end, :));
mn = mean(perm(2:end, :), 'omitnan');

% variance explained w/ null region
fh = figure; hold on

% plot filled in shape of min/max null permutations, skipping first row (the real data)
patch([1:Nkeep, Nkeep:-1:1], [ lb, ub ], ...
      [.1 .45 .95], 'FaceColor', [.9 .9 .9], 'EdgeColor', [ 0.5 0.5 0.5 ]);

% plot the average null line
plot(1:Nkeep, mn, 'black');

% plot the observed variability accounted for
plot(1:Nkeep, obs, 'blue', 'LineWidth', 1.5);

hold off

end

