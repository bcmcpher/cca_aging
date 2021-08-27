function [ mag, sd, fh ] = ccaProjectFactor(cca, ccf, age, Nperm)
%[ mag, sd, fh ] = ccaProjectFactor(grotU, grotV, ccf, age, Nperm);
%   Project the axes of the requested CC factor as a single axis vector.
%   The significance is determined through a bootstrap and plot.
%
%   The inputs of this will change when a CCA object is defined
%   
%   INPUTS:
%       cca   - contains CCA variables for plotting
%       ccf   - the CC factor to project
%       age   - vector of subject ages; holdout variable to determine prediction
%       Nperm - number of permutations for a bootstrap estimate of sd
%   OUTPUTS:
%       mag - the magnitude vector of the projected CC
%       sd  - the bootstrap estimated standard deviation for mag
%       fh  - figure handle for the plot
%
% Copyright (c) Brent McPherson (Indiana University), 2019. All rights reserved.
%

grotU = cca.dat1.factor;
grotV = cca.dat2.factor;

nsubj = size(grotU, 1);

% fit trendline through first cca axis
tln = fit(grotU(:, ccf), grotV(:, ccf), 'poly1'); % x = 0.8617
% x = 3 * 0.8617
% 3 is max value on the other axis, .8617 is the x value from linear fit

% project data onto axis of cca plot
mag = nan(nsubj, 1);
for ii = 1:nsubj
    v = [ 1, 1*tln.p1 ]; % estimated from tln
    mag(ii) = dot([grotU(ii, 1) grotV(ii, 1) ], v) * v / [ norm(v) norm(v) ]; 
end

% correlation that predicts age hold out
%corr(mag, age)
%mean(mag) % functionally zero, as it should be

% resample for error bars
rcor = nan(Nperm, 1);
for ii = 1:Nperm
    ridx = randsample(1:nsubj, nsubj, true); % bootstrap
    rcor(ii) = corr(mag(ridx), age(ridx));
end

% bootstrapped CI
sd = std(rcor);

fh = figure('Position', [ 380 645 800 450 ]); hold on;
for ii = 1:nsubj
    plot(mag(ii), age(ii), 'o', 'MarkerFaceColor', [ 0.1 0.45 0.95 ], ...
         'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'MarkerSize', 6);
end
title([ 'Predicted Correlation Between CCA and Age: ' num2str(corr(mag, age)) ]);
xlabel('CCA Vector');
ylabel('Age');

end
