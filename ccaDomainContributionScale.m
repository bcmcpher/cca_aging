function [ fh ] = ccaDomainContributionScale(scorr, svars)
%[ fh ] = ccaDomainContributionScale(scorr, svars);
%   Creates a plot of the domain names of the variables in the CCA are
%   scaled by the sum of the loadings of the individual variables.
%   
%   This will change when the CCA structure exists
%   This may not exist afterwards b/c the variable stems are hard coded and
%   will not be generalized for a tool.
%
%   INPUTS:
%       scorr - cca weight output (?)
%       svars - variable names w/ weight loadings for scaling
%   OUTPUTS:
%       fh - figure handle of the plot
%
% Copyright (c) Brent McPherson (Indiana University), 2019. All rights reserved.
%

% the logic to grab the relevant indices is in a nested cellfun call
% ~cellfun(@isempty, cellfun(@(x) regexp(x, '^PAT_'), svar', 'UniformOutput', false));
% first a cellfun is used to find all variables that begin w/ the regex pattern
% that call is wrapped in a cellfun call to isempty, finding cells w/o that pattern
% the opposite of that call is the resulting indices that are the domain

% indices of variables by domain

att_idx = logical(~cellfun(@isempty, cellfun(@(x) regexp(x, '^att_'), svars', 'UniformOutput', false)));
mtr_idx = logical(~cellfun(@isempty, cellfun(@(x) regexp(x, '^mtr_'), svars', 'UniformOutput', false)));
emt_idx = logical(~cellfun(@isempty, cellfun(@(x) regexp(x, '^emt_'), svars', 'UniformOutput', false)));
lng_idx = logical(~cellfun(@isempty, cellfun(@(x) regexp(x, '^lng_'), svars', 'UniformOutput', false)));
mem_idx = logical(~cellfun(@isempty, cellfun(@(x) regexp(x, '^mem_'), svars', 'UniformOutput', false)));
soc_idx = logical(~cellfun(@isempty, cellfun(@(x) regexp(x, '^soc_'), svars', 'UniformOutput', false)));
cog_idx = logical(~cellfun(@isempty, cellfun(@(x) regexp(x, '^co'), svars', 'UniformOutput', false))) | logical(~cellfun(@isempty, cellfun(@(x) regexp(x, '^hint_'), svars', 'UniformOutput', false)));

% make copies of positive / negative
pcorr = scorr;
ncorr = scorr;
acorr = scorr;

% zero out the opposite
pcorr(pcorr < 0) = NaN;
ncorr(ncorr > 0) = NaN;
acorr = abs(acorr);

% positive correlation
pscale = [ mean(pcorr(att_idx), 'omitnan'), mean(pcorr(mtr_idx), 'omitnan'), mean(pcorr(emt_idx), 'omitnan'), mean(pcorr(lng_idx), 'omitnan'), mean(pcorr(mem_idx), 'omitnan'), mean(pcorr(soc_idx), 'omitnan'), mean(pcorr(cog_idx), 'omitnan') ];
nscale = [ mean(ncorr(att_idx), 'omitnan'), mean(ncorr(mtr_idx), 'omitnan'), mean(ncorr(emt_idx), 'omitnan'), mean(ncorr(lng_idx), 'omitnan'), mean(ncorr(mem_idx), 'omitnan'), mean(ncorr(soc_idx), 'omitnan'), mean(ncorr(cog_idx), 'omitnan') ];
ascale = [ mean(acorr(att_idx), 'omitnan'), mean(acorr(mtr_idx), 'omitnan'), mean(acorr(emt_idx), 'omitnan'), mean(acorr(lng_idx), 'omitnan'), mean(acorr(mem_idx), 'omitnan'), mean(acorr(soc_idx), 'omitnan'), mean(acorr(cog_idx), 'omitnan') ];

% names for sorting
cornam = {'att', 'mtr', 'emt', 'lng', 'mem', 'soc', 'cog'};

% sort the most positive / most negative domains
[ ~, psi ] = sort(pscale, 'descend');
[ ~, nsi ] = sort(nscale);
[ ~, asi ] = sort(ascale, 'descend');

% plot the cirlces scaled off of each other

% positive
fh = figure; hold on
for ii = 1:7
    plot(ii, 1, 'o', 'MarkerFaceColor', [.1 .45 .95], ...
         'MarkerEdgeColor', 'k', 'LineWidth', 0.75, ...
         'MarkerSize', 100 * pscale(psi(ii)));
    text(ii, 1, cornam{psi(ii)});
end
set(gca, 'XLim', [ -1 8 ]);
hold off;

% negative
figure; hold on
for ii = 1:7
    plot(ii, 1, 'o', 'MarkerFaceColor', [.1 .45 .95], ...
         'MarkerEdgeColor', 'k', 'LineWidth', 0.75, ...
         'MarkerSize', 100 * abs(nscale(nsi(ii))));
    text(ii, 1, cornam{nsi(ii)});
end
set(gca, 'XLim', [ -1 8 ]);
hold off;

% absolute
figure; hold on
for ii = 1:7
    plot(ii, 1, 'o', 'MarkerFaceColor', [.1 .45 .95], ...
         'MarkerEdgeColor', 'k', 'LineWidth', 0.75, ...
         'MarkerSize', 100 * abs(ascale(asi(ii))));
    text(ii, 1, cornam{asi(ii)});
end
set(gca, 'XLim', [ -1 8 ]);
hold off;

end

