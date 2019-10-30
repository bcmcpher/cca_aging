function [ rr, rse, pval, null ] = ccaLinRegCorr(cca, ccf, age, Nperm)
%[ r2, pval, null ] = ccaLinRegCorr(cca, ccf, age, Nperm);
%   Estimate a multiway correlation and run a permutation test to determine
%   it's significance.
%   
%   INPUTS:
%       cca   - data structure that holds factors
%       ccf   - the canonical correlation to pull loadings from
%       age   - a vector of each subjects age (holdout variable)
%       Nperm - the number of permutations to run for pval
%   OUTPUTS:
%       rr   - resampled mean r value estimated for multiway correlation
%       rse  - resampled standard deviation of the estimated r value
%       pval - p-value from bootstrap test to determine r significance
%       null - the null values estimated for the bootstrap (debug)
%
% Copyright (c) Brent McPherson (Indiana University), 2019. All rights reserved.
%

% adjusted r2 outs - dropped b/c odd behavior: , ar2, apval, anull
%
% mean((x-mean(x)).*(y-mean(y)).*(age-mean(age))) / (std(x).*std(y).*std(age))
% Degree1: 0.3888; Degree1: 0.3505
% BTW1: 0.3106; BTW2: 0.1692
%

% load data
%load([ 'camcan_594_' stem '_cca.mat' ], 'grotU', 'grotV');

grotU = cca.dat1.factor;
grotV = cca.dat2.factor;

% load age
%load('/N/dc2/projects/lifebid/HCP/Brent/camcan/canoncorr_analysis_full_data.mat', 'age');

% grab data
x = grotU(:, ccf); 
y = grotV(:, ccf);

% build design matrix with cca axes
z = [ x, y ];

% fit simple linear model of age ~ cca factor #
out = fitlm(z, age);

% return multiway correlation
r2 = sqrt(out.Rsquared.Ordinary);
%ar2 = sqrt(out.Rsquared.Adjusted);

%% perform a permuatation test to determine R2 significance

% preallocate null tests
null = nan(Nperm, 2);
nsub = size(z, 1);

% for every permutation
for ii = 1:Nperm
    
    % if I resort age as well, I can get a resampled mean + sd for the R2
    % if I only resort the data, I get the null distribution and find the p-value
    
    % grab a random set of subjects w/ replacement
    rsub = randsample(1:nsub, nsub, 'true');    
    
    % randomly sample input w/ replacement
    rdat = z(rsub, :);
    
    % fitlm on resorted w/ preplacement z
    tout1 = fitlm(rdat, age);
    tout2 = fitlm(rdat, age(rsub));    
    
    % catch relevant output
    null(ii, 1) = sqrt(tout1.Rsquared.Ordinary); % this is R, not R2 ???
    null(ii, 2) = sqrt(tout2.Rsquared.Ordinary);
    
end

% determine is R2 is signifcant
pval = 1 - (sum(null(:, 1) < r2) / Nperm);
rr = mean(null(:, 2));
rse = std(null(:, 2));

end
