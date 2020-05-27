function [ rr, rse, pval, null ] = ccaLinRegCorr(cca, ccf, age, Nperm, full)
%[ r2, pval, null ] = ccaLinRegCorr(cca, ccf, age, Nperm);
%   Estimate a multiway correlation and run a permutation test to determine
%   it's significance.
%   
%   INPUTS:
%       cca   - data structure that holds factors
%       ccf   - the canonical correlation to pull loadings from
%       age   - a vector of each subjects age (holdout variable)
%       Nperm - the number of permutations to run for pval
%       full  - logical, use the full dataset CCA (not cross-validated)
%   OUTPUTS:
%       rr   - resampled mean r value estimated for multiway correlation
%       rse  - resampled standard deviation of the estimated r value
%       pval - p-value from bootstrap test to determine r significance
%       null - the null values estimated for the bootstrap (debug)
%
% Copyright (c) Brent McPherson (Indiana University), 2019. All rights reserved.
%

% adjusted r2 outs - dropped b/c odd behavior: , ar2, apval, anull

if(~exist('full', 'var') || isempty(full))
    full = false;
end

if full
    grotU = cca.full.grotU;
    grotV = cca.full.grotV;
else
    grotU = cca.dat1.factor;
    grotV = cca.dat2.factor;
end

% grab data
x = grotU(:, ccf); 
y = grotV(:, ccf);

% build design matrix with cca axes
z = [ x, y ];

% fit simple linear model of age ~ cca factor #
out = fitlm(z, age);

% return multiway correlation
r = sqrt(out.Rsquared.Ordinary);
%ar2 = sqrt(out.Rsquared.Adjusted);

%% perform a permuatation test to determine R2 significance

% preallocate null tests
null = nan(Nperm, 2);
nsub = size(z, 1);

% for every permutation
for ii = 1:Nperm
    
    % if I resort age as well, I can get a resampled mean + sd for the R
    % if I only resort the data, I get the null distribution and find the p-value
    
    % grab a random set of subjects w/ replacement
    rsub = randsample(1:nsub, nsub, 'true');    
    
    % randomly sample input w/ replacement
    rdat = z(rsub, :);
    
    % fitlm on resorted w/ replacement z
    tout1 = fitlm(rdat, age);
    tout2 = fitlm(rdat, age(rsub));    
    
    % catch relevant output
    null(ii, 1) = sqrt(tout1.Rsquared.Ordinary); % this is R, not R2
    null(ii, 2) = sqrt(tout2.Rsquared.Ordinary);
    
end

% determine if R is signifcant
pval = 1 - (sum(null(:, 1) < r) / Nperm);
rr = mean(null(:, 2));
rse = std(null(:, 2));

end
