function [ pval, tdiff, ndiff ] = ccaModuleNullTest(mat1, mat2, nperm)
%[ pval, tdiff, ndiff ] = ccaModuleNullTest(mat1, mat2, nperm);
%   Perform a non-parametric test between the upper diagonal of 2 matrices
%   to determine if they are significantly different from one another.
%
%   INPUTS:
%       mat1  - the first symmetric data matrix
%       mat2  - the second symmetric data matrix
%       nperm - the number of permutations to perform (default = 10000)
%
%   OUTPUTS:
%       pval  - the estimated p-value for the difference
%       tdiff - the observed true difference
%       ndiff - the distribution of null differences
%
% Copyright (c) Brent McPherson (Indiana University), 2020. All rights reserved.
%

% default values for the number of permutations
if(~exist('nperm', 'var') || isempty(nperm))
    nperm = 10000;
end

% grab the size of the inputs
size1 = size(mat1, 1);
size2 = size(mat2, 1);

% pull the upper diagonal of values
dat1 = mat1(logical(triu(ones(size1), 1)));
dat2 = mat2(logical(triu(ones(size2), 1)));

% build the combined null distribution
dat0 = [ dat1; dat2 ];

% preallocate null pulls
rand1 = zeros(size1, nperm);
rand2 = zeros(size2, nperm);
ndiff = zeros(nperm, 1);

% for every permutation
for ii = 1:nperm
    
    % randomly sample mean
    rand1(:, ii) = randsample(dat0, size1, true);
    rand2(:, ii) = randsample(dat0, size2, true);
    
    % null difference
    ndiff(ii) = mean(rand1(:, ii), 1) - mean(rand2(:, ii), 1);
    
end

% pull the true differnece
tdiff = mean(dat1) - mean(dat2);

% estimate p-value
pval = sum(abs(ndiff) > abs(tdiff)) / nperm;

% print results
disp([ 'Oberved diff: ' num2str(tdiff) ]);
disp([ 'Null diff:    ' num2str(mean(ndiff)) ' +/- ' num2str(std(ndiff)) ]);
disp([ 'p-value:      ' num2str(pval) ]);

end
