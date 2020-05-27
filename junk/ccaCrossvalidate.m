function [ hcor, se, pval ] = ccaCrossvalidate(dat, prop, Nperm)
%[ hcor, se, pval ] = ccaCrossvalidate(uu1, uu2, prop, Nperm);
%   Perform a crossvalidation of the CCA using a proportion split of
%   prepared data from the data preparation for a variable number of
%   permutations.
%   
%   The inputs will change when I figure out the structure for the CCA type
%
%   INTPUTS:
%       dat   - preprocessed data object from CCA
%       prop  - the proportion of data to split and validate with
%       Nperm - the number of permutations to run on during crossvalidation
%   OUTPUTS:
%       hcor - the crossvalidated canonical correlation
%       se   - the crossvalidated standard error of each canonical correlation
%       pval - the crossvalidated p-value of each canonical correlation
%
% Probably deprecated b/c of internal crossvalidation procedure
%
% Copyright (c) Brent McPherson (Indiana University), 2019. All rights reserved.
%

uu1 = dat.dat1.uu1;
uu2 = dat.dat2.uu2;

% preallocate outputs based on input size
nsubj = size(uu1, 1);
ncorr = size(uu1, 2);
hocca = nan(Nperm, ncorr);
hcor = nan(ncorr, 1);
pval = nan(ncorr, 1);
se = nan(ncorr, 1);

% randomly sample a proportion of subjects to train
ntrn = round(nsubj * prop);
ntst = nsubj - ntrn;

% for every permutation
for ii = 1:Nperm

    % sample w/ replacement from a random sample for training
    trn_sub = randsample(1:nsubj, ntrn, 'true');
    
    % sample w/ replacement from the subjects not included for testing
    tst_pop = setdiff(1:nsubj, trn_sub);
    tst_sub = randsample(tst_pop, ntst, 'true');
    
    % this computes CCA on the training data
    [ A, B ] = canoncorr(uu1(trn_sub, :), uu2(trn_sub, :));
    
    % projecting the test data
    Xtestpr = uu1(tst_sub, :) * A;
    Ytestpr = uu2(tst_sub, :) * B;
    
    % for every correlation 
    for jj = 1:ncorr
        
        % compute correlation based on holdout loadings
        hocca(ii, jj) = corr(Xtestpr(:, jj), Ytestpr(:, jj));
    
    end
    
    % % loop over all train samples
    % correct = 0;
    % for i=1:size(Xtestpr,1)
    %     % Using only the first CCA projection, find the sample in Y
    %     % closest to the one in X. Euclidean distance is used as a
    %     % similarity measure.
    %     [~, ind] = min(sum((Xtestpr(i,1) - Ytestpr(:,1)).^2, 2));
    %
    %     % if classified correctly
    %     if ind==i
    %         correct = correct+1;
    %     end
    % end
    %
    % % compute the probability that so many correct matchings could be obtained by chance
    % pval = 1 - binocdf(correct-1, size(Xtestpr,1), 1/size(Xtestpr,1));
    %
    % % compute confidence interval on correct matching rate
    % [ ~, cinf ] = binofit(correct, size(Xtestpr,1));
    
end

clear ii jj trn_sub tst_pop tst_sub A B r Xtestpr Ytestpr 

% get the true cca correlations
[ ~, ~, r ] = canoncorr(uu1, uu2);

% estimate hold out correlation, standard error, and p-value for each CC
for ii = 1:ncorr
    
    % find the mean hold out cca
    hcor(ii) = mean(hocca(:, ii));
    
    % create standard error estimate from permuted holdouts (?)
    se(ii) = std(hocca(:, ii));
        
    % get true cca value
    tcor = r(ii);
    
    % determine how often hold out is as high or higher than full data
    pval(ii) = 1 - (sum(hocca(:, ii) >= tcor) / Nperm);
    
end

end

