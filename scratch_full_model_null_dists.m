%% CCA permutation testing for parameters
    
Nsub = size(m1.dat.dat1.uu1, 1);
Nkeep = size(m1.dat.dat1.uu1, 2);
Nperm = 10000;

% preallocate output
m1R = zeros(Nkeep, Nperm);
m1U = zeros(Nsub, Nkeep, Nperm);
m1V = zeros(Nsub, Nkeep, Nperm);

m2R = zeros(Nkeep, Nperm);
m2U = zeros(Nsub, Nkeep, Nperm);
m2V = zeros(Nsub, Nkeep, Nperm);

%    grotAtr = zeros(Nnet, Nkeep, Nperm);
%    grotBtr = zeros(Nvar, Nkeep, Nperm);

% predefine permutation sets, add 1 to Nperm b/c the first one isn't randomized
PAPset = palm_quickperms([], ones(Nsub, 1), Nperm+1);

% for every permutation
for ii = 1:Nperm
    
    % get cross-validated loadings for each permutation
    % use ii+1 to skip first, original ordering of data, only use random sortings
    [ ~, ~, ~, m1U(:,:,ii) ] = canoncorr(m1.dat.dat1.uu1, m1.dat.dat2.uu2(PAPset(:, ii+1), :));
    [ ~, ~, ~, ~, m1V(:,:,ii) ] = canoncorr(m1.dat.dat1.uu1(PAPset(:, ii+1), :), m1.dat.dat2.uu2);
    
    [ ~, ~, ~, m2U(:,:,ii) ] = canoncorr(m2.dat.dat1.uu1, m2.dat.dat2.uu2(PAPset(:, ii+1), :));
    [ ~, ~, ~, ~, m2V(:,:,ii) ] = canoncorr(m2.dat.dat1.uu1(PAPset(:, ii+1), :), m2.dat.dat2.uu2);
    
    % for every cc
    for jj = 1:Nkeep

        % get the corr between the null splits
        m1R(jj, ii) = corr(m1U(:,jj,ii), m1V(:,jj,ii));
        m2R(jj, ii) = corr(m2U(:,jj,ii), m2V(:,jj,ii));

    end
end

% get the percentiles from here? are these the absolute values?
prctile(squeeze(m1R(1,:)), [ 0.05 0.01 0.001 0.0001 ])
prctile(squeeze(m2R(1,:)), [ 0.05 0.01 0.001 0.0001 ])

% the plots are normal around 0, the abs of prctile would be the tails...
figure; hist(squeeze(m1R(1,:)), 128);
figure; hist(squeeze(m2R(1,:)), 128);

% direct estimation of pval
m1pval = nan(Nkeep, 1);
m2pval = nan(Nkeep, 1);
for ii = 1:Nkeep
    m1pval(ii) = sum(m1R(ii,:) > m1.cca.full.grotR(ii)) / Nperm;
    m2pval(ii) = sum(m2R(ii,:) > m2.cca.full.grotR(ii)) / Nperm;
end

% show pvals for all axes of full models above cv models
[ m1pval'; m1.cca.cca.pval ]'
[ m2pval'; m2.cca.cca.pval ]'

%% compute 2:end average correlation w/ age

zzz = nan(Nkeep, 3);
for ii = 1:Nkeep
    [ zzz(, ~, zzz(ii) ] = ccaLinRegCorr(m2.cca, ii, age, 1000);
end
mean(zzz(2:end))

%% catch age and pvalue for both models

% preallocate full / cross-validated corr w/ age and pvalues for both models
flca = nan(Nkeep, 6);
cvca = nan(Nkeep, 6);

for ii = 1:Nkeep

    % full model
    [ flca(ii, 1), flca(ii, 2), flca(ii, 3) ] = ccaLinRegCorr(m1.cca, ii, age, 10000, true);
    [ flca(ii, 4), flca(ii, 5), flca(ii, 6) ] = ccaLinRegCorr(m2.cca, ii, age, 10000, true);

    % cross-validated model
    [ cvca(ii, 1), cvca(ii, 2), cvca(ii, 3) ] = ccaLinRegCorr(m1.cca, ii, age, 10000);
    [ cvca(ii, 4), cvca(ii, 5), cvca(ii, 6) ] = ccaLinRegCorr(m2.cca, ii, age, 10000);
    
end


