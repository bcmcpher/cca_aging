function [ dat, cca ] = ccaFullAnalysis(NET, vars, conf, netsNames, varsNames, confNames, varsLabel, Nkeep, Nperm)
%[ dat, cca ] = ccaFullAnalysis(NET, vars, conf, netsNames, varsNames, confNames, Nkeep, Nperm);
%   This function takes brain, behavior, and confound matrices of subj x vars
%   size, normalizes them, and performs canonical correlation analysis.
%   It returns simple structures for plotting / reporting different analyses.
%
%   INPUTS:
%       NET  - subj x net of brain network data (node measures or upper diagonal) 
%       vars - subj x var of behavior data
%       conf - subj x reg of the confound regressors to be removed from NET and vars
%       netsNames - net x 1 cell array of network variable labels
%       varsNames - var x 1 cell array of behavior variable labels
%       confNames - reg x 1 cell array of confound variable labels
%       Nkeep - number of PCA components / CCA dimensions to use
%       Nperm - number of permutations for bootstrapping null hypothesis
%   OUTPUTS:
%       dat - any normalized / PCA data stored for each data set
%       cca - all CCA component values of interest for reporting or plotting
%
% References:
% http://www.statsoft.com/Textbook/Canonical-Analysis
% https://stats.idre.ucla.edu/spss/output/canonical-correlation-analysis/
%
% Copyright (c) Brent McPherson (Indiana University), 2019. All rights reserved.
%

%% setup confounds matrix

disp('Normalizing confounds data...');

% normalize, impute w/ 0s, and add squared terms
conf = palm_inormal(conf);
conf(isnan(conf)) = 0;
conf = nets_normalise([ conf conf.^2 ]);

% create squared nuissance names
confNames = [ confNames cellfun(@(c)[ c '^2' ], confNames, 'uni', false) ];

% idenfity the badly specified confounds
badconf = zeros(size(conf, 2), 1);
for ii = 1:size(conf, 2)
    
    Y = conf(:, ii); 
    grotKEEP =~ isnan(Y);
    grot = (Y(grotKEEP) - median(Y(grotKEEP))).^2; 
    grot = max(grot / mean(grot));  % do we have extreme outliers?
    
    % if the variables are not poorly defined
    if (sum(grotKEEP) > 250) && (std(Y(grotKEEP)) > 0) && (max(sum(nets_class_vectomat(Y(grotKEEP)))) / length(Y(grotKEEP)) < 0.95) && (grot < 100)
        % do nothing
    else
        badconf(ii) = 1;
    end
end

badconf = logical(badconf);

clear ii Y grotKEEP grot

disp(['Removing ' num2str(sum(badconf)) ' poorly defined variables...']);

% drop bad nuisance regressors / labels
conf(:, badconf) = [];
confNames(badconf) = [];
% they should all be kept

% store confounds in the output
dat.conf.conf = conf;
dat.conf.confNames = confNames;

%% prepare main netmat matrix 

disp('Normalizing and deconfounding network data...');

NET1 = nets_demean(NET); % if pass log10 of degree
NET1 = NET1 / std(NET1(:)); % need to zero out inf vals here
amNET = abs(mean(NET)); % need to zero out inf vals here
NET3 = nets_demean(NET ./ repmat(amNET, size(NET, 1), 1)); % need to zero out NaN here
NET3(:, amNET < 0.1) = [];
NET3 = NET3 / std(NET3(:)); % norm by mean of columns, removing badly conditioned ones
grot = [ NET1 NET3 ];       % concat horizontally
NETd = nets_demean(grot - conf * (pinv(conf) * grot)); % deconfound and demean % need to zero nan

disp('Performing PCA on network data...');

% SVD reduction
[ uu1, ss1 ] = nets_svds(NETd, Nkeep);

% store network data in the output
dat.dat1.uu1 = uu1;
dat.dat1.ss1 = ss1;
dat.dat1.raw = NET;
dat.dat1.nrm = NETd;
dat.dat1.names = netsNames;

%% identify "bad" SMs - e.g. because of bad outliers or not enough distinct values

disp('Beginning normalization of behaior data...');

badvars = zeros(size(vars, 2), 1);
for ii = 1:size(vars, 2)
    
    Y = vars(:, ii); 
    grotKEEP =~ isnan(Y);
    grot = (Y(grotKEEP) - median(Y(grotKEEP))).^2; 
    grot = max(grot / mean(grot));  % do we have extreme outliers?
    
    % if the variables are not poorly defined
    if (sum(grotKEEP) > 250) && (std(Y(grotKEEP)) > 0) && (max(sum(nets_class_vectomat(Y(grotKEEP)))) / length(Y(grotKEEP)) < 0.95) && (grot < 100)
        % do nothing
    else
        badvars(ii) = 1;
    end
end

badvars = logical(badvars);

clear ii Y grotKEEP grot

disp(['Removing ' num2str(sum(badvars)) ' poorly defined variables...']);

% keep well defined variables
vars(:, badvars) = [];
varsNames(badvars) = [];
varsLabel(badvars) = [];

% "impute" missing vars data / normalize
varsd = palm_inormal(vars);

disp('Normalizing and deconfounding behavior data...');

% deconfound ignoring missing data
for ii = 1:size(varsd, 2)
    grot = (isnan(varsd(:, ii)) == 0); 
    grotconf = nets_demean(conf(grot, :)); 
    varsd(grot, ii) = normalise(varsd(grot, ii) - grotconf * (pinv(grotconf) * varsd(grot, ii)));
end

clear ii grot grotconf

% estimate "pairwise" covariance, ignoring missing data
varsdCOV = zeros(size(varsd, 1));
for ii = 1:size(varsd, 1)
    for jj = 1:size(varsd, 1)
        grot = varsd([ ii jj ], :); 
        grot = cov(grot(:, sum(isnan(grot)) == 0)'); 
        varsdCOV(ii, jj) = grot(1, 2);
    end
end

clear ii jj grot

% minor adjustment: project onto the nearest valid covariance matrix
varsdCOV2 = nearestSPD(varsdCOV);

disp('Performing PCA on network data...');

% SVD (eigs actually)
[ uu, dd ] = eigs(varsdCOV2, Nkeep);

% deconfound again just to be safe
uu2 = uu - conf * (pinv(conf) * uu);

% normalize data for factor creation
varsgrot = palm_inormal(vars);

% store behavior data in the output structure
dat.dat2.uu2 = uu2;
dat.dat2.ss2 = dd;
dat.dat2.raw = vars;
dat.dat2.nrm = varsgrot;
dat.dat2.names = varsNames;
dat.dat2.label = varsLabel;

disp('Storing preprocessed data...');

%% run the analysis

disp('Performing CCA...');

[ A, B, grotR, grotU, grotV, grotstats ] = canoncorr(uu1, uu2);

% save fields to output
cca.dat1.scale = A;
cca.dat1.factor = grotU;

cca.dat2.scale = B;
cca.dat2.factor = grotV;

cca.cca.corrs = grotR;
cca.cca.stats = grotstats;

%% extract parameters back to raw measures

disp('Converting CCA loadings back to raw network / behavior measures...');

% preallocate 
grotAAd = nan(size(NET, 2), Nkeep);
grotBBd = nan(size(varsgrot, 2), Nkeep);
grotAAv = nan(Nkeep, 1);
grotBBv = nan(Nkeep, 1);
grotAAr = nan(Nkeep, 1);
grotBBr = nan(Nkeep, 1);

for ii = 1:Nkeep
    
    % network weights after deconfounding
    grotAAd(:, ii) = corr(grotU(:, ii), NETd(:, 1:size(NET, 2)))';
    % only use the first columns, ignore the additional pca confounds

    % behavior weights after deconfounding
    grotBBd(:, ii) = corr(grotV(:, ii), varsgrot, 'rows', 'pairwise')';
    
    % variability extracted by each cc data set
    grotAAv(ii) = mean(grotAAd(:, ii).^2, 'omitnan');
    grotBBv(ii) = mean(grotBBd(:, ii).^2, 'omitnan');
    
    % redundancy of each cc component
    grotAAr(ii) = mean(grotAAd(:, ii).^2, 'omitnan') * grotR(ii)^2;
    grotBBr(ii) = mean(grotBBd(:, ii).^2, 'omitnan') * grotR(ii)^2;
    
end

clear ii

% proportion of variability accounted for in each CC
grotRv = grotR .^ 2;

% percent accounted for in each CC
grotPc = (grotRv ./ sum(grotRv)) * 100;

% save fields to output
cca.dat1.loading = grotAAd;
cca.dat1.variability = grotAAv;
cca.dat1.redundancy = grotAAr;

cca.dat2.loading = grotBBd;
cca.dat2.variability = grotBBv;
cca.dat2.reduncancy = grotBBr;

cca.cca.varExplained = grotRv;
cca.cca.percentExp = grotPc;

%% CCA permutation testing for parameters

disp('Performing a permutation test to determine significant number of CCs...');

% preallocate output
grotRp = zeros(Nperm, Nkeep, 2); 
grotRpval = nan(Nkeep, 2);

% preallocate null noise data
grotAtr = zeros(Nperm, Nkeep, size(NET, 2));
grotBtr = zeros(Nperm, Nkeep, size(varsgrot, 2));
grotAv = zeros(Nperm, Nkeep);
grotBv = zeros(Nperm, Nkeep);
grotAr = zeros(Nperm, Nkeep);
grotBr = zeros(Nperm, Nkeep);

% predefine permutation sets
PAPset = palm_quickperms([], ones(size(varsgrot, 1), 1), Nperm);

disp('Reconstructing loadings from permutations...');

% for every permutation
for ii = 1:Nperm
    
    % run cca w/ scrambled data set
    [ ~, ~, grotRp(ii, :, 1), grotUr, ~ ] = canoncorr(uu1, uu2(PAPset(:, ii), :));
    [ ~, ~, grotRp(ii, :, 2), ~, grotVr ] = canoncorr(uu1(PAPset(:, ii), :), uu2);
    
    % for every cc
    for jj = 1:Nkeep
        
        % network weights after deconfounding
        grotAtr(ii, jj, :) = corr(grotUr(:, jj), NETd(:, 1:size(NET, 2)))';
        
        % behavior weights after deconfounding
        grotBtr(ii, jj, :) = corr(grotVr(:, jj), varsgrot, 'rows', 'pairwise')';
        
        % store variability captured in A / B
        grotAv(ii, jj) = mean(grotAtr(ii, jj, :) .^ 2, 'omitnan');
        grotBv(ii, jj) = mean(grotBtr(ii, jj, :) .^ 2, 'omitnan');
        
        % store redundancy in A / B
        grotAr(ii, jj) = mean(grotAtr(ii, jj, :) .^ 2, 'omitnan') * grotRp(ii, jj);
        grotBr(ii, jj) = mean(grotBtr(ii, jj, :) .^ 2, 'omitnan') * grotRp(ii, jj);
        
    end
end

clear ii jj grotUr grotVr

for ii = 1:Nkeep  % get FWE-corrected pvalues
    grotRpval(ii, 1) = (1 + sum(grotRp(2:end, 1) >= grotR(ii))) / Nperm;
    grotRpval(ii, 2) = (1 + sum(grotRp(2:end, 2) >= grotR(ii))) / Nperm;
end

clear ii 

% count # of CCA factors that pass correction
Ncca1 = sum(grotRpval(:, 1) < 0.05); % number of FWE-significant CCA components
Ncca2 = sum(grotRpval(:, 2) < 0.05); % number of FWE-significant CCA components

% save fields to output
cca.dat1.var = grotAv;
cca.dat1.red = grotAr;
cca.dat1.rld = grotAtr;

cca.dat2.var = grotBv;
cca.dat2.red = grotBr;
cca.dat2.rld = grotBtr;

cca.cca.pval = grotRpval;
cca.cca.grotRp = grotRp;
cca.cca.ncca1 = Ncca1;
cca.cca.ncca2 = Ncca2;

end
