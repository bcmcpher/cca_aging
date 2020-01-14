function [ dat, cca, cc2 ] = ccaMapFullAnalysis(NET, vars, conf, netsNames, varsNames, confNames, varsLabel, Nkeep1, Nkeep2, Nperm, Nfold, Nrep)
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

%% set default values

% if(~exist('ct', 'var') || isempty(ct))
%     ct = 'median';
% end

if(~exist('Nperm', 'var') || isempty(Nperm))
    Nperm = 0;
end

if(~exist('Nfold', 'var') || isempty(Nfold))
    Nfold = 5;
end

if(~exist('Nrep', 'var') || isempty(Nrep))
    Nrep = 1000;
end

if Nrep < 1000
    warning('You need at least 1000 reps.');
    Nrep = 1000;
end

% if to few permutations are passed, skip them
if Nperm < 2
    disp('Bootstrap permutation testing disabled.');
end

% force the minimum number of PCA components to become canonical factors for global comparisons
Nkeep = min([ Nkeep1, Nkeep2 ]);

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
[ uu1, ss1 ] = nets_svds(NETd, Nkeep1);

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

disp('Performing PCA on behavior data...');

% SVD (eigs actually)
[ uu, dd ] = eigs(varsdCOV2, Nkeep2);

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

% perform the full analysis
[ A, B, grotR, grotU, grotV, grotstats ] = canoncorr(uu1, uu2);

disp('Cross-validating CCA estimates...');

% create repeated k-fold estimated loadings
[ cvU, cvV, trR, hoR, grotAAd, grotBBd, ...
  grotAAv, grotBBv, grotAAr, grotBBr ] = cvCCA(uu1, uu2, NETd, varsgrot, Nfold, Nrep);

% save full data cca fields
cca.full.A = A;
cca.full.B = B;
cca.full.grotR = grotR;
cca.full.grotU = grotU;
cca.full.grotV = grotV;
cca.full.stats = grotstats;

% save cross-validated fields to output
cca.dat1.factor = cvU.mn;
cca.dat1.factor_sd = cvU.std;
cca.dat1.loading = grotAAd.mn;
cca.dat1.loading_sd = grotAAd.std;
cca.dat1.variability = grotAAv.mn;
cca.dat1.variability_sd = grotAAv.std;
cca.dat1.redundancy = grotAAr.mn;
cca.dat1.redundancy_sd = grotAAr.std;

cca.dat2.factor = cvV.mn;
cca.dat2.factor_sd = cvV.std;
cca.dat2.loading = grotBBd.mn;
cca.dat2.loading_sd = grotBBd.std;
cca.dat2.variability = grotBBv.mn;
cca.dat2.variability_sd = grotBBv.std;
cca.dat2.redundancy = grotBBr.mn;
cca.dat2.redundancy_sd = grotBBr.std;

cca.cca.trcorrs = trR.mn;
cca.cca.trcorrs_sd = trR.std;
cca.cca.hocorrs = hoR.mn;
cca.cca.hocorrs_sd = hoR.std;

% proportion of variability accounted for in each CC / percent accounted for in each CC
cca.cca.trRv = mean(trR.mn, 'omitnan') .^ 2;
cca.cca.trPc = (cca.cca.trRv ./ sum(cca.cca.trRv)) * 100;
cca.cca.hoRv = mean(hoR.mn, 'omitnan') .^ 2;
cca.cca.hoPc = (cca.cca.hoRv ./ sum(cca.cca.hoRv)) * 100;

%% extract parameters back to raw measures
% this rebuilds loadings from ct of cv factors, not ct of loadings w/in cv
% should be equivalent: wasn't the first time, is now?
% if it is equivalent, drop this part

disp('Converting crossvalidated CCA loadings back to raw network / behavior measures...');

% preallocate 
cgrotAAd = nan(size(NET, 2), Nkeep);
cgrotBBd = nan(size(varsgrot, 2), Nkeep);
% cgrotAAv = nan(Nkeep, 1);
% cgrotBBv = nan(Nkeep, 1);
% cgrotAAr = nan(Nkeep, 1);
% cgrotBBr = nan(Nkeep, 1);
cgrotR = nan(Nkeep, 1);

% for every CC
for ii = 1:Nkeep
        
        % network weights after deconfounding
        cgrotAAd(:, ii) = corr(cca.dat1.factor(:, ii), NETd(:, 1:size(NET, 2)))';
        % only use the first normalized network columns, ignore the additional pca confounds
        
        % behavior weights after deconfounding
        cgrotBBd(:, ii) = corr(cca.dat2.factor(:, ii), varsgrot, 'rows', 'pairwise')';
        
        % pull observed correlation
        cgrotR(ii) = corr(cca.dat1.factor(:, ii), cca.dat2.factor(:, ii));
        
%         % variability extracted by each cc data set
%         cgrotAAv(ii) = mean(cgrotAAd(:, ii).^2, 'omitnan');
%         cgrotBBv(ii) = mean(cgrotBBd(:, ii).^2, 'omitnan');
%         
%         % redundancy of each cc component
%         cgrotAAr(ii) = mean(cgrotAAd(:, ii).^2, 'omitnan') * grotR(ii)^2;
%         cgrotBBr(ii) = mean(cgrotBBd(:, ii).^2, 'omitnan') * grotR(ii)^2;

end

% THE CV LOADINGS ARE NEARLY THE SAME AS THE REBUILT LOADINGS
% HOW DO WE DETERMINE WHICH ONES TO USE?

clear ii

% proportion of variability accounted for in each CC
cgrotRv = cgrotR .^ 2;

% percent accounted for in each CC
cgrotPc = (cgrotRv ./ sum(cgrotRv)) * 100;

% save fields to output
cc2.dat1.factor = cvU.mn;
cc2.dat1.factor_sd = cvU.std;
cc2.dat2.factor = cvV.mn;
cc2.dat2.factor_sd = cvV.std;

cc2.dat1.loading = cgrotAAd;
% cc2.dat1.variability = cgrotAAv;
% cc2.dat1.redundancy = cgrotAAr;

cc2.dat2.loading = cgrotBBd;
% cc2.dat2.variability = cgrotBBv;
% cc2.dat2.reduncancy = cgrotBBr;

cc2.cca.varExplained = cgrotRv;
cc2.cca.percentExp = cgrotPc;

% FIGURE OUT HOW TO TURN ON / OFF THE PERMUTATION TESTING
%% CCA permutation testing for parameters

if Nperm >= 2
    
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
    % add 1 to Nperm b/c the first one isn't random
    PAPset = palm_quickperms([], ones(size(varsgrot, 1), 1), Nperm+1);
    
    disp('Reconstructing loadings from permutations...');
    
    % for every permutation, skipping the first, non-random entry
    for ii = 2:Nperm+1
        
        % get cross-validated loadings for each permutation
        [ grotUr, ~, ~, thoR1 ] = cvCCA(uu1, uu2(PAPset(:, ii), :), NETd, varsgrot, Nfold, 5000);
        [ ~, grotVr, ~, thoR2 ] = cvCCA(uu1(PAPset(:, ii), :), uu2, NETd, varsgrot, Nfold, 5000);
        
        % catch the average across the holdouts - should this be min/max?
        grotRp(ii, :, 1) = mean(mean(thoR1.mn), 'omitnan');
        grotRp(ii, :, 2) = mean(mean(thoR2.mn), 'omitnan');
        
        % for every cc
        for jj = 1:Nkeep
            
            % network weights after deconfounding
            grotAtr(ii, jj, :) = corr(grotUr.mn(:, jj), NETd(:, 1:size(NET, 2)))';
            
            % behavior weights after deconfounding
            grotBtr(ii, jj, :) = corr(grotVr.mn(:, jj), varsgrot, 'rows', 'pairwise')';
            
            % store variability captured in A / B
            grotAv(ii, jj) = mean(grotAtr(ii, jj, :) .^ 2, 'omitnan');
            grotBv(ii, jj) = mean(grotBtr(ii, jj, :) .^ 2, 'omitnan');
            
            % store redundancy in A / B
            grotAr(ii, jj) = mean(grotAtr(ii, jj, :) .^ 2, 'omitnan') * grotRp(ii, jj);
            grotBr(ii, jj) = mean(grotBtr(ii, jj, :) .^ 2, 'omitnan') * grotRp(ii, jj);
            
        end
    end
    
    clear ii jj grotUr grotVr
    
    % count # of CCA factors that pass correction from the hold-out data
    Ncca1 = sum(grotRpval(:, 1) < 0.05); % number of FWE-significant CCA components
    Ncca2 = sum(grotRpval(:, 2) < 0.05); % number of FWE-significant CCA components
    
    % save fields to output
    cca.dat1.var = grotAv;
    cca.dat1.red = grotAr;
    cca.dat1.rld = grotAtr;
    
    cca.dat2.var = grotBv;
    cca.dat2.red = grotBr;
    cca.dat2.rld = grotBtr;
    
    cca.cca.ncca1 = Ncca1;
    cca.cca.ncca2 = Ncca2;
    
    cca.cca.pval = grotRpval;
    cca.cca.grotRp = grotRp;
    
else
    
    % fill in fields as empty
    cca.dat1.var = [];
    cca.dat1.red = [];
    cca.dat1.rld = [];
    
    cca.dat2.var = [];
    cca.dat2.red = [];
    cca.dat2.rld = [];
    
    cca.cca.ncca1 = [];
    cca.cca.ncca2 = [];
    
    cca.cca.pval = [];
    cca.cca.grotRp = [];
    
end

end

%% cross-validation function for compartmentalization
function [ cvU, cvV, trR, hoR, grotAAd, grotBBd, grotAAv, grotBBv, grotAAr, grotBBr ] = cvCCA(uu1, uu2, NETd, varsgrot, Nfold, Nrep)

% the number of repeats to do per mapped estimate - always at least 5000 reps passed
Nmap = 5000;

% determine the number of map estimates to make
Nsplit = Nrep / Nmap;

% always set loop to 1, even if the size is less
if Nsplit < 1
    Nsplit = 1;
end

% fake the number kept by size of input columns
Nkeep = min([ size(uu1, 2) size(uu2, 2) ]);

% preallocate reduced output strucutres
cvU = initCVstruct(size(uu1, 1), Nkeep);
cvV = initCVstruct(size(uu2, 1), Nkeep);
trR = initCVstruct(Nfold, Nkeep);
hoR = initCVstruct(Nfold, Nkeep);
grotAAd = initCVstruct((size(NETd, 2)/2), Nkeep);
grotBBd = initCVstruct(size(varsgrot, 2), Nkeep);
grotAAv = initCVstruct(Nkeep, 1);
grotBBv = initCVstruct(Nkeep, 1);
grotAAr = initCVstruct(Nkeep, 1);
grotBBr = initCVstruct(Nkeep, 1);

% for every part
for map = 1:Nsplit
    
    % preallocate split outputs
    mcvU = zeros(size(uu1, 1), Nkeep, Nmap);
    mcvV = zeros(size(uu2, 1), Nkeep, Nmap);
    mtrR = zeros(Nfold, Nkeep, Nmap);
    mhoR = zeros(Nfold, Nkeep, Nmap);
    mgrotAAd = zeros((size(NETd, 2)/2), Nkeep, Nmap);
    mgrotBBd = zeros(size(varsgrot, 2), Nkeep, Nmap);
    mgrotAAv = zeros(Nkeep, Nmap);
    mgrotBBv = zeros(Nkeep, Nmap);
    mgrotAAr = zeros(Nkeep, Nmap);
    mgrotBBr = zeros(Nkeep, Nmap);
    
    % for every repeated k-fold
    for rep = 1:Nsplit
        
        % create the requested number of folds across all subjects
        cvlab = crossvalind('Kfold', 1:size(uu1, 1), Nfold);
        
        % for every fold
        for fold = 1:Nfold
            
            % make a copy of full data for dropping values
            truu1 = uu1;
            truu2 = uu2;
            
            % pull the hold-out indices
            hoidx = cvlab == fold;
            
            % remove the hold-out from the training data
            truu1(hoidx, :) = [];
            truu2(hoidx, :) = [];
            
            % estimate the weights and corr from the training data
            [ cvA, cvB, mtrR(fold, :, rep) ] = canoncorr(truu1, truu2);
            
            % get the estimates for the hold-out subjects
            mcvU(hoidx, :, rep) = uu1(hoidx, :) * cvA;
            mcvV(hoidx, :, rep) = uu2(hoidx, :) * cvB;
            
            % grab the validated correlations b/w the hold-out and rebuilt loadings
            for cf = 1:Nkeep
                mhoR(fold, cf, rep) = corr(mcvU(hoidx, cf, rep), mcvV(hoidx, cf, rep));
                mgrotAAd(:, cf, rep) = corr(mcvU(:, cf, rep), NETd(:, 1:(size(NETd, 2)/2)))';
                mgrotBBd(:, cf, rep) = corr(mcvV(:, cf, rep), varsgrot, 'rows', 'pairwise')';
                mgrotAAv(cf, rep) = mean(mgrotAAd(:, cf, rep), 'omitnan').^2;
                mgrotBBv(cf, rep) = mean(mgrotBBd(:, cf, rep), 'omitnan').^2;
                mgrotAAr(cf, rep) = mean(mgrotAAd(:, cf, rep), 'omitnan').^2 * mhoR(fold, cf, rep).^2;
                mgrotBBr(cf, rep) = mean(mgrotBBd(:, cf, rep), 'omitnan').^2 * mhoR(fold, cf, rep).^2;
            end
            
        end
        
    end

% summarize map/reduced values to return objects
cvU = mapCVstruct(cvU, mcvU);
cvV = mapCVstruct(cvV, mcvV);
trR = mapCVstruct(trR, mtrR);
hoR = mapCVstruct(hoR, mhoR);
grotAAd = mapCVstruct(grotAAd, mgrotAAd);
grotBBd = mapCVstruct(grotBBd, mgrotBBd);
grotAAv = mapCVstruct(grotAAv, mgrotAAv);
grotBBv = mapCVstruct(grotBBv, mgrotBBv);
grotAAr = mapCVstruct(grotAAr, mgrotAAr);
grotBBr = mapCVstruct(grotBBr, mgrotBBr);

end

end

function [ out ] = initCVstruct(x, y, z)
% return a simple structure

% don't print fold stats by defualt
if(~exist('z', 'var') || isempty(z))
    z = 1;
end

out.n = 0;
out.mn = squeeze(zeros(x, y, z));
out.var = squeeze(zeros(x, y, z));
out.std = squeeze(zeros(x, y, z));

end

function [ out ] = mapCVstruct(cv, dat)
% combine a mapped data structure to the next raw subset

% find central tendency across highest dimension of dat
ndim = ndims(dat);

% take current (potentially empty) structure
a_n = cv.n;
a_mean = cv.mn;
a_var = cv.var;
%a_std = cv.std;

% summarize across the next batch
b_n = size(dat, ndim);
b_mean = mean(dat, ndim, 'omitnan');
b_var = var(dat, 1, ndim, 'omitnan');
%b_std = sqrt(b_var);

% compute and return the combination
out.n = a_n + b_n;
out.mn = ((a_mean .* a_n) + (b_mean .* b_n)) ./ out.n;
out.var = (((a_n .* a_var) + (b_n .* b_var)) ./ out.n) + ((a_n .* b_n) .* ((b_mean - a_mean) ./ out.n).^2);
out.std = sqrt(out.var);

end


% %% extimate central tendency / variability across 3rd dimension
% function [ cv, sd ] = ct3d(mat, val)
% % this is a simple wrapper to easily call mean / median values across
% % permutations w/o excessive if cases.
% 
% % grab mean / median based on passed argument (median default)
% switch val
%     case {'mean'}
%         cv = mean(mat, 3, 'omitnan');
%     case {'median'}
%         cv = median(mat, 3, 'omitnan');
%     otherwise
%         cv = median(mat, 3, 'omitnan');
% end
% 
% % always estimate central tendency
% sd = std(mat, [], 3, 'omitnan');
% 
% end

%% validation of sequential estimates of mean / var / sd
% 
% % split a long vector into n parts
% x1 = randi(100, 10, 1);
% x2 = randi(100, 10, 1);
% %x1 = [ 2; 4; 6; 8; 9; 10; 6; 2; 3; 7 ];
% %x2 = [ 1; 3; 5; 7; 4; 6; 4; 4; 9; 2 ];
% 
% % the means are the same
% full_mean = mean([ x1; x2 ]);
% %full_median = median([ x1; x2 ]);
% full_std = std([ x1; x2 ], 1);
% full_var = var([ x1; x2 ], 1);
% 
% % would initialize as empty
% a_n = size(x1, 1);
% a_mean = mean(x1, 'omitnan');
% %a_median = median(x1, 'omitnan');
% a_var = var(x1, 1, 'omitnan');
% 
% % would be the first batch
% b_n = size(x2, 1);
% b_mean = mean(x2, 'omitnan');
% %b_median = median(x2, 'omitnan');
% b_var = var(x2, 1, 'omitnan');
% 
% % compute the combination
% c_n = a_n + b_n;
% c_mean = ((a_mean .* a_n) + (b_mean .* b_n)) ./ c_n;
% %c_median = ((a_median .* a_n) + (b_median .* b_n)) ./ c_n;
% c_var = (((a_n .* a_var) + (b_n .* b_var)) ./ c_n) + ((a_n .* b_n) .* ((b_mean - a_mean) ./ c_n)^2);
% c_std = sqrt(c_var);
% 
% % not boolean precise, but display precision precise
% [ full_mean c_mean full_mean == c_mean ]
% %[ full_median c_median full_median == c_median ]
% [ full_var c_var full_var == c_var ]
% [ full_std c_std full_std == c_std ]
% 
% % update a to be c, next block becomes b 
% % ad nauseum
