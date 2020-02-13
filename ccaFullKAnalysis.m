function [ dat, cca ] = ccaMapNonLinearAnalysis(NET, vars, conf, netsNames, varsNames, confNames, varsLabel, Nkeep1, Nkeep2, Nperm, Nfold, Nrep, ct, kernel, kernelpar, reg)
%[ dat, cca ] = ccaFullAnalysis(NET, vars, conf, netsNames, varsNames, confNames, Nkeep, Nperm);
%   This function takes brain, behavior, and confound matrices of subj x vars
%   size, normalizes them, and performs either a linear canonical correlation analysis or a 
%   kernel based non-linear CCA analysis. If arguments for the non-linear
%   KCCA are passed, then that is the analysis run.
%
%   The KCCA method is implemented from the KMBox tools. The functions have
%   been modified to return internal fields so the cross-validation may be
%   applied to hold-out data.
%
%   INPUTS:
%       NET  - subj x net of brain network data (node measures or upper diagonal) 
%       vars - subj x var of behavior data
%       conf - subj x reg of the confound regressors to be removed from NET and vars
%       netsNames - net x 1 cell array of network variable labels
%       varsNames - var x 1 cell array of behavior variable labels
%       confNames - reg x 1 cell array of confound variable labels
%       varsLabel - var x 1 cell array of behavior label (full text for plots)
%       Nkeep1 - number of PCA components to use for brain
%       Nkeep2 - number of PCA components to use for behavior
%       Nperm - number of permutations for bootstrapping null hypothesis
%       Nfold - number of folds for the cross-validation procedure
%       Nrep - the number of times to repeat the cross-validation procedure
%       ct - central tendency across repeated measures. Either {'mn', 'mean'} or {'md', 'median'}
%       kernel - the type of kernel to apply to data for xform to non-linear space
%                'gauss' (default), 'gauss-diag', 'poly', 'linear'
%       kernelpar - kernel width ('gauss', 'gauss-diag') or 
%                   polynomial order and additive ('poly')
%       reg - regularization to be applied (default = 1E-5);
%   OUTPUTS:
%       dat - any normalized / PCA data stored for each data set
%       cca - all CCA component values of interest for reporting or plotting
%
% References:
% http://www.statsoft.com/Textbook/Canonical-Analysis
% https://stats.idre.ucla.edu/spss/output/canonical-correlation-analysis/
%
% Copyright (c) Brent McPherson (Indiana University), 2020. All rights reserved.
%

%% set default values

if(~exist('Nperm', 'var') || isempty(Nperm))
    Nperm = 0;
end

if(~exist('Nfold', 'var') || isempty(Nfold))
    Nfold = 5;
end

if(~exist('Nrep', 'var') || isempty(Nrep))
    Nrep = 50;
end

if(~exist('ct', 'var') || isempty(ct))
    ct = 'median';
end

% if nonlinear parameters are not passed
if (~exist('kernel', 'var') && ~exist('kernelpar', 'var') && ~exist('reg', 'var'))
    
    % perform a regular, linear CCA
    lcca = true;
    
    % set fields for empty so cvCCA fxn runs linear model too
    kernel = [];
    kernelpar = [];
    reg = [];
    
else
    
    % perform a non-linear KCCA, selecting defaults if they're not passed
    lcca = false;
    
    if(~exist('kernel', 'var') || isempty(kernel))
        kernel = 'gauss';
    end
    
    if(~exist('kernelpar', 'var') || isempty(kernelpar))
        kernelpar = 5;
    end
    
    if(~exist('reg', 'var') || isempty(reg))
        reg = 1E-5;
    end
    
end

% if Nrep < 1000
%     warning('You need at least 1000 reps.');
%     Nrep = 1000;
% end

% if to few permutations are passed, skip them
if Nperm < 2
    disp('Bootstrap permutation testing disabled.');
end

% force the minimum number of PCA components to become canonical factors for global comparisons
Nkeep = min([ Nkeep1, Nkeep2 ]);

% enforce that the input data has the same observations
if size(NET, 1) ~= size(vars, 1)
    error('Input datasets have a different number of observations.');
end

% pull the number of subjects / variables from inputs
Nsub = size(NET, 1);
Nnet = size(NET, 2);

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

% pull the number of variables from corrected data
Nvar = size(varsgrot, 2);

disp('Storing preprocessed data...');

clear amNET badconf badvars conf confNames dd NET NET1 NET3 netsNames ss1 uu vars varsd varsdCOV varsdCOV2 varsLabel varsNames

%% run the analysis

if lcca
    
    disp('Performing CCA on all observations...');
    [ A, B, grotR, grotU, grotV, grotstats ] = canoncorr(uu1, uu2);
   
else
    
    disp('Performing KCCA on all observations...');
    [ grotU, grotV, grotR, alpha1, alpha2, K1, K2 ] = km_kcca(uu1, uu2, kernel, kernelpar, reg, Nkeep);
    
end

% get the full analysis loadings, variability, sign, etc.
fgrotAAd = nan(Nnet, Nkeep);
fgrotBBd = nan(Nvar, Nkeep);
fgrotAAp = nan(Nnet, Nkeep);
fgrotBBp = nan(Nvar, Nkeep);

% for every CC
for ii = 1:Nkeep
        
        % network weights after deconfounding
        fgrotAAd(:, ii) = corr(grotU(:, ii), NETd(:, 1:Nnet))';
        fgrotAAp(:, ii) = fgrotAAd(:, ii) > 0;
        % only use the first normalized network columns, ignore the additional pca confounds
        
        % behavior weights after deconfounding
        fgrotBBd(:, ii) = corr(grotV(:, ii), varsgrot, 'rows', 'pairwise')';
        fgrotBBp(:, ii) = fgrotBBd(:, ii) > 0;

end

% removed measures that I never used and how to compute them
% variability = mean(fgrotAAd(:, ii).^2, 'omitnan');
% redundancy = mean(fgrotAAd(:, ii).^2, 'omitnan') * grotR(ii)^2;

clear ii

% save full data cca fields
cca.full.grotR = grotR;
cca.full.grotU = grotU;
cca.full.grotV = grotV;

% kernal transformed data and alpha scaling
if lcca
    
    cca.full.A = A;
    cca.full.B = B;
    cca.full.stats = grotstats;
    
else
    
    cca.full.alpha1 = alpha1;
    cca.full.alpha2 = alpha2;
    cca.full.K1 = K1;
    cca.full.K2 = K2;
    
end

% save full recovered values
cca.full.fgrotAAd = fgrotAAd;
cca.full.fgrotAAp = logical(fgrotAAp);

cca.full.fgrotBBd = fgrotBBd;
cca.full.fgrotBBp = logical(fgrotBBp);

clear fgrotAAd fgrotAAp fgrotBBd fgrotBBp A B grotR grotstats grotU grotV
clear alpha1 alpha2 K1 K2

%% cross-validation of CCA loadings

% create repeated k-fold estimated loadings
[ cvU, cvV, trR, hoR, grotAAd, grotBBd ] ...
    = cvCCA(uu1, uu2, NETd, varsgrot, Nfold, Nrep, kernel, kernelpar, reg, true);

% save cross-validated fields to output

% pull either mean or median from cv values
switch ct
    
    case {'mn', 'mean'}
    
        cca.dat1.factor = cvU.mn;
        cca.dat1.loading = grotAAd.mn;
       
        cca.dat2.factor = cvV.mn;
        cca.dat2.loading = grotBBd.mn;
        
        cca.cca.trcorrs = trR.mn;
        cca.cca.hocorrs = hoR.mn;
    
    otherwise
        
        cca.dat1.factor = cvU.md;
        cca.dat1.loading = grotAAd.md;
        
        cca.dat2.factor = cvV.md;
        cca.dat2.loading = grotBBd.md;
        
        cca.cca.trcorrs = trR.md;
        cca.cca.hocorrs = hoR.md;
        
end

% pull standard error from cv values
cca.dat1.factor_se = cvU.std ./ sqrt(Nrep);
cca.dat1.loading_se = grotAAd.std ./ sqrt(Nrep);

cca.dat2.factor_se = cvV.std ./ sqrt(Nrep);
cca.dat2.loading_se = grotBBd.std ./ sqrt(Nrep);

cca.cca.trcorrs_se = trR.std ./ sqrt(Nrep);
cca.cca.hocorrs_se = hoR.std ./ sqrt(Nrep);

% proportion of variability accounted for in each CC / percent accounted for in each CC
cca.cca.trRv = mean(trR.mn, 'omitnan') .^ 2;
cca.cca.trPc = (cca.cca.trRv ./ sum(cca.cca.trRv)) * 100;
cca.cca.hoRv = hoR.mn .^ 2;
cca.cca.hoPc = cca.cca.hoRv ./ sum(cca.cca.hoRv) * 100;

% grab the permutation parameters
cca.param.Nkeep1 = Nkeep1;
cca.param.Nkeep2 = Nkeep2;
cca.param.Nperm = Nperm;
cca.param.Nfold = Nfold;
cca.param.Nrep = Nrep;
cca.param.ct = ct;
cca.param.kernel = kernel;
cca.param.kernelpar = kernelpar;
cca.param.reg = reg;

clear cvU cvV trR hoR grotAAd grotBBd

%% CCA permutation testing for parameters

if Nperm >= 2
    
    disp('Performing a permutation test to determine significant number of CCs...');
    
    % preallocate output
    grotRr = zeros(2, Nkeep, Nperm);
    grotUtr = zeros(Nsub, Nkeep, Nperm);
    grotVtr = zeros(Nsub, Nkeep, Nperm);
    grotAtr = zeros(Nnet, Nkeep, Nperm);
    grotBtr = zeros(Nvar, Nkeep, Nperm);
    
    % predefine permutation sets 
    % add 1 to Nperm b/c the first one isn't randomized
    PAPset = palm_quickperms([], ones(Nsub, 1), Nperm+1);
    
    disp('Reconstructing null parameters from random permutations...');
    
    % for every permutation
    for ii = 1:Nperm
        
        % get cross-validated loadings for each permutation
        % use ii+1 to skip first, original ordering of data, only use random sortings
        % only run 5 repeats for null, noise is noise afterall... (?)
        [ grotUr, ~, ~, thoR1 ] = cvCCA(uu1, uu2(PAPset(:, ii+1), :), NETd, varsgrot, Nfold, 5, kernel, kernelpar, reg);
        [ ~, grotVr, ~, thoR2 ] = cvCCA(uu1(PAPset(:, ii+1), :), uu2, NETd, varsgrot, Nfold, 5, kernel, kernelpar, reg);
                
        % for every cc
        for jj = 1:Nkeep
            
            switch ct
                
                case {'mn', 'mean'}
                    
                    % pull null factor scores
                    grotUtr(:, :, ii) = grotUr.mn;
                    grotVtr(:, :, ii) = grotVr.mn;
                    
                    % null network / behavior loadings after deconfounding
                    grotAtr(:, jj, ii) = corr(grotUr.mn(:, jj), NETd(:, 1:size(NET, 2)))';
                    grotBtr(:, jj, ii) = corr(grotVr.mn(:, jj), varsgrot, 'rows', 'pairwise')';
                    
                    % catch the holdout correlation
                    grotRr(1, jj, ii) = thoR1.mn(jj);
                    grotRr(2, jj, ii) = thoR2.mn(jj);
        
                otherwise
                    
                    % pull null factor scores
                    grotUtr(:, :, ii) = grotUr.md;
                    grotVtr(:, :, ii) = grotVr.md;
                    
                    % null network / behavior loadings after deconfounding
                    grotAtr(:, jj, ii) = corr(grotUr.md(:, jj), NETd(:, 1:Nnet))';
                    grotBtr(:, jj, ii) = corr(grotVr.md(:, jj), varsgrot, 'rows', 'pairwise')';

                    % catch the holdout correlation
                    grotRr(1, jj, ii) = thoR1.md(jj);
                    grotRr(2, jj, ii) = thoR2.md(jj);
                    
            end            
        end
    end
    
    clear ii jj grotUr grotVr thoR1 thoR2
    
    % collapse null correlations to a single null observation
    grotRr = mean(grotRr, 1, 'omitnan');

    % count # of null correlations greater than trained hold-outs
    pval = nan(1, Nkeep);
    for ii = 1:Nkeep
        pval(ii) = sum(grotRr(1, ii, :) > cca.cca.hocorrs(ii)) / Nperm;
    end
    
    % determine number of significant FWE elements
    Ncca = sum(pval < 0.05); % number of FWE-significant CCA components
    
    % save fields to output
    cca.dat1.factor_null = grotUtr;
    cca.dat1.loading_null = grotAtr;
    
    cca.dat2.factor_null = grotVtr;
    cca.dat2.loading_null = grotBtr;
    
    cca.cca.ncca = Ncca;
    cca.cca.pval = pval;
    cca.cca.hocorrs_null = mean(grotRr, 1, 'omitnan');
    
else
    
    % fill in fields as empty
    cca.dat1.factor_null = [];
    cca.dat1.loading_null = [];
    
    cca.dat2.factor_null = [];
    cca.dat2.loading_null = [];
    
    cca.cca.ncca = [];
    cca.cca.pval = [];   
    cca.cca.hocorrs_null = [];
    
end

end

%% cross-validation function for compartmentalization
function [ cvU, cvV, trR, hoR, grotAAd, grotBBd ] = cvCCA(uu1, uu2, NETd, varsgrot, Nfold, Nrep, kernel, kernelpar, reg, verb)

if(~exist('verb', 'var') || isempty(verb))
    verb = false;
end

% if the kernel parameters are empty, cv linear cca
if (isempty(kernel) || isempty(kernelpar) || isempty(reg))
    if verb
        disp('Cross-validating linear CCA model...');
    end
    lcca = true; 
else
    if verb
        disp('Cross-validating non-linear KCCA model...');
    end
    lcca = false;
end

% fake the number kept by size of input columns
Nkeep = min([ size(uu1, 2) size(uu2, 2) ]);

% preallocate data for loop
rcvU = nan(size(uu1, 1), Nkeep, Nrep);
rcvV = nan(size(uu2, 1), Nkeep, Nrep);
rtrR = nan(Nfold, Nkeep, Nrep);
rhoR = nan(1, Nkeep, Nrep);
rgrotAAd = nan((size(NETd, 2)/2), Nkeep, Nrep);
rgrotBBd = nan(size(varsgrot, 2), Nkeep, Nrep);

% for every repeat requested
for rep = 1:Nrep
    
    % create the requested number of folds across all subjects
    cvlab = crossvalind('Kfold', 1:size(uu1, 1), Nfold);
    
    % for every fold within a repeat
    for fold = 1:Nfold
        
        % make a copy of full data for dropping values
        truu1 = uu1;
        truu2 = uu2;
        
        % pull the hold-out indices
        hoidx = cvlab == fold;
        
        % remove the hold-out from the training data
        truu1(hoidx, :) = [];
        truu2(hoidx, :) = [];
        
        if lcca
            
            % estimate the weights and corr from the training data
            [ cvA, cvB, rtrR(fold, :, rep) ] = canoncorr(truu1, truu2);
            
            % get the estimates for the hold-out subjects
            rcvU(hoidx, :, rep) = uu1(hoidx, :) * cvA;
            rcvV(hoidx, :, rep) = uu2(hoidx, :) * cvB;
            
        else
            
            % train model on a subset of observations
            [ ~, ~, rtrR(fold, :, rep), ta1, ta2, ~, ~, n0 ] = km_kcca(truu1, truu2, kernel, kernelpar, reg, Nkeep);
            
            % build test filter space - these are identical, based on subject number
            [ ~, ~, ~, ~, ~, ~, ~, tnx ] = km_kcca(uu1(hoidx, :), uu1(hoidx, :), kernel, kernelpar, reg, Nkeep);
            [ ~, ~, ~, ~, ~, ~, ~, tny ] = km_kcca(uu2(hoidx, :), uu2(hoidx, :), kernel, kernelpar, reg, Nkeep);
            
            % rebuild the estimates for the hold out observations
            rcvU(hoidx, :, rep) = (n0*km_kernel(truu1, uu1(hoidx, :), kernel, kernelpar)*tnx)'*ta1;
            rcvV(hoidx, :, rep) = (n0*km_kernel(truu2, uu2(hoidx, :), kernel, kernelpar)*tny)'*ta2;
            
        end
        
    end
    
    % grab the validated correlations b/w the hold-out and rebuilt loadings
    for cf = 1:Nkeep
        rhoR(1, cf, rep) = corr(rcvU(:, cf, rep), rcvV(:, cf, rep));
        rgrotAAd(:, cf, rep) = corr(rcvU(:, cf, rep), NETd(:, 1:(size(NETd, 2)/2)))';
        rgrotBBd(:, cf, rep) = corr(rcvV(:, cf, rep), varsgrot, 'rows', 'pairwise')';
    end    
    
end

% summarize the full values to return objects compatible w/ map-reduce loop

cvU.n = Nrep;
cvU.mn = mean(rcvU, 3, 'omitnan');
cvU.md = median(rcvU, 3, 'omitnan');
cvU.var = var(rcvU, 1, 3, 'omitnan');
cvU.std = sqrt(cvU.var);

cvV.n = Nrep;
cvV.mn = mean(rcvV, 3, 'omitnan');
cvV.md = median(rcvV, 3, 'omitnan');
cvV.var = var(rcvV, 1, 3, 'omitnan');
cvV.std = sqrt(cvV.var);

trR.n = Nrep;
trR.mn = mean(rtrR, 3, 'omitnan');
trR.md = median(rtrR, 3, 'omitnan');
trR.var = var(rtrR, 1, 3, 'omitnan');
trR.std = sqrt(trR.var);

hoR.n = Nrep;
hoR.mn = mean(rhoR, 3, 'omitnan');
hoR.md = median(rhoR, 3, 'omitnan');
hoR.var = var(rhoR, 1, 3, 'omitnan');
hoR.std = sqrt(hoR.var);

grotAAd.n = Nrep;
grotAAd.mn = mean(rgrotAAd, 3, 'omitnan');
grotAAd.md = median(rgrotAAd, 3, 'omitnan');
grotAAd.var = var(rgrotAAd, 1, 3, 'omitnan');
grotAAd.std = sqrt(grotAAd.var);

grotBBd.n = Nrep;
grotBBd.mn = mean(rgrotBBd, 3, 'omitnan');
grotBBd.md = median(rgrotBBd, 3, 'omitnan');
grotBBd.var = var(rgrotBBd, 1, 3, 'omitnan');
grotBBd.std = sqrt(grotBBd.var);

end

%% kcca development notes

% this is the right size, but can't be applied to hold out data
% I added alpha1/2, K1/K2, N0 output to try and reconstruct loadings
% but, their size is dependent on input - unclear if this generalized to
% subsets or could be used for train test at all.

% FIGURE OUT HOW TO APPLY FIT MODEL TO NEW OBSERVATIONS
%grotU = K1*alpha1;

% % split training / test
% rng(47405);
% trni = randi(594, round(594*.75), 1);
% tsti = setdiff(1:594, trni)';
% 
% % create training / test 
% trnx = uu1(trni, :);
% trny = uu2(trni, :);
% 
% tstx = uu1(tsti, :);
% tsty = uu2(tsti, :);
% 
% % train model on a subset of observations
% [ ~, ~, ~, ta1, ta2, tk1, tk2, n0 ] = km_kcca(trnx, trny, kernel, kernelpar, reg, Nkeep);
% 
% % build test filter space - these are identical, based on subject number
% [ ~, ~, ~, ~, ~, ~, ~, tnx ] = km_kcca(tstx, tstx, kernel, kernelpar, reg, Nkeep);
% [ ~, ~, ~, ~, ~, ~, ~, tny ] = km_kcca(tsty, tsty, kernel, kernelpar, reg, Nkeep);
% 
% % rebuild on hold out observations (?)
% % use training kernel to scale train to test, match scaling with kernel
% % from training set for common space (?). Scale by trained alphas -
% % these are the mappings to the filter.
% % Is this actually applying transform to a hold out?
% hox = n0*km_kernel(trnx, tstx, 'gauss', 5)*tnx;
% zzz1 = hox'*ta1; % not sure why this needs a transpose now
% 
% hoy = n0*km_kernel(trny, tsty, 'gauss', 5)*tny;
% zzz2 = hoy'*ta2; % not sure why this needs a transpose now
% 
% % inspect
% corr(zzz1(:, 1), zzz2(:, 1))
% plot(zzz1(:, 1), zzz2(:, 1), '.');

% transform data with kernel; data to itself for data in kernal space
% this is an exponential chi^2 - it doesn't work 
%[train_a_ker, omegaA] = kernel_expchi2(uu1, uu1);
%[train_b_ker, omegaB] = kernel_expchi2(uu2, uu2);

% % the number of subjects
% M = size(uu1, 1);
% 
% % build the empty filters
% Kx = diag(ones(M, 1));
% Ky = diag(ones(M, 1));
% 
% % apply x/y datasets - apply to self / training data, pass
% X = uu1;
% Y = uu2;
% 
% % create Gram matrix with polynomial filter
% p = 1;
% for i = 1:M
%     for j = 1:M
%         % do X and Y need to be crossed instead of together? Is that how
%         % train/test can be applied?
%         Kx(i, j) = power((X(i, :) * X(j, :)'), p);
%         Ky(i, j) = power((Y(i, :) * Y(j, :)'), p);
%     end
% end
% 
% % used to apply nonlinear space to hold out
% omega_Kx = mean(Kx(:));
% omega_Ky = mean(Ky(:));

% apply trained kernal to test data - take omega from train set and apply
% it to hold out w/ train data for correct transform
%[test_a_ker] = kernel_expchi2(test_a', train_a', omegaA);
%[test_b_ker] = kernel_expchi2(test_b', train_b', omegaB);

% apply kcca to transformed data
% [ nalpha, nbeta, r, Kx_c, Ky_c ] = kcanonca_reg_ver2(Kx, Ky, 1, 0.5, 0, 0);
% is always empty...

% project to test data?
%[train_a_ker,test_a_ker,train_b_ker,test_b_ker] = center_kcca(train_a_ker,test_a_ker,train_b_ker,test_b_ker);
%test_b_ker_proj = test_b_ker*Wx;
%test_a_ker_proj = test_a_ker*Wy;

% score for evaluating performance?
%score_kcca = pdist2(test_b_ker_proj,test_a_ker_proj,'cosine');
%score_nn = pdist2(test_b',test_a','euclidean');

% based on demo, these are reasonable params for non-linear kernal CCA
%[ W1, b1, W2, b2, P1, P2, m1, m2, D ]  = randKCCA(uu1, uu2, Nkeep, 5000, 5, 1E-5, 5000, 5, 1E-5);
%[ ~, ~, ~, ~, grotU, grotV, ~, ~, grotR ]  = randKCCA(uu1, uu2, Nkeep, 5000, 5, 1E-5, 5000, 5, 1E-5);
% this changes the size drastically, but can be applied to hold out data
