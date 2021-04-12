function [ out, dat, cca ] = pcaParamSweep(pca1, pca2, rep, dat, outpath)
%pcaParamSweep(pca1, pca2); 
%   simple function to parallelize iterations of pca pairs for parameter
%   tuning correlation b/w data, w/ age, and testing dissimilarity.
%
% Copyright (c) Brent McPherson (Indiana University), 2019. All rights reserved.
%

%% load everything

% % path to output directory
% outdir = '/N/dc2/projects/lifebid/HCP/Brent/camcan/param_sweep/';
% 
% % load the edge measure
% load('/N/dc2/projects/lifebid/HCP/Brent/camcan/canoncorr_analysis_full_data.mat', 'deg');
% 
% % load the rest of the behavior and labels
% load('/N/dc2/projects/lifebid/HCP/Brent/camcan/canoncorr_analysis_full_data.mat', 'age', 'vars', 'varsQconf', ...
%      'netNames', 'varsNames', 'confNames');
% 
% % replace varsNames w/ corresponding labels
% load('/N/dc2/projects/lifebid/HCP/Brent/camcan/camcan_vars_labels.mat', 'varsLabels');
% 
% % get the indices of intersection between the datasets
% vidx = contains(varsLabels(:, 1), varsNames);
% 
% % grab plaintext names in the right order
% varsLabel = varsLabels(vidx, 2);
% 
% clear vidx
% 
% % load the yeo labels
% load('/N/dc2/projects/lifebid/HCP/Brent/camcan/yeoLabs.mat', 'yeoLabs');

pca1 = str2num(pca1);
pca2 = str2num(pca2);
rep = str2num(rep);

% load the whole workspace
load(dat);

% randomly set the seed differently for each repeat
rng(seed(rep));

% merge age into confounds?
%varsQconf = [ age, varsQconf ];
%confNames = [ 'Age', confNames ];

% parse output name from pca sweeps
% stem = sprintf('pca_brain_%03d_behavior_%03d_250k_rep%04d.mat', pca1, pca2, rep);
stem = sprintf('pca_brain_%03d_behavior_%03d_15k.mat', pca1, pca2);

disp([ 'Writing to output file: ' stem ]);

%% run the iteration

% preallocate the simple data output
out = nan(2, 2);

% run the cca for the parameter vals
[ dat, cca ] = ccaFullKAnalysis(deg, vars, varsQconf, ...
                                netNames, varsNames, confNames, varsLabels, ...
                                pca1, pca2, 0, 5, 15000);

% boostrap mean / sd of the correlation with age
[ out(1, 1), out(1, 2) ] = ccaLinRegCorr(cca, 1, age, 1000);

% pull the mean / sd of the dataset correlation
out(2, 1) = cca.cca.hocorrs(1);
out(2, 2) = cca.cca.hocorrs_se(1);

% create dissimilarity b/w all variables in CCA
dmat = ccaDissimilarityMatrix(cca);

% get index of sorted behaviors
svar = regexprep(dat.dat2.names', '_.*', '');
svar = regexprep(svar, 'hint', 'comp');
svar{1} = 'comp'; % replace the first dumb label
[ ~, ~, ib ] = unique(svar);

mdDat = fnModuleDensity(dmat, [ yeoLabs.yeo7; ib+10 ], 'mean');

%% save the data down

save([ outpath stem ], 'dat', 'cca', 'out', 'dmat', 'mdDat', '-v7.3');

disp('Done.');

end
