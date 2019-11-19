function [ out, dat, cca ] = pcaParamSweep(pca1, pca2, dat, outpath)
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

% load the whole workspace
load(dat);

% merge age into confounds?
%varsQconf = [ age, varsQconf ];
%confNames = [ 'Age', confNames ];

% parse output name from pca sweeps
stem = sprintf('pca_brain_%03d_behavior_%03d_test.mat', pca1, pca2);

%% run the iteration

% preallocate the simple data output
out = nan(4, 2);

% run the cca for the parameter vals
[ dat, cca ] = ccaTestFullAnalysis(deg, vars, varsQconf, ...
                                   netNames, varsNames, confNames, varsLabel, ...
                                   pca1, pca2, 2, true, 'mean', 5, 1000);

% boostrap mean / sd of the correlation with age
[ out(1, 1), out(1, 2) ] = ccaLinRegCorr(cca, 1, age, 1000);

% pull the mean / sd of the dataset correlation
cdat = squeeze(cca.cca.hocorrs(:, 1, :));
out(2, 1) = mean(cdat(:));
out(2, 2) = std(cdat(:));

% pull the mean / sd of the loading standard deviation
out(3, 1) = mean([ cca.dat1.loading_sd(:, 1); cca.dat2.loading_sd(:, 1) ]);
out(3, 2) = std([ cca.dat1.loading_sd(:, 1); cca.dat2.loading_sd(:, 1) ]);

% pull the mean / sd of the loadings
out(4, 1) = mean([ cca.dat1.loading(:, 1); cca.dat2.loading(:, 1) ]);
out(4, 2) = std([ cca.dat1.loading(:, 1); cca.dat2.loading(:, 1) ]);

% create dissimilarity b/w all variables in CCA
dmat = ccaDissimilarityMatrix(cca);

% sort the brain regions
%[ y7_lab, y7_brn ] = sort(yeoLabs.yeo7);

% get index of sorted behaviors
svar = regexprep(dat.dat2.names', '_.*', '');
svar = regexprep(svar, 'hint', 'comp');
svar{1} = 'comp'; % replace the first dumb label
[ ~, ~, ib ] = unique(svar);
%[ si, sv ] = sort(ib);

%y7_axis = [ y7_brn; sv+376 ];

% lines splitting domains
%dlin = [ find(diff(y7_lab)); find(diff(si))+376 ];
%dlin = sort([ dlin; 376 ]);

mdDat = fnModuleDensity(dmat, [ yeoLabs.yeo7; ib+10 ], 'mean');

% figure; imagesc(mdDat);
% axis square; axis equal; axis tight; colorbar; caxis([ 0 1 ]);
% title('Dissimiliarity Between Brain and Task Domains');
% mdLabs = [ yeoLabs.yeo7Names'; S ];
% set(gca, 'XTick', 1:size(mdDat, 1), 'XTickLabels', mdLabs, 'XTickLabelRotation', 45, ...
%     'YTick', 1:size(mdDat, 1), 'YTickLabels', mdLabs);
% set(gca, 'XLim', [ 0.5 10.5 ], 'YLim', [ 10.5 17.5 ]); caxis([ 0.5 1 ]);
% print([ outdir 'figs/' stem '_dissimilarity.eps' ], '-painters', '-depsc');
% close all; clc

%% save the data down

%save([ outdir stem '.mat' ], 'dat', 'cca', 'out', 'dmat', 'mdDat');
save([ outpath stem ], 'dat', 'cca', 'out', 'dmat', 'mdDat');

disp('Done.');

end

