%% just run simple example all the way through

% define the number of domains and permutations
Nkeep = 50;
Nperm = 2;

% path to output directory
outdir = '/N/dc2/projects/lifebid/HCP/Brent/camcan/figs/panels_20191114/';

% load the edge measure
load('canoncorr_analysis_full_data.mat', 'deg');

% load the rest of the behavior and labels
load('canoncorr_analysis_full_data.mat', 'age', 'vars', 'varsQconf', ...
     'netNames', 'varsNames', 'confNames');

% replace varsNames w/ corresponding labels
load('camcan_vars_labels.mat', 'varsLabels');

% get the indices of intersection between the datasets
vidx = contains(varsLabels(:, 1), varsNames);

% grab plaintext names in the right order
varsLabel = varsLabels(vidx, 2);

clear vidx

% load the yeo labels
load('~/hcp_mmp_vol/working/yeoLabs.mat');
 
% merge age into behavior / confounds?
%varsQconf = [ age, varsQconf ];
%confNames = [ 'Age', confNames ];

stem = 'ho_age_';
%stem = 'rg_age_';

%% run the cv cca

% 38, 40; % max corr w/ age
% 57, 48; % min fano factor

[ dat, cca ] = ccaTestFullAnalysis(deg, vars, varsQconf, ...
                                   netNames, varsNames, confNames, varsLabel, ...
                                   38, 40, Nperm, true, 'median', 5, 1000);

%% get corr with age

[ mcrR1, mcrS1, mcpval1 ] = ccaLinRegCorr(cca, 1, age, 1000);
[ mcrR2, mcrS2, mcpval2 ] = ccaLinRegCorr(cca, 2, age, 1000);

%% plot of the main axis

ccaPlotAxisCon(cca, 1, age, parula(88), true);
%print([ outdir stem 'test_cca1_finding.eps' ], '-painters', '-depsc');

ccaPlotAxisCon(cca, 2, age, parula(88));
%print([ outdir stem 'test_cca2_finding.eps' ], '-painters', '-depsc');

%% plot the ranked loadings to see if they're consistent

ccaPlotRankedTrends(dat, cca, age, 'brain', 'load', 1, 30);
%print([ outdir stem 'test_brain_loading.eps' ], '-painters', '-depsc');

ccaPlotRankedTrends(dat, cca, age, 'behavior', 'load', 1, 30);
%print([ outdir stem 'test_behavior_loading.eps' ], '-painters', '-depsc');

%% create dissimilarity

% create dissimilarity b/w all variables in CCA
dmat = ccaDissimilarityMatrix(cca);

% sort the brain regions
[ y7_lab, y7_brn ] = sort(yeoLabs.yeo7);

% get index of sorted behaviors
svar = regexprep(dat.dat2.names', '_.*', '');
svar = regexprep(svar, 'hint', 'comp');
svar{1} = 'comp'; % replace the first dumb label
[ S, ~, ib ] = unique(svar);
[ si, sv ] = sort(ib);

y7_axis = [ y7_brn; sv+376 ];

% lines splitting domains
dlin = [ find(diff(y7_lab)); find(diff(si))+376 ];
dlin = sort([ dlin; 376 ]);

% visualize all - plotAdjacencyMatrix zeros out negatives %[231 102 1187 959]
figure; imagesc(dmat(y7_axis, y7_axis)); hold on;
axis square; axis equal; axis tight; colorbar;
set(gca, 'XTick', [], 'YTick', []); caxis([ 0 1 ]);
title('Dissimilarity between all variables in CCA');
for ii = 1:size(dlin, 1)
    line([ dlin(ii)+.5 dlin(ii)+.5 ], [ .5 size(dmat, 1)+.5 ], 'color', 'black', 'LineWidth', .25);
    line([ .5 size(dmat, 1)+.5 ], [ dlin(ii)+.5 dlin(ii)+.5 ], 'color', 'black', 'LineWidth', .25);
end
hold off;
%print([ outdir stem 'cca_dissimilarity_sorted.eps' ], '-painters', '-depsc');
%close all;

% plot the averaged together values
[ mdDat, ~, mdBtw ] = fnModuleDensity(dmat, [ yeoLabs.yeo7; ib+10 ], 'mean');
figure; imagesc(mdDat);
axis square; axis equal; axis tight; colorbar; caxis([ 0 1 ]);
title('Dissimiliarity Between Brain and Task Domains');
mdLabs = [ yeoLabs.yeo7Names'; S ];
set(gca, 'XTick', 1:size(mdDat, 1), 'XTickLabels', mdLabs, 'XTickLabelRotation', 45, ...
    'YTick', 1:size(mdDat, 1), 'YTickLabels', mdLabs);
%print([ outdir stem 'cca_dissimilarity_modules.eps' ], '-painters', '-depsc');

set(gca, 'XLim', [ 0.5 10.5 ], 'YLim', [ 10.5 17.5 ]); caxis([ 0.5 1 ]);
%print([ outdir stem 'cca_dissimilarity_brnxbeh_modules.eps' ], '-painters', '-depsc');
%close all;
