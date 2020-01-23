%% load fit data and make the plots

% path to output directory
outdir = '/N/dc2/projects/lifebid/HCP/Brent/camcan/figs/panels_20200122/';

% load the rest of the behavior and labels
load('canoncorr_analysis_full_data.mat', 'age');

% load the yeo labels
load('~/hcp_mmp_vol/working/yeoLabs.mat');

% % hold out age
% stem = 'ho_age';
% load('pca_brain_038_behavior_040_250k_ho_age.mat', 'dat', 'cca', 'dmat1', 'mdDat1');

% regress out age
stem = 'rg_age';
load('pca_brain_038_behavior_040_250k_rg_age.mat', 'dat', 'cca', 'dmat1', 'mdDat1');

% simple reassign of outputs
dmat = dmat1;
mdDat = mdDat1;
clear dmat1 mdDat1
                               
%% get corr with age

[ mcrR1, mcrS1, mcpval1 ] = ccaLinRegCorr(cca, 1, age, 10000);
[ mcrR2, mcrS2, mcpval2 ] = ccaLinRegCorr(cca, 2, age, 10000);

figure; hold on

plot([ 1 1 ], [ (mcrR1 - mcrS1) (mcrR1 + mcrS1) ], 'color', 'black');
plot(1, mcrR1, 'o', 'MarkerFaceColor', [.1 .45 .95], ...
     'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'MarkerSize', 12);

plot([ 2 2 ], [ (mcrR2 - mcrS2) (mcrR2 + mcrS2) ], 'color', 'black');
plot(2, mcrR2, 'o', 'MarkerFaceColor', [.1 .45 .95], ...
     'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'MarkerSize', 12);

hold off

set(gca, 'XLim', [ 0 3 ], 'XTick', 1:2, 'XTickLabel', {'CC1', 'CC2'});
title('Correlation with Age for CCA Components 1 & 2');
print([ outdir stem '_corr_with_age.eps' ], '-painters', '-depsc'); close all

%% plot of the main axis

ccaPlotAxisCon(cca, 1, age, parula(88), true);
print([ outdir stem '_cca1_finding.eps' ], '-painters', '-depsc'); close all

ccaPlotAxisCon(cca, 2, age, parula(88), true);
print([ outdir stem '_cca2_finding.eps' ], '-painters', '-depsc'); close all

%% plot the ranked loadings to see if they're consistent

ccaPlotRankedTrends(dat, cca, age, 'brain', 'load', 1, 'line', 30);
set(gcf, 'Position', [ 150 150 1300 875 ]);
print([ outdir stem '_cca1_brain_loading.eps' ], '-painters', '-depsc'); close all

ccaPlotRankedTrends(dat, cca, age, 'behavior', 'load', 1, 'line', 30);
set(gcf, 'Position', [ 150 150 1300 875 ]);
print([ outdir stem '_cca1_behavior_loading.eps' ], '-painters', '-depsc'); close all

%% plot the points / error bars of all variable loadings

ccaPlotRankedTrends(dat, cca, age, 'brain', 'load', 1, 'points');
set(gcf, 'Position', [ 150 150 1300 875 ]);
print([ outdir stem '_cca1_brain_loading_points.eps' ], '-painters', '-depsc'); close all

ccaPlotRankedTrends(dat, cca, age, 'behavior', 'load', 1, 'points');
set(gcf, 'Position', [ 150 150 1300 875 ]);
print([ outdir stem '_cca1_behavior_loading_points.eps' ], '-painters', '-depsc'); close all

%% create dissimilarity

% create dissimilarity b/w all variables in CCA
%dmat = ccaDissimilarityMatrix(cca);

% sort the brain regions
[ y7_lab, y7_brn ] = sort(yeoLabs.yeo7);

% get index of sorted behaviors
svar = regexprep(dat.dat2.names', '_.*', '');
svar = regexprep(svar, 'hint', 'comp');
svar{1} = 'comp'; % replace the first dumb label
[ S, ~, ib ] = unique(svar);
[ si, sv ] = sort(ib);

% pull label names
mdLabs = [ yeoLabs.yeo7Names'; S ];

% build single axis index
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
print([ outdir stem '_cca_dissimilarity_sorted.eps' ], '-painters', '-depsc');
close all;

% plot the averaged together values
%[ mdDat, ~, mdBtw ] = fnModuleDensity(dmat, [ yeoLabs.yeo7; ib+10 ], 'mean');
figure; imagesc(mdDat);
axis square; axis equal; axis tight; colorbar; caxis([ 0 1 ]);
title('Dissimiliarity Between Brain and Task Domains');
set(gca, 'XTick', 1:size(mdDat, 1), 'XTickLabels', mdLabs, 'XTickLabelRotation', 45, ...
    'YTick', 1:size(mdDat, 1), 'YTickLabels', mdLabs);
print([ outdir stem '_cca_dissimilarity_modules.eps' ], '-painters', '-depsc');

set(gca, 'XLim', [ 0.5 10.5 ], 'YLim', [ 10.5 17.5 ]); caxis([ 0 .70 ]);
print([ outdir stem '_cca_dissimilarity_brnxbeh_modules.eps' ], '-painters', '-depsc');
close all;

%% misc extras

% 'word cloud' size contribution
ccaScaledContribution(cca, 'behavior', 1, si, S)
set(gcf, 'Position', [ 250 125 1300 775 ]);
print([ outdir stem '_cca1_behavior_domain_contribution.eps' ], '-painters', '-depsc'); close all

% 'word cloud' size contribution
ccaScaledContribution(cca, 'brain', 1, yeoLabs.yeo7, yeoLabs.yeo7Names');
set(gcf, 'Position', [ 250 125 1300 775 ]);
print([ outdir stem '_cca1_brain_domain_contribution.eps' ], '-painters', '-depsc'); close all

% null distribution of CCA
null = mean(squeeze(cca.cca.grotRp(:, 1, :))');
figure; hold on;
hist(null, 64);
print([ outdir stem '_cca_corr_null_distribution.eps' ], '-painters', '-depsc'); close all;
%plot([ cca.full.grotR(1) cca.full.grotR(1) ], [ 0 600 ], 'color', 'red');

%% NEED A FULL MODULE CHORD PLOT
% useful to describe change of within brain / between behavior when age is removed

% assign placeholder name
mat = mdDat;
thr = 0.75;

% assign semi-transparent colors for all edges, color by domain
cmap = [ repmat([.1 .45 .95 .5], 10, 1); repmat([.95 .1 .25 .5 ], 7, 1) ];

% readable labels
clabs = {'Visual', 'Somatomotor', 'Dorsal Attention', 'Ventral Attention', ...
         'Limbic', 'Frontoparietal', 'DMN', 'Subcortical', 'Hippocampus', 'Amygdala', ...
         'Attention', 'Clinical', 'Emotional', 'Language', 'Memory', 'Motor', 'Social'};

% pull the min / max
mn = min(mat(:));
mx = max(mat(:));

% create the output
out = (mat - mn) / (mx - mn);

% threshold
out(out > thr) = 0;
     
% plot the labeled graph
circularGraph(out, 'Colormap', cmap, 'label', clabs);

%% module significance?

% brain / behavior group labels for modules
%grp = [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2 ];

% run simple significance fxn
%[ Mden, fh, Mprm ] = fnModuleSignificance(dmat, [ yeoLabs.yeo7; ib+10 ], grp, 10000);

% use bootstrap test to draw significance
%[ Mden, fh, Mprm ] = fnModulePermutationSignificance(dmat, ???, [ yeoLabs.yeo7; ib+10 ], grp);
% second input should be 3d mat of null dissimilarity matix?
