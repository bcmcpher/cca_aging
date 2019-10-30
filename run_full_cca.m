%% run the analysis and make all the plots
%

%% set up the analysis

% define the number of domains and permutations
Nkeep = 100;
Nperm = 10001;

% path to output directory
outdir = '/N/dc2/projects/lifebid/HCP/Brent/camcan/figs/panels_20190906/';

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
varsQconf = [ age, varsQconf ];
confNames = [ 'Age', confNames ];

% determine the stem of this output
%stem = 'ho_age_deg_';
%stem = 'bh_age_deg_';
stem = 'rg_age_deg_';

% is this a full set of edges? 1 = yes; 0 = no
edge = 0;

%% run the analysis

[ dat, cca ] = ccaFullAnalysis(deg, vars, varsQconf, ...
                               netNames, varsNames, confNames, ...
                               varsLabel, Nkeep, Nperm);

%% get all the data / analyses up to this point

% crossvalidate the cca
[ crossval.Rmean, crossval.Rstd, crossval.pval ] = ccaCrossvalidate(dat, .8, Nperm);

% estimate the multiway correlation of age with CCA axes
mcrR = nan(100, 1);
mcrS = nan(100, 1);
mcpval = nan(100, 1);
for ii = 1:100
    [ mcrR(ii), mcrS(ii), mcpval(ii) ] = ccaLinRegCorr(cca, ii, age, Nperm);
end

corrAge.rmean = mcrR;
corrAge.rstd = mcrS;
corrAge.pval = mcpval;

clear ii

%% make plots showing explained variance in PCAs

fh = ccaPcaVarianceExplained(dat.dat1.ss1, 1, 'Brain');
print([ outdir stem 'brain_pca_variance_explained.eps' ], '-painters', '-depsc'); 
close all; clear fh

fh = ccaPcaVarianceExplained(dat.dat2.ss2, 1, 'Behavior');
print([ outdir stem 'behavior_pca_variance_explained.eps' ], '-painters', '-depsc'); 
close all; clear fh

%% determine the significance of the cc

fh = ccaPermSigTestPlot(cca, 1);
print([ outdir stem 'cca1_sig_test.eps' ], '-painters', '-depsc'); 
close all; clear fh

fh = ccaPermSigTestPlot(cca, 2);
print([ outdir stem 'cca2_sig_test.eps' ], '-painters', '-depsc'); 
close all; clear fh

%% plot the main cca finding

fh = ccaPlotAxisCon(cca, 1, age, parula(88));
print([ outdir stem 'cca1_finding.eps' ], '-painters', '-depsc'); 
close all; clear fh

fh = ccaPlotAxisCon(cca, 2, age, parula(88));
print([ outdir stem 'cca2_finding.eps' ], '-painters', '-depsc'); 
close all; clear fh

%% all the ranked trends

[ brain.loading, brain.ld_bd ] = ccaPlotRankedTrends(dat, cca, age, 'brain', 'load', 1, 30);
set(gcf, 'Position', [ 150 175 1200 850 ]);
print([ outdir stem 'cca1_brain_loading_sorted.eps' ], '-painters', '-depsc'); 
close all; 

[ brain.slope, brain.sp_bd ] = ccaPlotRankedTrends(dat, cca, age, 'brain', 'slope', 1, 30);
set(gcf, 'Position', [ 150 175 1200 850 ]);
print([ outdir stem 'cca1_brain_slope_sorted.eps' ], '-painters', '-depsc'); 
close all; 

[ behavior.loading, behavior.ld_bd ] = ccaPlotRankedTrends(dat, cca, age, 'behavior', 'load', 1, 30);
set(gcf, 'Position', [ 150 175 1200 850 ]);
print([ outdir stem 'cca1_behavior_loading_sorted.eps' ], '-painters', '-depsc'); 
close all; 

[ behavior.slope, behavior.sp_ld ] = ccaPlotRankedTrends(dat, cca, age, 'behavior', 'slope', 1, 30);
set(gcf, 'Position', [ 150 175 1200 850 ]);
print([ outdir stem 'cca1_behavior_slope_sorted.eps' ], '-painters', '-depsc'); 
close all; 

%% project data to cca trend axis

[ projection.vals, projection.sd, fh ] = ccaProjectFactor(cca, 1, age, Nperm);
print([ outdir stem 'cca1_projected_axis.eps' ], '-painters', '-depsc'); 
close all; clear fh

%% get variance explained in cca axes

plotVarianceExplained(cca, 'dat1');
print([ outdir stem 'brain_cca_variance_explained.eps' ], '-painters', '-depsc'); 
close all; 

plotVarianceExplained(cca, 'dat2');
print([ outdir stem 'behavior_cca_variance_explained.eps' ], '-painters', '-depsc'); 
close all;

%% get the mean / std of the R2 for the individual data sets

% get the trend through the data
[ trends.brn_lcoeff, trends.brn_lr2 ] = ccaTrendByDecile(dat.dat1.raw, age);
[ trends.beh_lcoeff, trends.beh_lr2 ] = ccaTrendByDecile(dat.dat2.raw, age);

% reasampled central tendency
rr2_brn = nan(Nperm, 2);
rr2_beh = nan(Nperm, 2);

for ii = 1:Nperm
    
    % pull a random sample
    tmp_brn = randsample(trends.brn_lr2, size(trends.brn_lr2, 1), 'true');
    tmp_beh = randsample(trends.beh_lr2, size(trends.beh_lr2, 1), 'true');
    
    % get mean / std of each sample
    rr2_brn(ii, 1) = mean(tmp_brn);
    rr2_brn(ii, 2) = std(tmp_brn);
    
    rr2_beh(ii, 1) = mean(tmp_beh);
    rr2_beh(ii, 2) = std(tmp_beh);
    
end

clear ii tmp_brn tmp_beh

% the resampled average / std of the R2 for each individual dataset
rr2.brain = mean(rr2_brn);
rr2.behavior = mean(rr2_beh);

%% optional / depends stuff

if edge
    
    bmat = ccaEdgeWeights(dat.dat1.loading, 1);
    figure; imagesc(bmat); hold on;
    axis square; axis equal; axis tight; colorbar;
    set(gca, 'XTick', [], 'YTick', []); caxis([ -1 1 ]); 
    title('Edge-wise loading of all network edges in CCA');
    print([ outdir stem 'brain_cca_edge_loadings.eps' ], '-painters', '-depsc'); 
    close all;    
    
    edges = bmat;
    
else

    % create dissimilarity b/w all variables in CCA
    [ dissimilarity.dmat, dissimilarity.cmat ] = ccaDissimilarityMatrix(cca);
    dmat = dissimilarity.dmat;
    
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
    
    dissimilarity.y7_axis = y7_axis;
    dissimilarity.bounds = dlin;
    
    figure; imagesc(dissimilarity.dmat); hold on;
    axis square; axis equal; axis tight; colorbar;
    set(gca, 'XTick', [], 'YTick', []); caxis([ 0 1 ]); 
    title('Dissimilarity between all variables in CCA');
    print([ outdir stem 'cca_dissimilarity_unsorted.eps' ], '-painters', '-depsc'); 
    close all;    
    
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
    
    % rescale for brain only
    set(gca, 'XLim', [ 1 376 ], 'YLim', [ 1 376 ]);
    print([ outdir stem 'cca_dissimilarity_brain_sorted.eps' ], '-painters', '-depsc'); 
    dissimilarity.rbrn = dissimilarity.dmat(1:376, 1:376);
    
    % rescale for behavior only
    set(gca, 'XLim', [ 377 710 ], 'YLim', [ 377 710 ]);
    print([ outdir stem 'cca_dissimilarity_behavior_sorted.eps' ], '-painters', '-depsc'); 
    dissimilarity.rbeh = dissimilarity.dmat(377:710, 377:710);
    
    % rescale for brain x behavior
    set(gca, 'XLim', [ 1 376 ], 'YLim', [ 377 710 ]);
    print([ outdir stem 'cca_dissimilarity_brnxbeh_sorted.eps' ], '-painters', '-depsc'); 
    close all;    
    dissimilarity.rbtw = dissimilarity.dmat(377:710, 1:376);
       
    % plot the averaged together values
    [ dissimilarity.mdDat, ~, dissimilarity.mdBtw ] = fnModuleDensity(dissimilarity.dmat, [ yeoLabs.yeo7; ib+10 ], 'mean');
    figure; imagesc(dissimilarity.mdDat);
    axis square; axis equal; axis tight; colorbar; caxis([ 0 1 ]);
    title('Dissimiliarity Between Brain and Task Domains');
    mdLabs = [ yeoLabs.yeo7Names'; S ];
    set(gca, 'XTick', 1:size(dissimilarity.mdDat, 1), 'XTickLabels', mdLabs, 'XTickLabelRotation', 45, ...
        'YTick', 1:size(dissimilarity.mdDat, 1), 'YTickLabels', mdLabs);
    
    print([ outdir stem 'cca_dissimilarity_modules.eps' ], '-painters', '-depsc'); 
    
    set(gca, 'XLim', [ 0.5 10.5 ], 'YLim', [ 0.5 10.5 ]); caxis([ 0 1 ]);
    print([ outdir stem 'cca_dissimilarity_brn_modules.eps' ], '-painters', '-depsc'); 
    
    set(gca, 'XLim', [ 10.5 17.5 ], 'YLim', [ 10.5 17.5 ]); caxis([ 0 1 ]);
    print([ outdir stem 'cca_dissimilarity_beh_modules.eps' ], '-painters', '-depsc'); 
    
    set(gca, 'XLim', [ 0.5 10.5 ], 'YLim', [ 10.5 17.5 ]); caxis([ 0.7 1 ]);
    print([ outdir stem 'cca_dissimilarity_brnxbeh_modules.eps' ], '-painters', '-depsc'); 
    close all;    
    
    dissimilarity.mbrn = dissimilarity.mdDat(1:10, 1:10);
    dissimilarity.mbeh = dissimilarity.mdDat(11:17, 11:17);
    dissimilarity.mbtw = dissimilarity.mdDat(11:17, 1:10);
    
end

%% save the final output

save([ stem 'data.mat' ], 'dat', 'cca', 'brain', 'behavior', ... 
     'corrAge', 'crossval', 'dissimilarity', 'projection', 'rr2');
