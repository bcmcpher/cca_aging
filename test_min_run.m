%% look at a specific subset of options

% define the number of domains and permutations
Nperm = 2;

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

%% run the cross validated analysis testign for the number of PCA components

% the number of PCA components to test
vals = [ 5 10:20 50 75 100 ];
Nkeep = length(vals);

% preallocate output for parameter sweep
outm = nan(Nkeep, Nkeep, 6);
maxC = 0;
bdat = [];
bcca = [];

% for every combination of values up to Nkeep
for ii = 2:Nkeep
    for jj = 2:Nkeep
        
        % run the cross-validated analysis for every pca pair
        [ dat, cca ] = ccaTestFullAnalysis(deg, vars, varsQconf, ...
                                           netNames, varsNames, confNames, varsLabel, ...
                                           vals(ii), vals(jj), Nperm, true, 5, 1000);

        % catch the first crossvalidated loadings cc correlation
        outm(ii, jj, 1) = corr(cca.dat1.factor(:, 1), cca.dat2.factor(:, 1));
        
        % catch the mean / sd crossvalidated correlations across all holdouts for the first cc
        tcorr = cca.cca.hocorrs(:, 1, :);
        outm(ii, jj, 2) = mean(tcorr(:), 'omitnan');
        outm(ii, jj, 3) = std(tcorr(:), 'omitnan');
        
        % estimate correlation / sd / significance with age
        [ outm(ii, jj, 4), outm(ii, jj, 5), outm(ii, jj, 6) ] = ccaLinRegCorr(cca, 1, age, 1000);
        
        % main plot of figure
        if outm(ii, jj, 4) > maxC
        
            % catch the best dat / cca objects
            bdat = dat;
            bcca = cca;
        
            % plot the results
            ccaPlotAxisCon(cca, 1, age, parula(88));
            title([ 'Best CCA - Brain: ' num2str(vals(ii)) ' Behavior: ' num2str(vals(jj)) ]);
            print('pca_best_sweep.eps', '-painters', '-depsc'); 
            close all;
            
        end
    end
end

% save the output b/c this takes a long time
save('pca_param_min_sweep.mat', 'outm', 'bdat', 'bcca');

%% inspect correlation b/w data sets

figure;
subplot(1, 2, 1);
imagesc(outm(:,:,2)); colorbar;
axis equal; axis square; axis tight
title('CC1 Mean Correlation of Data by PCA Components');
xlabel('Brain PCA Components'); ylabel('Behavior PCA Components');
set(gca, 'XTick', 1:Nkeep, 'YTick', 1:Nkeep, 'XTickLabel', vals, 'YTickLabel', vals);

subplot(1, 2, 2);
imagesc(outm(:,:,3)); colorbar;
axis equal; axis square; axis tight
title('CC1 SD Correlation of Data by PCA Components');
xlabel('Brain PCA Components'); ylabel('Behavior PCA Components');
set(gca, 'XTick', 1:Nkeep, 'YTick', 1:Nkeep, 'XTickLabel', vals, 'YTickLabel', vals);

%% inspect correlation w/ age

figure;
subplot(1, 2, 1);
imagesc(outm(:,:,4)); colorbar;
axis equal; axis square; axis tight
title('CC1 Mean Correlation with Age by PCA Components');
xlabel('Brain PCA Components'); ylabel('Behavior PCA Components');
set(gca, 'XTick', 1:Nkeep, 'YTick', 1:Nkeep, 'XTickLabel', vals, 'YTickLabel', vals);

subplot(1, 2, 2);
imagesc(outm(:,:,5)); colorbar;
axis equal; axis square; axis tight
title('CC1 SD Correlation with Age by PCA Components');
xlabel('Brain PCA Components'); ylabel('Behavior PCA Components');
set(gca, 'XTick', 1:Nkeep, 'YTick', 1:Nkeep, 'XTickLabel', vals, 'YTickLabel', vals);

%% quickly check dissimilarity

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

% fit into modules
[ mdDat, ~, mdBtw ] = fnModuleDensity(dmat, [ yeoLabs.yeo7; ib+10 ], 'mean');
figure; imagesc(mdDat);
axis square; axis equal; axis tight; colorbar; caxis([ 0 1 ]);
title('Dissimiliarity Between Brain and Task Domains');
mdLabs = [ yeoLabs.yeo7Names'; S ];
set(gca, 'XTick', 1:size(mdDat, 1), 'XTickLabels', mdLabs, 'XTickLabelRotation', 45, ...
    'YTick', 1:size(mdDat, 1), 'YTickLabels', mdLabs);
