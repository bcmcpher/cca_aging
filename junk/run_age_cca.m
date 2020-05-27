%% compare CCA run w/o and w/ age regressed out
%
%

%% prep data

% load the edge measure
load('canoncorr_analysis_full_data.mat', 'deg');

% load the rest of the behavior and labels
load('canoncorr_analysis_full_data.mat', 'age', 'vars', 'varsQconf', ...
     'netNames', 'varsNames', 'confNames');

% load the yeo labels
load('~/hcp_mmp_vol/working/yeoLabs.mat');

% sort the brain regions
[ y7_lab, y7_brn ] = sort(yeoLabs.yeo7);

% define the number of domains and permutations
Nkeep = 100;
Nperm = 101;

%% create confounds w/ and w/o age

% confounds w/o
conf = varsQconf;

% confounds w/ age
acon = [ age, varsQconf ];
aconNames = [ 'age', confNames ];

%% run the CCAs

% original
[ odat, occa ] = ccaFullAnalysis(deg, vars, conf, ...
                                 netNames, varsNames, confNames, ...
                                 Nkeep, Nperm);

% regress out age                             
[ adat, acca ] = ccaFullAnalysis(deg, vars, acon, ...
                                 netNames, varsNames, aconNames, ...
                                 Nkeep, Nperm);

%% plot the cca axes to see correlation w/o and w/ age regressed out
               
ccaPlotAxisCon(occa, 2, age, parula(88));

ccaPlotAxisCon(acca, 2, age, parula(88));

%% finish axes labels

% get index of sorted behaviors
svar = regexprep(odat.dat2.names', '_.*', '');
svar = regexprep(svar, 'hint', 'comp');
svar{1} = 'comp'; % replace the first dumb label
[ S, ~, ib ] = unique(svar);
[ si, sv ] = sort(ib);

% create sorted axis labels
y7_axis = [ y7_brn; sv+376 ];

%% create modules

% create dissimilarity b/w all variables in CCA
[ odmat, ocmat ] = ccaDissimilarityMatrix(occa);
[ admat, acmat ] = ccaDissimilarityMatrix(acca);

% plot the averaged together values
[ omdDat, ~, omdBtw ] = fnModuleDensity(odmat, [ yeoLabs.yeo7; ib+10 ], 'mean');
[ amdDat, ~, amdBtw ] = fnModuleDensity(admat, [ yeoLabs.yeo7; ib+10 ], 'mean');

%% plot the dissimilarity

% lines splitting domains
dlin = [ find(diff(y7_lab)); find(diff(si))+376 ];
dlin = sort([ dlin; 376 ]);
nval = size(odmat, 1);

figure; 
imagesc(odmat(y7_axis, y7_axis)); hold on;
axis square; axis equal; axis tight; colorbar;
set(gca, 'XTick', [], 'YTick', []);
caxis([ 0 1 ]); %colormap(redbluecmap);
title('Dissimilarity between all variables in CCA');
for ii = 1:size(dlin, 1)
    line([ dlin(ii)+.5 dlin(ii)+.5 ], [ .5 nval+.5 ], 'color', 'black', 'LineWidth', .25);
    line([ .5 nval+.5 ], [ dlin(ii)+.5 dlin(ii)+.5 ], 'color', 'black', 'LineWidth', .25);
end
hold off;

%% plot the modules

figure('Position', [ 300 100 1500 775 ]); 

subplot(1, 2, 1);
imagesc(omdDat);
axis square; axis equal; axis tight; colorbar; caxis([ 0 1 ]);
title('Dissimiliarity Between Brain and Task Domains w/o Age');
mdLabs = [ yeoLabs.yeo7Names'; S ];
set(gca, 'XTick', 1:size(omdDat, 1), 'XTickLabels', mdLabs, 'XTickLabelRotation', 45, ...
    'YTick', 1:size(omdDat, 1), 'YTickLabels', mdLabs);

subplot(1, 2, 2);
imagesc(amdDat);
axis square; axis equal; axis tight; colorbar; caxis([ 0 1 ]);
title('Dissimiliarity Between Brain and Task Domains w/ Age Regressed');
mdLabs = [ yeoLabs.yeo7Names'; S ];
set(gca, 'XTick', 1:size(amdDat, 1), 'XTickLabels', mdLabs, 'XTickLabelRotation', 45, ...
    'YTick', 1:size(amdDat, 1), 'YTickLabels', mdLabs);

% plot the difference
plotDifferenceMatrix(amdDat, omdDat, [ -0.3 0.3 ]);
title('Change in Dissimiliarity Between Brain and Task Domains w/ Age Regressed');
set(gca, 'XTick', 1:size(amdDat, 1), 'XTickLabels', mdLabs, 'XTickLabelRotation', 45, ...
    'YTick', 1:size(amdDat, 1), 'YTickLabels', mdLabs);
