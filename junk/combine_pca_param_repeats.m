%% compile together the pca parameter sweep

% pick the folder to load param sweep data from
%vers = 'osg';
%vers = 'osg-mean';
vers = 'osg-repeat';

% define the number of permutations
Nperm = 10001;

% build the path to the data files
datadir = [ '/N/dc2/projects/lifebid/HCP/Brent/camcan/param_sweep/' vers ];

% build the list of files
files = dir([ datadir '/*.mat' ]);
files = {files.name}';

% get the n of files to merge
nrep = size(files, 1);

% load first subject for object sizes and prepped data
tdat = load(fullfile(datadir, files{1}));
Nkeep = min(size(tdat.cca.dat1.factor, 2), size(tdat.cca.dat2.factor, 2));
uu1 = tdat.dat.dat1.uu1;
uu2 = tdat.dat.dat2.uu2;
NET = tdat.dat.dat1.raw;
NETd = tdat.dat.dat1.nrm;
varsgrot = tdat.dat.dat2.nrm;
dat = tdat.dat;
clear tdat

% preallocate output fields
fl_A = nan(38, 38, nrep);
fl_B = nan(40, 38, nrep);
fl_grotR = nan(1, 38, nrep);
fl_grotU = nan(594, 38, nrep);
fl_grotV = nan(594, 38, nrep);

d1_factor = nan(594, 38, nrep);
d1_loading = nan(376, 38, nrep);
d1_variability = nan(38, 1, nrep);
d1_redundancy = nan(38, 1, nrep);

d2_factor = nan(594, 38, nrep);
d2_loading = nan(334, 38, nrep);
d2_variability = nan(38, 1, nrep);
d2_redundancy = nan(38, 1, nrep);

cc_trcorr = nan(5, 38, nrep);
cc_hocorr = nan(5, 38, nrep);
cc_trRv = nan(1, 38, nrep);
cc_trPc = nan(1, 38, nrep);
cc_hoRv = nan(1, 38, nrep);
cc_hoPc = nan(1, 38, nrep);

% for every repeat
for ii = 1:nrep
    
    % grab the file
    file = files{ii};
    
    % load the file
    tdat = load(fullfile(datadir, file));
    
    % fill in the preallocated outputs for each repeat
    fl_A(:, :, ii) = tdat.cca.full.A;
    fl_B(:, :, ii) = tdat.cca.full.B;
    fl_grotR(:, :, ii) = tdat.cca.full.grotR;
    fl_grotU(:, :, ii) = tdat.cca.full.grotU;
    fl_grotV(:, :, ii) = tdat.cca.full.grotV;
    
    d1_factor(:, :, ii) = tdat.cca.dat1.factor;
    d1_loading(:, :, ii) = tdat.cca.dat1.loading;
    d1_variability(:, :, ii) = tdat.cca.dat1.variability;
    d1_redundancy(:, :, ii) = tdat.cca.dat1.redundancy;
    
    d2_factor(:, :, ii) = tdat.cca.dat2.factor;
    d2_loading(:, :, ii) = tdat.cca.dat2.loading;
    d2_variability(:, :, ii) = tdat.cca.dat2.variability;
    d2_redundancy(:, :, ii) = tdat.cca.dat2.redundancy;
    
    cc_trcorr(:, :, ii) = tdat.cca.cca.trcorrs;
    cc_hocorr(:, :, ii) = tdat.cca.cca.hocorrs;
    cc_trRv(:, :, ii) = tdat.cca.cca.trRv;
    cc_trPc(:, :, ii) = tdat.cca.cca.trPc;
    cc_hoRv(:, :, ii) = tdat.cca.cca.hoRv;
    cc_hoPc(:, :, ii) = tdat.cca.cca.hoPc;
    
end

clear ii file tdat

% create the average / std object across all repeats
cca.full.A = mean(fl_A, 3);
cca.full.A_sd = std(fl_A, 1, 3);
cca.full.B = mean(fl_B, 3);
cca.full.B_sd = std(fl_B, 1, 3);
cca.full.grotR = mean(fl_grotR, 3);
cca.full.grotR_sd = std(fl_grotR, 1, 3);
cca.full.grotU = mean(fl_grotU, 3);
cca.full.grotU_sd = std(fl_grotU, 1, 3);
cca.full.grotV = mean(fl_grotV, 3);
cca.full.grotV_sd = std(fl_grotV, 1, 3);
cca.full.stats = nan;

cca.dat1.factor = mean(d1_factor, 3);
cca.dat1.factor_sd = std(d1_factor, 1, 3);
cca.dat1.loading = mean(d1_loading, 3);
cca.dat1.loading_sd = std(d1_loading, 1, 3);
cca.dat1.variability = mean(d1_variability, 3);
cca.dat1.variability_sd = std(d1_variability, 1, 3);
cca.dat1.redundancy = mean(d1_redundancy, 3);
cca.dat1.reduncancy_sd = std(d1_redundancy, 1, 3);

cca.dat2.factor = mean(d2_factor, 3);
cca.dat2.factor_sd = std(d2_factor, 1, 3);
cca.dat2.loading = mean(d2_loading, 3);
cca.dat2.loading_sd = std(d2_loading, 1, 3);
cca.dat2.variability = mean(d2_variability, 3);
cca.dat2.variability_sd = std(d2_variability, 1, 3);
cca.dat2.redundancy = mean(d2_redundancy, 3);
cca.dat2.reduncancy_sd = std(d2_redundancy, 1, 3);

cca.cca.trcorrs = mean(cc_trcorr, 3);
cca.cca.trcorrs_sd = std(cc_trcorr, 1, 3);
cca.cca.hocorrs = mean(cc_hocorr, 3);
cca.cca.hocorrs_sd = std(cc_hocorr, 1, 3);
cca.cca.trRv = mean(cc_trRv, 3);
cca.cca.trRv_sd = std(cc_trRv, 1, 3);
cca.cca.trPc = mean(cc_trPc, 3);
cca.cca.trPc_sd = std(cc_trPc, 1, 3);
cca.cca.hoRv = mean(cc_hoRv, 3);
cca.cca.hoRv_sd = std(cc_hoRv, 1, 3);
cca.cca.hoPc = mean(cc_hoPc, 3);
cca.cca.hoPc_sd = std(cc_hoPc, 1, 3);

%% add permuatation test to merged results

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
PAPset = palm_quickperms([], ones(594, 1), Nperm);

% for every permutation
for ii = 2:Nperm
    
    % get cross-validated loadings for each permutation
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
        grotAr(ii, jj) = mean(grotAtr(ii, jj, :) .^ 2, 'omitnan') * grotRp(ii, jj, 1);
        grotBr(ii, jj) = mean(grotBtr(ii, jj, :) .^ 2, 'omitnan') * grotRp(ii, jj, 2);
        
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

%% clear utilized / redundant fields

clear cc_hocorr cc_hoPc cc_hoRv cc_trcorr cc_trPc cc_trRv 
clear d1_factor d1_loading d1_redundancy d1_variability
clear d2_factor d2_loading d2_redundancy d2_variability
clear fl_A fl_B fl_grotR fl_grotU fl_grotV
clear grotAr grotAtr grotAv grotBr grotBtr grotBv grotRp grotRpval Ncca1 Ncca2
clear uu1 uu2 varsgrot PAPset NET NETd

%% save the merged result

save([ '/N/dc2/projects/lifebid/HCP/Brent/camcan/param_sweep/combined_' vers '_n' num2str(nrep) '.mat' ], 'cca', 'dat', '-v7.3');

%% plot the results

% correlation with age
figure('Position', [ 100 500 1300 500 ]);
subplot(1, 2, 1);
imagesc(dat(:, :, 1)); colorbar;
axis equal; axis square; axis tight;
title('Correlation with Age by PCA Component');
xlabel('Brain PCA Components');
ylabel('Behavior PCA Components');
subplot(1, 2, 2);
imagesc(dat(:, :, 2)); colorbar;
axis equal; axis square; axis tight;
title('SE of Correlation with Age by PCA Component');
xlabel('Brain PCA Components');
ylabel('Behavior PCA Components');

% correlation between datasets
figure('Position', [ 100 500 1300 500 ]);
subplot(1, 2, 1);
imagesc(dat(:, :, 3)); colorbar;
axis equal; axis square; axis tight;
title('Correlation between Brain and Behavior by PCA Component');
xlabel('Brain PCA Components');
ylabel('Behavior PCA Components');
subplot(1, 2, 2);
imagesc(dat(:, :, 4)); colorbar;
axis equal; axis square; axis tight;
title('SE of Correlation Between Brain and Behavior by PCA Component');
xlabel('Brain PCA Components');
ylabel('Behavior PCA Components');

% variable loading variability
figure('Position', [ 100 500 1300 500 ]);
subplot(1, 2, 1);
imagesc(dat(:, :, 5)); colorbar;
axis equal; axis square; axis tight;
title('Variability of Loadings by PCA Component');
xlabel('Brain PCA Components');
ylabel('Behavior PCA Components');
subplot(1, 2, 2);
imagesc(dat(:, :, 6)); colorbar;
axis equal; axis square; axis tight;
title('SE of Variability of Loadings by PCA Component');
xlabel('Brain PCA Components');
ylabel('Behavior PCA Components');

figure('Position', [ 100 500 1300 500 ]);
subplot(1, 2, 1);
imagesc(dat(:, :, 7)); colorbar;
axis equal; axis square; axis tight;
title('Mean of Loadings by PCA Component');
xlabel('Brain PCA Components');
ylabel('Behavior PCA Components');
subplot(1, 2, 2);
imagesc(dat(:, :, 8)); colorbar;
axis equal; axis square; axis tight;
title('SE of Mean of Loadings by PCA Component');
xlabel('Brain PCA Components');
ylabel('Behavior PCA Components');

% plot the surfaces
figure; 
ax1 = subplot(1, 3, 1);
surf(dat(:,:,1));
ax2 = subplot(1, 3, 2);
surf(dat(:,:,3));
ax3 = subplot(1, 3, 3);
surf(dat(:,:,5));
Link = linkprop([ ax1, ax2, ax3 ],{'CameraUpVector', 'CameraPosition', 'CameraTarget', ...
                'XLim', 'YLim', 'ZLim'});
setappdata(gcf, 'StoreTheLink', Link);

% find the ideal parameters?
ii = 38;
jj = 40;

[ dat(ii, jj, 1), dat(ii, jj, 3), dat(ii, jj, 5); dat(ii, jj, 2), dat(ii, jj, 4), dat(ii, jj, 6)]

ii = 57;
jj = 48;

[ dat(ii, jj, 1), dat(ii, jj, 3), dat(ii, jj, 5); dat(ii, jj, 2), dat(ii, jj, 4), dat(ii, jj, 6)]

%% load the repeats at 100k / 250k

cdat = [];

figure('Position', [ 50 550 1500 400 ]);
for ii = 1:5
    
    subplot(1, 5, ii);
    tfile = sprintf('param_sweep/local-01/pca_brain_038_behavior_040_250k_rep%02d.mat', ii);
    tdat = load(tfile, 'mdDat1', 'dat');
    imagesc(tdat.mdDat1); colorbar; caxis([ 0 1 ]);
    title([ 'Repeat: ' num2str(ii) ]); ndim = size(tdat.mdDat1, 1);
    axis square; axis equal; axis tight;
    set(gca, 'XTick', [], 'YTick', []);
    %set(gca, 'XTick', 1:ndim, 'XTickLabels', mdLabs, 'XTickLabelRotation', 45, ...
    %    'YTick', 1:ndim, 'YTickLabels', mdLabs);
    set(gca, 'XLim', [ 0.5 10.5 ], 'YLim', [ 10.5 17.5 ]);
    cdat = cat(3, cdat, tdat.mdDat1);
    
end

clear ii
dat = tdat.dat;

% get index of sorted behaviors
svar = regexprep(dat.dat2.names', '_.*', '');
svar = regexprep(svar, 'hint', 'comp');
svar{1} = 'comp'; % replace the first dumb label
[ S, ~, ib ] = unique(svar);

% load the yeo labels
load('~/hcp_mmp_vol/working/yeoLabs.mat');

mdLabs = [ yeoLabs.yeo7Names'; S ];

figure;
subplot(1, 2, 1)
imagesc(mean(cdat, 3)); colorbar; caxis([ 0 1 ]);
title('Mean of 10 Repeated Calls');
axis square; axis equal; axis tight;
set(gca, 'XTick', 1:size(cdat, 1), 'XTickLabels', mdLabs, 'XTickLabelRotation', 45, ...
    'YTick', 1:size(cdat, 1), 'YTickLabels', mdLabs);
set(gca, 'XLim', [ 0.5 10.5 ], 'YLim', [ 10.5 17.5 ]);

subplot(1, 2, 2)
imagesc(std(cdat, [], 3)); colorbar; caxis([ 0 0.2 ]);
title('SD of 10 Repeated Calls');
axis square; axis equal; axis tight;
set(gca, 'XTick', 1:size(cdat, 1), 'XTickLabels', mdLabs, 'XTickLabelRotation', 45, ...
    'YTick', 1:size(cdat, 1), 'YTickLabels', mdLabs);
set(gca, 'XLim', [ 0.5 10.5 ], 'YLim', [ 10.5 17.5 ]);

%% load and plot average modules

% preallocate data for output
amdDat1 = nan(17, 17, 1000);
amdDat2 = nan(17, 17, 1000);

for ii = 1:1000
    
    file = sprintf('pca_brain_038_behavior_040_10k_rep%04d.mat', ii);
    
    % if the repeat does not exist
    if exist(fullfile(datadir, file), 'file')
        
        tdat = load(fullfile(datadir, file));
        amdDat1(:,:,ii) = tdat.mdDat1;
        amdDat2(:,:,ii) = tdat.mdDat2;
                
    end
    
end

% cross-validated loadings
figure;
subplot(1, 2, 1)
imagesc(mean(amdDat1, 3, 'omitnan')); colorbar; caxis([ 0 0.50 ]);
title('Mean of 971 Repeated Calls');
axis square; axis equal; axis tight;
set(gca, 'XTick', 1:size(amdDat1, 1), 'XTickLabels', mdLabs, 'XTickLabelRotation', 45, ...
    'YTick', 1:size(amdDat1, 1), 'YTickLabels', mdLabs);
set(gca, 'XLim', [ 0.5 10.5 ], 'YLim', [ 10.5 17.5 ]);

subplot(1, 2, 2)
imagesc(std(amdDat1, 1, 3, 'omitnan')); colorbar; caxis([ 0 0.05 ]);
title('SD of 971 Repeated Calls');
axis square; axis equal; axis tight;
set(gca, 'XTick', 1:size(amdDat1, 1), 'XTickLabels', mdLabs, 'XTickLabelRotation', 45, ...
    'YTick', 1:size(amdDat1, 1), 'YTickLabels', mdLabs);
set(gca, 'XLim', [ 0.5 10.5 ], 'YLim', [ 10.5 17.5 ]);

% cross-validated factors
figure;
subplot(1, 2, 1)
imagesc(mean(amdDat2, 3, 'omitnan')); colorbar; caxis([ 0 1 ]);
title('Mean of 971 Repeated Calls');
axis square; axis equal; axis tight;
set(gca, 'XTick', 1:size(amdDat2, 1), 'XTickLabels', mdLabs, 'XTickLabelRotation', 45, ...
    'YTick', 1:size(amdDat2, 1), 'YTickLabels', mdLabs);
set(gca, 'XLim', [ 0.5 10.5 ], 'YLim', [ 10.5 17.5 ]);

subplot(1, 2, 2)
imagesc(std(amdDat2, 1, 3, 'omitnan')); colorbar; caxis([ 0 0.07 ]);
title('SD of 971 Repeated Calls');
axis square; axis equal; axis tight;
set(gca, 'XTick', 1:size(amdDat2, 1), 'XTickLabels', mdLabs, 'XTickLabelRotation', 45, ...
    'YTick', 1:size(amdDat2, 1), 'YTickLabels', mdLabs);
set(gca, 'XLim', [ 0.5 10.5 ], 'YLim', [ 10.5 17.5 ]);

%% look at a bunch of these

figure;
count = 1;
start = 1;
for ii = 1:5
    for jj = 1:5
        subplot(5, 5, count)
        imagesc(amdDat2(:, :, start)); colorbar; caxis([ 0 1 ]);
        title([ 'Repeat ' num2str(start) ]);
        axis square; axis equal; axis tight;
        set(gca, 'XTick', [], 'YTick', []);
        %set(gca, 'XTick', 1:size(amdDat2, 1), 'XTickLabels', mdLabs, 'XTickLabelRotation', 45, ...
        %    'YTick', 1:size(amdDat2, 1), 'YTickLabels', mdLabs);
        set(gca, 'XLim', [ 0.5 10.5 ], 'YLim', [ 10.5 17.5 ]);
        count = count + 1;
        start = start + 1;
        
    end
end


%% find the empty indices and rebuild a subset job section for osg

fileID = fopen('param_sweep/osg-repeat_redo3.txt','w');
for ii = 1:1000
    
    file = sprintf('pca_brain_038_behavior_040_10k_rep%04d.mat', ii);
    
    % if the repeat does not exist
    if ~exist(fullfile(datadir, file), 'file')
        
        disp([ 'Re-doing Repeat: ' num2str(ii) ]);
        fprintf(fileID, '\narguments = %04d pca_raw_data.mat\noutput = logs/rep_%04d_stdout.log\nerror = logs/rep_%04d_stderr.log\nlog = logs/run_%04d_condor.out\nqueue\n', ii, ii, ii, ii);
        
    end
    
end
fclose(fileID);

clear ii jj fileID file

%% run the two that are left 

for ii = 1:1000
    
    file = sprintf('pca_brain_038_behavior_040_10k_rep%04d.mat', ii);
    
    % if the repeat does not exist
    if ~exist(fullfile(datadir, file), 'file')
        
        disp([ 'Re-doing Repeat: ' num2str(ii) ]);
        pcaParamSweep('38', '40', num2str(ii), 'pca_raw_data.mat', 'param_sweep/osg-repeat/');
        
    end
    
end
