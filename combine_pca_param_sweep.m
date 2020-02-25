%% compile together the pca parameter sweep

% pick the folder to load param sweep data from
%vers = 'osg';
vers = 'osg-mean-perm';

datadir = [ '/N/dc2/projects/lifebid/HCP/Brent/camcan/param_sweep/' vers ];

% if requested output exists 
if exist([ '/N/dc2/projects/lifebid/HCP/Brent/camcan/param_sweep/combined_' vers '.mat' ], 'file')
    % load it into the workspace
    load([ '/N/dc2/projects/lifebid/HCP/Brent/camcan/param_sweep/combined_' vers '.mat' ]);
else
    % preallocate the data output
    dat = nan(100, 100, 2);
end

% for every pca pair
for ii = 2:100
    for jj = 2:100
        
        % build the file
        file = fullfile(datadir, sprintf('pca_brain_%03d_behavior_%03d_15k.mat', ii, jj));
        
        % if the file exists and the data isn't filled in
        if exist(file, 'file') && isnan(dat(ii, jj, 1))
            
            % try to load the data - skip if it fails
            try
                out = load(file, 'out');
            catch
                disp([ 'Bad file - brain ' num2str(ii) ' behavior ' num2str(jj) ]);
                continue
            end
            
            % fill in the relevant fields
            dat(ii, jj, 1) = out.out(1, 1);
            dat(ii, jj, 2) = out.out(1, 2);
            dat(ii, jj, 3) = out.out(2, 1);
            dat(ii, jj, 4) = out.out(2, 2);
%             dat(ii, jj, 5) = out.out(3, 1);
%             dat(ii, jj, 6) = out.out(3, 2);
%             dat(ii, jj, 7) = out.out(4, 1);
%             dat(ii, jj, 8) = out.out(4, 2);
            clear out
            
        end
        
    end
end

clear ii jj

% append data back output
if exist([ '/N/dc2/projects/lifebid/HCP/Brent/camcan/param_sweep/combined_' vers '.mat' ], 'file')
    save([ '/N/dc2/projects/lifebid/HCP/Brent/camcan/param_sweep/combined_' vers '.mat' ], 'dat', '-append');
else
    save([ '/N/dc2/projects/lifebid/HCP/Brent/camcan/param_sweep/combined_' vers '.mat' ], 'dat');
end

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

%% load the repeats at 10k

cdat = [];

figure('Position', [ 50 550 1500 400 ]);
for ii = 1:10
    subplot(2, 5, ii);
    tfile = sprintf('param_sweep/local/pca_brain_038_behavior_040_10k_rep%02d.mat', ii);
    try
        tdat = load(tfile, 'mdDat', 'dat');
        imagesc(tdat.mdDat); colorbar; caxis([ 0.5 1 ]);
        title([ 'Repeat: ' num2str(ii) ]); ndim = size(tdat.mdDat, 1);
        axis square; axis equal; axis tight;
        set(gca, 'XTick', [], 'YTick', []);
        %set(gca, 'XTick', 1:ndim, 'XTickLabels', mdLabs, 'XTickLabelRotation', 45, ...
        %    'YTick', 1:ndim, 'YTickLabels', mdLabs);
        set(gca, 'XLim', [ 0.5 10.5 ], 'YLim', [ 10.5 17.5 ]);
        cdat = cat(3, cdat, tdat.mdDat);
    catch
        continue
    end
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
imagesc(mean(cdat, 3)); colorbar; caxis([ 0.5 1 ]);
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

%% find the empty indices and rebuild a subset job section for osg

fileID = fopen('param_sweep/osg_redo.txt','w');
redo = [];
for ii = 2:100
    for jj = 2:100
        
        % if the file exists and the data isn't filled in
        if isnan(dat(ii, jj, 1))
            
            %disp([ 'Brain: ' num2str(ii) ' - Behavior: ' num2str(jj) ]);
            redo = [ redo; [ ii, jj ] ];
            
            fprintf(fileID, '\narguments = %03d %03d pca_raw_data.mat\noutput = logs/run_%03d_%03d_stdout.log\nerror = logs/run_%03d_%03d_stderr.log\nlog = logs/run_%03d_%03d_condor.out\nqueue\n', ii, jj, ii, jj, ii, jj, ii, jj)

        end
        
    end
end
fclose(fileID);

clear ii jj fileID

%% check dissimilarity

% build the indices to sort the modules in the final panel
[ y7_lab, y7_brn ] = sort(yeoLabs.yeo7);
svar = regexprep(dat.dat2.names', '_.*', '');
svar = regexprep(svar, 'hint', 'comp');
svar{1} = 'comp'; % replace the first dumb label
[ S, ~, ib ] = unique(svar);
mdLabs = [ yeoLabs.yeo7Names'; S ];

% run the analysis
[ dat, cca, cc2 ] = ccaMapFullAnalysis(deg, vars, varsQconf, netNames, varsNames, confNames, varsLabel, 38, 40, 0, 5, 10000);

% build the dissimilarity
dmat1 = ccaDissimilarityMatrix(cca);
dmat2 = ccaDissimilarityMatrix(cc2);

% create the modules
mdDat1 = fnModuleDensity(dmat1, [ yeoLabs.yeo7; ib+10 ], 'mean');
mdDat2 = fnModuleDensity(dmat2, [ yeoLabs.yeo7; ib+10 ], 'mean');

figure('Position', [ 525 400 900 700 ]);

subplot(2, 2, 1);
imagesc(dmat1); axis equal; axis square; axis tight;
colorbar; caxis([ 0 .1 ]); 
set(gca, 'XTick', [], 'YTick', []);
title({'Crossvalidated Variable Loadings', '(Has error bars)'});

subplot(2, 2, 2);
imagesc(dmat2); axis equal; axis square; axis tight;
colorbar; caxis([ 0 1 ]); 
set(gca, 'XTick', [], 'YTick', []);
title({'Loadings from Crossvalidated Factors', '(Has no error bars)'});

subplot(2, 2, 3); imagesc(mdDat1);
axis square; axis equal; axis tight; colorbar; caxis([ 0 1 ]);
title({'Dissimiliarity Between Brain and Task Domains', '(Has error bars)'});
set(gca, 'XTick', 1:size(mdDat1, 1), 'XTickLabels', mdLabs, 'XTickLabelRotation', 45, ...
    'YTick', 1:size(mdDat1, 1), 'YTickLabels', mdLabs);
set(gca, 'XLim', [ 0.5 10.5 ], 'YLim', [ 10.5 17.5 ]); caxis([ 0 0.5 ]);

subplot(2, 2, 4); imagesc(mdDat2);
axis square; axis equal; axis tight; colorbar; caxis([ 0 1 ]);
title({'Dissimiliarity Between Brain and Task Domains', '(Has no error bars)'});
set(gca, 'XTick', 1:size(mdDat2, 1), 'XTickLabels', mdLabs, 'XTickLabelRotation', 45, ...
    'YTick', 1:size(mdDat2, 1), 'YTickLabels', mdLabs);
set(gca, 'XLim', [ 0.5 10.5 ], 'YLim', [ 10.5 17.5 ]); caxis([ 0.5 1 ]);
