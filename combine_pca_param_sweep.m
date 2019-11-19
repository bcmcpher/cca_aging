%% compile together the pca parameter sweep

% pick the folder to load param sweep data from
%vers = 'local';
vers = 'osg';

datadir = [ '/N/dc2/projects/lifebid/HCP/Brent/camcan/param_sweep/' vers ];

% if requested output exists 
if exist([ '/N/dc2/projects/lifebid/HCP/Brent/camcan/param_sweep/combined_' vers '.mat' ], 'file')
    % load it into the workspace
    load([ '/N/dc2/projects/lifebid/HCP/Brent/camcan/param_sweep/combined_' vers '.mat' ]);
else
    % preallocate the data output
    dat = nan(100, 100, 8);
end

% for every pca pair
for ii = 2:100
    for jj = 2:100
        
        % build the file
        file = fullfile(datadir, sprintf('pca_brain_%03d_behavior_%03d_test.mat', ii, jj));
        
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
            dat(ii, jj, 5) = out.out(3, 1);
            dat(ii, jj, 6) = out.out(3, 2);
            dat(ii, jj, 7) = out.out(4, 1);
            dat(ii, jj, 8) = out.out(4, 2);
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
