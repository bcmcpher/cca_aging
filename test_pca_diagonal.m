%% parameter tuning on the diagonal of the PCA space

%% load everything

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
load('yeoLabs.mat');
 
% merge age into behavior / confounds?
%varsQconf = [ age, varsQconf ];
%confNames = [ 'Age', confNames ];

%% run the parameter sweep

% preallocate values
params = 5:5:100;
out = nan(length(params), 4, 2);

% for every parameter
for ii = 1:length(params)
    
    % run the cca for the parameter vals
    [ dat, cca ] = ccaTestFullAnalysis(deg, vars, varsQconf, ...
                                       netNames, varsNames, confNames, varsLabel, ...
                                       params(ii), params(ii), 2, true, 'mean', 5, 1000);
    
    % boostrap mean / sd of the correlation with age
    [ out(ii, 1, 1), out(ii, 1, 2) ] = ccaLinRegCorr(cca, 1, age, 1000);
    
    % pull the mean / sd of the dataset correlation
    cdat = squeeze(cca.cca.hocorrs(:, 1, :));
    out(ii, 2, 1) = mean(cdat(:));
    out(ii, 2, 2) = std(cdat(:));
    
    % pull the mean / sd of the loading standard deviation
    out(ii, 3, 1) = mean([ cca.dat1.loading_sd(:, 1); cca.dat2.loading_sd(:, 1) ]);
    out(ii, 3, 2) = std([ cca.dat1.loading_sd(:, 1); cca.dat2.loading_sd(:, 1) ]);
    
    % pull the mean / sd of the loadings
    out(ii, 4, 1) = mean([ cca.dat1.loading(:, 1); cca.dat2.loading(:, 1) ]);
    out(ii, 4, 2) = std([ cca.dat1.loading(:, 1); cca.dat2.loading(:, 1) ]);
    
end

%% plot the data from this run

fh = figure; 

% corr with age
subplot(3, 1, 1); hold on;

% plot the line
plot(params, squeeze(out(:, 1, 1)));

% plot the points with errorbars
for ii = 1:size(out, 1)
    
    % pull values
    xv = params(ii);
    pt = out(ii, 1, 1);
    sd = out(ii, 1, 2);
    
    % plot points and error bars
    plot([ xv xv ], [ pt pt ], '.', 'color', 'black');
    plot([ xv xv ], [ (pt + 2*sd) (pt - 2*sd) ], 'color', 'black');
     
end

% label and format subplot
title('Correlation with Age across # of Prinicpal Components');
xlabel('N PCs Used'); 
ylabel('Correrlation with Age');
hold off;

% corr with data
subplot(3, 1, 2); hold on;

% plot the line
plot(params, squeeze(out(:, 2, 1)));

% plot the points with errorbars
for ii = 1:size(out, 1)
    
    % pull values
    xv = params(ii);
    pt = out(ii, 2, 1);
    sd = out(ii, 2, 2);
    
    % plot points and error bars
    plot([ xv xv ], [ pt pt ], '.', 'color', 'black');
    plot([ xv xv ], [ (pt + 2*sd) (pt - 2*sd) ], 'color', 'black');
     
end

% label and format subplot
title('Correlation Between Brain/Behavior across # of Prinicpal Components');
xlabel('N PCs Used'); 
ylabel('Correrlation Between Domains');
hold off;

% sd of loadings
subplot(3, 1, 3); hold on;

% plot the line
plot(params, squeeze(out(:, 3, 1)));

% plot the points with errorbars
for ii = 1:size(out, 1)
    
    % pull values
    xv = params(ii);
    pt = out(ii, 3, 1);
    sd = out(ii, 3, 2);
    
    % plot points and error bars
    plot([ xv xv ], [ pt pt ], '.', 'color', 'black');
    plot([ xv xv ], [ (pt + 2*sd) (pt - 2*sd) ], 'color', 'black');
     
end

% label and format subplot
title('Loading Variability across # of Prinicpal Components');
xlabel('N PCs Used'); 
ylabel('Loading Variability');
hold off;
