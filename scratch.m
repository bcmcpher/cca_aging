%% scrach space for testing things
%

%% ridge regression of individual data sets to get an R2 value from them

load('canoncorr_analysis_full_data.mat', 'age', 'deg', 'vars');
vars(isnan(vars)) = 0;

% use ridge regression to estimate the betas w/ l2 norm
% this does not mean center the data so estimates are on the right scale
b0_brain = ridge(age, deg, 1, 0);
b0_behavior = ridge(age, vars, 1, 0);
b0_both = ridge(age, [ deg vars ], 1, 0);

% create estimated from ridge regression
yh_brain = [ ones(size(deg, 1), 1) deg ] * b0_brain;
yh_behavior = [ ones(size(vars, 1), 1) vars ] * b0_behavior;
yh_both = [ ones(size([ deg vars ], 1), 1) [ deg vars ] ] * b0_both;

% find R2 b/c matlab ridge fxns don't estimate it
Rsq_brain = 1 - sum((age - yh_brain).^2)/sum((age - mean(age)).^2);
Rsq_behavior = 1 - sum((age - yh_behavior).^2)/sum((age - mean(age)).^2);
Rsq_both = 1 - sum((age - yh_both).^2)/sum((age - mean(age)).^2);

% Brain R2 = 0.9142
% Behavior R2 = 0.9270
% Both R2 = 0.9986

mean(abs(age - yh_brain))    % +/- 4.2 years
mean(abs(age - yh_behavior)) % +/- 3.8 years
mean(abs(age - yh_both))     % +/- 0.5 years

%% permute the cca

vars(isnan(vars)) = 0;

r = nan(Nperm, 2);
r2 = nan(Nperm, 2);
for ii = 1:Nperm
    
    % randomly split data sets in half
    rd1 = unique(randsample(1:size(deg, 2), size(deg, 2), 'true'));
    rd2 = setdiff(1:size(deg, 2), rd1);
    rb1 = unique(randsample(1:size(vars, 2), size(vars, 2), 'true'));
    rb2 = setdiff(1:size(vars, 2), rb1);
    
    % permuted canoncorr
    [ ~, ~, R1, x1, y1 ] = canoncorr(deg(:, rd1), deg(:, rd2));
    [ ~, ~, R2, x2, y2 ] = canoncorr(vars(:, rb1), vars(:, rb2));
    
    r(ii, 1) = R1(1);
    r(ii, 2) = R2(1);
    
    % build design matrix with cca axes
    z1 = [ x1, y1 ];
    z2 = [ x2, y2 ];
    
    % fit simple linear model of age ~ cca factor #
    out1 = fitlm(z1, age);
    out2 = fitlm(z2, age);
    
    % return multiway correlation
    r2(ii, 1) = sqrt(out1.Rsquared.Ordinary);
    r2(ii, 2) = sqrt(out2.Rsquared.Ordinary);
    
end

clear ii rd1 rd2 rb1 rb2 x1 x2 y1 y2 z1 z2 out1 out2

cmap = parula(88);

figure; hold on;
for ii = 1:size(x1, 1)
    plot(x1(ii, 1), y1(ii, 1), '.', 'color', cmap(age(ii), :));
end
hold off;

figure; hold on;
for ii = 1:size(x2, 1)
    plot(x2(ii, 1), y2(ii, 1), '.', 'color', cmap(age(ii), :));
end
hold off;

%% null slope test

pdat = nan(Nperm, 1);
for ii = 1:Nperm
    tdat = randsample(beh_lcoeff, 334, 'true');
    pdat(ii) = mean(tdat);
end

sum(beh_lcoeff > mean(pdat))

ub = mean(pdat) + 2*(std(pdat));
lb = mean(pdat) - 2*(std(pdat));

prcntile(pdat, [ 5 95 ]);

size(beh_lcoeff < lb, 1)
mean(beh_lcoeff(beh_lcoeff < lb))
std(beh_lcoeff(beh_lcoeff < lb))

sum(beh_lcoeff > ub)
mean(beh_lcoeff(beh_lcoeff > ub))
std(beh_lcoeff(beh_lcoeff > ub))

size((brn_lcoeff > lb) & (brn_lcoeff < ub), 1)
mean(brn_lcoeff(((brn_lcoeff > lb) & (brn_lcoeff < ub))))
std(brn_lcoeff(((brn_lcoeff > lb) & (brn_lcoeff < ub))))

%% plot cv factors / error bars

% pull loadings
x1 = cca.dat1.loading(:, 1);
y1 = cca.dat1.loading_sd(:, 1);
x2 = cca.dat2.loading(:, 1);
y2 = cca.dat2.loading_sd(:, 1);

% plot
figure;
subplot(2, 1, 1);
for ii = 1:size(x1, 1)
    plot(1:size(x1, 1), x1(1), '.');
end
subplot(2, 1, 2);
for ii = 1:size(x1, 1)
    plot(1:size(x2, 1), x2(1), '.');
end

%% plot param sweep plots

% plot everything stored
figure;
for ii = 1:8
    
    subplot(2, 4, ii); 
    imagesc(dat(:,:,ii)); 
    axis equal; axis square; axis tight;
    
end

% just pull age
cage = dat(:,:,1);
sage = dat(:,:,2);

figure; 
% subplot(1, 2, 1);
imagesc(cage); colorbar;
axis equal; axis square; axis tight;
set(gca, 'XTick', [ 2, 50, 100 ], 'YTick', [ 2, 50, 100 ]);
% subplot(1, 2, 2);
% imagesc(sage); colorbar;
% axis equal; axis square; axis tight;
% set(gca, 'XTick', [ 2, 50, 100 ], 'YTick', [ 2, 50, 100 ]);

% max corr w/ age?
[ x, y ] = find(cage == max(cage(:))); % not what I ran. Why?

%% load an individual param from sweep and make basic plots

% pick brain / behavior number of PCs
npca = [ 38, 40; 50, 50; 100, 100 ];
datadir = '/N/dc2/projects/lifebid/HCP/Brent/camcan/figs/param_search_fig';

for ii = 2:size(npca, 1)
    
    % create the stem for the output
    stem = fullfile(datadir, sprintf('pca_brain_%03d_behavior_%03d', npca(ii, 1), npca(ii, 2)));
    
    % run the cca for the parameter vals
    [ dat, cca, cc2 ] = ccaMapFullAnalysis(deg, vars, varsQconf, ...
                                           netNames, varsNames, confNames, varsLabel, ...
                                           npca(ii, 1), npca(ii, 2), 0, 5, 250000);
    
    % plot the basics of them
    ccaPlotAxisCon(cca, 1, age, parula(88));
    print([ stem '_main.eps' ], '-painters', '-depsc'); close all;
    
    ccaPlotRankedTrends(dat, cca, age, 'brain', 'load', 1, 'lines', 10);
    print([ stem '_brain.eps' ], '-painters', '-depsc'); close all;
    
    ccaPlotRankedTrends(dat, cca, age, 'behavior', 'load', 1, 'lines', 10);
    print([ stem '_behavior.eps' ], '-painters', '-depsc'); close all;
    
end

%% test reproducibility of fixed cv repeat loop

% get 2 runs
[ dat, cca1 ] = ccaFullKAnalysis(deg, vars, varsQconf, netNames, varsNames, confNames, varsLabels, 38, 40, 0, 5, 2500);
[ ~, cca2 ] = ccaFullKAnalysis(deg, vars, varsQconf, netNames, varsNames, confNames, varsLabels, 38, 40, 0, 5, 2500);

% compare main plots
ccaPlotAxisCon(cca1, 1, age, parula(88), true);
ccaPlotAxisCon(cca2, 1, age, parula(88), true);

% build dissimilarity
dmat1 = ccaDissimilarityMatrix(cca1);
dmat2 = ccaDissimilarityMatrix(cca2);

% plot the modules
mdDat1 = fnModuleDensity(dmat1, [ yeoLabs.yeo7; ib+10 ], 'mean');
figure; imagesc(mdDat1);
axis square; axis equal; axis tight; colorbar; caxis([ 0 1 ]);
title('Dissimiliarity Between Brain and Task Domains');
mdLabs = [ yeoLabs.yeo7Names'; S ];
set(gca, 'XTick', 1:size(mdDat1, 1), 'XTickLabels', mdLabs, 'XTickLabelRotation', 45, ...
    'YTick', 1:size(mdDat1, 1), 'YTickLabels', mdLabs);

mdDat2 = fnModuleDensity(dmat2, [ yeoLabs.yeo7; ib+10 ], 'mean');
figure; imagesc(mdDat2);
axis square; axis equal; axis tight; colorbar; caxis([ 0 1 ]);
title('Dissimiliarity Between Brain and Task Domains');
mdLabs = [ yeoLabs.yeo7Names'; S ];
set(gca, 'XTick', 1:size(mdDat2, 1), 'XTickLabels', mdLabs, 'XTickLabelRotation', 45, ...
    'YTick', 1:size(mdDat2, 1), 'YTickLabels', mdLabs);
