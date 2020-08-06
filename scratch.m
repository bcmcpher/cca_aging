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

sqrt(Rsq_brain)
sqrt(Rsq_behavior)
sqrt(Rsq_both)

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

%% plot the parameter search space

z = load('param_sweep/combined_osg-mean-perm.mat');

figure;
for ii = 1:4
    subplot(2, 2, ii);
    imagesc(z.dat(:,:,ii)); colorbar;
    axis equal; axis square; axis tight;
end



%% error estimate by axis

% % R code
%
% # standardize data
% > Xs <- scale(X)
% > Ys <- scale(Y)
%
% # canonical correlations of correlation (standardized data)
% > ccas <- cancor(Xs, Ys)
%
% # cca (the normal way)
% > Sx <- cov(Xs)
% > Sy <- cov(Ys)
% > Sxy <- cov(Xs,Ys)
%
% # error of approximation matrices (with r=2)
% > Ainv <- solve(Ahat)
% > Binv <- solve(Bhat)%
% > r <- 2
%
% > Ex <- Sx - crossprod(Ainv[1:r,])
% > Ey <- Sy - crossprod(Binv[1:r,])
% > Exy <- Sxy - crossprod(diag(rho[1:r]) %*% Ainv[1:r,], Binv[1:r,])
%
% # get norms of error matrices
% > sqrt(mean(Ex^2))[1] 0.2432351 
% > sqrt(mean(Ey^2))[1] 0.2296716
% > sqrt(mean(Exy^2))[1] 0.07458264
% 

% run from debugger
% must happen after PCA - need same # of columns
Sx = cov(uu1);
Sy = cov(uu2(:, 1:38));
Sxy = cov(uu1, uu2(:, 1:38));

Ainv = inv(A);
Binv = inv(B(1:38, 1:38));

Ex = Sx - (Ainv * Ainv');
Ey = Sy - (Binv * Binv');
Exy = Sxy - (Ainv * Binv);

%% the count of variables

size(find(contains(dat.dat2.names', 'RT')), 1) + ...
size(find(contains(dat.dat2.names', 'rt')), 1) + ...
size(find(contains(dat.dat2.names', 'reITRT')), 1) + ...
size(find(contains(dat.dat2.names', 'invRT')), 1) + ...
size(find(contains(dat.dat2.names', 'ReactionTime')), 1) + ...
size(find(contains(dat.dat2.names', 'MovementTime')), 1)

% the accuracy measures
size(find(contains(dat.dat2.names', 'ACC')), 1) + ...
size(find(contains(dat.dat2.names', 'pValid')), 1) + ...
size(find(contains(dat.dat2.names', 'pYes')), 1) + ...
size(find(contains(dat.dat2.names', 'FAM')), 1) + ...
size(find(contains(dat.dat2.names', 'fam')), 1) + ...
size(find(contains(dat.dat2.names', '_mean')), 1) + ...
size(find(contains(dat.dat2.names', 'emm')), 1) + ...
size(find(contains(dat.dat2.names', 'ekm')), 1) 

%% plot between data and loading

% data modules
zdat

% cca modules
mdDat

% unique indices + diagonal
xyi = nchoosek(1:17, 2);
xyi = [ xyi; [ 1:17; 1:17 ]' ];

% create figure
fh = figure;

% make 3 panels
for sp = 1:3
    
    subplot(1, 3, sp); hold on;
        
    % for every point
    for ii = 1:size(xyi, 1)
        
        % pull the points
        ptx = xyi(ii, 1);
        pty = xyi(ii, 2);
        
        % pull the values
        cval = 1-mcDat(ptx, pty);
        dval = 1-mdDat(ptx, pty);
        
        % plot different colors if x/y are within/between brain/behavior
        
        % check x index
        if ptx > 10
            xc = 'behavior';
        else
            xc = 'brain';
        end
        
        % check y index
        if pty > 10
            yc = 'behavior';
        else
            yc = 'brain';
        end
        
        if (strcmp(xc, 'brain') && strcmp(yc, 'brain')) && sp == 2
            ttl = 'Brain-Brain Interactions';
            xlim = [ 0 1 ]; ylim = [ 0 1 ];
            plot(cval, dval, 'o', 'MarkerFaceColor', [ 0.1 0.3 0.7 ], ...
                'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'MarkerSize', 6); % plot within brain
        end
        
        if (strcmp(xc, 'behavior') && strcmp(yc, 'behavior')) && sp == 3
            ttl = 'Behavior-Behavior Interactions';
            xlim = [ 0 1 ]; ylim = [ 0 1 ];
            plot(cval, dval, 'o', 'MarkerFaceColor', [ 0.3 0.7 0.1 ], ...
                'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'MarkerSize', 6); % plot within behavior
        end
        
        if (strcmp(xc, 'brain') && strcmp(yc, 'behavior')) && sp == 1
            ttl = 'Brain-Behavior Interactions';
            xlim = [ 0 1 ]; ylim = [ 0 1 ];
            plot(cval, dval, 'o', 'MarkerFaceColor', [ 0.7 0.3 0.1 ], ...
                'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'MarkerSize', 6); % plot brain x behavior
        end
        
    end
    
    hold off;
    title({ttl, 'Comparison Between Data and CCA'});
    xlabel('Data Similarity');
    ylabel('CCA Similarity');
    axis square; axis equal; axis tight
    set(gca, 'XLim', xlim, 'YLim', ylim);
    
end

%% display yeo labels as a convenient table

% 
[ {yeoLabs.yeo7Names{yeoLabs.yeo7}}' ];

%% plot rank order of modules from chord plot w/ readable names / numbers



%% get new parameters of trends across variables

% pull normalized behavior data
beh = dat.dat1.raw; % dat1 brain, dat2 behavior
beh_names = dat.dat1.names;

% or pull network statistics
z = load('canoncorr_analysis_full_workspace.mat', 'node', 'glob');
for ii = 1:size(z.glob, 1)
    beh(ii, 1) = z.glob{ii}.density;
    beh(ii, 2) = z.glob{ii}.efficiency;
    beh(ii, 3) = max(z.node{ii}.degree);
end

% preallocate output
out = nan(size(beh, 2), 4);

% for every variable
for ii = 1:size(beh, 2)
    
    % create the trend
    tmp = plotVarByDeciles(beh(:, ii), age, 'quadratic', beh_names{ii});
    
    % pull the values
    out(ii, 1) = tmp.lcoeff;
    out(ii, 2) = tmp.R2;
    out(ii, 3) = tmp.AIC;
    out(ii, 4) = tmp.AICc;
    
    % close the figure
    close all

end

% something is bad b/c complex, but just drop it
out = real(out);

% find positive / negative trends
sp = out(:, 1) > 0;

sum(sp)

% positive parameters
mean(out(sp, :))
std(out(sp, :), [], 1)

% negative parameters
mean(out(~sp, :))
std(out(~sp, :), [], 1)

%% Fig2 plot showing corrs w/ and w/o age for full and xval models

m1 = load('sub_pca_brain_38_behavior_40_15k_10k_boot.mat');
m2 = load('sub_pca_brain_38_behavior_40_15k_10k_boot_m2.mat');

figure('Position', [ 680 325 200 600]); hold on;

plot(1, m1.cca.full.grotR(1), 'o', 'MarkerFaceColor', 'red', ...
     'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'MarkerSize', 10);

plot([ 2, 2 ], [ 3*m1.cca.cca.hocorrs(1)+m1.cca.cca.hocorrs_se(1) 3*m1.cca.cca.hocorrs(1)-m1.cca.cca.hocorrs_se(1) ], 'color', 'black');
plot(2, m2.cca.full.grotR(1), 'o', 'MarkerFaceColor', 'blue', ...
     'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'MarkerSize', 10);

plot(3, m1.cca.cca.hocorrs(1), 'o', 'MarkerFaceColor', 'red', ...
     'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'MarkerSize', 10);

plot([ 4, 4 ], [ 3*m2.cca.cca.hocorrs(1)+m2.cca.cca.hocorrs_se(1) 3*m2.cca.cca.hocorrs(1)-m2.cca.cca.hocorrs_se(1) ], 'color', 'black');
plot(4, m2.cca.cca.hocorrs(1), 'o', 'MarkerFaceColor', 'blue', ...
     'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'MarkerSize', 10);

set(gca, 'Xlim', [ 0.5 4.5 ], 'Ylim', [ 0 1 ]);
print('figs/panels_20200728/m1_vs_m2_ca1_corr.eps', '-painters', '-depsc');

%% Fig2 plot showing corrs w/ and w/o age for full and xval models

% load the models
m1 = load('sub_pca_brain_38_behavior_40_15k_10k_boot.mat');
m2 = load('sub_pca_brain_38_behavior_40_15k_10k_boot_m2.mat');

% compute the corr +/- with age
[ p1, se1 ] = ccaLinRegCorr(m1.cca, 1, age, 1000);
[ p2, se2 ] = ccaLinRegCorr(m2.cca, 1, age, 1000, true);
[ p3, se3 ] = ccaLinRegCorr(m1.cca, 1, age, 1000);
[ p4, se4 ] = ccaLinRegCorr(m2.cca, 1, age, 1000, true);

figure('Position', [ 680 325 200 600]); hold on;

plot([ 1, 1 ], [ p1+(3*se1) p1-(3*se1) ], 'color', 'black');
plot(1, p1, 'o', 'MarkerFaceColor', 'red', ...
     'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'MarkerSize', 10);

plot([ 2, 2 ], [ p2+(3*se2) p2-(3*se2) ], 'color', 'black');
plot(2, p2, 'o', 'MarkerFaceColor', 'blue', ...
     'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'MarkerSize', 10);

plot([ 3, 3 ], [ p3+(3*se3) p3-(3*se3) ], 'color', 'black');
plot(3, p3, 'o', 'MarkerFaceColor', 'red', ...
     'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'MarkerSize', 10);

plot([ 4, 4 ], [ p4+(3*se4) p4-(3*se4) ], 'color', 'black');
plot(4, p4, 'o', 'MarkerFaceColor', 'blue', ...
     'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'MarkerSize', 10);

set(gca, 'Xlim', [ 0.5 4.5 ], 'Ylim', [ 0 1 ]);
print('figs/panels_20200728/m1_vs_m2_ca1_age.eps', '-painters', '-depsc');

%% plot null bootstrap around corr

% preallocate null shape
lbub = nan(size(cca.cca.hocorrs_null, 2), 2);

% pull the 5th / 95th percentiles
for ii = 1:size(cca.cca.hocorrs_null, 2)
    lbub(ii,:) = prctile(squeeze(cca.cca.hocorrs_null(1,ii,:)), [ 5 95 ]);
end

figure; hold on

% plot null background
fill([ 1:38, fliplr(1:38) ], [ lbub(:,1)', lbub(:,2)' ], [ .75 .75 .75 ]);

% for every point
for ii = 1:size(cca.cca.hocorrs, 2)
    
    % plot error bars
    plot([ ii ii ], ...
         [ (cca.cca.hocorrs(ii) + 2*cca.cca.hocorrs_se(ii)) ...
         (cca.cca.hocorrs(ii) - 2*cca.cca.hocorrs_se(ii)) ], 'k');
    
    % plot points
    plot(ii, cca.cca.hocorrs(ii), 'o', 'MarkerEdgeColor', [ 0 0 0 ], ...
        'MarkerFaceColor', [ 0 0 .8 ], 'MarkerSize', 8);
    
end

hold off

%% mds color plot - not used

% pull the raw data
zzz = [ dat.dat1.nrm(:, 1:376) dat.dat2.nrm ];
zzz(isnan(zzz)) = 0;

% pull the dissimilarity of the data
zz0 = squareform(pdist(zzz));

% run mds
[ zz1, zz2, zz3, zz4, zz5 ] = ccaMDScale(zz0, 10);

idx = [ yeoLabs.yeo7; ib+10 ];

fh = figure; hold on;

% by variable
for ii = 1:size(idx, 1)
   
    val = idx(ii);
    
    switch val
        case {1,2,3,4,5,6,7,8,9,10}
            clr = 'red';
        case {11,12,13,14,15,16,17}
            clr = 'blue';
        otherwise
            clr = 'green';
    end
            
    plot(zz1(ii, 1), zz1(ii, 2), 'o', 'MarkerFaceColor', clr, ...
         'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'MarkerSize', 6);
    
    %plot3(zz1(ii, 1), zz1(ii, 2), zz1(ii, 3), '.', 'color', clr);
            
end

hold off

cmap = parula(88);

% by subject
for ii = 1:size(zz1, 1)
   
    clr = cmap(age(ii), :);
            
    plot(zz1(ii, 1), zz1(ii, 2), 'o', 'MarkerFaceColor', clr, ...
         'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'MarkerSize', 6);
    
    %plot3(zz1(ii, 1), zz1(ii, 2), zz1(ii, 3), '.', 'color', clr);
            
end

hold off
