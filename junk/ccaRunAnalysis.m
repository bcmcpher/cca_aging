function [ grotAv, grotAAd, grotAAv, grotBv, grotBBd, grotBBv, grotU, grotV, grotRp, grotRv, grotPc, grotstats, grotRpval, Ncca1, Ncca2 ] = ccaRunAnalysis(uu1, uu2, NETd, NET, varsgrot, Nkeep, Nperm)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%
% uu1: nets
% uu2: behavior
%
% References:
% http://www.statsoft.com/Textbook/Canonical-Analysis
% https://stats.idre.ucla.edu/spss/output/canonical-correlation-analysis/
%

%% CCA

disp('Performing CCA...');

[ A, B, grotR, grotU, grotV, grotstats ] = canoncorr(uu1, uu2);

%% extract parameters back to raw measures

disp('Converting CCA loadings back to raw network / behavior measures...');

% preallocate 
grotAAd = nan(size(NET, 2), Nkeep);
grotBBd = nan(size(varsgrot, 2), Nkeep);
grotAAv = nan(Nkeep, 1);
grotBBv = nan(Nkeep, 1);
grotAAr = nan(Nkeep, 1);
grotBBr = nan(Nkeep, 1);

for ii = 1:Nkeep
    
    % network edge weights for CCA mode 1
    %grotAA = corr(grotU(:, ccf), NET)';
    
    % weights after deconfounding
    grotAAd(:, ii) = corr(grotU(:, ii), NETd(:, 1:size(NET, 2)))';
    
    % demo weights for CCA mode 1
    %grotBB = corr(grotV(:, ccf), vars, 'rows', 'pairwise');
    
    % weights after deconfounding
    grotBBd(:, ii) = corr(grotV(:, ii), varsgrot, 'rows', 'pairwise')';
    
    % proportion of variance explained in each cc
    %grotRv(ii) = grotR(ii)^2;
    
    % variance extracted by each cc data set
    grotAAv(ii) = mean(grotAAd(:, ii).^2, 'omitnan');
    grotBBv(ii) = mean(grotBBd(:, ii).^2, 'omitnan');
    
    % redundancy of each cc component
    grotAAr(ii) = mean(grotAAd(:, ii).^2, 'omitnan') * grotR(ii)^2;
    grotBBr(ii) = mean(grotBBd(:, ii).^2, 'omitnan') * grotR(ii)^2;
    
end

clear ii

% proportion of variance accounted for in each CC
grotRv = grotR .^ 2;

% percent accounted for in each CC
grotPc = (grotRv ./ sum(grotRv)) * 100;

% grotR.^2 : eigenvalues of CCA, the proportion of variance accounted for (NOT EXPLAINED)
% [ mean(grotAAd.^2); mean(grotBBd.^2) ] is the "variance explained" in each data set
% redundancy is mean(grotXX .^ 2) * grotR(X)^2

%% create mapping from variables in uu1 to uu2

% create scaled data in original scale x subject space 

% A/B is:             pca components x cca coeffs
% grotU/grotV is:     subj x cca coeffs
% grotAAd/grotBBd is: raw input x cca coeffs
%
% grot??d? * A' scales the observation values to the CCA
%     x? is a raw inputs scaled to each cca component
% x? * grot??d transforms the cca scaled measure/subj values into a measure x subj contribution in cca space?
%     z? is a scaled contribution for every input from each subject

x1 = grotAAd * A'; % scale original measures (nodes) to CC space
z1 = x1 * grotU'; % scaled contribution of nodes by subject
%z1 = grotAAd * grotU'; % doesn't scale vars to CCA space, right?

x2 = grotBBd * B'; % scale original measures (behaviors) to CC space
z2 = x2 * grotV'; % scaled contribution of behaviors by subject
%z2 = grotBBd * grotV';

% combine into one big [variable x subject] matrix
z = [ z1; z2 ];
% this is a variable x subj matrix that has the scaled constribution of
% each field stored in it. The goal is to see how the vars correlate w/
% each other using the unique subject contributions
nval = size(z, 1);

% for every unique value in the merged matrix
xy = nchoosek(1:size(z, 1), 2);
omat = zeros(size(z, 1));

% find the correlation b/w subjects scaled values b/w each other variable
% in CCA from both domains. 
for ii = 1:size(xy, 1) % 1 - corr() for a dissimilarity matrix
    %omat(xy(ii, 1), xy(ii, 2)) = corr(z(xy(ii, 1), :)', z(xy(ii, 2), :)');
    %omat(xy(ii, 2), xy(ii, 1)) = corr(z(xy(ii, 1), :)', z(xy(ii, 2), :)');
    omat(xy(ii, 1), xy(ii, 2)) = 1 - abs(corr(z(xy(ii, 1), :)', z(xy(ii, 2), :)'));
    omat(xy(ii, 2), xy(ii, 1)) = 1 - abs(corr(z(xy(ii, 1), :)', z(xy(ii, 2), :)'));
end

% hierarchical clusters
%mc_brn = hierClust(z1);
%mc_beh = hierClust(z2);

% simple louvain detection
%lc_brn = community_louvain(omat(1:376, 1:376), 1, [], 'negative_sym');
%lc_beh = community_louvain(omat(377:end, 377:end), 1, [], 'negative_sym');

% get index of edges into communities
%[ ~, lc_brn ] = sort(lc_brn);
%[ ~, lc_beh ] = sort(lc_beh);

%load('~/hcp_mmp_vol/working/yeoLabs.mat');
[ y7_lab, y7_brn ] = sort(yeoLabs.yeo7);
%[ ~, y17_brn ] = sort(yeoLabs.yeo17);

% get index of sorted behaviors
%load('camcan_594_deg_cca.mat', 'varsNames');
svar = regexprep(varsNames', '_.*', '');
svar = regexprep(svar, 'hint', 'comp');
svar{1} = 'comp'; % replace the first dumb label
[ S, ~, ib ] = unique(svar);
[ si, sv ] = sort(ib);
%[ bh_lab, sd_beh ] = sort(svar); 

% merge axes
%lc_axis = [ lc_brn; (lc_beh + 376) ];
%mc_axis = [ mc_brn; (mc_beh + 376) ];

y7_axis = [ y7_brn; sv+376 ];
%y17_axis = [ y17_brn; sd_beh+376 ];

% lines splitting domains
dlin = [ find(diff(y7_lab)); find(diff(si))+376 ];
dlin = sort([ dlin; 376 ]);

% visualize all - plotAdjacencyMatrix zeros out negatives
figure; imagesc(omat(y7_axis, y7_axis)); hold on;
axis square; axis equal; axis tight; colorbar;
set(gca, 'XTick', [], 'YTick', []);
%set(gca, 'XTick', 1:710, 'YTick', 1:710);
%set(gca, 'XTickLabel', [ netNames; varsNames' ], 'YTickLabel', [ netNames; varsNames' ]);
caxis([ 0 1 ]); %colormap(redbluecmap);
title('Dissimilarity between all variables in CCA');
for ii = 1:size(dlin, 1)
    line([ dlin(ii)+.5 dlin(ii)+.5 ], [ .5 nval+.5 ], 'color', 'black', 'LineWidth', .25);
    line([ .5 nval+.5 ], [ dlin(ii)+.5 dlin(ii)+.5 ], 'color', 'black', 'LineWidth', .25);
end
hold off;

% rescale for brain only
set(gca, 'XLim', [ 1 376 ], 'YLim', [ 1 376 ]);

% rescale for behavior only
set(gca, 'XLim', [ 376 710 ], 'YLim', [ 376 710 ]);

% rescale for brain x behavior
set(gca, 'XLim', [ 1 376 ], 'YLim', [ 376 710 ]);

% % visualize brain x behavior only
% int_mat = omat(sv+376, y7_brn);
% figure; imagesc(int_mat); axis square; axis equal; axis tight; 
% set(gca, 'XTick', [], 'YTick', []); colorbar;
% caxis([ 0 1 ]); %colormap(redbluecmap);
% title('Dissimilarity between brain regions and tasks');
% xlabel('Brain'); ylabel('Behavior');

% plot the averaged together values
[ mdDat, ~, mdBtw ] = fnModuleDensity(omat, [ yeoLabs.yeo7; ib+10 ], 'mean');
figure; 
imagesc(mdDat);
%imagesc(mdDat - mdBtw); mean center b/w weights?
axis square; axis equal; axis tight; colorbar;
caxis([ 0 1 ]); 
title('Dissimiliarity Between Brain and Task Domains');
mdLabs = [ yeoLabs.yeo7Names'; S ];
set(gca, 'XTick', 1:size(mdDat, 1), 'XTickLabels', mdLabs, 'XTickLabelRotation', 45, ...
    'YTick', 1:size(mdDat, 1), 'YTickLabels', mdLabs);

% plot value for each module
for ii = 1:size(mdDat, 1)
    for jj = 1:size(mdDat, 2)
        text(ii-.25, jj, num2str(round(mdDat(ii, jj), 2)));
    end
end

% subset and scale to off diagonal
set(gca, 'XLim', [ 0.5 10.5 ], 'YLim', [ 10.5 17.5 ]); caxis([ 0.85 1 ]);

%% CCA permutation testing for parameters

disp('Performing a permutation test to determine significant number of CCs...');

% preallocate output
grotRp = zeros(Nperm, Nkeep, 2); 
grotRpval = nan(Nkeep, 2);

% preallocate null noise data
grotAv = zeros(Nperm, Nkeep);
grotBv = zeros(Nperm, Nkeep);
grotAr = zeros(Nperm, Nkeep);
grotBr = zeros(Nperm, Nkeep);

% predefine permutation sets
PAPset = palm_quickperms([], ones(size(varsgrot, 1), 1), Nperm);

% for every permutation
for ii = 1:Nperm
    
    % run cca w/ scrambled data set
    %[ ~, ~, grotRp(ii, :), grotUr, grotVr ] = canoncorr(uu1, uu2(PAPset(:, ii), :));
    [ ~, ~, grotRp(ii, :, 1), grotUr, ~ ] = canoncorr(uu1, uu2(PAPset(:, ii), :));
    [ ~, ~, grotRp(ii, :, 2), ~, grotVr ] = canoncorr(uu1(PAPset(:, ii), :), uu2);
    
    % for every cc
    for jj = 1:Nkeep
        
        % edge weights after deconfounding
        grotAtr = corr(grotUr(:, jj), NETd(:, 1:size(NET, 2)))';
        
        % behavior weights after deconfounding
        grotBtr = corr(grotVr(:, jj), varsgrot, 'rows', 'pairwise')';
        
        % store variance captured in A / B
        grotAv(ii, jj) = mean(grotAtr .^ 2, 'omitnan');
        grotBv(ii, jj) = mean(grotBtr .^ 2, 'omitnan');
        
        % store redundancy in A / B
        grotAr(ii, jj) = mean(grotAtr .^ 2, 'omitnan') * grotRp(ii, jj);
        grotBr(ii, jj) = mean(grotBtr .^ 2, 'omitnan') * grotRp(ii, jj);
        
    end
end

clear ii jj grotAtr grotBtr grotUr grotVr

for ii = 1:Nkeep  % get FWE-corrected pvalues
    %grotRpval(ii) = (1 + sum(grotRp(2:end, 1) >= grotR(ii))) / Nperm;
    grotRpval(ii, 1) = (1 + sum(grotRp(2:end, 1) >= grotR(ii))) / Nperm;
    grotRpval(ii, 2) = (1 + sum(grotRp(2:end, 2) >= grotR(ii))) / Nperm;
end

clear ii 

% count # of CCA factors that pass correction
%Ncca = sum(grotRpval < 0.05); % number of FWE-significant CCA components
Ncca1 = sum(grotRpval(:, 1) < 0.05); % number of FWE-significant CCA components
Ncca2 = sum(grotRpval(:, 2) < 0.05); % number of FWE-significant CCA components

end
