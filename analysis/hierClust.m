function [ T, coph, copd, inco, fh ] = hierClust(data, cval, dist, lmeth, depth)
%[ T, coph, copd, inco, fh ] = hierClust(data, cval, dist, lmeth, depth);
%   Estimate and plot a useful hierarchical clustering based on MATLAB built in matrix operations.
%
%   THIS WAS NOT USED IN THE PAPER. IT IS INCOMPLETE.
%

% NOT WORKING, BUT USEFUL CODE FOR HIERARCHICAL CLUSTERS

% data is loadings or factors (?)

% dist: 'euclidean', 'squaredeuclidean', 'seuclidean', 'cityblock',
% 'minkowski', 'chebychev', 'mahalanobis', 'cosine', 'correlation',
% 'spearman', 'hamming', 'jaccard'

% are all the same, but on differenc scales:
% 'euclidean', 'squaredeuclidean', 'seuclidean', 'cityblock', 'minkowski', 'chebychev'
% OR: 'cosine', 'correlation', 'spearman'
% don't work / are all the same:
% 'mahalanobis', hamming', 'jaccard'

% lmeth: 'single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward'
% are all slightly different
% warning, don't use: 'centroid', 'median'

% matlab hierarchical clustering based on distance between all variables
% split into brain / behavior to ID clusters separately...
% https://www.mathworks.com/help/stats/hierarchical-clustering.html

% parse defualt arguments
if(~exist('cval', 'var') || isempty(cval))
    cval = .8;
end

if(~exist('dist', 'var') || isempty(dist))
    dist = 'euclidean';
end

if(~exist('lmeth', 'var') || isempty(lmeth))
    lmeth = 'single';
end

if(~exist('depth', 'var') || isempty(depth))
    depth = 2;
end

% determine size of data
nvar = size(data, 1);

% distance between each data point based on subject vectors
diss = pdist(data, dist); % the same as below, but w/ different scale

% estimate the linkage b/w 
link = linkage(diss, lmeth);

% compute the inconsistency - distance b/w cluster megerd at a node
inco = inconsistent(link, depth);
% inco(i,1) = mean(link(Si, 3)); % the mean height of nodes in S_i
% inco(i,2) = std(link(Si, 3)); % the standard deviation of node heights in S_i
% inco(i,3) = length(Si); % the number of nodes in S_i
% inco(i,4) = (link(i, 3) - inco(i, 1)) / inco(i, 2); % the inconsistent value
       
% correlation b/w observed distances and estimated link strength
[ coph, copd ] = cophenet(link, diss);
disp([ 'Cophenetic Correlation Coefficient: ' num2str(coph) ]);

% cluster the data
ci = cluster(link, 'Cutoff', cval, 'Depth', depth);
%ci = cluster(link, 'Cutoff', cval, 'Criterion', 'distance');
%ci = cluster(link, 'MaxClust', cval);
% multiple arguments to tweak

% grab the total number of clusters found to print
nclust = size(unique(ci), 1);
disp([ 'Total clusters found: ' num2str(nclust) ]);

% try and detect the number of clusters? if set to nvar, will find nvar
E = evalclusters(data, 'linkage', 'Silhouette', 'Distance', dist, 'klist', 1:nvar);
disp([ 'Silhouette Plot Ideal N Clusters: ' num2str(E.OptimalK) ]);

% create sorted indices in original space
[ sf, T ] = sort(ci);

% create indices to show communities
cidx = nchoosek(1:nvar, 2);

% determine communities of edges share nodes
cmat = zeros(nvar);

% for every edge
for ii = 1:size(cidx, 1)
    % if the nodes are shared, fill in the label
    if sf(cidx(ii, 1)) == sf(cidx(ii, 2))
        cmat(cidx(ii, 1), cidx(ii, 2)) = sf(cidx(ii, 1));
        cmat(cidx(ii, 2), cidx(ii, 1)) = sf(cidx(ii, 1));
    % otherwise they're zero
    else
        cmat(cidx(ii, 1), cidx(ii, 2)) = 0;
        cmat(cidx(ii, 2), cidx(ii, 1)) = 0;
    end
end

% create the matrix
pmat = squareform(diss);

%% plot the data

% plot all the hierarchical clustering
fh = figure('Position', [ 175 75 1300 800 ]); 

haxis(1) = subplot(3, 3, 2);
imagesc(pmat(T, T)); colorbar;
title('Sorted Distance Matrix');
axis square; axis equal; axis tight
set(gca, 'XTick', [], 'YTick', []);

haxis(2) = subplot(3, 3, 3);
imagesc(cmat); colorbar; 
title([ 'Communities - N = ' num2str(nclust) ]);
axis square; axis equal; axis tight
set(gca, 'XTick', [], 'YTick', []);

haxis(3) = subplot(3, 3, 5:6);
dendrogram(link, nclust, 'ColorThreshold', 'default'); % nclust was nvar
title('Dendrogram Clusters of All Variables');
set(gca, 'XTickLabelRotation', 90);

haxis(4) = subplot(3, 3, 8:9);
plot(E); set(gca, 'XLim', [ 0 nvar+1 ]);
title([ 'Silhouette Plot Ideal N Clusters: ' num2str(E.OptimalK) ]);

% show how clusters want to group?
haxis(5) = subplot(3, 3, [ 1; 4; 7 ]);
silhouette(data, ci, dist);
title('Silhouette Cluster of Nodes');

%% set all axes to be larger

% https://www.mathworks.com/matlabcentral/answers/233818-how-to-create-subplots-with-little-vertical-spacing

%pos1 = get(haxis(1), 'Position');
pos1 = [ .35 .67 .30 .30 ];
set(haxis(1), 'Position', pos1);

%pos2 = get(haxis(2), 'Position');
pos2 = [ .67 .67 .30 .30 ];
set(haxis(2), 'Position', pos2);

%pos3 = get(haxis(3), 'Position');
pos3 = [ .40 .375 .59 .25 ];
set(haxis(3), 'Position', pos3);

%pos4 = get(haxis(4), 'Position');
pos4 = [ .40 .05 .59 .25 ];
set(haxis(4), 'Position', pos4);

%pos5 = get(haxis(5), 'Position');
pos5 = [ .02 .05 .33 .90 ];
set(haxis(5), 'Position', pos5);

end
