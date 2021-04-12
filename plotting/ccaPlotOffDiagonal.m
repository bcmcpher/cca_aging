function [ fh ] = ccaPlotOffDiagonal(cca)
%[ fh ] = ccaPlotOffDiagonal(dat1, dat2, offd);
%   Plot a hierarchical cluster on each dataset and compare the
%   relationship between them within a CCA.
%
%   Can this be built into the CCA structure?
%
%   INPUTS:
%       cca - output from the CCA for plotting the data
%       %   dat1 - dissimilarity / distance between vars in first data set
%       %   dat2 - dissimilarity / distance between vars in second data set
%       %   offd - dissimilarity / distance matrix b/w dat1 and dat2 (dat1 x dat2)
%   OUTPUTS:
%       fh - figure handle of the plot
%
% Copyright (c) Brent McPherson (Indiana University), 2019. All rights reserved.
%

ndat = size(cca.dat1.loading, 1);
dmat = cca.dist.dmat;
dat1 = dmat(1:ndat, 1:ndat);
dat2 = dmat(ndat+1:end, ndat+1:end);
offd = dmat(ndat+1:end, 1:ndat);

dist = 'euclidean';
lmeth = 'average';
depth = 2;
cval = 0.8;

% distance between each data point based on subject vectors
dis1 = pdist(dat1, dist); % the same as below, but w/ different scale
dis2 = pdist(dat2, dist);

% estimate the linkage b/w 
link1 = linkage(dis1, lmeth);
link2 = linkage(dis2, lmeth);

% correlation b/w observed distances and estimated link strength
[ coph1, copd1 ] = cophenet(link1, dis1);
[ coph2, copd2 ] = cophenet(link2, dis2);

disp([ 'Brain Cophenetic Correlation Coefficient: ' num2str(coph1) ]);
disp([ 'Behavior Cophenetic Correlation Coefficient: ' num2str(coph2) ]);

% cluster the data
ci1 = cluster(link1, 'Cutoff', cval, 'Depth', depth);
ci2 = cluster(link2, 'Cutoff', cval, 'Depth', depth);

%% plot the data

fh = figure('Position', [ 160 -10 1170 880 ]); hold on

haxis(1) = subplot(3, 3, 2:3);
[ ~, dg1 ] = dendrogram(link1, size(dat1, 1), 'ColorThreshold', 'default', ...
                        'labels', repmat(' ', size(dat1, 1))); 
title('Brain');

haxis(2) = subplot(3, 3, [ 4; 7 ]);
[ ~, dg2 ] = dendrogram(link2, size(dat2, 1), 'ColorThreshold', 'default', ...
                        'Orientation', 'left', 'labels', repmat(' ', size(dat2, 1)));
bt = title('Behavior');
set(gca, 'YTickLabelRotation', 90);

haxis(3) = subplot(3, 3, [ 5, 6, 8, 9 ]);
imagesc(offd(dg2, dg1)); colorbar; axis tight;
set(gca, 'XTick', [], 'YTick', []); colorbar;
caxis([ 0 1 ]); %colormap(redbluecmap);
xlabel('Dissimilarity Between Brain and Behavior');

% add text descriptions
haxis(4) = subplot(3, 3, 1);
set(gca, 'XTick', [], 'YTick', []);
axis off;

text(0, .95, [ 'N Subjects: ' num2str(size(dat1, 2)) ]);
text(0, .85, [ 'N Brain Regions: ' num2str(size(dat1, 1)) ]);
text(0, .75, [ 'N Behaviors: ' num2str(size(dat2, 1)) ]);
text(0, .65, ' ');
text(0, .55, [ 'Brain - N Clusters: ' num2str(size(unique(ci1), 1)) ]);
text(0, .45, [ 'Brain Cophenetic: ' num2str(coph1) ]);
text(0, .35, ' ');
text(0, .25, [ 'Behavior - N Clusters: ' num2str(size(unique(ci2), 1)) ]);
text(0, .15, [ 'Behavior Cophenetic: ' num2str(coph2) ]);

%% fix axes

%pos1 = get(haxis(1), 'Position');
pos1 = [ 0.30    0.70    0.60    0.25 ];
set(haxis(1), 'Position', pos1);

%pos2 = get(haxis(2), 'Position');
pos2 = [ 0.05    0.10    0.25    0.60 ];
set(haxis(2), 'Position', pos2);

%set(bt, 'Position', [ 1050 170 0 ], 'Rotation', 90);
set(bt, 'Position', [ 5 175 0 ], 'Rotation', 90);

%pos3 = get(haxis(3), 'Position');
pos3 = [ .30 .10 .60 .60 ];
set(haxis(3), 'Position', pos3);

%pos4 = get(haxis(4), 'Position');
pos4 = [ .05 .70 .25 .25 ];
set(haxis(4), 'Position', pos4);

end

