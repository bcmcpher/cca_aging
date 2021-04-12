function [ mds, cmd, mderror, eigvals, neigval, fh ] = ccaMDScale(mat, ndim)
%[ mds, cmd, mderror, eigvals, neigvals, fh ] = ccaMDScale(mat, ndim);
%   This function performs a classic and non-metric multidimensional
%   scaling on an input matrix to a user specified number of dimensions.
%
% THIS IS NOT USED IN THE PAPER - IT WAS AN EXPLORATORY IDEA AND EASY ENOUGH TO MAKE A FXN
%   
%   INPUTS:
%       mat  - nxn square input matrix of distances / dissimilarities b/w all values
%       ndim - the number of input dimensions for the non-metric scaling
%   OUTPUTS:
%       mds     - n x ndim coordinate matrix of the non-metric MDS
%       cmd     - n x dim coordinate matrix of the classical MDS
%       mderror - the estimations error for each subsequent embedding dimension
%       eigvals - the eigenvalues for each classical MDS dimension
%       neigval - the normalized eigenvalues for each classical MDS dimension
%       fh      - figure handle for the qc plot
%
% https://www.mathworks.com/help/stats/examples/non-classical-multidimensional-scaling.html
%
% Copyright (c) Brent McPherson (Indiana University), 2019. All rights reserved.
%

%% nonmectic multidimensional scaling

% vectorize the upper diagonal from the input matrix
dissim = squareform(mat, 'tovector');

disp([ 'Compute multidimensional scaling on input data to find ' num2str(ndim) ' factors...' ]);
[ mds, stress, disparities ] = mdscale(dissim, ndim, 'criterion', 'metricstress');

% print estimated stress from mdscale - lower is better
disp([ 'MD Scaling Minimized Stress: ' num2str(stress) ]);

% convert multidimensional embedding into distances
distances = pdist(mds);

%% classical multidimensional scaling

% estimate classical multidimensional scaling
[ cmd, eigvals ] = cmdscale(dissim);

% grab the estimated number of dimensions
ncdim = size(cmd, 2);

% normalize the eigenvalues
neigval =  eigvals / max(abs(eigvals));

% find the highest value in the matrix
mxval = max(mat(:));

% preallocate reconstruction error from n dimensions
mderror = nan(ncdim, 1);

% for every md axis
for ii = 1:ncdim
    
    % estimate the reconstruction error for only the first N axes
    mderror(ii) = max(abs(dissim - pdist(cmd(:, 1:ii))));
    
end

clear ii

%% plot the data

% autodetect the highest value and scale appropriately
lim = max([ dissim distances ]);
lim = lim + (lim * 0.15);

fh = figure('Position', [ 60 400 1435 400 ]);

subplot(1, 3, 1); hold on;
plot(1:ncdim, mderror, '.-');
plot([ 1 ncdim ], [ mxval mxval ], 'r-');
axis square;
title({'Classical Multidimensional Scaling', [ num2str(ncdim) ' Dimensions' ]});
xlabel('N dimensions');
ylabel('Max Reconstruction Error');
set(gca, 'XLim', [ 0 ncdim ]);
hold off;

subplot(1, 3, 2);
plot(dissim', distances', 'bo', [ 0 lim ], [ 0 lim ], 'k--');
axis square;
title([ 'Shepard Plot; Stress = ' num2str(stress) ]);
xlabel('Dissimilarities');
ylabel('Distances');

subplot(1, 3, 3);

% sort the rows of the inputs to draw the path throught the point cloud
% the smoother the line is to the points, the better the fit
[ ~, ord ] = sortrows([ disparities' dissim' ]);
plot(dissim, distances, 'bo', ...
     dissim(ord), disparities(ord), 'r.-');
axis square;
title('Nonmetric Shepard Plot');
xlabel('Dissimilarities')
ylabel('Distances/Disparities')
legend({'Distances' 'Disparities'}, 'Location','NorthWest');
set(gca, 'Xlim', [ 0 lim ], 'Ylim', [ 0 lim ]);

end

