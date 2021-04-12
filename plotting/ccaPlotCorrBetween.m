function [ fh ] = ccaPlotCorrBetween(mat1, mat2, ci, niter)
%[ fh ] = ccaPlotCorrBetwen(mat1, mat2, ci);
%   This will plot the values between 2 dissimilarity matrices and compute
%   a boostrapped correlation between them to determine their relationship.
%
%   INPUTS:
%       mat1  - the first dissimilarity matrix
%       mat2  - the second dissimilarity matrix
%       ci    - the index to divide modules by (default = 10)
%       niter - the number of iterations for correlation bootstrap test
%               (default = 10000)
%
%   OUTPUTS:
%       fh - the figure handle of the plot
%
% Copyright (c) Brent McPherson (Indiana University), 2020. All rights reserved.
%

if(~exist('ci', 'var') || isempty(ci))
    ci = 10;
end

if(~exist('niter', 'var') || isempty(niter))
    niter = 10000;
end

% make sure matrices are the same size
if size(mat1) ~= size(mat2)
    error('Matrices have to have the same size.');
end

% pull the number of module lables
nmod = size(mat1, 1);

% unique indices + diagonal
xyi = nchoosek(1:nmod, 2);
xyi = [ xyi; [ 1:nmod; 1:nmod ]' ];

% create figure
fh = figure;

% make 3 panels
for sp = 1:3
    
    subplot(1, 3, sp); hold on;
    
    % preallocate empty matrices for correlation estimate
    data = [];
        
    % for every point
    for ii = 1:size(xyi, 1)
        
        % pull the points
        ptx = xyi(ii, 1);
        pty = xyi(ii, 2);
        
        % pull the values
        cval = mat1(ptx, pty);
        dval = mat2(ptx, pty);
        
        % plot different colors if x/y are within/between brain/behavior
        
        % check x index
        if ptx > ci
            xc = 'behavior';
        else
            xc = 'brain';
        end
        
        % check y index
        if pty > ci
            yc = 'behavior';
        else
            yc = 'brain';
        end
        
        if (strcmp(xc, 'brain') && strcmp(yc, 'brain')) && sp == 2
            ttl = 'Brain-Brain Interactions';
            %xlim = [ 0 1 ]; ylim = [ 0 1 ];
            plot(cval, dval, 'o', 'MarkerFaceColor', [ 0.1 0.3 0.7 ], ...
                'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'MarkerSize', 6); % plot within brain
            data = [ data; cval, dval ];
        end
        
        if (strcmp(xc, 'behavior') && strcmp(yc, 'behavior')) && sp == 3
            ttl = 'Behavior-Behavior Interactions';
            %xlim = [ 0 1 ]; ylim = [ 0 1 ];
            plot(cval, dval, 'o', 'MarkerFaceColor', [ 0.3 0.7 0.1 ], ...
                'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'MarkerSize', 6); % plot within behavior
            data = [ data; cval, dval ];
        end
        
        if (strcmp(xc, 'brain') && strcmp(yc, 'behavior')) && sp == 1
            ttl = 'Brain-Behavior Interactions';
            %xlim = [ 0 1 ]; ylim = [ 0 1 ];
            plot(cval, dval, 'o', 'MarkerFaceColor', [ 0.7 0.3 0.1 ], ...
                'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'MarkerSize', 6); % plot brain x behavior
            data = [ data; cval, dval ];
        end
        
    end
    
    % pull the size of data for correlation tests
    ndat = size(data, 1);
    
    % store the resampling
    rdat = nan(niter, 2);
    
    % for every random iteration
    for iter = 1:niter
        
        % create a random index w/ replacement
        ri = randsample(1:ndat, ndat, true);
        
        % estimate the resampled corr
        rdat(iter, 1) = corr(data(ri, 1), data(ri, 2));
    
        % estimate the random null corr
        rdat(iter, 2) = corr(data(ri, 1), data(:, 2));
        
    end
        
    hold off;
    title({ttl, [num2str(mean(rdat(:, 1))) ' +/- ' num2str(std(rdat(:, 1))) ' (p = ' num2str(1-(sum(rdat(:, 2) > mean(rdat(:, 1)))/niter)) ')' ]});
    xlabel('mat1 Dissimilarity');
    ylabel('mat2 Dissimilarity');
    axis square; axis equal; axis tight
    set(gca, 'XLim', [ .90*min(data(:, 1)) 1.1*max(data(:, 1)) ], ...
        'YLim', [ .90*min(data(:, 2)) 1.1*max(data(:, 2)) ]);
    
end

end

