function [ fh, out ] = ccaScaledContribution(cca, type, ccf, si, labels, summary, show)
%[ fh ] = ccaScaledContribution(dat, cca, type, si);
%   Take a CCA output and provide indices for how to sort variables from a 
%   specific dataset into a scaled contribution from each domain in the datasets.
%
%   INPUTS:
%       cca     - fit cca object
%       type    - 'brain' or 'behavior' - identify the data to summarize
%       ccf     - identify the axis from the CCA to inspect
%       si      - a [1 x var] vector assigning variables to domains for summary
%       labels  - a [1 x labels] cell arrary of names for modules. Must correspond 
%                 to the ascending order of the values provided in si.
%       summary - the summary operation to perform on the data. Default 'mean'.
%                 Also allows: 'sum', 'std', 'median'  
%       show    - optionally plot a point behind the scaled text for the
%                 values as well.
%
%   OUTPUTS:
%       fh  - the figure handle for the requested plot
%       out - the data summarized in the figure
%
% Copyright (c) Brent McPherson (Indiana University), 2020. All rights reserved.
%

% parse arguments to pull requested data
switch type
    case 'brain'
        dat = cca.dat1.loading(:, ccf);
    case 'behavior'
        dat = cca.dat2.loading(:, ccf);
    otherwise
        error('The requested type doesn''t exist. Pick ''brain'' or ''behavior''.');
end

% sanity checks
if size(si, 1) ~= size(dat, 1)
    error('The provided sorting index doesn''t match the size of the data.');
end

% pull the unique indices
ui  = unique(si);

if size(ui, 1) ~= length(labels)
    error('There are a different number of labels in ''si'' than in ''labels''.');
end

if(~exist('summary', 'var') || isempty(summary))
    summary = 'mean';
end

if(~exist('show', 'var') || isempty(show))
    showPoints = false;
end

% get the size of the labels
nlab = size(labels, 1);

% split +/- dat
pdat = dat(dat > 0);
ndat = dat(dat < 0);
adat = abs(dat);

% pull the size
nvar = size(dat, 1);
ndmn = size(ui, 1);

% grab the positive / negative count
%p_nvar = size(pdat, 1);
%n_nvar = size(ndat, 1);
%a_nvar = nvar;

% preallocate variables
vprop = size(nlab, 1);
pv = size(nlab, 1);
nv = size(nlab, 1);
av = size(nlab, 1);

% build the pos/neg/abs proportion for each domain
for ii = 1:size(ui, 1)
    
    % pull the absolute count of each domain
    vprop(ii) = sum(si == ii);
    
    % pull the indices for each domain in each subset
    pi = si(dat > 0) == ii;
    ni = si(dat < 0) == ii;
    ai = si == ii;
    
    % pull the mean?
    switch summary
        case 'mean'
            pv(ii) = mean(pdat(pi), 'omitnan');
            nv(ii) = mean(ndat(ni), 'omitnan');
            av(ii) = mean(adat(ai), 'omitnan');
            
        case 'std' % doesn't make much sense?
            pv(ii) = std(pdat(pi), 'omitnan');
            nv(ii) = std(ndat(ni), 'omitnan');
            av(ii) = std(adat(ai), 'omitnan');
            % should this be the size of the circle somehow?
            
        case 'sum'
            pv(ii) = sum(pdat(pi), 'omitnan');
            nv(ii) = sum(ndat(ni), 'omitnan');
            av(ii) = sum(adat(ai), 'omitnan');
        
        case 'median'
            pv(ii) = median(pdat(pi), 'omitnan');
            nv(ii) = median(ndat(ni), 'omitnan');
            av(ii) = median(adat(ai), 'omitnan');
            
        otherwise
            warning('Invalid summary requested. Defaulting to ''mean''.');
            pv(ii) = mean(pdat(pi), 'omitnan');
            nv(ii) = mean(ndat(ni), 'omitnan');
            av(ii) = mean(adat(ai), 'omitnan');
            
    end 
    
end
   
clear ii pi ni ai

% pull the proportion observed
vprop = vprop / nvar;

% replace nan w/ 0 (no pos/neg w/in domain)
pv(isnan(pv)) = min(pv);
nv(isnan(nv)) = max(nv);
av(isnan(av)) = min(av);

% replace any zeros with a small number
pv(pv == 0) = 0.0001;
nv(nv == 0) = 0.0001;
av(av == 0) = 0.0001;

% scale the values
pv = pv ./ sum(pv);
nv = nv ./ sum(nv);
av = av ./ sum(av);

% start figure
fh = figure; hold on;

% for every domain, plot the scaled points
for ii = 1:size(ui, 1)
    
    if showPoints == true
        
        % plot the proportion of questions at 1
        plot(1, ii, 'o', 'MarkerFaceColor', [.1 .45 .95], ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.75, ...
            'MarkerSize', 100 * vprop(ii));
        
        % plot the proportion of positive loadings at 2
        plot(2, ii, 'o', 'MarkerFaceColor', [.1 .45 .95], ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.75, ...
            'MarkerSize', 100 * pv(ii));
        
        % plot the proportion of negative loadings at 3
        plot(3, ii, 'o', 'MarkerFaceColor', [.1 .45 .95], ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.75, ...
            'MarkerSize', 100 * abs(nv(ii)));
        
        % plot the proportion of absolute loadings at 4
        plot(4, ii, 'o', 'MarkerFaceColor', [.1 .45 .95], ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.75, ...
            'MarkerSize', 100 * av(ii));
        
    end
    
    text(1, ii, labels(ii), 'HorizontalAlignment', 'center', 'FontSize', (100 * vprop(ii))+6);
    text(2, ii, labels(ii), 'HorizontalAlignment', 'center', 'FontSize', (100 * pv(ii))+6);
    text(3, ii, labels(ii), 'HorizontalAlignment', 'center', 'FontSize', (100 * nv(ii))+6);
    text(4, ii, labels(ii), 'HorizontalAlignment', 'center', 'FontSize', (100 * av(ii))+6);
    
    
end

set(gca, 'XLim', [ 0 5 ], 'YLim', [ 0.5 ndmn+.5 ], ...
    'XTick', 1:4, 'YTick', [], ...
    'XTickLabel', {'Proportion', 'Positive', 'Negative', 'Absolute'});

hold off;

out = [ vprop' pv' nv' av' ];

end

% % find unique domains for behavior
% zz = cellfun(@(x) x(1:3), svar', 'UniformOutput', false);
% 
% % find unique values, counts, and 
% [ labels, zz2, si ] = unique(zz);
% 
% clear zz zz1
% 
% % not quite right, but close enough for testing
%

% labels for brain just have to be loaded
