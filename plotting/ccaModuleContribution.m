function [ out, fh ] = ccaModuleContribution(cca, type, ccf, si, labels)
%[ out, fh ] = ccaModuleContribution(dat, cca, type, si);
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
%
%   OUTPUTS:
%       out - the data summarized in the figure
%       fh  - the figure handle for the requested plot
%
% Copyright (c) Brent McPherson (Indiana University), 2020. All rights reserved.
%

% the number of null permutation to determine module significance
nperm = 10000;

% parse arguments to pull requested data
switch type
    case 'brain'
        
        % pull data
        dat = cca.dat1.loading(:, ccf);
        ste = cca.dat1.loading_se(:, ccf);
        
        % find significance threshold
        lb = dat - (2*ste);
        ub = dat + (2*ste);
        
        % compute pvalue based on CI
        pvl = ones(size(cca.dat1.loading, 1), 1);
        % if CI is all below zero or all above zero, it's significant
        pvl((lb < 0) & (ub < 0) | ((lb > 0) & (ub > 0))) = 0;
        
    case 'behavior'
        
        % pull data
        dat = cca.dat2.loading(:, ccf);
        ste = cca.dat2.loading_se(:, ccf);
        
        % find significance threshold
        lb = dat - (2*ste);
        ub = dat + (2*ste);
        
        % compute pvalue based on CI
        pvl = ones(size(cca.dat1.loading, 1), 1);
        % if CI is all below zero or all above zero, it's significant
        pvl((lb < 0) & (ub < 0) | ((lb > 0) & (ub > 0))) = 0;
        
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

% get the size of the labels
nlab = size(labels, 2);

% split +/- dat
pdat = dat(dat > 0);
pste = ste(dat > 0);
ppvl = pvl(dat > 0);

ndat = dat(dat < 0);
nste = ste(dat < 0);
npvl = pvl(dat < 0);

% pull the size
nvar = size(dat, 1);
ndmn = size(ui, 1);

% grab the positive / negative count
%p_nvar = size(pdat, 1);
%n_nvar = size(ndat, 1);
%a_nvar = nvar;

% preallocate variables
vprop = nan(nlab, 1);
pv = nan(nlab, 3);
nv = nan(nlab, 3);

% build the pos/neg/abs proportion for each domain
for ii = 1:size(ui, 1)
    
    % pull the absolute count of each domain
    vprop(ii) = sum(si == ii);
    
    % pull the indices for each domain in each subset
    pi = si(dat > 0) == ii;
    ni = si(dat < 0) == ii;
    
    % pull the mean of the observations
    pv(ii, 1) = mean(pdat(pi), 'omitnan');
    nv(ii, 1) = mean(ndat(ni), 'omitnan');
            
    % pull the mean of the estimated standard error
    % we want the average variation of the measures within each domain, 
    % not the variation across the values of the domain.
    pv(ii, 2) = mean(pste(pi), 'omitnan');
    nv(ii, 2) = mean(nste(ni), 'omitnan');
    
    % run a permutation test to determine if the module is significant
    pvi = nan(nperm, 1);
    nvi = nan(nperm, 1);
    for perm = 1:nperm
        
        % if randsample fails it's because the module is empty
        
        try
            pvi(perm) = mean(randsample(ppvl(pi), size(pv, 1), true));
        catch
            pvi(perm) = nan;
        end
        
        try
            nvi(perm) = mean(randsample(npvl(ni), size(nv, 1), true));        
        catch
            nvi(perm) = nan;
        end
    end
    
    % if the simple majority across the permutations are significant, 
    % the module is labeled as significant.
    if all(isnan(pvi))
        pv(ii, 3) = nan;
    else
        if mean(pvi, 'omitnan') > 0.50
            pv(ii, 3) = 1;
        else
            pv(ii, 3) = 0;
        end
    end
    
    if all(isnan(nvi))
        nv(ii, 3) = nan;
    else
        
        if mean(nvi, 'omitnan') > 0.50
            nv(ii, 3) = 1;
        else
            nv(ii, 3) = 0;
        end
    end
    
end
   
clear ii pi ni ai

% pull the proportion observed
vprop = vprop / nvar;

% scale across +/- for scaling word cloud output
wcsz = [ pv(:,1); abs(nv(:,1)) ];
wcsz = wcsz ./ sum(wcsz, 'omitnan');

% clear nan from data, plot as not significant zeros
pv(isnan(pv(:,1)),3) = 1;
pv(isnan(pv)) = 0;
nv(isnan(nv(:,1)),3) = 1;
nv(isnan(nv)) = 0;

% basic plotting options
offset = 0.15;
marksz = 6;

% start figure
fh{1} = figure; hold on;

% for every domain, plot the scaled points
for ii = 1:ndmn
    
    % determine the color
    if pv(ii, 3) == 0
        pcolor = 'blue';
    else
        pcolor = 'cyan';
    end
    
    if nv(ii, 3) == 0
        ncolor = 'red';
    else
        ncolor = 'magenta';
    end
    
    % estimate upper / lower bounds for axes
    pub = pv(ii,1) + pv(ii,2);
    plb = pv(ii,1) - pv(ii,2);
    nub = nv(ii,1) + nv(ii,2);
    nlb = nv(ii,1) - nv(ii,2);
    
    % plot errorbar / point for positive contribution
    plot([ ii-offset ii-offset ], [ pub plb ], 'k');
    plot(ii-offset, pv(ii,1), 'o', 'MarkerFaceColor', pcolor, ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'MarkerSize', marksz);
    
    % plot errorbar / point for negative contribution
    plot([ ii+offset ii+offset ], [ nub nlb ], 'k');
    plot(ii+offset, nv(ii,1), 'o', 'MarkerFaceColor', ncolor, ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'MarkerSize', marksz);
    
end

set(gca, 'XLim', [ 0 ndmn+1 ], 'XTick', 1:ndmn, 'XTickLabel', labels, 'XTickLabelRotation', 90);
hold off;

% start figure
fh{2} = figure; hold on;

% for every domain, plot the scaled points
for ii = 1:size(ui, 1)
        
    % determine the color
    if pv(ii, 3) == 0
        pcolor = 'blue';
    else
        pcolor = 'cyan';
    end
    
    if nv(ii, 3) == 0
        ncolor = 'red';
    else
        ncolor = 'magenta';
    end
    
    % plot the scaled text for the word cloud if the values are non-zero
    
    %if pv(ii,1) ~= 0
    if ~isnan(wcsz(ii))
        text(1, ii, labels(ii), 'HorizontalAlignment', 'center', ...
             'color', pcolor, 'FontSize', (100 * wcsz(ii)));
    end
    
    %if nv(ii,1) ~= 0
    if ~isnan(wcsz(ii+ndmn))
        text(2, ii, labels(ii), 'HorizontalAlignment', 'center', ...
             'color', ncolor, 'FontSize', (100 * wcsz(ii+ndmn)));
    end
        
end

set(gca, 'XLim', [ 0 3 ], 'YLim', [ 0.5 ndmn+.5 ], ...
    'XTick', 1:2, 'YTick', [], ...
    'XTickLabel', {'Positive', 'Negative'});

hold off;

% simple debugging output
try % lazy fix b/w brain and behavior
    out = [ vprop pv nv ];
catch
    out = [ vprop' pv nv ];
end

end
