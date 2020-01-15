function [ wght, mnmxpct, fh ] = ccaPlotRankedTrends(dat, cca, age, mod, type, ccf, nshow)
%[ wght, fh ] = ccaPlotRankedTrends(dat, cca, age, mod, type, ccf, nshow);
%   Estimate and plot the ranked trends of the requested data sets as a
%   line plot with the labels, indicating the highest / lowest loading
%   values.
%
%   INPUTS:
%       dat   - the preprocessed data from the CCA
%       cca   - the findings from the CCA
%       age   - a 1 x subj vector of everyone's age
%       mod   - 'brain' or 'beahavior'; selects the data to plot
%       type  - 'slope' or 'load'; rank the CCA loadings or slopes from the data
%       ccf   - the cc to plot data from
%       nshow - the number of loadings to show on both top and bottom.
%   OUTPUTS:
%       wght - the vector of slopes or weights estimated for the plot
%       mnmxpct - the upper / lower 2.5 perctile bounds from a resampled mean
%       fh - the figure handle of the created images
%
% Copyright (c) Brent McPherson (Indiana University), 2019. All rights reserved.
%

%stem = 'deg';
%mod = 'brain';
%mod = 'behavior';
%type = 'slope';
%type = 'load';
%ccf = 1;
%nshow = 20;

%% parse inputs

if strcmp(mod, 'brain')
    
    if strcmp(type, 'load')
        varLabs = dat.dat1.names;
        vals = cca.dat1.loading;
    end
    
    if strcmp(type, 'slope')
        varLabs = dat.dat1.names;
        vals = ccaTrendByDecile(dat.dat1.raw, age);
    end
    
end

if strcmp(mod, 'behavior')
    
    % load behavior data
    if strcmp(type, 'load')
        vals = cca.dat2.loading;
        varLabs = dat.dat2.label;
    end
    
    if strcmp(type, 'slope')
        vals = ccaTrendByDecile(dat.dat2.nrm, age);
        varLabs = dat.dat2.label';
    end       
    
end

%% extract and sort weights

% grab and extract the variable weights
wght = vals(:, ccf);

% this is fast enough, run a bunch
Nperm = 100000;

% null distribtuion of average weight to pull from
ndist = nan(Nperm, 1);
for ii = 1:Nperm
    rvals = randsample(wght, size(wght, 1), 'true');
    ndist(ii) = mean(rvals);
end

% pull p=0.05 two tail test
mnmxpct = prctile(ndist, [ 2.5 97.5 ]);

% figure; hold on;
% hist(ndist, 64);
% line([ mnmxpct(1) mnmxpct(1) ], [ 0 6000 ]);
% line([ mnmxpct(2) mnmxpct(2) ], [ 0 6000 ]);
% hold off;

% logical to id weights above / below threshold
sig = wght < mnmxpct(1) | wght > mnmxpct(2);

% autodetect the low / high points to plot text
mnmx = minmax(wght');
mnb = mnmx(1)*1.15;
%mxb = mnmx(1)*1.15;

disp([ 'The lower / upper bounds of the values are: ' num2str(mnmx) ]);

% sort the weights
[ ws, wi ] = sort(wght);

% sort the variable names
sortNames = varLabs(wi);
ssig = sig(wi);

% grab the nshow top / bottom to show
showValues = ws([ 1:nshow, end-(nshow-1):end ]);
showNames = sortNames([ 1:nshow, end-(nshow-1):end ]);

%% plot the lines

fh = figure('Position', [ 150 600 1300 420 ]); hold on;

% for positive / negative
for ii = 1:(2*nshow)
    
    % for the positive values
    if showValues(ii) > 0
        if ssig(ii) == 1
            color = 'blue';
        else
            color = 'cyan';
        end
    else % for the negative values
        if ssig(ii) == 1
            color = 'red';
        else
            color = 'magenta';
        end
    end
    plot([ 0 showValues(ii) ], [ ii ii ], '-', 'color', color);
    text(mnb, ii, showNames{ii}, 'HorizontalAlignment', 'right', 'Interpreter', 'none');
end

set(gca, 'XTick', [], 'YTick', []);
title([ 'Top and Bottom ' num2str(nshow) ' Contributing Weights to CCA' ]);

hold off

% return original order of weights w/ null sig test
wght = [ wght sig ];

end

