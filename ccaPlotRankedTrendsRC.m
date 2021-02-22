function [ out, fh ] = ccaPlotRankedTrendsRC(dat, cca, age, mod, type, ccf, plotType, nshow, rc)
%[ out, fh ] = ccaPlotRankedTrends(dat, cca, age, mod, type, ccf, nshow);
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
%       plotType - either 'lines' or 'points'; determines how data is displayed
%       nshow - the number of loadings to show on both top and bottom.
%   OUTPUTS:
%       out - the vector of slopes/weights  and 1/0 significance estimated for the plot
%       fh  - the figure handle of the created images
%
% Copyright (c) Brent McPherson (Indiana University), 2019. All rights reserved.
%

%% parse inputs

% determine default plot type
if(~exist('plotType', 'var') || isempty(plotType))
    plotType = 'points';
end

% determine default nshow
if(~exist('nshow', 'var') || isempty(nshow))
    nshow = -1;
end

% pull the right data from inputs
if strcmp(mod, 'brain')
    
    % grab names / postitive logical of var w/ axis
    varLabs = dat.dat1.names;
    valp = cca.full.fgrotAAp;
    
    % grab data / variability
    if strcmp(type, 'load')
        vals = cca.dat1.loading;
        vars = cca.dat1.loading_se;
    end
    
    if strcmp(type, 'slope')
        [ vals, ~, vars ] = ccaTrendByDecile(dat.dat1.raw, age);
    end
    
end

% if strcmp(mod, 'behavior')
%     
%     % grab names / postitive logical of var w/ axis
%     varLabs = dat.dat2.label;
%     valp = cca.full.fgrotBBp;
%     
%     % grab data / variability
%     if strcmp(type, 'load')
%         vals = cca.dat2.loading;
%         vars = cca.dat2.loading_se;
%     end
%     
%     if strcmp(type, 'slope')
%         [ vals, ~, vars ] = ccaTrendByDecile(dat.dat2.nrm, age);
%     end       
%     
% end

%% extract and sort weights

% grab and extract the variable weights, standard error, and positive
wght = vals(:, ccf);
sder = vars(:, ccf);
valp = valp(:, ccf);

% create sig, logical vector of vars that indicate if 0 is w/in 2 SE
sig = nan(size(wght));

% for every weight
for ii = 1:size(wght, 1)
    
    % build upper / lower bounds
    lb = wght(ii) - (2 * sder(ii));
    ub = wght(ii) + (2 * sder(ii));
    
    % if CI is fully above or below zero
    if ((lb < 0) && (ub < 0)) || ((lb > 0) && (ub > 0))
        sig(ii) = 1;
    % if CI is not fully above or below zero    
    else
        sig(ii) = 0;
    end
    
    % debug significance logic
    % disp([ num2str(lb) ' | ' num2str(ub) '; ' num2str(sig(ii)) ]);
    
end

% autodetect the low / high points to plot text
mnmx = minmax(wght');
mnb = mnmx(1)*1.25;
%mxb = mnmx(1)*1.25;

disp([ 'The lower / upper bounds of the values are: ' num2str(mnmx) ]);

% sort the weights
[ sWgh, wi ] = sort(wght);

% sort the variable names
sNam = varLabs(wi);
sErr = sder(wi);
sSig = sig(wi);
sPos = valp(wi);
sRC = rc(wi);

% if nshow is empty, show all of them
if nshow < 0
    showVal = sWgh;
    showErr = sErr;
    showPos = sPos;
    showNam = sNam;
    showSig = sSig;
    showRC = sRC;
    nplot = size(sWgh, 1);
else
    % grab the nshow top / bottom to plot
    showVal = sWgh([ 1:nshow, end-(nshow-1):end ]);
    showErr = sErr([ 1:nshow, end-(nshow-1):end ]);
    showPos = sPos([ 1:nshow, end-(nshow-1):end ]);
    showNam = sNam([ 1:nshow, end-(nshow-1):end ]);
    showSig = sSig([ 1:nshow, end-(nshow-1):end ]);
    showRC = sRC([ 1:nshow, end-(nshow-1):end ]);
    nplot = 2*nshow;
end

%% plot the lines

fh = figure('Position', [ 150 600 1300 420 ]); hold on;

% for values to plot
for ii = 1:nplot
    
%     % determine the color
%     if showVal(ii) > 0
%         if showSig(ii) == 1
%             color = 'blue';
%         else
%             color = 'cyan';
%         end
%     else % for the negative values
%         if showSig(ii) == 1
%             color = 'red';
%         else
%             color = 'magenta';
%         end
%     end
    
    if showRC(ii) == 2
        color = [ 0.1 0.2 0.7 ];
    else
        color = [ 0.7 0.7 0.7 ];
    end
    
    % assign +/- sign to label for correlation
    if showPos(ii) == 1
        pNam = strcat(showNam(ii), ' (+)');
    else
        pNam = strcat(showNam(ii), ' (-)');
    end
    
    % either plot points or lines of scaled contributions
    switch plotType
        case 'points'
            
            plot([ (showVal(ii) - showErr(ii)) (showVal(ii) + showErr(ii)) ], [ ii ii ], 'color', 'black');
            plot(showVal(ii), ii, 'o', 'MarkerFaceColor', color, ...
                 'MarkerEdgeColor', 'none', 'LineWidth', 0.75, 'MarkerSize', 6);
            text(mnb, ii, pNam, 'HorizontalAlignment', 'right', 'Interpreter', 'none');
            
        otherwise
            
            % plot the line and text
            plot([ 0 showVal(ii) ], [ ii ii ], '-', 'color', color);
            text(mnb, ii, pNam, 'HorizontalAlignment', 'right', 'Interpreter', 'none');
            
    end
        
end

% vertical line at zero
%plot([ 0 0 ], [ 0 nplot ], 'black');

% format final axes and add title
set(gca, 'YTick', []);        
title([ 'Top and Bottom ' num2str(nshow) ' Contributing Weights to CCA' ]);

hold off

% return sorted order of weights w/ null sig test
out = flipud([ sWgh sErr sSig sRC ]);

end
