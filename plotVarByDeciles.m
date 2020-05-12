function [ out, fh ] = plotVarByDeciles(x, age, model, vname)
%[ lcoeff, fh ] = plotVarByDeciles(x, age, vname);
%   Plot a point for each decile (decade) of age in the dataset.
%
%   INPUTS:
%       x - subj x 1 vector of data to plot
%       age - subj x 1 vector of age for splitting data into deciles
%       model - string indicating the type of model to fit
%       vname - text field of variable for plot title
%   OUTPUTS:
%       out - relevant fit coefficients based on call
%       fh - figure handle for the generated plot
%
% Copyright (c) Brent McPherson (Indiana University), 2019. All rights reserved.
%

% create 7 decile categories
abin = 18:10:89;
abin(8) = abin(8) + 1;
xlabs = {'18-27', '28-37', '38-47', '48-57', '58-67', '68-77', '78-88', ''};

% preallocate plot data
data = nan(7, 2);
nsub = nan(7, 1);

% for every age bin
for ii = 1:size(abin, 2)-1 % indexes 1 ahead, so short 1 here
    
    % pull the possible ages between the categories
    ages = eval([num2str(abin(ii)) ':' num2str(abin(ii+1)-1)])';
    
    % grab index of every age in that category 
    aidx = ismember(age, ages); 
    
    % compute the mean, se, and total subject count of the measure
    data(ii, 1) = mean(x(aidx), 'omitnan');
    data(ii, 2) = std(x(aidx), 'omitnan') ./ sqrt(sum(aidx));
    nsub(ii) = sum(aidx);
    
end

%% fit trend line through binned points

% create data from binned points
xd = abin(1:end-1)'+5;
yd = data(:, 1);

% user defined optimization
switch model
    
    case {'exp', 'exponential'}
        
        % exponential function to optimize
        fun = @(p,x) p(1).*x.^p(2) + p(3);
        
        % initial linear fit for starting params
        L = polyfit(xd, yd, 1);
        
        % turn off verbose output; set up and estimate exponential function
        options = optimset('Display','off');
        fitv = lsqcurvefit(fun, [ L(1) 1 L(2) ], xd, yd, [], [], options);
        
        % get predicted trends
        fitp = fun(fitv, xd);
        
        % define term / value for plot
        pname = 'Exponential';
        lcoeff = fitv(2);
        
    case 'linear'
        fitv = polyfit(xd, yd, 1);
        fitp = polyval(fitv, xd);
        pname = 'Linear';
        lcoeff = fitv(1);
        
    case 'quadratic'
        fitv = polyfit(xd, yd, 2);
        fitp = polyval(fitv, xd);
        pname = 'Quadratic';
        lcoeff = fitv(1);
        
    otherwise
        error('Not a valid model request.');
        
end

%% estimate R2

% pull parameters from data
nobs = size(xd, 1);
nparam = size(fitv, 2);

% calculate R2
R2 = 1 - ((nobs-1/nobs-nparam)*(sum((yd-fitp).^2)/sum(yd.^2)));

% estimate AIC and AICc
AIC = (2*nparam) - (2*log(R2));
AICc = AIC * (2*nparam.^2 + 2*nparam)/(nobs - nparam - 1);

%% plot the data and the requested trend

% create figure
fh = figure; hold on;

% plot fit trend line behind points
plot(xd, fitp, 'b-');

% plot points with error bars
for ii = 1:size(xd, 1)
    plot([ xd(ii) xd(ii) ], ...
         [ (data(ii, 1) + 2*data(ii, 2)) (data(ii, 1) - 2*data(ii, 2)) ], 'k');
    plot(xd(ii), data(ii, 1), 'o', 'MarkerEdgeColor', [ 0 0 0 ], ...
         'MarkerFaceColor', [ .8 0 0 ], 'MarkerSize', 8);
end
hold off;
set(gca, 'Xlim', [ 16 90 ], 'XTick', xd, 'XTickLabels', xlabs);

xlabel('Age Category');
ylabel(vname, 'Interpreter', 'none');
title({vname, [ pname ' Trendline - p = ' num2str(lcoeff) ]});

%% create out

out.fit = pname;
out.param = fitv;
out.lcoeff = lcoeff;
out.R2 = R2;
out.AIC = AIC;
out.AICc = AICc;

end
