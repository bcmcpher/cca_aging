function [ lcoeff, fh ] = plotVarByDeciles(x, age, vname)
%[ lcoeff, fh ] = plotVarByDeciles(x, age, vname);
%   Plot a point for each decile (decade) of age in the dataset.
%
%   INPUTS:
%       x - subj x 1 vector of data to plot
%       age - subj x 1 vector of age for splitting data into deciles
%       vname - text field of variable for plot title
%   OUTPUTS:
%       lcoeff - coefficient of the trend line through the plot of data
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

% plot the mean point and se bars for each decile bin
fh = figure; hold on;
for ii = 1:size(data, 1)
    plot([ ii ii ], ...
         [ (data(ii, 1) + 2*data(ii, 2)) (data(ii, 1) - 2*data(ii, 2)) ], 'k');
    plot(ii, data(ii, 1), 'o', 'MarkerEdgeColor', [ 0 0 0 ], ...
         'MarkerFaceColor', [ .8 0 0 ], 'MarkerSize', 8);
end
set(gca, 'Xlim', [ 0.5, 7.5 ], 'XTickLabels', xlabs);

% fit trend through average points
z = polyfit([1:7]', data(:, 1), 1);

% extract the coefficients
tl = polyval(z, [1:7]');

% plot the trend line
plot(1:7, tl, '-');

lcoeff = z(1);

xlabel('Age Category');
ylabel(vname, 'Interpreter', 'none');
title([ vname ' by Age' ], 'Interpreter', 'none');
hold off;

end

