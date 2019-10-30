function [ lcoeff, fh ] = plotVarByAge(x, age, vname)
%[ fh ] = plotVarByAge(x, age, vname);
%   Create a plot to show [ age x score ] for identifying trends in a
%   dataset.
%
%   INPUTS:
%       x   - a subj x 1 vector of scores
%       age - a subj x 1 vector of ages
%       vname - a string to plot as the title of the plot
%   OUTPUTS:
%       lcoeff - slope of the trend line through the points
%       fh - a figure handle of the plot
%
% Copyright (c) Brent McPherson (Indiana University), 2019. All rights reserved.
%

%% pull data

% handle random missing values
y = age;
y(isnan(x)) = [];
x(isnan(x)) = [];

% fit trend line
z = fit(y, x, 'poly1');
%z = fit(y, x, 'poly2');
%z = fit(y, x, 'exp1');

% plot data w/ trend line
fh = figure; hold on;
set(gca, 'Xlim', [ 15, 95 ]);
tl = plot(z);
set(tl, 'LineWidth', 3, 'Color', [ .8 0 0 ]);
plot(y, x, 'o', 'MarkerEdgeColor', [ 0 0 0 ], ...
    'MarkerFaceColor', [ .8 0 0 ]);
set(gca, 'XLim', [ 15 90 ]);
xlabel('Age', 'Interpreter', 'none');
ylabel(vname, 'Interpreter', 'none');
title([ 'Age by ' vname ], 'Interpreter', 'none');
lg = gca; 
legend(lg,'off');

% store slope output
lcoeff = z.p1;

hold off;


end

