function [ fh ] = plotViolinByDeciles(x, age, out)
%[ lcoeff, fh ] = plotViolinByDeciles(x, age, vname);
%   Plot a point for each decile (decade) of age in the dataset while showing all observations.
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

% ADD OPTION TO PREDICT ON AGE DIRECTLY, NOT JUST BINS

% create 7 decile categories
abin = 18:10:89;
abin(8) = abin(8) + 1;
xlabs = {'18-27', '28-37', '38-47', '48-57', '58-67', '68-77', '78-88', '88+'};

% make age a category
agec = age;
agec(agec<=27 & agec>=18) = 1;
agec(agec<=37 & agec>=28) = 2;
agec(agec<=47 & agec>=38) = 3;
agec(agec<=57 & agec>=48) = 4;
agec(agec<=67 & agec>=58) = 5;
agec(agec<=77 & agec>=68) = 6;
agec(agec<=87 & agec>=78) = 7;
agec(agec>=88) = 8;

% plot
fh = figure;
violinplot(x, agec, 'ViolinColor', [ 0.1 0.2 0.7 ]);
set(gca, 'TickDir', 'out', 'Xlim', [ -0.5 8.5 ], 'XTickLabel', xlabs);
set(gcf, 'renderer', 'painters');

saveas(fh, fullfile('figs', 'violin', out), 'epsc');

end
