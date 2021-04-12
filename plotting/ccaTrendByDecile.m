function [ lcoeff, lr2, lse ] = ccaTrendByDecile(data, age)
%[ lcoeff, lr2 ] = ccaTrendByDecile(data, age);
%   Plot a trend line through the decade split of the data. The function
%   will return a slope and R^2 value for each variable. The function will
%   internally split into decades.
%
%   INPUTS:
%       data - dat.dat?.raw(:, ?); subj x 1 to estimate the trend line through
%       age  - subj x 1 age of each subject
%   OUTPUTS:
%       lcoeff - vector of slopes for each variable
%       lr2    - vector of R^2 values for each variable
%       lse    - vector of slope standard errors for each variable
%
% Copyright (c) Brent McPherson (Indiana University), 2019. All rights reserved.
%

% load data
%load([ 'camcan_594_' stem '_cca.mat' ], 'varsgrot');
%load('canoncorr_analysis_full_data.mat', 'age');

% create 7 decile categories
abin = 18:10:89;
abin(8) = abin(8) + 1;

%% this will be a loop over every variable?

lcoeff = nan(size(data, 2), 1);
lr2 = nan(size(data, 2), 1);
lse = nan(size(data, 2), 1);

for ii = 1:size(data, 2)
    
    % preallocate plot data
    tdat = nan(7, 2);
    nsub = nan(7, 1);
    
    % grab a variable
    x = data(:, ii);
    
    % for every age bin
    for jj = 1:size(abin, 2)-1 % indexes 1 ahead, so short 1 here
        
        % pull the possible ages between the categories
        ages = eval([num2str(abin(jj)) ':' num2str(abin(jj+1)-1)])';
        
        % grab index of every age in that category
        aidx = ismember(age, ages);
        
        % compute the mean, se, and total subject count of the measure
        tdat(jj, 1) = mean(x(aidx), 'omitnan');
        tdat(jj, 2) = std(x(aidx), 'omitnan') ./ sqrt(sum(aidx));
        nsub(jj) = sum(aidx);
        
    end
    
    % fit trend through average points
    z = polyfit([1:7]', tdat(:, 1), 1);
    lm = fitlm([1:7]', tdat(:, 1));

    
    % grab the trend
    lcoeff(ii) = z(1);
    lr2(ii) = lm.Rsquared.Ordinary;
    lse(ii) = lm.Coefficients.SE(2);

end

end
