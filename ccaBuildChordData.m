function [ out ] = ccaBuildChordData(mat, fname)
%[ out ] = ccaBuildChordData(mat);
%   This will take the m x n off diagonal matrix mat and create a scaled
%   version for creating an axe / bat / directed chord plot in d3.
%
%   INPUTS:
%       mat   - the off diagonal input matrix
%       fname - the name of the output .csv to load into d3
%   OUTPUTS:
%       out - the created matrix with the correct space / scaling added
%
% Copyright (c) Brent McPherson (Indiana University), 2019. All rights reserved.
%

% grab the size of each dimension, get the final data size
n1 = size(mat, 1);
n2 = size(mat, 2);
dm = n1 + n2 + 2;

%% scale the values better for the plot

% convert from dissimilarity to similarity
dat = 1 - mat;

% round and flip so more similar == higher number
%val = round(dat*1000);          % linearly scale
val = round(dat*1000).^2;       % linearly scale and square - I think this one
%val = log10(round(dat*100000)); % linearly scale and log
%val = zscore(dat);               % z-score data - need to figure out negative axes

% threshold by column?
%vmn  = mean(val(:));
%vsd = std(val(:));
%val(val < (vmn + vsd)) = 0;

% get the total count for spacing between ends
nsubj = sum(val(:));

% preallocate output
out = zeros(dm, dm);

%% fill the data in output

out(end-n1:end-1, 1:n2) = val;
out(1:n2, (n2+2):end-1) = val';
out(n2+1, dm) = nsubj;
out(dm, n2+1) = nsubj;

%% write out the csv

dlmwrite([ fname '.csv' ], out, 'delimiter', ',');

end
