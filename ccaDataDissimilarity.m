function [ dmat, cmat ] = ccaDataDissimilarity(dat, val)
%[ dmat, cmat ] = ccaDataDissimilarity(dat);
%   Create a dissimilarity matrix between all variables in both datasets.
%
%   INPUTS:
%       dat - a fit data structure
%       val - which dataset to request: 'both' (default), 'brain', 'behavior'
%   OUTPUTS:
%       dmat - the dissimilarity matrix between all variables in the CCA
%       cmat - the correlation matrix between all variables in the CCA
%
% Copyright (c) Brent McPherson (Indiana University), 2019. All rights reserved.
%

if(~exist('val', 'var') || isempty(val))
    val = 'both';
end

% pull the data from the structure
A = dat.dat1.raw';
B = dat.dat2.raw';

% dumb, but I don't know how else to clean this
B(isnan(B)) = 0;

%% create mapping from variables in uu1 to uu2

disp('Inspecting dissimilarity between all variables in raw space...');

% create scaled data in original scale x subject space 

% combine into one big [variable x subject] matrix
%z = [ z1; z2 ];

% select / merge the requested data
switch val
    case 'brain'
        disp('Only using brain data...');
        z = A;
    case 'behavior'
        disp('Only using behavior data...');
        z = B;
    otherwise
        disp('Using all data...');
        z = [ A; B ];
end

% this is a variable x subj matrix that has the scaled constribution of
% each field stored in it. The goal is to see how the vars correlate w/
% each other using the unique subject contributions in CCA space.

% for every unique value in the merged matrix
xy = nchoosek(1:size(z, 1), 2);

% preallocate output
cmat = zeros(size(z, 1));
dmat = zeros(size(z, 1));

% find the correlation b/w subjects scaled values b/w each other variable
% in CCA from both domains. 
for ii = 1:size(xy, 1)
    
    % find values - correlation and dissimilarity; others?
    cval = corr(z(xy(ii, 1), :)', z(xy(ii, 2), :)');
    dval = 1 - abs(corr(z(xy(ii, 1), :)', z(xy(ii, 2), :)'));
    
    % store in matrices
    cmat(xy(ii, 1), xy(ii, 2)) = cval;
    cmat(xy(ii, 2), xy(ii, 1)) = cval;
    dmat(xy(ii, 1), xy(ii, 2)) = dval;
    dmat(xy(ii, 2), xy(ii, 1)) = dval;
    
end

end

