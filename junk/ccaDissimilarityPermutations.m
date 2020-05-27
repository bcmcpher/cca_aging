function [ dmat, dnull ] = ccaDissimilarityPermutations(cca)
%[ dmat, cmat ] = ccaDissimilarityMatrix(cca);
%   Create a dissimilarity matrix between all variables in both datasets
%   used in a CCA. Make additional plots based on stored permutations.
%
%   INPUTS:
%       cca - a fit CCA structure
%   OUTPUTS:
%       dmat - the dissimilarity matrix between all variables in the CCA
%       cmat - the correlation matrix between all variables in the CCA
%
% Copyright (c) Brent McPherson (Indiana University), 2019. All rights reserved.
%

%A = cca.dat1.scale;
%grotU = cca.dat1.factor;
grotAAd = cca.dat1.loading;
grotAAr = cca.dat1.rld;

%B = cca.dat2.scale;
%grotV = cca.dat2.factor;
grotBBd = cca.dat2.loading;
grotBBr = cca.dat2.rld;

% catch permutation size
nperm = size(grotAAr, 1);

%% create mapping from variables in uu1 to uu2

disp('Inspecting dissimilarity between all variables in CCA space across permutations...');

% create scaled data in original scale x subject space 

% A/B is:             pca components x cca coeffs
% grotU/grotV is:     subj x cca coeffs
% grotAAd/grotBBd is: raw input x cca coeffs
%
% grot??d? * A' scales the observation values to the CCA
%     x? is a raw inputs scaled to each cca component
% x? * grot??d transforms the cca scaled measure/subj values into a measure x subj contribution in cca space?
%     z? is a scaled contribution for every input from each subject

%x1 = grotAAd * A'; % scale original measures (nodes) to CC space
%z1 = x1 * grotU'; % scaled contribution of nodes by subject

%x2 = grotBBd * B'; % scale original measures (behaviors) to CC space
%z2 = x2 * grotV'; % scaled contribution of behaviors by subject

% combine into one big [variable x subject] matrix
%z = [ z1; z2 ];
z = [ grotAAd; grotBBd ];
zr = cat(3, grotAAr, grotBBr);
% this is a variable x subj matrix that has the scaled constribution of
% each field stored in it. The goal is to see how the vars correlate w/
% each other using the unique subject contributions in CCA space.

% for every unique value in the merged matrix
xy = nchoosek(1:size(z, 1), 2);

% preallocate output
%cmat = zeros(size(z, 1));
dmat = zeros(size(z, 1));

%cnull = zeros(size(z, 1), size(z, 1), size(zr, 1));
dnull = zeros(size(z, 1), size(z, 1), size(zr, 1));

% find the correlation b/w subjects scaled values b/w each other variable
% in CCA from both domains. 
for ii = 1:size(xy, 1)
    
    % find values - correlation and dissimilarity; others?
    dval = 1 - abs(corr(z(xy(ii, 1), :)', z(xy(ii, 2), :)'));
    
    % store in matrices
    dmat(xy(ii, 1), xy(ii, 2)) = dval;
    dmat(xy(ii, 2), xy(ii, 1)) = dval;
    
    for jj = 1:nperm
        
        % find values - correlation and dissimilarity; others?
        dval = 1 - abs(corr(zr(jj, :, xy(ii, 1))', zr(jj, :, xy(ii, 2))'));
        
        % store in matrices
        dnull(xy(ii, 1), xy(ii, 2), jj) = dval;
        dnull(xy(ii, 2), xy(ii, 1), jj) = dval;
        
    end
    
end

end


for ii = 1:100
    for jj = 1:100
        if ii == jj
            continue
        end
        figure; hold on;
        title([ 'Dissimilarity b/w ' num2str(ii) ' and ' num2str(jj) ]);
        hist(squeeze(y(ii, jj, :)), 32);
        line([ x(ii, jj) x(ii, jj) ], [ 0 30 ], 'color', 'red');
        hold off;
        pause; close all        
    end
end
        
