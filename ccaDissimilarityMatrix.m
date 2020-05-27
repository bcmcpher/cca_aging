function [ dmat, cmat ] = ccaDissimilarityMatrix(cca)
%[ dmat, cmat ] = ccaDissimilarityMatrix(cca);
%   Create a dissimilarity matrix between all variables in both datasets
%   used in a CCA.
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

%B = cca.dat2.scale;
%grotV = cca.dat2.factor;
grotBBd = cca.dat2.loading;

%% create mapping from variables in uu1 to uu2

disp('Inspecting dissimilarity between all variables in CCA space...');

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
    dval = 1 - abs(cval);
    
    % store in matrices
    cmat(xy(ii, 1), xy(ii, 2)) = cval;
    cmat(xy(ii, 2), xy(ii, 1)) = cval;
    dmat(xy(ii, 1), xy(ii, 2)) = dval;
    dmat(xy(ii, 2), xy(ii, 1)) = dval;
    
end

end

% for ii = 1:size(xy, 1)
%     
%     % find values - correlation and dissimilarity; others?
%     dval1 = 1 - abs(corr(z1(xy(ii, 1), :)', z1(xy(ii, 2), :)'));
%     dval2 = 1 - abs(corr(z2(xy(ii, 1), :)', z2(xy(ii, 2), :)'));
%     
%     % store in matrices
%     dmat1(xy(ii, 1), xy(ii, 2)) = dval1;
%     dmat1(xy(ii, 2), xy(ii, 1)) = dval1;
%     dmat2(xy(ii, 1), xy(ii, 2)) = dval2;
%     dmat2(xy(ii, 2), xy(ii, 1)) = dval2;
%     
% end
% 
% figure; 
% subplot(1, 2, 1); imagesc(dmat1); axis equal; axis square; axis tight; colorbar; caxis([ 0 1 ]);
% subplot(1, 2, 2); imagesc(dmat2); axis equal; axis square; axis tight; colorbar; caxis([ 0 1 ]);
