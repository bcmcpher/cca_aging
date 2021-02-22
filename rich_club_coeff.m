function [ out ] = rich_club_coeff(mat, k)
%[ out ] = rich_club_coeff(mat, k);
%   Actually estimate the Rich Club Coefficient b/c the BCT does not
%

niter = 10;

%% pull the observed rich club coefficient

% the "true" rich club coefficient
rco = rcf(mat, k);

%% perform null test

% preallocate null rich club vector
nrc = nan(niter,1);

% for a fixed number of null networks
for iter = 1:niter
    
    % create a random, preserved network
    nmat = null_model_und_sign(mat);
    %nmat = randomio_und(mat, 5);
    
    % estimate the rich clube on the null
    nrc(iter) = rcf(nmat, k);
    
end

%% estimate corrected rich club coefficient

% estimate the corrected rich-club
out = rco / mean(nrc, 'omitnan');

% print simple interpretation of result
if out > 1
    disp('Rich-club structure detected.');
else
    disp('No Rich-club structure detected.');
end

end

% function to estimate rich club coefficient
function phi = rcf(mat, k)

% estimate node degree
deg = degrees_und(mat);

% find the number of nodes
Nk = sum(deg > k);

% find the number of edges in this subset network
[ kmat, ksz ] = kcore_bu(mat, k);

% get the count of edges in this subset network
Ek = sum(kmat(logical(triu(ones(ksz), 1))) > 0);

% estimate the rich club coefficient
phi = (2 * Ek) / (Nk * (Nk-1));

end