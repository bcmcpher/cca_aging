function [] = ccaMakeContributionCifti(grotAAd)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

%% THERE ARE PROBLEMS HERE, BUT WE'RE NOT FIXING THESE NOW

%% run w/ the input 

% if a full matrix is passed
if size(grotAAd, 2) > 377
    
    % reform matrix
    x = nchoosek(1:376, 2);
    y = x(:, 2); x = x(:, 1);
    mcorr = nan(376);
    for ii = 1:size(x, 1)
        mcorr(x(ii), y(ii)) = ecorr(ii);
        mcorr(y(ii), x(ii)) = ecorr(ii);
    end
    % CHANGE ecorr TO mcorr IF NETWORK EDGES ARE USED INSTEAD OF DEGREE
    
    % positive correlations
    pcorr = mcorr;
    pcorr(pcorr < 0) = NaN;
    pwght = nanmean(pcorr, 2);

    % negative correlations
    ncorr = mcorr;
    ncorr(ncorr > 0) = NaN;
    nwght = nanmean(ncorr, 2);
    
    % ecorr?
    
else
    
    % does nothing, just pass ecorr
%     % positive correlations
%     pcorr = ecorr;
%     pcorr(pcorr < 0) = NaN;
%     pwght = nanmean(pcorr, 2);
%     
%     % negative correlations
%     ncorr = ecorr;
%     ncorr(ncorr > 0) = NaN;
%     nwght = nanmean(ncorr, 2);
    
end

% load template image 
data = ft_read_cifti('glasser_2016/mmp.32k_fs_LR.dlabel.nii');

% assume the labels are ascending like they appear to be
% file is ordered:    ??? cortical rh/lh ???
% labels are ordered: cortical lh/rh, subcortical lh/rh

% for the whole network
if size(grotAAd, 2) > 377
    
    % make positive / negative / combined cifti files
    for ii = 1:3
        
        switch ii
            case 1
                tcorr = pcorr;
                tstem = 'mmp.pos.32k_fs_LR.dlabel.nii';
            case 2
                tcorr = ncorr;
                tstem = 'mmp.neg.32k_fs_LR.dlabel.nii';
            case 3
                tcorr = ecorr;
                tstem = 'mmp.all.32k_fs_LR.dlabel.nii';
        end
        
        % copy data
        cdat = data;
        
        % for every ascending labels (that's not what Josh added)
        for jj = 1:376
            
            % grab the indices of this label
            indx = cdat.indexmax == jj;
            
            % reassign the correlation
            cdat.indexmax(indx) = tcorr(jj);
            
        end
        
        ft_write_cifti([ stem tstem ], cdat, ...
            'parameter', 'indexmax');
        
    end
    
else
    
    % copy data
    cdat = data;
    
    % for every ascending labels (that's not what Josh added)
    for jj = 1:376
        
        % grab the indices of this label
        indx = cdat.indexmax == jj;
        
        % reassign the correlation
        cdat.indexmax(indx) = ecorr(jj);
        
    end
    
    ft_write_cifti([ stem 'mmp.all.32k_fs_LR.dlabel.nii' ], cdat, ...
        'parameter', 'indexmax');
    
end

clear ii jj data cdat tcorr tstem 
    
end

