function [] = ccaMakeContributionCifti(cca, ccf)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

%% OLD

% THERE ARE PROBLEMS HERE, BUT WE'RE NOT FIXING THESE NOW
% run w/ the input 

% if a full matrix is passed
% if size(grotAAd, 2) > 377
%     
%     % reform matrix
%     x = nchoosek(1:376, 2);
%     y = x(:, 2); x = x(:, 1);
%     mcorr = nan(376);
%     for ii = 1:size(x, 1)
%         mcorr(x(ii), y(ii)) = ecorr(ii);
%         mcorr(y(ii), x(ii)) = ecorr(ii);
%     end
%     % CHANGE ecorr TO mcorr IF NETWORK EDGES ARE USED INSTEAD OF DEGREE
%     
%     % positive correlations
%     pcorr = mcorr;
%     pcorr(pcorr < 0) = NaN;
%     pwght = nanmean(pcorr, 2);
% 
%     % negative correlations
%     ncorr = mcorr;
%     ncorr(ncorr > 0) = NaN;
%     nwght = nanmean(ncorr, 2);
%     
%     % ecorr?
%     
% else
%     
%     % does nothing, just pass ecorr
% %     % positive correlations
% %     pcorr = ecorr;
% %     pcorr(pcorr < 0) = NaN;
% %     pwght = nanmean(pcorr, 2);
% %     
% %     % negative correlations
% %     ncorr = ecorr;
% %     ncorr(ncorr > 0) = NaN;
% %     nwght = nanmean(ncorr, 2);
%     
% end

%% parse cca structure

% pull requesed cca loadings
cval = cca.dat1.loading(:, ccf);

% load template image 
data = ft_read_cifti('glasser_2016/mmp.32k_fs_LR.dlabel.nii');

% assume the labels are ascending like they appear to be
% file is ordered:    ??? cortical rh/lh ???
% labels are ordered: cortical lh/rh, subcortical lh/rh

% check the size of the input model
if size(cval, 2) > 377
    error('Too many regions in model - This is hardcoded to the HCP-MMP parcellation.');
end

% %% This loop makes a +/-/combined map of cca correlations?
%
% % make positive / negative / combined cifti files
% for ii = 1:3
%     
%     switch ii
%         case 1
%             tcorr = pcorr;
%             tstem = 'mmp.pos.32k_fs_LR.dlabel.nii';
%         case 2
%             tcorr = ncorr;
%             tstem = 'mmp.neg.32k_fs_LR.dlabel.nii';
%         case 3
%             tcorr = ecorr;
%             tstem = 'mmp.all.32k_fs_LR.dlabel.nii';
%     end
%     
%     % copy data
%     cdat = data;
%     
%     % for every ascending labels (that's not what Josh added)
%     for jj = 1:376
%         
%         % grab the indices of this label
%         indx = cdat.indexmax == jj;
%         
%         % reassign the correlation
%         cdat.indexmax(indx) = tcorr(jj);
%         
%     end
%     
%     ft_write_cifti([ stem tstem ], cdat, ...
%         'parameter', 'indexmax');
%     
% end
    
%% create final output label
    
% copy data for output
cdat = data;

% for every ascending labels (that's not what Josh added)
for jj = 1:376
    
    % grab the indices of this label
    indx = cdat.indexmax == jj;
    
    % reassign the correlation
    cdat.indexmax(indx) = cval(jj);
    
    % simple debug out
    if jj == 363
        disp('Begin subcortical export...');
    end
    
    switch jj
        case 363
            sdat = gifti('glasser_2016/subcortical/lh.thalamus.shape.gii');
            sdat.cdata = repmat(cval(jj), size(sdat.cdata));
            save(sdat, 'glasser_2016/subcortical/lh.thalamus.cca.shape.gii', 'Base64Binary');
        case 364
            sdat = gifti('glasser_2016/subcortical/lh.caudate.shape.gii');
            sdat.cdata = repmat(cval(jj), size(sdat.cdata));
            save(sdat, 'glasser_2016/subcortical/lh.caudate.cca.shape.gii', 'Base64Binary');
        case 365
            sdat = gifti('glasser_2016/subcortical/lh.putamen.shape.gii');
            sdat.cdata = repmat(cval(jj), size(sdat.cdata));
            save(sdat, 'glasser_2016/subcortical/lh.putamen.cca.shape.gii', 'Base64Binary');
        case 366
            sdat = gifti('glasser_2016/subcortical/lh.pallidum.shape.gii');
            sdat.cdata = repmat(cval(jj), size(sdat.cdata));
            save(sdat, 'glasser_2016/subcortical/lh.pallidum.cca.shape.gii', 'Base64Binary');
        case 367
            sdat = gifti('glasser_2016/subcortical/lh.hippocampus.shape.gii');
            sdat.cdata = repmat(cval(jj), size(sdat.cdata));
            save(sdat, 'glasser_2016/subcortical/lh.hippocampus.cca.shape.gii', 'Base64Binary');
        case 368
            sdat = gifti('glasser_2016/subcortical/lh.amygdala.shape.gii');
            sdat.cdata = repmat(cval(jj), size(sdat.cdata));
            save(sdat, 'glasser_2016/subcortical/lh.amygdala.cca.shape.gii', 'Base64Binary');
        case 369
            sdat = gifti('glasser_2016/subcortical/lh.accumbens.shape.gii');
            sdat.cdata = repmat(cval(jj), size(sdat.cdata));
            save(sdat, 'glasser_2016/subcortical/lh.accumbens.cca.shape.gii', 'Base64Binary'); 
        case 370
            sdat = gifti('glasser_2016/subcortical/rh.thalamus.shape.gii');
            sdat.cdata = repmat(cval(jj), size(sdat.cdata));
            save(sdat, 'glasser_2016/subcortical/rh.thalamus.cca.shape.gii', 'Base64Binary');
        case 371
            sdat = gifti('glasser_2016/subcortical/rh.caudate.shape.gii');
            sdat.cdata = repmat(cval(jj), size(sdat.cdata));
            save(sdat, 'glasser_2016/subcortical/rh.caudate.cca.shape.gii', 'Base64Binary');
        case 372
            sdat = gifti('glasser_2016/subcortical/rh.putamen.shape.gii');
            sdat.cdata = repmat(cval(jj), size(sdat.cdata));
            save(sdat, 'glasser_2016/subcortical/rh.putamen.cca.shape.gii', 'Base64Binary');
        case 373
            sdat = gifti('glasser_2016/subcortical/rh.pallidum.shape.gii');
            sdat.cdata = repmat(cval(jj), size(sdat.cdata));
            save(sdat, 'glasser_2016/subcortical/rh.pallidum.cca.shape.gii', 'Base64Binary');
        case 374
            sdat = gifti('glasser_2016/subcortical/rh.hippocampus.shape.gii');
            sdat.cdata = repmat(cval(jj), size(sdat.cdata));
            save(sdat, 'glasser_2016/subcortical/rh.hippocampus.cca.shape.gii', 'Base64Binary');
        case 375
            sdat = gifti('glasser_2016/subcortical/rh.amygdala.shape.gii');
            sdat.cdata = repmat(cval(jj), size(sdat.cdata));
            save(sdat, 'glasser_2016/subcortical/rh.amygdala.cca.shape.gii', 'Base64Binary');
        case 376
            sdat = gifti('glasser_2016/subcortical/rh.accumbens.shape.gii');
            sdat.cdata = repmat(cval(jj), size(sdat.cdata));
            save(sdat, 'glasser_2016/subcortical/rh.accumbens.cca.shape.gii', 'Base64Binary');
        otherwise
            warning('Somehow got here.');
    end
            
% 363, lh.thalamus
% 364, lh.caudate
% 365, lh.putamen
% 366, lh.globuspallidus
% 367, lh.hippocampus
% 368, lh.amygdala
% 369, lh.accubens
% 370, rh.thalamus
% 371, rh.caudate
% 372, rh.putamen
% 373, rh.globuspallidus
% 374, rh.hippocampus
% 375, rh.amygdala
% 376, rh.accubens

end

%% ADD SERPATE OUTPUT / ASSIGNMENT FOR EACH SUBCORTICAL LABEL

% output path
stem = 'glasser_2016';

% write the relabeled output
ft_write_cifti([ stem 'mmp.cca_axis' num2str(ccf) '.32k_fs_LR.dlabel.nii' ], cdat, ...
    'parameter', 'indexmax');

end
