function [ grotAAd, grotAAv ] = plotTaskBrainRegions(vars, NETd, conf, outfile)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% hard coded values - make optional?
Nkeep = 100;
ccf = 1;

% see if domain specific behaviors / regions come out of selected variables

% behaviors
% 334-8:334 % visual memory (VSMT)

% nodes
% motor strip: ?h.1, 2, 3a, 3b, 4
% mtr = [ 9, 10, 52:54, ];
% netNames([ mtr mtr+181 ])

% % decide what to keep (either task or networks
% vKEEP = 334-8:334; % visual memory task
% nKEEP = 0; % network category
% 
% % drop the variables
% NETs = NET(:, [ mtr mtr+181 ]);
% pvars = vars(:, vKEEP);
% 
% % normalize reduced networks
% NET1 = nets_demean(NETs);
% NET1 = NET1 / std(NET1(:));
% amNET = abs(mean(NET));
% NET3 = nets_demean(NET ./ repmat(amNET, size(NET, 1), 1));
% NET3(:, amNET < 0.1) = [];
% NET3 = NET3 / std(NET3(:)); % norm by mean of columns, removing badly conditioned ones
% grot = [ NET1 NET3 ];       % concat horizontally
% NETd = nets_demean(grot - conf * (pinv(conf) * grot)); % deconfound and demean
% 

% SVD reduction
uu1 = nets_svds(NETd, Nkeep); % can't run if reduced too far

%% normalize variables

% normalize reduced variables
pvarsd = palm_inormal(vars);
pvarsd(isnan(pvarsd)) = 0; % just fill in zeros - too little avaiable for other missing corrections

% deconfound ignoring missing data
for ii = 1:size(pvarsd, 2)
    grot = (isnan(pvarsd(:, ii)) == 0); 
    grotconf = nets_demean(conf(grot, :)); 
    pvarsd(grot, ii) = normalise(pvarsd(grot, ii) - grotconf * (pinv(grotconf) * pvarsd(grot, ii)));
end

clear ii grot grotconf

% % estimate "pairwise" covariance, ignoring missing data
% pvarsdCOV = zeros(size(pvarsd, 1));
% for ii = 1:size(pvarsd, 1)
%     for jj = 1:size(pvarsd, 1)
%         grot = pvarsd([ ii jj ], :); 
%         grot = cov(grot(:, sum(isnan(grot)) == 0)'); 
%         pvarsdCOV(ii, jj) = grot(1, 2);
%     end
% end
% 
% clear ii jj grot
% 
% % minor adjustment: project onto the nearest valid covariance matrix
% pvarsdCOV2 = nearestSPD(pvarsdCOV);
% 
% % SVD (eigs actually)
% uu = eigs(pvarsdCOV2, Nkeep); % can't run if this is reduced too much
% 
% % deconfound again just to be safe
% uu2 = uu - conf * (pinv(conf) * uu);

uu2 = pvarsd;

% normalize data for factor creation?
%pvarsgrot = palm_inormal(pvars);

%% CCA w/ reduced variables

[ ~, ~, pgrotR, pgrotU, pgrotV ] = canoncorr(uu1, uu2);
%[ ~, ~, pgrotR, pgrotU, pgrotV ] = canoncorr(NETs, puu2);

%% rebuilt weights

% network weights after deconfounding
grotAAd = corr(pgrotU(:, ccf), NETd(:, 1:376))';
%[ ~, pNETs ] = sort(pgrotAAd, 'descend');
%netNames(pNETs) % see what's on top

% behavior weights after deconfounding
%pgrotBBd = corr(pgrotV(:, ccf), pvarsgrot, 'rows', 'pairwise')';
%[ ~, pvarsS ] = sort(pgrotBBd, 'descend');
%varsNames(pvarsS) % see what's on top

% variance explained by each cc data set
grotAAv = mean(grotAAd(:, ccf).^2);
%pgrotBBv = mean(grotBBd(:, ccf).^2);

% redundancy of each cc component
%pgrotAAr = mean(grotAAd(:, ccf).^2) * pgrotR(ccf)^2;
%pgrotBBr = mean(grotBBd(:, ccf).^2) * pgrotR(ccf)^2;

%% map the ranks onto the surface

% load template image 
data = ft_read_cifti('Glasser_et_al_2016_HCP_MMP1.0_RVVG/mmp.32k_fs_LR.dlabel.nii');

% copy data
cdat = data;

% same order, but L/R flipped. Thanks josh...

% for every ascending label on the surface
for jj = 1:376
    
    % if the index is > 181, it's on the left
    if jj >= 181
        % grab the indices of this label
        indx = cdat.indexmax == jj - 180;
    else
        indx = cdat.indexmax == jj + 181;
    end
    
    % reassign the correlation
    cdat.indexmax(indx) = grotAAd(jj);
    
end

% save cifti mapping out the 
ft_write_cifti(outfile, cdat, 'parameter', 'indexmax');

end

