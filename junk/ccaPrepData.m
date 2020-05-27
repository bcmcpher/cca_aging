function [ uu1, uu2, ss1, dd, NETd, varsgrot, varsNames, conf, confNames ] = ccaPrepData(NET, vars, varsQconf, varsNames, confNames, Nkeep)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Nkeep = 100;
%Nperm = 1001;

%% setup confounds matrix

disp('Normalizing confounds data...');

% clean up confounds data
%varsQconf(isnan(varsQconf)) = 0;
%varsQconf(isinf(varsQconf)) = 0;

% normalize, impute w/ 0s, and add squared terms
conf = palm_inormal(varsQconf);
conf(isnan(conf)) = 0;
conf = nets_normalise([ conf conf.^2 ]);

% create squared nuissance names
confNames = [ confNames cellfun(@(c)[ c '^2' ], confNames, 'uni', false) ];

% idenfity the badly specified confounds
badconf = [];
for ii = 1:size(conf, 2)
    
    Y = conf(:, ii); 
    grotKEEP =~ isnan(Y);
    grot = (Y(grotKEEP) - median(Y(grotKEEP))).^2; 
    grot = max(grot / mean(grot));  % do we have extreme outliers?
    
    % if the variables are not poorly defined
    if (sum(grotKEEP) > 250) & (std(Y(grotKEEP)) > 0) & (max(sum(nets_class_vectomat(Y(grotKEEP)))) / length(Y(grotKEEP)) < 0.95) & (grot < 100)
        
        % the 3rd thing above is:  is the size of the largest equal-values-group too large?
        % ii = ii; % do nothing
                
    else
        
        % print a debug message and add values to bad variables 
        %[ii sum(grotKEEP) std(Y(grotKEEP)) max(sum(nets_class_vectomat(Y(grotKEEP))))/length(Y(grotKEEP)) grot]
        badconf = [ badconf ii ];
    
    end
end

clear ii Y grotKEEP grot

disp(['Removing ' num2str(length(badconf)) ' poorly defined variables...']);

% drop bad nuisance regressors / labels
conf(:, badconf) = [];
confNames(badconf) = [];
% they should all be kept

%% prepare main netmat matrix 

disp('Normalizing and deconfounding network data...');

NET1 = nets_demean(NET); % if pass log10 of degree
NET1 = NET1 / std(NET1(:)); % need to zero out inf vals here
amNET = abs(mean(NET)); % need to zero out inf vals here
NET3 = nets_demean(NET ./ repmat(amNET, size(NET, 1), 1)); % need to zero out NaN here
NET3(:, amNET < 0.1) = [];
NET3 = NET3 / std(NET3(:)); % norm by mean of columns, removing badly conditioned ones
grot = [ NET1 NET3 ];       % concat horizontally
NETd = nets_demean(grot - conf * (pinv(conf) * grot)); % deconfound and demean % need to zero nan

disp('Performing PCA on network data...');

% SVD reduction
[ uu1, ss1, vv1 ] = nets_svds(NETd, Nkeep);

%% identify "bad" SMs - e.g. because of bad outliers or not enough distinct values

disp('Beginning normalization of behaior data...');

badvars = [];
for ii = 1:size(vars, 2)
    
    Y = vars(:, ii); 
    grotKEEP =~ isnan(Y);
    grot = (Y(grotKEEP) - median(Y(grotKEEP))).^2; 
    grot = max(grot / mean(grot));  % do we have extreme outliers?
    
    % if the variables are not poorly defined
    if (sum(grotKEEP) > 250) & (std(Y(grotKEEP)) > 0) & (max(sum(nets_class_vectomat(Y(grotKEEP)))) / length(Y(grotKEEP)) < 0.95) & (grot < 100)
        
        % the 3rd thing above is:  is the size of the largest equal-values-group too large?
        % ii = ii; % do nothing
        
    else
        
        % print a debug message and add values to bad variables 
        %[ii sum(grotKEEP) std(Y(grotKEEP)) max(sum(nets_class_vectomat(Y(grotKEEP))))/length(Y(grotKEEP)) grot]
        badvars = [ badvars ii ];
    
    end
end

clear ii Y grotKEEP grot

% print badly conditioned variables
%varsNames(badvars)'
% make sure they are ok to drop w/ next line

disp(['Removing ' num2str(size(badvars, 1)) ' poorly defined variables...']);

% keep well defined variables
vars(:, badvars) = [];
varsNames(:, badvars) = [];

% "impute" missing vars data / normalize
varsd = palm_inormal(vars);

disp('Normalizing and deconfounding behavior data...');

% deconfound ignoring missing data
for ii = 1:size(varsd, 2)
    grot = (isnan(varsd(:, ii)) == 0); 
    grotconf = nets_demean(conf(grot, :)); 
    varsd(grot, ii) = normalise(varsd(grot, ii) - grotconf * (pinv(grotconf) * varsd(grot, ii)));
end

clear ii grot grotconf

% estimate "pairwise" covariance, ignoring missing data
varsdCOV = zeros(size(varsd, 1));
for ii = 1:size(varsd, 1)
    for jj = 1:size(varsd, 1)
        grot = varsd([ ii jj ], :); 
        grot = cov(grot(:, sum(isnan(grot)) == 0)'); 
        varsdCOV(ii, jj) = grot(1, 2);
    end
end

clear ii jj grot

% minor adjustment: project onto the nearest valid covariance matrix
varsdCOV2 = nearestSPD(varsdCOV);

disp('Performing PCA on network data...');

% SVD (eigs actually)
[ uu, dd ] = eigs(varsdCOV2, Nkeep);

% deconfound again just to be safe
uu2 = uu - conf * (pinv(conf) * uu);

% normalize data for factor creation
varsgrot = palm_inormal(vars);

disp('Returning preprocessed data.');

end
