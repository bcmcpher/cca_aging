function [ fh ] = ccaPcaParameterSearch(NET, vars, conf, Nvals)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

disp('Normalizing confounds data...');

% normalize, impute w/ 0s, and add squared terms
conf = palm_inormal(conf);
conf(isnan(conf)) = 0;
conf = nets_normalise([ conf conf.^2 ]);

% idenfity the badly specified confounds
badconf = zeros(size(conf, 2), 1);
for ii = 1:size(conf, 2)
    
    Y = conf(:, ii); 
    grotKEEP =~ isnan(Y);
    grot = (Y(grotKEEP) - median(Y(grotKEEP))).^2; 
    grot = max(grot / mean(grot));  % do we have extreme outliers?
    
    % if the variables are not poorly defined
    if (sum(grotKEEP) > 250) && (std(Y(grotKEEP)) > 0) && (max(sum(nets_class_vectomat(Y(grotKEEP)))) / length(Y(grotKEEP)) < 0.95) && (grot < 100)
        % do nothing
    else
        badconf(ii) = 1;
    end
end

badconf = logical(badconf);

clear ii Y grotKEEP grot

disp(['Removing ' num2str(sum(badconf)) ' poorly defined variables...']);

% drop bad nuisance regressors / labels
conf(:, badconf) = [];

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

%% identify "bad" SMs - e.g. because of bad outliers or not enough distinct values

disp('Beginning normalization of behaior data...');

badvars = zeros(size(vars, 2), 1);
for ii = 1:size(vars, 2)
    
    Y = vars(:, ii); 
    grotKEEP =~ isnan(Y);
    grot = (Y(grotKEEP) - median(Y(grotKEEP))).^2; 
    grot = max(grot / mean(grot));  % do we have extreme outliers?
    
    % if the variables are not poorly defined
    if (sum(grotKEEP) > 250) && (std(Y(grotKEEP)) > 0) && (max(sum(nets_class_vectomat(Y(grotKEEP)))) / length(Y(grotKEEP)) < 0.95) && (grot < 100)
        % do nothing
    else
        badvars(ii) = 1;
    end
end

badvars = logical(badvars);

clear ii Y grotKEEP grot

disp(['Removing ' num2str(sum(badvars)) ' poorly defined variables...']);

% keep well defined variables
vars(:, badvars) = [];

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

%% run PCAs

Nkeep = length(Nvals);

% preallocate PCA findings
ss1 = cell(Nkeep, 1);
ss2 = cell(Nkeep, 1);

% for every parameter, keep the variability explained
for ii = 1:Nkeep
    
    % network SVD reduction
    [ ~, ss1{ii} ] = nets_svds(NETd, Nvals(ii));
    
    % behavior SVD (eigs actually)
    [ ~, ss2{ii} ] = eigs(varsdCOV2, Nvals(ii));

end

clear ii uu

%% plot the data

fh = figure('Position', [ 350 175 750 925 ]);
panel = 1;

% for every PCA size
for ii = 1:Nkeep

    % pull the variability explained in PCA
    dat1 = normalize(diag(ss1{ii}));
    dat2 = normalize(diag(ss2{ii}));
    
    % pull the size of each factor
    n1 = length(dat1);
    n2 = length(dat2);
    
    % compute the variability of each component (?)
    tve1 = dat1.^2/sum(dat1.^2);
    tve2 = dat2.^2/sum(dat2.^2);
    
    % plot brain
    subplot(Nkeep, 2, panel); hold on;
    plot(1:n1, dat1);
    plot(1:n1, dat1, 'o', 'MarkerEdgeColor', [ 0 0 0 ], 'MarkerFaceColor', [ 0 0 1 ]);
    set(gca, 'YLim', [ 0 1 ], 'XLim', [ 0 n1 ]);
    title(['Brain PCA Variance Explained: ' num2str(sum(tve1)) ]);
    %title([ 'Brain PCA Variance Explained: ' num2str(tve1(1)) ]);
    panel = panel + 1;
    hold off;
    
    % plot brain
    subplot(Nkeep, 2, panel); hold on;
    plot(1:n2, dat2);
    plot(1:n2, dat2, 'o', 'MarkerEdgeColor', [ 0 0 0 ], 'MarkerFaceColor', [ 0 0 1 ]);
    set(gca, 'YLim', [ 0 1 ], 'XLim', [ 0 n2 ]);
    title(['Behavior PCA Variance Explained: ' num2str(sum(tve2)) ]);
    %title(['Behavior PCA Variance Explained: ' num2str(tve2(1)) ]);
    panel = panel + 1;
    hold off;
    
end

end
