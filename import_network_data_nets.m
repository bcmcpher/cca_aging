%% import and fill in empty nodes based on labels

% index data
files = dir('data/*density.csv');
label = dir('data/*labels_dilated.nii.gz');

% get file names
mats = {files.name}';
labs = {label.name}';

subjs = strrep(mats, '_density.csv', '');
%filePh = fopen('ccid.txt', 'w');
%fprintf(filePh, '%s\n', nsubj{:});
%fclose(filePh);

% preallocate network data
data = cell(size(mats));
prob = cell(size(mats));

% loop over and load / fix networks
for ii = 1:size(mats, 1)
    
    % load temporary data
    tmat = dlmread([ 'data/' mats{ii} ]);
    tlab = niftiRead([ 'data/' labs{ii} ]);

    % grab unique index of labels
    ulab = unique(tlab.data(:));
            
    % fix all bad - catch runs, too
    while ~all(diff(ulab) == 1)
    
        %disp(['N labels: ' num2str(length(ulab) - 1) ' - File: ' labs{ii} ' - N nodes:  ' num2str(size(tmat, 1))]);
        disp(['Iter: ' num2str(ii) ' - Subj: ' num2str(subjs{ii}) ' - N nodes:  ' num2str(size(tmat, 1))]);
        
        % find indices where the difference is too large
        dlab = diff(ulab);
        
        % get the bad indices - indx + 1 are rows / columns that are empty
        indx = ulab(dlab > 1);
        
        % for every busted index
        for jj = 1:size(indx, 1)
            
            %disp([ 'Missing ROI index: ' num2str(indx(jj)) ]);
            
            % add the empty row and column to the matrix
            tmat = [ tmat(1:indx(jj), :); zeros(1, size(tmat, 1)); tmat(indx(jj)+1:end, :) ];
            tmat = [ tmat(:, 1:indx(jj)) zeros(size(tmat, 2)+1, 1) tmat(:, indx(jj)+1:end) ];
            
            % add index to problem count
            prob{ii} = [ prob{ii} indx(jj) ];
            
            % add index to included labels
            ulab = sort([ ulab; indx(jj)+1 ]);
            
            %disp([ 'N labels: ' num2str(size(ulab, 1)) ]);
            
        end
    
        disp([ 'Fixed size: ' num2str(size(tmat, 1)) ]);
    
    end
    
    if size(tmat, 1) ~= 376
        keyboard
    end
    
    % store output
    data{ii} = tmat;
    prob{ii} = sort(prob{ii}); % what does this do?
    
end

% Subjects missing 1 or more nodes:
% CC210051
% CC310086
% CC310414
% CC320448
% CC520624
% CC620821
% CC620935
% CC710462
% CC710494
% CC711244
% CC721244
% CC721729
% CC723197

%sum(cellfun(@(x) size(x, 1), data) ~= 376)

% combine into a single array
nets = cat(3, data{:});

%save('camcan_netstats.mat', 'nsubj', 'data', 'glob', 'node', 'nets', '-v7.3'); 

clear ii jj indx tmat tlab ulab dlab

%% mask the networks

% grab the size and set the threshold of edges to keep
nrep = size(nets, 3);
prop = 0.50;

% convert to logical
lnet = nets > 0;

% get a sum of how many edges exist across subjects
smat = sum(lnet, 3);

% define the edges w/ above proprotion threshold
lmat = smat >= (nrep*prop);

%plotAdjacencyMatrix(smat);

% preallocate output
omat = zeros(size(nets));

% loop and threshold each subject
for ii = 1:size(nets, 3)
    omat(:, :, ii) = squeeze(nets(:, :, ii)) .* lmat;
end

% clear out any impossible values
omat(isinf(omat)) = 0;
omat(isnan(omat)) = 0;

clear drop lmat lnet 

%plotAdjacencyMatrix(mean(omat, 3));

%% create the network measures
% % network measures to write out
% % edges: raw, communicability_wei
% % nodes: degrees_und, betweenness_wei

% grab indices for the upper diagonal
upd = logical(triu(ones(376), 1));

% % raw - edges

% preallocate edge measure
raw = nan(size(omat, 3), 70500);

% for every subject, put raw data in cca matrix
for ii = 1:size(omat, 3)
    tmat = squeeze(omat(:,:,ii));
    raw(ii, :) = tmat(upd);
end

clear ii tmat

% % communicability - edges

% preallocate edge measure
com = nan(size(omat, 3), 70500);

% preallocate communicability matrix
cmat = nan(size(omat));

% estimate communicability
for ii = 1:size(omat, 3)
    tmat = squeeze(omat(:,:,ii));
    cmat(:, :, ii) = communicability_wei(tmat);
end

clear ii tmat

% for every subject, put communicability data in cca matrix
for ii = 1:size(cmat, 3)
    tmat = squeeze(cmat(:,:,ii));
    com(ii, :) = tmat(upd);
end

clear ii tmat

% % degree - nodes

% preallocate node measure
deg = nan(size(omat, 3), 376);

% estimate node degree
for ii = 1:size(omat, 3)
    tmat = squeeze(omat(:,:,ii));
    deg(ii, :) = degrees_und(tmat);
end

clear ii tmat

% % betweeness centrality - node

% preallocate node measure
btw = nan(size(omat, 3), 376);

% estimate node degree
for ii = 1:size(omat, 3)
    tmat = squeeze(omat(:,:,ii));
    btw(ii, :) = betweenness_wei(tmat);
end

clear ii tmat

% % also get Average Controllability, Modal Controlloability, and Boundary Controllability
% x = averMeasTieredValsDirected(net); ave_control()
% y = moduMeasTieredVals(net); modal_control()
% z = getTieredControlStatesSH({net}, [.9], [1]); % takes all networks, need to determine threshold (.9)

% controllability measures

% preallocate node measures
ave_control_rank = nan(size(omat, 3), 376);
ave_control_vals = nan(size(omat, 3), 376);
mod_control_rank = nan(size(omat, 3), 376);
mod_control_vals = nan(size(omat, 3), 376);
bnd_control_rank = nan(size(omat, 3), 376);

% estimate controllability measures
for ii = 1:size(omat, 3)
    tmat = squeeze(omat(:,:,ii));
    ave_control_rank(ii, :) = averMeasTieredValsDirected(tmat);
    ave_control_vals(ii, :) = ave_control(tmat);
    mod_control_rank(ii, :) = moduMeasTieredVals(tmat);
    mod_control_vals(ii, :) = modal_control(tmat);
    tbd = getTieredControlStatesSH({tmat}, 1, 1);
    bnd_control_rank(ii, :) = tbd{1};
end

clear ii tmat tbd

%save('camcan_594_net_control.mat', 'ave_control_vals', 'ave_control_rank', 'mod_control_vals', 'mod_control_rank', 'bnd_control_rank');

%% call all network measures from FiNE

glob = cell(size(omat, 3), 1);
node = cell(size(omat, 3), 1);
netw = cell(size(omat, 3), 1);

% estimate controllability measures
for ii = 1:size(omat, 3)
    tmat = squeeze(omat(:,:,ii));
    [ glob{ii}, node{ii}, netw{ii} ] = fnNetworkStats(tmat);
end

%% save data as plain text

% write out as plain text
dlmwrite('camcan_594_net_raw.csv', raw, 'delimiter', ',');
dlmwrite('camcan_594_net_com.csv', com, 'delimiter', ',');
dlmwrite('camcan_594_net_deg.csv', deg, 'delimiter', ',');
dlmwrite('camcan_594_net_btw.csv', btw, 'delimiter', ',');

%% read in demographic data

% load demographic data
ddemo = dlmread('camcan_raw_data.csv');

% load subject ID's for demo
filePh = fopen('camcan_raw_id.csv', 'r');
C = textscan(filePh, '%s\n'); dsubj = C{1};
fclose(filePh); clear filePh C

% load demo variable names
filePh = fopen('camcan_raw_var_names.csv', 'r');
C = textscan(filePh, '%s\n'); varsNames = C{1}';
fclose(filePh); clear filePh C

%% read in regressors

% load demographic data
rdemo = dlmread('camcan_reg_data.csv');

% load subject ID's for demo
filePh = fopen('camcan_reg_id.csv', 'r');
C = textscan(filePh, '%s\n'); rsubj = C{1};
fclose(filePh); clear filePh C

% load demo variable names
filePh = fopen('camcan_reg_var_names.csv', 'r');
C = textscan(filePh, '%s\n'); confNames = C{1}';
fclose(filePh); clear filePh C

%% subset data based on existing subjects

% line up the subject ID's
[ subj, dind, nind ] = intersect(dsubj, subjs); % warning but works?
    
% subset the behavior / nuissance data so they match
vars = ddemo(dind, :);
varsQconf = rdemo(dind, :);

% subset the network values (probably not necessary)
%dens = ndens(nind, :);
raw = raw(nind, :);
com = com(nind, :);
deg = deg(nind, :);
btw = btw(nind, :);

% remove bad values
raw(isnan(raw)) = 0;
com(isnan(com)) = 0;
deg(isnan(deg)) = 0;
btw(isnan(btw)) = 0;

% remove bad values
raw(isinf(raw)) = 0;
com(isinf(com)) = 0;
deg(isinf(deg)) = 0;
btw(isinf(btw)) = 0;

% decided to drop: parent alive/age/age_died
vars(:, 39:44) = [];
varsNames(39:44) = [];

% extract mmse summary score
mmse = vars(:, 1);

% extract age
age = varsQconf(:, 1); 
varsQconf(:, 1) = [];
confNames(1) = [];

% create color map for age, color each plotted point by age
cmap = parula(size(18:88, 2));
nage = age - 17;

% pull gender for 2 cca factor
male = varsQconf(:, 2);

% drop down to most relevant
varsQconf = varsQconf(:, 1:15);
confNames = confNames(1:15);

% read in parcellation labels
filePh = fopen('hcp_mmp_key.txt', 'r');
C = textscan(filePh, '%d, %s\n');
netIndex = C{1};
netNames = C{2};
fclose(filePh); clear filePh C

%% write out data workspace for analysis

save('canoncorr_analysis_full_data.mat', 'omat', 'vars', 'varsQconf', ...
     'raw', 'com', 'deg', 'btw', 'subj', 'varsNames', 'confNames', ...
     'male', 'mmse', 'age', 'nage', 'cmap', 'netIndex', 'netNames', ...
     'glob', 'node', 'netw', 'ave_control_vals', 'ave_control_rank', ...
     'mod_control_vals', 'mod_control_rank', 'bnd_control_rank');
