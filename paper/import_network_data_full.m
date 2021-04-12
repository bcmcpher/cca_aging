%% import and fill in empty nodes based on labels

% index data
files = dir('data/*density.csv');
label = dir('data/*labels_dilated.nii.gz');

% get file names
mats = {files.name}';
labs = {label.name}';

nsubj = strrep(mats, '_density.csv', '');
filePh = fopen('ccid.txt', 'w');
fprintf(filePh, '%s\n', nsubj{:});
fclose(filePh);

% preallocate network data
data = cell(size(mats));
prob = cell(size(mats));

% preallocate network summary for alternate
degr = nan(size(mats, 1), 376);
btwc = nan(size(mats, 1), 376);
pgrk = nan(size(mats, 1), 376);

% loop over and load / fix networks
for ii = 1:size(mats, 1)
    
    % load temporary data
    tmat = dlmread([ 'data/' mats{ii} ]);
    tlab = niftiRead([ 'data/' labs{ii} ]);

    % grab unique index of labels
    ulab = unique(tlab.data(:));
        
    disp(['N labels: ' num2str(length(ulab) - 1) ' - Label File ID: ' labs{ii} ]);
        
    % fix all bad - catch runs, too
    while ~all(diff(ulab) == 1)
        
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
        
    end
    
    % store output
    data{ii} = tmat;
    prob{ii} = sort(prob{ii});
    
    % compute network statistics for alternate CCA
    %degr(ii, :) = degrees_und(tmat);
    %btwc(ii, :) = betweenness_wei(tmat);
    %pgrk(ii, :) = pagerank_centrality(tmat, 0.85);
    
    % compute the full set
    %[ glob{ii}, node{ii}, nets{ii} ] = fnNetworkStats(tmat);
    
end

save('camcan_netstats.mat', 'nsubj', 'data', 'glob', 'node', 'nets', '-v7.3'); 

clear ii jj indx tmat tlab ulab dlab

%% linearize the upper diagonal 

% grab indices for the upper diagonal
upd = logical(triu(ones(376), 1));

% preallocate output
out = nan(size(data, 1), 70500);

% for every subject, put data in matrix
for ii = 1:size(data, 1)
    out(ii, :) = data{ii}(upd);
end

clear ii

% write out as plain text
%dlmwrite('camcan_593_density.csv', out, 'delimiter', ',');

% catch network data in a useful structure
ndens = out;

clear out upd filePh labs label data mats files

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
[ subj, dind, nind ] = intersect(dsubj, nsubj); % warning but works?
    
% subset the behavior / nuissance data so they match
vars = ddemo(dind, :);
varsQconf = rdemo(dind, :);

% subset the networks (probably not necessary)
dens = ndens(nind, :);

%% write out data workspace for analysis

save('canoncorr_analysis_full_data.mat', 'vars', 'varsQconf', 'dens', 'subj', 'varsNames', 'confNames', 'degr', 'btwc', 'pgrk');
