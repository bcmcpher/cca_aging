%% analyses for the paper

%% load the data

% load the edge measure
load('canoncorr_analysis_full_data.mat', 'deg');

% load the rest of the behavior and labels
load('canoncorr_analysis_full_data.mat', 'age', 'vars', 'varsQconf', ...
     'netNames', 'varsNames', 'confNames');

% replace varsNames w/ corresponding labels
load('camcan_vars_labels.mat', 'varsLabels');

% get the indices of intersection between the datasets
vidx = contains(varsLabels(:, 1), varsNames);

% grab plaintext names in the right order
varsLabel = varsLabels(vidx, 2);

clear vidx

% load the yeo labels
load('~/hcp_mmp_vol/working/yeoLabs.mat');

% load the average network for RC estimate
load('camcan_average_network.mat');

% path to the output
outpath = '/N/slate/bcmcpher/camcan/cb_data';

%% The figures in order

% 1b - an example input network
% copied over CC110045, the subject used

%% 1e - MMSE binned by decade

fig1e = plotVarByDeciles(vars(:,1), age, 'linear', 'MMSE');
writematrix(fig1e.data, fullfile(outpath, 'fig1e.csv'), 'delimiter', ',');

%% 1f - Max node degree binned by decade

mdeg = nan(594, 1);
for ii = 1:594
    mdeg(ii) = max(deg(ii,:));
end

fig1f = plotVarByDeciles(mdeg, age, 'linear', 'MMSE');
writematrix(fig1f.data, fullfile(outpath, 'fig1f.csv'), 'delimiter', ',');

%% 2a - CCA axes

fig2a = [ cca.dat1.factor(:,1), cca.dat1.factor_se(:,1), cca.dat2.factor(:,1), cca.dat2.factor_se(:,1), age ];
writematrix(fig2a, fullfile(outpath, 'fig2a.csv'), 'delimiter', ',');

%% 2b - r_{age} - Full and CV 

% this call produces the values - permutations make in slighlty different
% run to run, so I'm just writing out the values from the manuscript

% fig2b = nan(2,2);
% [ fig2b(1,1), fig2b(1,2) ] = ccaLinRegCorr(cca, 1, age, 10000, true);
% [ fig2b(2,1), fig2b(2,2) ] = ccaLinRegCorr(cca, 1, age, 10000);

% values from the paper
fig2b = [ 0.6270, 0.0220; 0.6110, 0.0230 ];
writematrix(fig2b, fullfile(outpath, 'fig2b.csv'), 'delimiter', ',');

%% 2c - CA1 - Full and CV

fig2c = [ cca.cca.hocorrs(1), cca.cca.hocorrs_se(1); cca.full.grotR(1), 0 ];
writematrix(fig2c, fullfile(outpath, 'fig2c.csv'), 'delimiter', ',');

%% 3a - brian loadings

% pull the data
load3a = cca.dat1.loading(:,1);
ldse3a = cca.dat1.loading_se(:,1);
pval3a = nan(376,1);

% pull pval
for ii = 1:376
    pval3a(ii) = sum(cca.dat1.loading_null(ii,1,:) > cca.dat1.loading(ii,1))/10000; 
end

% sort like plot did
[ ~, l3ai ] = sort(load3a, 'descend');

% sort the names
fig3aNames = dat.dat1.names(l3ai);

% build the output
fig3a = [ cca.dat1.loading(l3ai), cca.dat1.loading_se(l3ai), pval3a(l3ai) ];

fileID = fopen('cb_data/fig3a.csv','w');
for ii = 1:376
    fprintf(fileID,'%s,%d,%d,%d\n',fig3aNames{ii}, fig3a(ii,1), fig3a(ii,2), fig3a(ii,3));
end
fclose(fileID);

%% 3b - behavior loadings

% pull the data
load3b = cca.dat2.loading(:,1);
ldse3b = cca.dat2.loading_se(:,1);
pval3b = nan(334,1);

% pull pval
for ii = 1:334
    pval3b(ii) = sum(cca.dat2.loading_null(ii,1,:) > cca.dat2.loading(ii,1))/10000; 
end

% sort like plot did
[ ~, l3bi ] = sort(load3b, 'descend');

% sort the names
fig3bNames = dat.dat2.names(l3bi);

% build the output
fig3b = [ cca.dat2.loading(l3bi), cca.dat2.loading_se(l3bi), pval3b(l3bi) ];

fileID = fopen('cb_data/fig3b.csv','w');
for ii = 1:334
    fprintf(fileID,'%s,%d,%d,%d\n',fig3bNames{ii}, fig3b(ii,1), fig3b(ii,2), fig3b(ii,3));
end
fclose(fileID);

%% 3c - brain networks

fig3c = ccaModuleContribution(cca, 'brain', 1, yeoLabs.yeo7, yeoLabs.yeo7Names);

fileID = fopen('cb_data/fig3c.csv','w');
for ii = 1:10
    fprintf(fileID,'%s,%d,%d,%d,%d,%d,%d,%d\n',yeoLabs.yeo7Names{ii}, fig3c(ii,1), fig3c(ii,2), fig3c(ii,3), fig3c(ii,4), fig3c(ii,5), fig3c(ii,6), fig3c(ii,7));
end
fclose(fileID);

%% 3d - behavior domains

fig3d = ccaModuleContribution(cca, 'behavior', 1, si, S);

fileID = fopen('cb_data/fig3d.csv','w');
for ii = 1:7
    fprintf(fileID,'%s,%d,%d,%d,%d,%d,%d,%d\n', S{ii}, fig3d(ii,1), fig3d(ii,2), fig3d(ii,3), fig3d(ii,4), fig3d(ii,5), fig3d(ii,6), fig3d(ii,7));
end
fclose(fileID);

%% 4b - average +/- loadings in RC/periphery

ndeg = degrees_und(mat);

% normalize by minmax deg / str
zdeg = (ndeg-6)/(298-6);

% create a set of labels for non-RC (1) and RC (2) based on relative cutoff
% to reference
degi = ones(376,1);
degi(zdeg > 0.43) = 2;

fig4b = ccaModuleContribution(cca, 'brain', 1, degi, {'Periphery', 'RC'});

labs = {'Periphery', 'RC'};

fileID = fopen('cb_data/fig4d.csv','w');
for ii = 1:2
    fprintf(fileID,'%s,%d,%d,%d,%d,%d,%d,%d\n', labs{ii}, fig3d(ii,1), fig3d(ii,2), fig3d(ii,3), fig3d(ii,4), fig3d(ii,5), fig3d(ii,6), fig3d(ii,7));
end
fclose(fileID);

%% 5a - the module values

dmat = ccaDissimilarityMatrix(cca);
[ mdDat, ~, mdBtw ] = fnModuleDensity(dmat, [ yeoLabs.yeo7; ib+10 ], 'mean');

% get index of sorted behaviors
svar = regexprep(dat.dat2.names', '_.*', '');
svar = regexprep(svar, 'hint', 'comp');
svar{1} = 'comp'; % replace the first dumb label
[ S, ~, ib ] = unique(svar);
[ si, sv ] = sort(ib);

y7_axis = [ y7_brn; sv+376 ];

fig5a = fnModuleDensity(dmat, [ yeoLabs.yeo7; ib+10 ], 'mean');
writematrix(fig5a, fullfile(outpath, 'fig5a.csv'), 'delimiter', ',');

%% 5b - the off diagonal

fig5b = fig5a(11:end,1:10);
writematrix(fig5b, fullfile(outpath, 'fig5b.csv'), 'delimiter', ',');

%% 5c is the same values in 5b
