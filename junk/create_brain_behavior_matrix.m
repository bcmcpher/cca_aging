%% create data for batplot

% load data
load('ho_age_deg_data.mat', 'dissimilarity');

% pull the module dissimilarity
x = dissimilarity.mdDat;

% round and flip so more similar == higher number
round((1-x(11:end,1:10))*1000)

% transpose for off diagonal
round((1-x(11:end,1:10))*1000)'

% the total size N for determining the blank space in the plot
sum(sum(round((1-x(11:end,1:10))*1000)))
