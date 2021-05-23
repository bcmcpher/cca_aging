%% estimate the rich club of the average network

% load the averaged data
load('camcan_average_network.mat');

% load the labels for the brain regions
load('canoncorr_analysis_full_data.mat', 'netIndex', 'netNames');
netNames = strrep(netNames, '???', 'CC'); % correct bad label

% load yeo label assignments for predefined communtities
load('yeoLabs.mat');

%% step through van den Heuval 2011

% including the (node-specific) degree, clustering coefficient, 
% characteristic path length, betweennesscentrality, normalized clustering 
% coefficient and normalized path length (both normalized relative to a set 
% of 100 comparable random graphs), global efficiency, assortativity, and modularity. 

% true observations of all these parameters
deg = degrees_und(mat);
str = strengths_und(mat);
ccf = clustering_coef_wu(mat);
cpl = charpath(mat);
btc = betweenness_wei(mat);
bte = edge_betweenness_wei(mat);
glb = efficiency_wei(mat);
ast = assortativity_wei(mat, 0);
mod = modularity(mat);

% these are the true vales / average values observed from 100 random permutations
% probably not important...
%ncc = 
%npl = 

%% create the rich club coefficient / k-core on binarized graph

% binarize intput
bat = weight_conversion(mat, 'binarize');

% is this the cutoff they used - 1 SD above mean? (used for Fig 2)
mean(deg) + std(deg)

% binarized rich club curve
[ rcb, rcn, rce ] = rich_club_bu(bat);

% run k-core on binarized network at arbitrary degree (this keeps every node)
[ kmat, ksz, kpo, kpl ] = kcore_bu(bat);
% how do I prune this down to a rich club?

%% create the rich club coefficient / k-core on binarized graph

% weighted rich club curve
rcw = rich_club_wu(mat);

% preallocate null rich club vector
niter = 100;
nrc = nan(niter, size(rcw,2));

% for a fixed number of null networks
for iter = 1:niter
    
    % create a random, preserved network
    nmat = null_model_und_sign(mat);
    %nmat = randomio_und(mat, 5);
    
    % estimate the rich clube on the null
    nrc(iter,:) = rich_club_wu(nmat);
    
end

% estimate the corrected rich-club
rcc = rcw ./ mean(nrc, 'omitnan');

% Figure 3
figure; hold on;
plot(1:length(rcw), rcw, '-o', 'color', 'black', 'MarkerFaceColor', 'white', 'MarkerEdgeColor', 'black', 'LineWidth', 1, 'MarkerSize', 2);
plot(1:length(nrc), mean(nrc), '-o', 'color', [.7 .7 .7], 'MarkerFaceColor', 'white', 'MarkerEdgeColor', [.7 .7 .7], 'LineWidth', 1, 'MarkerSize', 2);
plot(1:length(rcc), rcc, '-o', 'color', 'red', 'MarkerFaceColor', 'white', 'MarkerEdgeColor', 'red', 'LineWidth', 1, 'MarkerSize', 2);
title('Corrected Rich-Club Coefficient for Average Network');
xlabel('k (degree)');
ylabel('RCC \phi(k)');
legend('Observed', 'Null', 'Corrected');
hold off;

% use this to identify rich club somehow...?

%% use s-core

[ scr, ssz ] = score_wu(mat, 0.55712);
% this is very sensitive and is either very large or 0

%% just use normalized degree / strength to threshold a rich club?

% normalize by minmax deg / str
z = (deg-6)/(298-6);
%z = (str-0.0930)/(2.9555-0.0930);

% create a set of labels for non-RC (1) and RC (2)
zi = ones(376,1);
zi(z > 0.43) = 2;

% search through the threshold to find proportionally the same number
% of regions in the RC
netNames(z > 0.43)
%     {'lh.CC'            }
%     {'lh.4'             }
%     {'lh.3b'            }
%     {'lh.RSC'           }
%     {'lh.POS2'          }
%     {'lh.PCV'           }
%     {'lh.5mv'           }
%     {'lh.23c'           }
%     {'lh.7Am'           }
%     {'lh.2'             }
%     {'lh.3a'            }
%     {'lh.6mp'           }
%     {'lh.8Av'           }
%     {'lh.8C'            }
%     {'lh.6a'            }
%     {'lh.TE1p'          }
%     {'lh.PF'            }
%     {'lh.PFm'           }
%     {'lh.PGi'           }
%     {'lh.PGs'           }
%     {'lh.PoI1'          }
%     {'rh.CC'            }
%     {'rh.4'             }
%     {'rh.3b'            }
%     {'rh.RSC'           }
%     {'rh.POS2'          }
%     {'rh.PCV'           }
%     {'rh.5mv'           }
%     {'rh.23c'           }
%     {'rh.7Am'           }
%     {'rh.2'             }
%     {'rh.3a'            }
%     {'rh.8Av'           }
%     {'rh.8Ad'           }
%     {'rh.8C'            }
%     {'rh.6r'            }
%     {'rh.6a'            }
%     {'rh.RI'            }
%     {'rh.AIP'           }
%     {'rh.PF'            }
%     {'rh.PFm'           }
%     {'rh.PGi'           }
%     {'rh.PGs'           }
%     {'rh.PoI1'          }
%     {'lh.thalamus'      }
%     {'lh.caudate'       }
%     {'lh.putamen'       }
%     {'lh.globuspallidus'}
%     {'lh.hippocampus'   }
%     {'rh.thalamus'      }
%     {'rh.caudate'       }
%     {'rh.putamen'       }
%     {'rh.globuspallidus'}
%     {'rh.hippocampus'   }

% plot the point plot labeling RC vs. periphery
ccaPlotRankedTrendsRC(dat,cca,age,'brain','load',1,'points',-1,zi);
% the RC nodes go quite a ways down the figure...

% get the average bins for RC
fig4b = ccaModuleContribution(cca, 'brain', 1, zi, {'Periphery', 'RC'});

% spearman rank correlation
corr(deg', cca.dat1.loading(:,1),'type','spearman')

corr(tdeg(zi == 1), cca.dat1.loading(zi == 1,1),'type','spearman')
corr(tdeg(zi == 2), cca.dat1.loading(zi == 2,1),'type','spearman')

% permutation test for correlation
nrep = 10000;
sprc = nan(nrep,1);
sprn = nan(nrep,1);
sppp = nan(nrep,2);
sprr = nan(nrep,2);
tdeg = deg';
pdeg = tdeg(zi == 1);
rdeg = tdeg(zi == 2);
pcca = cca.dat1.loading(zi == 1,1);
rcca = cca.dat1.loading(zi == 2,1);

for ii = 1:nrep
    
    % pull random sort
    ridx = randsample(1:size(tdeg,1), size(tdeg,1), true);
    
    % compute the resampled mean / se / pval for full corr b/w rc and cca
    sprc(ii,1) = corr(tdeg(ridx), cca.dat1.loading(ridx,1),'type','spearman');
    sprn(ii,1) = corr(tdeg, cca.dat1.loading(ridx,1),'type','spearman');
    
    % random sort periphery (p) rich club (r) separately
    pidx = randsample(1:size(pdeg,1), size(pdeg,1), true);
    ridx = randsample(1:size(rdeg,1), size(rdeg,1), true);

    % periphery to cca only
    sppp(ii,1) = corr(pdeg(pidx), pcca(pidx),'type','spearman');
    sppp(ii,2) = corr(pdeg, pcca(pidx),'type','spearman');
    
    % rc to cca only
    sprr(ii,1) = corr(rdeg(ridx), rcca(ridx),'type','spearman');
    sprr(ii,2) = corr(rdeg, rcca(ridx),'type','spearman');

end

% answer
[ num2str(mean(sprc)) ' +/- ' num2str(std(sprc)) ' (p = ' num2str(sum(sprn > corr(deg',cca.dat1.loading(:,1),'type','spearman'))/nrep) ')' ]
[ num2str(mean(sppp(:,1))) ' +/- ' num2str(std(sppp(:,1))) ' (p = ' num2str(sum(sppp(:,2) > corr(pdeg,cca.dat1.loading(zi==1,1),'type','spearman'))/nrep) ')' ]
[ num2str(mean(sprr(:,1))) ' +/- ' num2str(std(sprr(:,1))) ' (p = ' num2str(sum(sprr(:,2) > corr(rdeg,cca.dat1.loading(zi==2,1),'type','spearman'))/nrep) ')' ]

% test the difference between the average cca
tdif = nan(nrep,1);
for rep = 1:nrep

    % random sort periphery (p) rich club (r) separately
    pidx = randsample(1:size(pdeg,1), size(pdeg,1), true);
    ridx = randsample(1:size(rdeg,1), size(rdeg,1), true);

    % take the reampled difference between the rc / periphery
    tdif(rep) = mean(rcca(ridx)) - mean(pcca(pidx));

end

% compare the explicit null test to the observed difference
sum(tdif > (0.045 - 0.032)) / nrep

%% evaluate what the RC nodes do

% participation coefficient in Yeo networks
%pcf = participation_coef(mat, yeoLabs.yeo7);
pcf = participation_coef(mat, zi);
% should be the RC assignments...

% connector hubs
pcf > 0.50

% provincial hubs
pcf <= 0.50

% threshold of 0.40 works...

%% plot the scatter plot of loading x participation coefficient, assigning RC

figure; hold on;

% plot error bars first
for ii = 1:length(pcf)
    ct = cca.dat1.loading(ii,1);
    sd = 2*cca.dat1.loading_se(ii,1);
    plot([ pcf(ii) pcf(ii) ], [ ct-sd ct+sd ], 'color', 'black');
end
clear ct sd

% plot the periphery nodes
for ii = 1:length(pcf)
    if zi(ii) == 1
        plot(pcf(ii), cca.dat1.loading(ii,1), 'o', 'MarkerFaceColor', [ 0.7 0.7 0.7 ], ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'MarkerSize', 6);
    end
end

% plot the RC nodes on top
for ii = 1:length(pcf)
    if zi(ii) == 2
        plot(pcf(ii), cca.dat1.loading(ii,1), 'o', 'MarkerFaceColor', [ 0.1 0.2 0.7 ], ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'MarkerSize', 6);
    end
end

title('Rich Club Participation x CCA Loading');
xlabel('Rich Club Participation');
ylabel('CCA Loading');
hold off;

%% plot the dots and lines

% load points to plot
centers = dlmread('data/CC110045_centers.csv');

% find points to keep
rcenter = centers(zi == 2,:);
dcenter = centers(zi ~= 2,:);

% pull the subset of node properties
rdeg = deg(zi == 2);
rstr = str(zi == 2);
rcca = cca.dat1.loading(zi == 2,1);
ncca = (rcca-min(rcca))/((max(rcca)-min(rcca)))+0.001;

ddeg = deg(zi ~= 2);
dcca = cca.dat1.loading(zi ~= 2,1);
ddcc = (dcca-min(dcca))/((max(dcca)-min(dcca)))+0.001;

% subset the graph
rmat = mat(zi == 2, zi == 2);
zmat = dmat(zi == 2, zi == 2);

% normalize graph edge weights beween 0 and 1
rmat = weight_conversion(rmat, 'normalize');
zmat = weight_conversion(zmat, 'normalize');

% build xy indices for plotting edges
[ xy ] = nchoosek(1:size(rmat,1), 2);

% build the default sphere
NumSphFaces = 15;
[ SX, SY, SZ ] = sphere(NumSphFaces);

% node size scale for plotting
nsz = 5;

figure; hold on;

% plot the periphery nodes
for node = 1:size(dcenter)
    
    %msz = (ddeg(node)/max(ddeg))*5; % scale size by node degree
    msz = ddcc(node)*nsz; % scale by cca contribution

    % scale / align the canonical sphere to the node center, scale size by input vector
    surf(SX*msz + dcenter(node, 1), SY*msz + dcenter(node, 2), SZ*msz + dcenter(node, 3), ...
         'EdgeColor', 'none', 'FaceColor', [ 0.7 0.7 0.7 ], 'FaceAlpha', 0.50);
    
end

% plot the edges between the points
for edge = 1:size(xy,1)
    
    % grab the indices of the edge
    xi = xy(edge,1);
    yi = xy(edge,2);
    
    % plot the edge if it exists above a minimum threshold
    % rmat > 0.05, *10 - zmat > 0.75, *3
    if rmat(xi, yi) > 0.05
        plot3([ rcenter(xi, 1) rcenter(yi, 1) ], [ rcenter(xi, 2) rcenter(yi, 2) ], [ rcenter(xi, 3) rcenter(yi, 3) ], ...
            'color', 'black', 'LineWidth', rmat(xi,yi)*10);
    end
        
end

% plot each node as a point scaled by degree over the lines
for node = 1:size(rcenter)
    
    %msz = (rdeg(node)/max(rdeg))*5; % scale size by node degree
    msz = ncca(node)*nsz; % scale size by cca contribution
        
    % scale / align the canonical sphere to the node center, scale size by input vector
    surf(SX*msz + rcenter(node, 1), SY*msz + rcenter(node, 2), SZ*msz + rcenter(node, 3), ...
         'EdgeColor', 'none', 'FaceColor', [ 0.1 0.2 0.7 ], 'FaceAlpha', 1.0);
    
end

% global plot configuration
title('Rich Club');
set(gca, 'xlim', [ -65 65 ], 'ylim', [ -90 70 ], 'zlim', [ -40 70 ]);
axis equal

% set lighting
material metal
lighting phong

% hide axes
set(gca, 'Visible', 'off');
hold off;

% save the figure out
view(0, 90);
camlight('right');
print('figs/20210111/cca_bs_axial.eps', '-painters', '-depsc', '-r300');

view(90,0); % right sagittal
camlight('right');
print('figs/20210111/cca_bs_rsag.eps', '-painters', '-depsc', '-r300');

view(-90,0); % left sagittal
camlight('right');
print('figs/20210111/cca_bs_lsag.eps', '-painters', '-depsc', '-r300');

view(0,0); % coronal
camlight('right');
print('figs/20210111/cca_bs_coronal.eps', '-painters', '-depsc', '-r300');

%% render loadings on the surface




