%% loop over different numbers of CV iterations to pick the number of repeat CVs

% all the relevant objects need to be loaded already

% pick the number of repeated CVs
nCV = 250000;
nRp = 5;

% preallocate the model fits
Rdat = cell(nRp, 1);
Rcca = cell(nRp, 1);
Rdmat = cell(nRp, 1);
RmdDat = cell(nRp, 1);

% repeat the CV 10 times
for ii = 1:nRp
    
    % run the analysis repeatedly
    [ Rdat{ii}, Rcca{ii} ] = ccaMapFullAnalysis(deg, vars, varsQconf, ...
                                                netNames, varsNames, confNames, ...
                                                varsLabel, 38, 40, 0, 5, nCV);
    
    % build the dissimilarity matrix
    Rdmat{ii} = ccaDissimilarityMatrix(Rcca{ii});
    
    % build the modules
    RmdDat{ii} = fnModuleDensity(Rdmat{ii}, [ yeoLabs.yeo7; ib+10 ], 'mean');
    
end

clear ii

fh = figure('Position', [ 375, 450 1300 350 ]);
for ii = 1:nRp
    
    subplot(1, nRp, ii); 
    imagesc(RmdDat{ii});
    axis square; axis equal; axis tight; caxis([ 0 1 ]);
    set(gca, 'XTick', [], 'YTick', [], ...
        'XLim', [ 0.5 10.5 ], 'YLim', [ 10.5 17.5 ]); 
    
end
