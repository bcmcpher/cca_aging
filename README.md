# CCA Aging

The analysis code for performing the analyses in [this paper](https://www.nature.com/articles/s42003-021-02451-0) distributed as a toolbox of MATLAB functions for reproducibility.

McPherson, B. C., & Pestilli, F. (2021). A single mode of population covariation associates brain networks structure and behavior and predicts individual subjects’ age. Communications biology, 4(1), 1-16.


## Basic usage

The domain datasets (domain1, domain2, confounds, holdout) need to have correspondence between the rows (same subjects in the same order) of observations for the analysis to work.

```
% add tool
addpath(genpath('/path/to/repo/cca_aging'));

% estimate the cross-validated CCA
[ dat, cca ] = ccaFullKAnalysis(domain1, domain2, confounds, ...
                                dom1Names, dom2Names, confNames, labels, ...
                                25, 25, 100, 5, 15000, 'median');

% estimate the fit of a holdout variable to the primary canonical axis
[ R, S, pval ] = ccaLinRegCorr(cca, 1, holdout, 1000);

%% example plots of panels from the paper

% Figure 2a - Main Finding
ccaPlotAxisCon(cca, 1, holdout, parula(88), false, true);

% Figure 3a/b
ccaPlotRankedTrends(dat, cca, holdout, 'brain', 'load', 1, 'lines', 30);
ccaPlotRankedTrends(dat, cca, holdout, 'behavior', 'load', 1, 'lines', 30);

```
