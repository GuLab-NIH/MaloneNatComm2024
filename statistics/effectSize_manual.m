%%

clear; clc

% data1 = [];
% data2 = [];

data1 = {};
data2 = {};


%% Perform ANOVA

testName = 'two-way repeated measures ANOVA';

% nUnits = 'cells';
nUnits = 'cells from 11,4 mice';
% nUnits = 'FOV from 11,4 mice';
% nUnits = 'mice';
% nUnits = 'cells from 4,2 mice';
% nUnits = 'sessions from 5 mice';
% nUnits = 'FOV from 4 mice';


outStats = anovaEffectSize(data1,data2,testName,nUnits);
% outStats = anovaEffectSize(data1,data2,testName,nUnits,[1 1]);



%% Perform Pearson Correlation

testName = 'two-tailed Pearson linear correlation';

% nUnits = 'time points';
% nUnits = 'time points from 15 mice';
% nUnits = 'time points from 11 mice';
% nUnits = 'time points from 660 cells from 11 mice';
% nUnits = 'time points from 366 cells from 11 mice';
% nUnits = 'time points from 35 cells from 11 mice';
% nUnits = 'time points from 21 FOV from 11 mice';
% nUnits = 'time points from 8 FOV from 4 mice';

% nUnits = 'time points from 4 mice';
% nUnits = 'time points from 23 cells from 4 mice';
% nUnits = 'time points from 527 cells from 4 mice';
% nUnits = 'time points from 284 cells from 4 mice';
% nUnits = 'time points from 23 cells from 4 mice';
% nUnits = 'time points from 9 FOV from 4 mice';

% nUnits = 'sessions from 15 mice';

% nUnits = 'mice';
nUnits = 'cell pairs';
% nUnits = 'sessions from 6 mice';


outStats = corrEffectSize(data1,data2,testName,nUnits);
% outStats = corrEffectSize(data1',data2',testName,nUnits);


%% Perform Student's t-test

% testName = 'two-tailed unpaired Students t-test';
testName = 'two-tailed unpaired Students t-test, no correction for multiple comparisons';
% testName = 'two-tailed paired Students t-test';
% testName = 'two-tailed paired Students t-test, no correction for multiple comparisons';
% testName = 'two-tailed unpaired Students t-test, Bonferroni-Holm correction for mutiple comparsions';
% testName = 'two-tailed paired Students t-test, Bonferroni-Holm correction for mutiple comparsions';
% testName = 'two-sided Mann-Whitney U test, no correction for multiple comparisons';
% testName = 'one-tailed unpaired Students t-test, no correction for multiple comparisons';


% nUnits = 'cells';
% nUnits = 'cells from 11 mice';
% nUnits = 'cells from 4 mice';
% nUnits = 'cells from 11,4 mice';
% nUnits = 'FOV from 11,4 mice';
% nUnits = 'FOV from 4 mice';
% nUnits = 'mice';
nUnits = 'bins';
% nUnits = 'trials';
% nUnits = 'cell pairs from 11 mice';
% nUnits = 'cell pairs from 4 mice';
% nUnits = 'cells from 6,4 mice';
% nUnits = 'slices from 16,10 mice';
% nUnits = 'sessions from 8,5 mice';


% pair = 1;
pair = 0;

MC = 0;

limitP = 0;
% limitP = 1;

outStats = ttestEffectSize(data1,data2,testName,nUnits,pair,MC,limitP);
% outStats = ttestEffectSize(data1',data2',testName,nUnits,pair,MC,limitP);
% outStats = ttestEffectSize(data1,data2,testName,nUnits,pair,MC,limitP,'ranksum');


%% Perform pairwise Student's t-test

% testName = 'two-tailed unpaired Students t-test';
testName = 'two-tailed unpaired Students t-test, no correction for multiple comparisons';
% testName = 'two-tailed paired Students t-test';
% testName = 'two-tailed paired Students t-test, no correction for multiple comparisons';
% testName = 'two-tailed unpaired Students t-test, Bonferroni-Holm correction for mutiple comparsions';
% testName = 'two-tailed paired Students t-test, Bonferroni-Holm correction for mutiple comparsions';

% nUnits = 'cells';
% nUnits = 'cells from 11 mice';
% nUnits = 'cells from 4 mice';
% nUnits = 'mice';
% nUnits = 'sessions from 8 mice';
% nUnits = 'sessions from 5 mice';
% nUnits = 'FOV from 11,4 mice';
% nUnits = 'FOV from 15 mice';
nUnits = 'bins from 10 sessions';


% pair = 1;
pair = 0;

% MC = 1;
MC = 0;

testIdxs = [1 2;1 3;2 3];
% testIdxs = [1 4;2 5;3 6;2 3;5 6];
% testIdxs = [1 2;3 4;5 6;7 8];
% testIdxs = [];

outStats = ttestAllPairsEffectSize(data1,testName,nUnits,pair,MC,testIdxs);


%% Get n's


