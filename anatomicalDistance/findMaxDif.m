%%

clear; close all; clc

load('distsFOVFx_bad_global.mat')

%%


dDif = distsCom{2}./distsCom{1};


difMean = mean(dDif,1);
[difMx,difIdx] = max(dDif(1:end,:));

difIdx = difIdx;

