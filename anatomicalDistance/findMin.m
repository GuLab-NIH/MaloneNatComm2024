%%

clear all; close all; clc

d = dir('data_*');
nFile = length(d);

minGlob = NaN;

%%

for f = 1:nFile
    load(d(f).name)
    
    minCur = min(size(allCen.common,1),size(allCen.unique,1));
    
    minGlob = min(minGlob,minCur,'omitnan');
end

save('minGlob.mat','minGlob')

