%in the "GridPlot.m", the field was actually number of in field bins, here
%let's calculate the real fields
cd ..\
load('nFieldReal.mat');
cd('gridCellsNewCueThresh');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\indicesAllNewCueThresh.mat');
i=indicesAllNewCueThresh.GridIdx;
nFieldRealGrid=nFieldReal(i,:);
save('nFieldRealGrid.mat','nFieldRealGrid');
figure,plot(nanmean(nFieldRealGrid,1))
title('nFields grid');
saveas(gcf,'nFieldsRealGrid.fig');