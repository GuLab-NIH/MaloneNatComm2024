
%To make the first day in new env and last day in old env consistent with
%all later days, we use M2 data for all
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\nCellsAllPerM2.mat');
GM=nCellsAllPerM2;

figure
subplot(141)
pie(mean(GM,1));
colormap(summer);

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\cellTypesNewAllFractinEachDayAllFOVs.mat');
G=cellTypesNewAllFractinEachDayAllFOVs;
idx=[1 4 5 6 7 8 9 10 11 12 13 14];
%since the two days (2 and 3) were removed, re normalize percentile 
A=G{1};
ANormFold=1-sum(A(:,2:3),2);
AN=[];
for n=1:size(A,1);
    AN(n,:)=A(n,idx)/ANormFold(n);
end

subplot(142)
pie(nanmean(AN));
colormap(summer);

A=G{2};
ANormFold=1-sum(A(:,2:3),2);
AN=[];
for n=1:size(A,1);
    AN(n,:)=A(n,idx)/ANormFold(n);
end

subplot(143)
pie(nanmean(AN));
colormap(summer);

A=G{3};
ANormFold=1-sum(A(:,2:3),2);
AN=[];
for n=1:size(A,1);
    AN(n,:)=A(n,idx)/ANormFold(n);
end

subplot(144)
pie(nanmean(AN));
colormap(summer);

saveas(gcf,'goodLearnersAllGridCue.fig');
%%

%To make the first day in new env and last day in old env consistent with
%all later days, we use M2 data for all
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\commonUniqueCells\nCellsAllPerM2.mat');
BM=nCellsAllPerM2;

figure
subplot(151)
pie(mean(BM,1));
colormap(spring);

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\commonUniqueCells\cellTypesNewAllFractinEachDayAllFOVs.mat');
B=cellTypesNewAllFractinEachDayAllFOVs;
idx=[1 4 5 6 7 8 9 10 11 12 13 14];
%since the two days (2 and 3) were removed, re normalize percentile 
A=B{1};
ANormFold=1-sum(A(:,2:3),2);
AN=[];
for n=1:size(A,1);
    AN(n,:)=A(n,idx)/ANormFold(n);
end

subplot(152)
pie(nanmean(AN));
colormap(spring);

A=B{2};
ANormFold=1-sum(A(:,2:3),2);
AN=[];
for n=1:size(A,1);
    AN(n,:)=A(n,idx)/ANormFold(n);
end

subplot(153)
pie(nanmean(AN));
colormap(spring);

A=B{3};
ANormFold=1-sum(A(:,2:3),2);
AN=[];
for n=1:size(A,1);
    AN(n,:)=A(n,idx)/ANormFold(n);
end

subplot(154)
pie(nanmean(AN));
colormap(spring);

saveas(gcf,'badLearnersAllGridCue.fig');

%% days tracked
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\daysAll.mat');

G=daysAll;
daysMeanG=[];
daysSEMG=[];
daysSTDG=[];
for n=1:length(daysAll);
    daysMeanG(n)=mean(daysAll{n});
    daysSEMG(n)=std(daysAll{n})/sqrt(length(daysAll{n}));
    daysSTDG(n)=std(daysAll{n});
end

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\commonUniqueCells\daysAll.mat');

B=daysAll;
daysMeanB=[];
daysSEMB=[];
daysSTDB=[];
for n=1:length(daysAll);
    daysMeanB(n)=mean(daysAll{n});
    daysSEMB(n)=std(daysAll{n})/sqrt(length(daysAll{n}));
    daysSTDB(n)=std(daysAll{n});
end
sig=[];
for n=1:length(G);
    [~,sig(n)]=ttest2(G{n},B{n});
end
figure,
subplot(141)
idx=[1 4 5 6:1:14];
plot([1:1:length(idx)],daysMeanG(idx),'g');
hold on
errorbar([1:1:length(idx)],daysMeanG(idx),daysSEMG(idx),'g','MarkerSize',15)
hold on
plot([1:1:length(idx)],daysMeanB(idx),'m');
hold on
errorbar([1:1:length(idx)],daysMeanB(idx),daysSEMB(idx),'m','MarkerSize',15)
xlim([0 13]);
ylim([0 11]);
title('all cells')

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\cellTypeDayTrackM.mat');
GM=cellTypeDayTrackM;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvbADMoreMice\commonUniqueCells\cellTypeDayTrackM.mat');
BM=cellTypeDayTrackM;

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\cellTypeDayTrackSEM.mat');
GSEM=cellTypeDayTrackSEM;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvbADMoreMice\commonUniqueCells\cellTypeDayTrackSEM.mat');
BSEM=cellTypeDayTrackSEM;


idx=[1 5 6 7:1:15];
n=1;
  subplot(142);
  plot([1:1:length(idx)],GM(n,idx),'g');
  hold on
    errorbar([1:1:12],GM(n,idx),GSEM(n,idx),'g','MarkerSize',15)
    hold on
     plot([1:1:length(idx)],BM(n,idx),'m');
    hold on
      errorbar([1:1:12],BM(n,idx),BSEM(n,idx),'m','MarkerSize',15)
      xlim([0 13]);
      ylim([0 11]);
title('grid')


n=2;
  subplot(143);
    plot([1:1:length(idx)],GM(n,idx),'g');
    hold on
    errorbar([1:1:12],GM(n,idx),GSEM(n,idx),'g','MarkerSize',15)
     hold on
     plot([1:1:length(idx)],BM(n,idx),'m');
    hold on
      errorbar([1:1:12],BM(n,idx),BSEM(n,idx),'m','MarkerSize',15)
      xlim([0 13]);
      ylim([0 11]);
title('cueR')

n=3;
  subplot(144);
    plot([1:1:length(idx)],GM(n,idx),'g');
    hold on
    errorbar([1:1:12],GM(n,idx),GSEM(n,idx),'g','MarkerSize',15)
     hold on
     plot([1:1:length(idx)],BM(n,idx),'m');
    hold on
      errorbar([1:1:12],BM(n,idx),BSEM(n,idx),'m','MarkerSize',15)
      xlim([0 13]);
      ylim([0 11]);
title('cueL')

saveas(gcf,'daysTrackedAllGridCueCells.fig');

%% send to taylor for statstics
pieData={};
pieData{1}={};
pieData{2}={};

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\nCellsAllPerM2.mat');
pieData{1}{1}=nCellsAllPerM2;

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\cellTypesNewAllFractinEachDayAllFOVs.mat');
G=cellTypesNewAllFractinEachDayAllFOVs;
idx=[1 4 5 6 7 8 9 10 11 12 13 14];
%since the two days (2 and 3) were removed, re normalize percentile 
A=G{1};
ANormFold=1-sum(A(:,2:3),2);
AN=[];
for n=1:size(A,1);
    AN(n,:)=A(n,idx)/ANormFold(n);
end
pieData{1}{2}=AN;

A=G{2};
ANormFold=1-sum(A(:,2:3),2);
AN=[];
for n=1:size(A,1);
    AN(n,:)=A(n,idx)/ANormFold(n);
end
pieData{1}{3}=AN;

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\commonUniqueCells\nCellsAllPerM2.mat');
pieData{2}{1}=nCellsAllPerM2;

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\commonUniqueCells\cellTypesNewAllFractinEachDayAllFOVs.mat');
B=cellTypesNewAllFractinEachDayAllFOVs;
idx=[1 4 5 6 7 8 9 10 11 12 13 14];
%since the two days (2 and 3) were removed, re normalize percentile 
A=B{1};
ANormFold=1-sum(A(:,2:3),2);
AN=[];
for n=1:size(A,1);
    AN(n,:)=A(n,idx)/ANormFold(n);
end
pieData{2}{2}=AN;

A=B{2};
ANormFold=1-sum(A(:,2:3),2);
AN=[];
for n=1:size(A,1);
    AN(n,:)=A(n,idx)/ANormFold(n);
end
pieData{2}{3}=AN;

save('pieData.mat','pieData');


