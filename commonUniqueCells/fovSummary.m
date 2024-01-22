% cd ..\
% load('foldersAll.mat');


%% Since some of the above FOVs are different from foldersAll, remake the folder
foldersAllN={};
foldersAllN{1}='E:\ID20210207\findCells\loc1\subsetFOVOldNO11';
foldersAllN{2}='E:\ID20210207\findCells\loc2\subsetFOVWithOldNo11';
foldersAllN{3}='E:\ID20210519_1\findCells\loc1\subsetFOVOldNew1_10';
foldersAllN{4}='E:\ID20210519_1\findCells\loc2\subsetFOVOld1-10';
foldersAllN{5}='E:\ID20210519_2\findCells\loc1\subsetFOVsOld1_10';
foldersAllN{6}='E:\ID20210519_2\findCells\loc2\subsetFOVOld1-10';
foldersAllN{7}='E:\ID20210519_2\findCells\loc3\subsetFOVOld1_10';
foldersAllN{8}='E:\ID20210811A\findCells\loc1\oldNew1_10';
foldersAllN{9}='E:\ID20210811A\findCells\loc3\oldNew1_10';

save('foldersAllN.mat','foldersAllN');

%% 
p=pwd;
cOldTrackAll={};
uOldTrackAll={};
uNewTrackAll={};
trackLaterAll={};
nCellsAll=[];
nCellsAllM2=[];
for n=1:11;
    trackLaterAll{n}={};
end

for n=1:length(foldersAllN);
    disp(n)
    cd(foldersAllN{n});
    load('trackInfo.mat');
    load('trackLater.mat');
    load('nCells.mat');
    load('nCellsM2.mat');
    cOldTrackAll{n}=cOldTrack;
    uOldTrackAll{n}=uOldTrack;
    uNewTrackAll{n}=uNewTrack;
    for m=1:length(trackLater);
        trackLaterAll{m}{end+1}=trackLater{m};
    end
    nCellsAll(n,:)=nCells;
    nCellsAllM2(n,:)=nCellsM2;
end
cd(p)

%nCells order: c old, uOld, uNew, uOld calculated with second method, uNew
%calculated with second method, u 3 to 11 (day 10 in new) calculated using
%the second method

daysAll={};
daysAll{1}=cell2mat(cOldTrackAll);%number of days common old cell last
daysAll{2}=cell2mat(uOldTrackAll);%number of days unique old cell last
daysAll{3}=cell2mat(uNewTrackAll);%number of days unique old cell last
for n=1:length(trackLaterAll);
    daysAll{3+n}=cell2mat(trackLaterAll{n});
end


daysMean=[];
daysSEM=[];
daysSTD=[];
for n=1:length(daysAll);
    daysMean(n)=mean(daysAll{n});
    daysSEM(n)=std(daysAll{n})/sqrt(length(daysAll{n}));
    daysSTD(n)=std(daysAll{n});
end
    
% daysToNormalize=[11 11 10 11 10 9 8 7 6 5 4 3 2];%not including the ones come up the last day
% 
% daysAllNorm={};
% for n=1:length(daysAll)-1;
%     daysAllNorm{n}=daysAll{n}/daysToNormalize(n);
% end
% 
% daysMeanNorm=[];
% daysSEMNorm=[];
% daysSTDNorm=[];
% for n=1:length(daysAllNorm);
%     daysMeanNorm(n)=mean(daysAllNorm{n});
%     daysSEMNorm(n)=std(daysAllNorm{n})/sqrt(length(daysAllNorm{n}));
%     daysSTDNorm(n)=std(daysAllNorm{n});
% end

nCellsAllPer=[];
for n=1:size(nCellsAll,1);
nCellsAllPer(n,:)=nCellsAll(n,:)/sum(nCellsAll(n,:));
end

nCellsAllPerM2=[];
for n=1:size(nCellsAllM2,1);
nCellsAllPerM2(n,:)=nCellsAllM2(n,:)/sum(nCellsAllM2(n,:));
end

save('nCellsAll.mat','nCellsAll');
save('nCellsAllPer.mat','nCellsAllPer');
save('nCellsAllM2.mat','nCellsAllM2');
save('nCellsAllPerM2.mat','nCellsAllPerM2');
save('daysAll.mat','daysAll');
% save('daysAllNorm.mat','daysAllNorm');
save('cOldTrackAll.mat','cOldTrackAll');
save('uOldTrackAll.mat','uOldTrackAll');
save('uNewTrackAll.mat','uNewTrackAll');


%%
figure
subplot(231)
c=flipud(gray(length(trackLaterAll)+2));
for n=1:length(trackLaterAll);
    [p,x]=ksdensity(cell2mat(trackLaterAll{n}));
    hold on
    plot(x,p/max(p),'Color',c(n+2,:));
end
hold on
 [p,x]=ksdensity(cell2mat(cOldTrackAll));
  plot(x,p/max(p),'r');
hold on
[p,x]=ksdensity(cell2mat(uOldTrackAll));
  plot(x,p/max(p),'g');
hold on
[p,x]=ksdensity(cell2mat(uNewTrackAll));
  plot(x,p/max(p),'b');
  title(['later days showup']);

  subplot(234)
  errorbar([1:1:length(daysSEM)],daysMean,daysSTD,'.','MarkerSize',8)
hold on
line([5 14],[10 1])
title('days tracked');

%   subplot(246)
%   
%   errorbar([1:1:length(daysMeanNorm)],daysMeanNorm,daysSTDNorm,'.','MarkerSize',8)
% 
%   title('days tracked norm');

  subplot(232)
pie(mean(nCellsAllPer,1));
title('cell fraction')

subplot(235)
bar([1:1:size(nCellsAllPer,2)],mean(nCellsAllPer,1));
hold on
errorbar([1:1:size(nCellsAllPer,2)],mean(nCellsAllPer,1),std(nCellsAllPer,1)/sqrt(size(nCellsAllPer,1)),'.');
title('cell fraction')

  subplot(233)
pie(mean(nCellsAllPerM2,1));
title('cell fraction M2')

subplot(236)
bar([1:1:size(nCellsAllPerM2,2)],mean(nCellsAllPerM2,1));
hold on
errorbar([1:1:size(nCellsAllPerM2,2)],mean(nCellsAllPerM2,1),std(nCellsAllPerM2,1)/sqrt(size(nCellsAllPerM2,1)),'.');
title('cell fraction M2')


saveas(gcf,'trackingCells.fig');


%% look at the activity correlation in all categories
p=pwd;
corrAllAllS={};

for n=1:14;
corrAllAllS{n}=[];
end

for n=1:length(foldersAllN);
    disp(n)
    cd(foldersAllN{n});
    load('corrAllAll.mat');
    for m=1:length(corrAllAll);
        corrAllAllS{m}(end+1:end+size(corrAllAll{m},1),:)=corrAllAll{m};
    end
end
        
cd(p)

c=gray(length(corrAllAllS)-2);

figure
A=corrAllAllS{1};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M,S,'r');
hold on
A=corrAllAllS{2};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M,S,'g');
hold on
A=corrAllAllS{3};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([2:1:length(M)],M(2:end),S(2:end),'b');

for n=1:length(corrAllAllS)-3;
    
        A=corrAllAllS{n+3};
M=nanmean(A,1);
S=nansem(A,1);
hold on
errorbar([n:1:length(M)],M(n:end),S(n:end),'Color',c(n,:));
   
end

save('corrAllAllS.mat');

c=gray(length(corrAllAllS)-2);

figure

subplot(141)
A=corrAllAllS{1};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M,S,'r');
hold on
A=corrAllAllS{2};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M,S,'g');
hold on
A=corrAllAllS{3};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([2:1:length(M)],M(2:end),S(2:end),'b');

for n=3:length(corrAllAllS)-3;
   
        A=corrAllAllS{n+3};
M=nanmean(A,1);
S=nansem(A,1);
hold on
errorbar([n:1:length(M)],M(n:end),S(n:end),'Color',c(n,:));
   
end

title('correlation');

subplot(142)
idx=[1 2 3 6 7 8 9 10 11 12 13 14];
A=corrAllAllS(idx);
M=[];
S=[];
for n=1:length(idx);
    B=A{n};
    B=reshape(B,[1,size(B,1)*size(B,2)]);
    M(n)=nanmean(B);
    S(n)=nansem(B);
end

bar([1:1:length(idx)],M);
hold on
errorbar([1:1:length(idx)],M,S,'.');
title('correlation mean days');


subplot(143)
A=corrAllAllS{1};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M,S,'r');
hold on
A=corrAllAllS{4};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M,S,'g');
hold on
A=corrAllAllS{5};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([2:1:length(M)],M(2:end),S(2:end),'b');

for n=3:length(corrAllAllS)-3;
   
        A=corrAllAllS{n+3};
M=nanmean(A,1);
S=nansem(A,1);
hold on
errorbar([n:1:length(M)],M(n:end),S(n:end),'Color',c(n,:));
   
end
title('correlation M2');

subplot(144)
idx=[1 4 5 6 7 8 9 10 11 12 13 14];
A=corrAllAllS(idx);
M=[];
S=[];
for n=1:length(idx);
    B=A{n};
    B=reshape(B,[1,size(B,1)*size(B,2)]);
    M(n)=nanmean(B);
    S(n)=nansem(B);
end

bar([1:1:length(idx)],M);
hold on
errorbar([1:1:length(idx)],M,S,'.');

title('correlation mean days M2');


saveas(gcf,'correlation.fig');

%% look at the activity correlation in all categories
p=pwd;
corrMeanActAllAllS={};
load('foldersAllN.mat');
for n=1:14;
corrMeanActAllAllS{n}=[];
end

for n=1:length(foldersAllN);
    disp(n)
    cd(foldersAllN{n});
    load('corrMeanActAllAll.mat');
    for m=1:length(corrMeanActAllAll);
        corrMeanActAllAllS{m}(end+1:end+size(corrMeanActAllAll{m},1),:)=corrMeanActAllAll{m};
    end
end
        
cd(p)

save('corrMeanActAllAllS.mat','corrMeanActAllAllS');
c=gray(length(corrMeanActAllAllS)-2);

figure

subplot(141)
A=corrMeanActAllAllS{1};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M,S,'r');
hold on
A=corrMeanActAllAllS{2};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M,S,'g');
hold on
A=corrMeanActAllAllS{3};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([2:1:length(M)],M(2:end),S(2:end),'b');

for n=3:length(corrMeanActAllAllS)-3;
   
        A=corrMeanActAllAllS{n+3};
M=nanmean(A,1);
S=nansem(A,1);
hold on
errorbar([n:1:length(M)],M(n:end),S(n:end),'Color',c(n,:));
   
end

title('correlation');

subplot(142)
idx=[1 2 3 6 7 8 9 10 11 12 13 14];
A=corrMeanActAllAllS(idx);
M=[];
S=[];
for n=1:length(idx);
    B=A{n}(:,2:end);%remove the env switch
    B=reshape(B,[1,size(B,1)*size(B,2)]);
    M(n)=nanmean(B);
    S(n)=nansem(B);
end

bar([1:1:length(idx)],M);
hold on
errorbar([1:1:length(idx)],M,S,'.');
title('correlation mean days');


subplot(143)
A=corrMeanActAllAllS{1};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M,S,'r');
hold on
A=corrMeanActAllAllS{4};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M,S,'g');
hold on
A=corrMeanActAllAllS{5};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([2:1:length(M)],M(2:end),S(2:end),'b');

for n=3:length(corrMeanActAllAllS)-3;
   
        A=corrMeanActAllAllS{n+3};
M=nanmean(A,1);
S=nansem(A,1);
hold on
errorbar([n:1:length(M)],M(n:end),S(n:end),'Color',c(n,:));
   
end
title('correlation M2');

subplot(144)
idx=[1 4 5 6 7 8 9 10 11 12 13 14];
A=corrMeanActAllAllS(idx);
M=[];
S=[];
for n=1:length(idx);
    B=A{n}(:,2:end);%remove the env switch
    B=reshape(B,[1,size(B,1)*size(B,2)]);
    M(n)=nanmean(B);
    S(n)=nansem(B);
end

bar([1:1:length(idx)],M);
hold on
errorbar([1:1:length(idx)],M,S,'.');

title('correlation mean days M2');


saveas(gcf,'correlationMeanAct.fig');

%% cell types

cellTypesPerAll{1}=[];%grid
cellTypesPerAll{2}=[];%cellR
cellTypesPerAll{3}=[];%cellL
cellTypesPerAll{4}=[];%speedP
cellTypesPerAll{5}=[];%speedN

for n=1:length(foldersAllN);
    load([foldersAllN{n} '\cellTypesPer.mat']);
    for  m=1:size(cellTypesPer,1);
        cellTypesPerAll{m}(n,:)=cellTypesPer(m,:);
    end
end
    
figure
for n=1:length(cellTypesPerAll);
    subplot(1,5,n)
   A=cellTypesPerAll{n};
   M=nanmean(A,1);
   S=nansem(A,1);
   
   errorbar([1:1:length(M)],M,S);
    if n==1;
        title('grid');
    elseif n==2;
        title('cueR');
          elseif n==3;
        title('cueL');
          elseif n==4;
        title('SpeedP');
    else 
        title('SpeedN');
    end
              
end
save('cellTypesPerAll.mat','cellTypesPerAll');
saveas(gcf,'cellTypesPer.fig')

%% cue score amd speed scores
p=pwd;

%max cue score on the first appear day
load('foldersAllN.mat');
maxCueScoresAppearDayAll={};
for n=1:15;
    maxCueScoresAppearDayAll{n}=[];
end

for n=1:length(foldersAllN);
    cd(foldersAllN{n})
    load('maxCueScoresAppearDay.mat');
    for m=1:length(maxCueScoresAppearDay);
       maxCueScoresAppearDayAll{m}(end+1:end+length(maxCueScoresAppearDay{m}))=maxCueScoresAppearDay{m};
    end
end
cd(p);

M=[];
S=[];
for n=1:length(maxCueScoresAppearDayAll);
    A=maxCueScoresAppearDayAll{n};
    if ~isempty(A);
    M(n)=nanmean(A);
    S(n)=nansem(A);
    else
        M(n)=nan;
        S(n)=nan;
    end
end

%cue score throughout days
cueScoresRAllDaysAll={};
cueScoresLAllDaysAll={};
speedScoresAllDaysAll={};
for n=1:15;
    cueScoresRAllDaysAll{n}=[];
    cueScoresLAllDaysAll{n}=[];
 speedScoresAllDaysAll{n}=[];
end

for n=1:length(foldersAllN);
    cd(foldersAllN{n})
     load('cueScoresRAllDays.mat');
load('cueScoresLAllDays.mat');
load('speedScoresAllDays.mat');
for m=1:15;
    cueScoresRAllDaysAll{m}(end+1:end+size(cueScoresRAllDays{m},1),:)=cueScoresRAllDays{m};
    cueScoresLAllDaysAll{m}(end+1:end+size(cueScoresLAllDays{m},1),:)=cueScoresLAllDays{m};
    speedScoresAllDaysAll{m}(end+1:end+size(speedScoresAllDays{m},1),:)=speedScoresAllDays{m};
end
end

cd(p)

save('maxCueScoresAppearDayAll.mat','maxCueScoresAppearDayAll')
save('cueScoresRAllDaysAll.mat','cueScoresRAllDaysAll');
save('cueScoresLAllDaysAll.mat','cueScoresLAllDaysAll');
save('speedScoresAllDaysAll.mat','speedScoresAllDaysAll');


figure,
subplot(241)
errorbar([1:1:length(M)],M,S,'r');
title('max cue score on appear day');

c=gray(length(cueScoresRAllDaysAll)-2);

subplot(242)
D=cueScoresRAllDaysAll;
A=D{1};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M,S,'r');
hold on
A=D{2};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M,S,'g');
hold on
A=D{4};%uNew
M=nanmean(A,1);
S=nansem(A,1);
errorbar([2:1:length(M)],M(2:end),S(2:end),'b');

% for n=1:length(D)-3;   
for n=3;% uNewM2 (6th data)
    if ~isempty(D{n+3});
        A=D{n+3};
M=nanmean(A,1);
S=nansem(A,1);
hold on
errorbar([n:1:length(M)],M(n:end),S(n:end),'Color',c(n,:));
    end
end

title('cue score R change');

subplot(246)
A=cueScoresRAllDaysAll;
M=[];
S=[];
for n=1:length(idx);
    B=A{n}(:,2:end);%remove the env switch
    B=reshape(B,[1,size(B,1)*size(B,2)]);
    M(n)=nanmean(B);
    S(n)=nansem(B);
end

bar([1:1:length(M)],M);
hold on
errorbar([1:1:length(M)],M,S,'.');


subplot(243)
D=cueScoresLAllDaysAll;
A=D{1};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M,S,'r');
hold on
A=D{2};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M,S,'g');
hold on
A=D{4};%uNew
M=nanmean(A,1);
S=nansem(A,1);
errorbar([2:1:length(M)],M(2:end),S(2:end),'b');

% for n=1:length(D)-3;   
for n=3;% uNewM2 (6th data)
    if ~isempty(D{n+3});
        A=D{n+3};
M=nanmean(A,1);
S=nansem(A,1);
hold on
errorbar([n:1:length(M)],M(n:end),S(n:end),'Color',c(n,:));
    end
end

title('cue score L change');

subplot(247)
A=cueScoresLAllDaysAll;
M=[];
S=[];
for n=1:length(idx);
    B=A{n}(:,2:end);%remove the env switch
    B=reshape(B,[1,size(B,1)*size(B,2)]);
    M(n)=nanmean(B);
    S(n)=nansem(B);
end

bar([1:1:length(M)],M);
hold on
errorbar([1:1:length(M)],M,S,'.');



subplot(244)
D=speedScoresAllDaysAll;%cOld
A=D{1};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M,S,'r');
hold on
A=D{2};%cNew
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M,S,'g');
hold on
A=D{4};%uNew
M=nanmean(A,1);
S=nansem(A,1);
errorbar([2:1:length(M)],M(2:end),S(2:end),'b');

% for n=1:length(D)-3;   
for n=3;% uNewM2 (6th data)
    if ~isempty(D{n+3});
        A=D{n+3};
M=nanmean(A,1);
S=nansem(A,1);
hold on
errorbar([n:1:length(M)],M(n:end),S(n:end),'Color',c(n,:));
    end
end

title('speed Score change');

subplot(248)
A=speedScoresAllDaysAll;
M=[];
S=[];
for n=1:length(idx);
    B=A{n}(:,2:end);%remove the env switch
    B=reshape(B,[1,size(B,1)*size(B,2)]);
    M(n)=nanmean(B);
    S(n)=nansem(B);
end

bar([1:1:length(M)],M);
hold on
errorbar([1:1:length(M)],M,S,'.');

saveas(gcf,'cueSpeed.fig');

%% 
load('foldersAllN.mat');
cellTypesNewAllFinalPercentAllFOVs={};
cellTypesNewAllFinalPercentAllFOVs{1}=[];%grid cells
cellTypesNewAllFinalPercentAllFOVs{2}=[];%cueR cells
cellTypesNewAllFinalPercentAllFOVs{3}=[];%cueL cells
cellTypesNewAllFinalPercentAllFOVs{4}=[];%speedP cells
cellTypesNewAllFinalPercentAllFOVs{5}=[];%speedN cells
p=pwd;
for n=1:length(foldersAllN);
    cd(foldersAllN{n});
    load('cellTypesNewAllFinalPercent.mat');
    cellTypesNewAllFinalPercentAllFOVs{1}(n,:)=cellTypesNewAllFinalPercent(1,:);
    cellTypesNewAllFinalPercentAllFOVs{2}(n,:)=cellTypesNewAllFinalPercent(2,:);
    cellTypesNewAllFinalPercentAllFOVs{3}(n,:)=cellTypesNewAllFinalPercent(3,:);
    cellTypesNewAllFinalPercentAllFOVs{4}(n,:)=cellTypesNewAllFinalPercent(4,:);
    cellTypesNewAllFinalPercentAllFOVs{5}(n,:)=cellTypesNewAllFinalPercent(5,:);
end
cd(p);

figure,

for n=1:length(cellTypesNewAllFinalPercentAllFOVs);
    A=cellTypesNewAllFinalPercentAllFOVs{n};
    subplot(1,5,n);
    errorbar([1:1:15],nanmean(A,1),nansem(A,1));
       if n==1;
        title('grid');
    elseif n==2;
        title('cueR');
    elseif n==3;
        title('cueL');
         elseif n==4;
        title('speedP');
    else
        title('speedN');
    end 
end
save('cellTypesNewAllFinalPercentAllFOVs.mat','cellTypesNewAllFinalPercentAllFOVs')

saveas(gcf,'cellTypesNewAllFinalPercentAllFOVs.fig');
%% here we count cells as a certain cell type in more than half of the track days, and see what fractions of cells appreaed on each day
load('foldersAllN.mat');
cellTypesNewAllFractinEachDayAllFOVs={};
cellTypesNewAllFractinEachDayAllFOVs{1}=[];%grid cells
cellTypesNewAllFractinEachDayAllFOVs{2}=[];%cueR cells
cellTypesNewAllFractinEachDayAllFOVs{3}=[];%cueL cells
cellTypesNewAllFractinEachDayAllFOVs{4}=[];%speedP cells
cellTypesNewAllFractinEachDayAllFOVs{5}=[];%speedN cells
p=pwd;
for n=1:length(foldersAllN);
    cd(foldersAllN{n});
    load('cellTypeFractionEachDayInAll.mat');
    cellTypesNewAllFractinEachDayAllFOVs{1}(n,:)=cellTypeFractionEachDayInAll(1,:);
    cellTypesNewAllFractinEachDayAllFOVs{2}(n,:)=cellTypeFractionEachDayInAll(2,:);
    cellTypesNewAllFractinEachDayAllFOVs{3}(n,:)=cellTypeFractionEachDayInAll(3,:);
    cellTypesNewAllFractinEachDayAllFOVs{4}(n,:)=cellTypeFractionEachDayInAll(4,:);
    cellTypesNewAllFractinEachDayAllFOVs{5}(n,:)=cellTypeFractionEachDayInAll(5,:);
end
cd(p);

figure,

for n=1:length(cellTypesNewAllFractinEachDayAllFOVs);
    A=cellTypesNewAllFractinEachDayAllFOVs{n};
    subplot(1,5,n);
    errorbar([1:1:14],nanmean(A,1),nansem(A,1));
       if n==1;
        title('grid');
    elseif n==2;
        title('cueR');
    elseif n==3;
        title('cueL');
         elseif n==4;
        title('speedP');
    else
        title('speedN');
    end 
end
save('cellTypesNewAllFractinEachDayAllFOVs.mat','cellTypesNewAllFractinEachDayAllFOVs')
saveas(gcf,'cellTypesNewAllFractinEachDayAllFOVs.fig');
%% here we count cells as a certain cell type in more than half of the track days, and see what fractions of cells appreaed on each day
load('foldersAllN.mat');
cellTypesNewAll_NDayCellType_passThreshAll={};
cellTypesNewAll_NDayCellType_passThreshAll{1}={};%grid cells
cellTypesNewAll_NDayCellType_passThreshAll{2}={};%cueR cells
cellTypesNewAll_NDayCellType_passThreshAll{3}={};%cueL cells
cellTypesNewAll_NDayCellType_passThreshAll{4}={};%speedP cells
cellTypesNewAll_NDayCellType_passThreshAll{5}={};%speedN cells
p=pwd;

for i=1:5;
    for m=1:15;
         cellTypesNewAll_NDayCellType_passThreshAll{i}{m}=[];
        for n=1:length(foldersAllN);
    cd(foldersAllN{n});
    load('cellTypesNewAll_NDayCellType_passThresh.mat');
       
        a=cellTypesNewAll_NDayCellType_passThresh{i}{m};
            cellTypesNewAll_NDayCellType_passThreshAll{i}{m}(end+1:end+length(a))=a;
    end
    end
end
cd(p);



cellTypeDayM=[];
cellTypeDaySEM=[];
for n=1:5;
    for m=1:15;
        if ~isempty(cellTypesNewAll_NDayCellType_passThreshAll{n}{m});
        cellTypeDayM(n,m)=nanmean(cellTypesNewAll_NDayCellType_passThreshAll{n}{m});
        cellTypeDaySEM(n,m)=nansem(cellTypesNewAll_NDayCellType_passThreshAll{n}{m},2);
        else
            cellTypeDayM(n,m)=nan;
        cellTypeDaySEM(n,m)=nan;
            
    end
    end
end

figure
for n=1:size(cellTypeDayM,1);
    
    subplot(1,5,n);
    errorbar([1:1:15],cellTypeDayM(n,:),cellTypeDaySEM(n,:));
       if n==1;
        title('grid');
    elseif n==2;
        title('cueR');
    elseif n==3;
        title('cueL');
         elseif n==4;
        title('speedP');
    else
        title('speedN');
    end 
end
save('cellTypesNewAll_NDayCellType_passThreshAll.mat','cellTypesNewAll_NDayCellType_passThreshAll')
saveas(gcf,'cellTypesNewAll_NDayCellType_passThreshAll.fig');
save('cellTypeDayM.mat','cellTypeDayM')
save('cellTypeDaySEM.mat','cellTypeDaySEM')

%% here we count cells as a certain cell type in more than half of the track days, and see what fractions of cells appreaed on each day
load('foldersAllN.mat');
cellTypesNewAll_NDayTrack_passThreshAll={};
cellTypesNewAll_NDayTrack_passThreshAll{1}={};%grid cells
cellTypesNewAll_NDayTrack_passThreshAll{2}={};%cueR cells
cellTypesNewAll_NDayTrack_passThreshAll{3}={};%cueL cells
cellTypesNewAll_NDayTrack_passThreshAll{4}={};%speedP cells
cellTypesNewAll_NDayTrack_passThreshAll{5}={};%speedN cells
p=pwd;

for i=1:5;
    for m=1:15;
         cellTypesNewAll_NDayTrack_passThreshAll{i}{m}=[];
        for n=1:length(foldersAllN);
    cd(foldersAllN{n});
    load('cellTypesNewAll_NDayTracked_passThresh.mat');
       
        a=cellTypesNewAll_NDayTracked_passThresh{i}{m};
            cellTypesNewAll_NDayTrack_passThreshAll{i}{m}(end+1:end+length(a))=a;
    end
    end
end
cd(p);



cellTypeDayTrackM=[];
cellTypeDayTrackSEM=[];
for n=1:5;
    for m=1:15;
        if ~isempty(cellTypesNewAll_NDayTrack_passThreshAll{n}{m});
        cellTypeDayTrackM(n,m)=nanmean(cellTypesNewAll_NDayTrack_passThreshAll{n}{m});
        cellTypeDayTrackSEM(n,m)=nansem(cellTypesNewAll_NDayTrack_passThreshAll{n}{m},2);
        else
            cellTypeDayTrackM(n,m)=nan;
        cellTypeDayTrackSEM(n,m)=nan;
            
    end
    end
end

figure
for n=1:size(cellTypeDayTrackM,1);
    
    subplot(1,5,n);
    errorbar([1:1:15],cellTypeDayTrackM(n,:),cellTypeDayTrackSEM(n,:));
       if n==1;
        title('grid');
    elseif n==2;
        title('cueR');
    elseif n==3;
        title('cueL');
         elseif n==4;
        title('speedP');
    else
        title('speedN');
    end 
end
save('cellTypesNewAll_NDayTrack_passThreshAll.mat','cellTypesNewAll_NDayTrack_passThreshAll')
saveas(gcf,'cellTypesNewAll_NDayTrack_passThreshAll.fig');

save('cellTypeDayTrackM.mat','cellTypeDayTrackM')
save('cellTypeDayTrackSEM.mat','cellTypeDayTrackSEM')

%% look at the activity correlation to template: RL, R, L
p=pwd;
corrAllAllSTempRL={};
corrAllAllSTempR={};
corrAllAllSTempL={};

load('foldersAllN.mat');
for n=1:15;
corrAllAllSTempRL{n}=[];
corrAllAllSTempR{n}=[];
corrAllAllSTempL{n}=[];
end

for n=1:length(foldersAllN);
    disp(n)
    cd(foldersAllN{n});
    load('corrAllAllTempRL.mat');
    load('corrAllAllTempR.mat');
    load('corrAllAllTempL.mat');
    for m=1:length(corrAllAllTempRL);
        corrAllAllSTempRL{m}(end+1:end+size(corrAllAllTempRL{m},1),:)=corrAllAllTempRL{m};
         corrAllAllSTempR{m}(end+1:end+size(corrAllAllTempR{m},1),:)=corrAllAllTempR{m};
          corrAllAllSTempL{m}(end+1:end+size(corrAllAllTempL{m},1),:)=corrAllAllTempL{m};
    end
end
        
cd(p)

save('corrAllAllSTempRL.mat','corrAllAllSTempRL');
save('corrAllAllSTempR.mat','corrAllAllSTempR');
save('corrAllAllSTempL.mat','corrAllAllSTempL');

c=gray(length(corrAllAllSTempRL)-1);

figure

subplot(341)
A=corrAllAllSTempRL{1};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([2:1:length(M)],M(2:end),S(2:end),'r');
hold on
A=corrAllAllSTempRL{2};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([2:1:length(M)],M(2:end),S(2:end),'r'); %new env temp only

hold on
A=corrAllAllSTempRL{3};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([2:1:length(M)],M(2:end),S(2:end),'g'); %new env temp only
hold on
A=corrAllAllSTempRL{4};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([2:1:length(M)],M(2:end),S(2:end),'b'); %new env temp only

for n=3:length(corrAllAllSTempRL)-4;
   
        A=corrAllAllSTempRL{n+4};
M=nanmean(A,1);
S=nansem(A,1);
hold on
errorbar([n:1:length(M)],M(n:end),S(n:end),'Color',c(n,:));
   
end

title('correlationNewTempRL');

subplot(342)
idx=[1 2 3 4 7 8 9 10 11 12 13 14 15];
A=corrAllAllSTempRL(idx);
M=[];
S=[];
for n=1:length(idx);
    B=A{n}(:,2:end); %new env temp only
    B=reshape(B,[1,size(B,1)*size(B,2)]);
    M(n)=nanmean(B);
    S(n)=nansem(B);
end

bar([1:1:length(idx)],M);
hold on
errorbar([1:1:length(idx)],M,S,'.');
title('correlationNewTempRL mean days');


subplot(343)
A=corrAllAllSTempRL{1};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([2:1:length(M)],M(2:end),S(2:end),'r'); %new env temp only
hold on
A=corrAllAllSTempRL{2};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([2:1:length(M)],M(2:end),S(2:end),'r'); %new env temp only

hold on
A=corrAllAllSTempRL{5};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([2:1:length(M)],M(2:end),S(2:end),'g'); %new env temp only
hold on
A=corrAllAllSTempRL{6};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([2:1:length(M)],M(2:end),S(2:end),'b'); %new env temp only

for n=3:length(corrAllAllSTempRL)-4;
   
        A=corrAllAllSTempRL{n+4};
M=nanmean(A,1);
S=nansem(A,1);
hold on
errorbar([n:1:length(M)],M(n:end),S(n:end),'Color',c(n,:));
   
end
title('correlationNewTempRL M2');

subplot(344)
idx=[1 2 5 6 7 8 9 10 11 12 13 14 15];
A=corrAllAllSTempRL(idx);
M=[];
S=[];
for n=1:length(idx);
    B=A{n}(:,2:end); %new env temp only
    B=reshape(B,[1,size(B,1)*size(B,2)]);
    M(n)=nanmean(B);
    S(n)=nansem(B);
end

bar([1:1:length(idx)],M);
hold on
errorbar([1:1:length(idx)],M,S,'.');

title('correlationNewTempRL mean days M2');

subplot(345)
A=corrAllAllSTempR{1};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([2:1:length(M)],M(2:end),S(2:end),'r');
hold on
A=corrAllAllSTempR{2};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([2:1:length(M)],M(2:end),S(2:end),'r'); %new env temp only

hold on
A=corrAllAllSTempR{3};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([2:1:length(M)],M(2:end),S(2:end),'g'); %new env temp only
hold on
A=corrAllAllSTempR{4};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([2:1:length(M)],M(2:end),S(2:end),'b'); %new env temp only

for n=3:length(corrAllAllSTempR)-4;
   
        A=corrAllAllSTempR{n+4};
M=nanmean(A,1);
S=nansem(A,1);
hold on
errorbar([n:1:length(M)],M(n:end),S(n:end),'Color',c(n,:));
   
end

title('correlationNewTempR');

subplot(346)
idx=[1 2 3 4 7 8 9 10 11 12 13 14 15];
A=corrAllAllSTempR(idx);
M=[];
S=[];
for n=1:length(idx);
    B=A{n}(:,2:end); %new env temp only
    B=reshape(B,[1,size(B,1)*size(B,2)]);
    M(n)=nanmean(B);
    S(n)=nansem(B);
end

bar([1:1:length(idx)],M);
hold on
errorbar([1:1:length(idx)],M,S,'.');
title('correlationNewTempR mean days');


subplot(347)
A=corrAllAllSTempR{1};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([2:1:length(M)],M(2:end),S(2:end),'r'); %new env temp only
hold on
A=corrAllAllSTempR{2};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([2:1:length(M)],M(2:end),S(2:end),'r'); %new env temp only

hold on
A=corrAllAllSTempR{5};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([2:1:length(M)],M(2:end),S(2:end),'g'); %new env temp only
hold on
A=corrAllAllSTempR{6};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([2:1:length(M)],M(2:end),S(2:end),'b'); %new env temp only

for n=3:length(corrAllAllSTempR)-4;
   
        A=corrAllAllSTempR{n+4};
M=nanmean(A,1);
S=nansem(A,1);
hold on
errorbar([n:1:length(M)],M(n:end),S(n:end),'Color',c(n,:));
   
end
title('correlationNewTempR M2');

subplot(348)
idx=[1 2 5 6 7 8 9 10 11 12 13 14 15];
A=corrAllAllSTempR(idx);
M=[];
S=[];
for n=1:length(idx);
    B=A{n}(:,2:end); %new env temp only
    B=reshape(B,[1,size(B,1)*size(B,2)]);
    M(n)=nanmean(B);
    S(n)=nansem(B);
end

bar([1:1:length(idx)],M);
hold on
errorbar([1:1:length(idx)],M,S,'.');

title('correlationNewTempR mean days M2');

subplot(349)
A=corrAllAllSTempL{1};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([2:1:length(M)],M(2:end),S(2:end),'r');
hold on
A=corrAllAllSTempL{2};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([2:1:length(M)],M(2:end),S(2:end),'r'); %new env temp only

hold on
A=corrAllAllSTempL{3};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([2:1:length(M)],M(2:end),S(2:end),'g'); %new env temp only
hold on
A=corrAllAllSTempL{4};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([2:1:length(M)],M(2:end),S(2:end),'b'); %new env temp only

for n=3:length(corrAllAllSTempL)-4;
   
        A=corrAllAllSTempL{n+4};
M=nanmean(A,1);
S=nansem(A,1);
hold on
errorbar([n:1:length(M)],M(n:end),S(n:end),'Color',c(n,:));
   
end

title('correlationNewTempL');

subplot(3,4,10)
idx=[1 2 3 4 7 8 9 10 11 12 13 14 15];
A=corrAllAllSTempL(idx);
M=[];
S=[];
for n=1:length(idx);
    B=A{n}(:,2:end); %new env temp only
    B=reshape(B,[1,size(B,1)*size(B,2)]);
    M(n)=nanmean(B);
    S(n)=nansem(B);
end

bar([1:1:length(idx)],M);
hold on
errorbar([1:1:length(idx)],M,S,'.');
title('correlationNewTempL mean days');

subplot(3,4,11)
A=corrAllAllSTempL{1};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([2:1:length(M)],M(2:end),S(2:end),'r'); %new env temp only
hold on
A=corrAllAllSTempL{2};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([2:1:length(M)],M(2:end),S(2:end),'r'); %new env temp only

hold on
A=corrAllAllSTempL{5};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([2:1:length(M)],M(2:end),S(2:end),'g'); %new env temp only
hold on
A=corrAllAllSTempL{6};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([2:1:length(M)],M(2:end),S(2:end),'b'); %new env temp only

for n=3:length(corrAllAllSTempL)-4;
   
        A=corrAllAllSTempL{n+4};
M=nanmean(A,1);
S=nansem(A,1);
hold on
errorbar([n:1:length(M)],M(n:end),S(n:end),'Color',c(n,:));
   
end
title('correlationNewTempL M2');

subplot(3,4,12)
idx=[1 2 5 6 7 8 9 10 11 12 13 14 15];
A=corrAllAllSTempL(idx);
M=[];
S=[];
for n=1:length(idx);
    B=A{n}(:,2:end); %new env temp only
    B=reshape(B,[1,size(B,1)*size(B,2)]);
    M(n)=nanmean(B);
    S(n)=nansem(B);
end

bar([1:1:length(idx)],M);
hold on
errorbar([1:1:length(idx)],M,S,'.');

title('correlationNewTempL mean days M2');

saveas(gcf,'correlationToNewTemp.fig');
%%
load('E:\ID20201118\20201226\2ndloc1\pcaica\cueAnalysisNew_sig\tempRL.mat');
load('E:\ID20201118\20201226\2ndloc1\pcaica\cueAnalysisNew_sig\tempR.mat');
load('E:\ID20201118\20201226\2ndloc1\pcaica\cueAnalysisNew_sig\tempL.mat');

load('foldersAllN.mat');
fieldDistriCommon=[];%field only activity on day1 for common cells between old and new day1
fieldDistriUnique=[];%field only activity on day1 for Unique cells on new day1
p=pwd;
corrRLNewDay1=[];
corrRNewDay1=[];
corrLNewDay1=[];
for n=1:length(foldersAllN);
    cd(foldersAllN{n});
    load('trackIDsAll.mat');
    load('useFolders.mat');
    ic=trackIDsAll{2};
    iu=trackIDsAll{4};
    load([useFolders{2},'\PValueClassifier_KY2_6_sig\allCellsCorrected.mat']);
ica=allCellsCorrected.dfofaveragesmoothFields(:,ic);
ica(ica>0)=1;
fieldDistriCommon(:,n)=mean(ica,2);
iua=allCellsCorrected.dfofaveragesmoothFields(:,iu);
iua(iua>0)=1;
fieldDistriUnique(:,n)=mean(iua,2);

corrRLNewDay1(n,1)=corr(tempRL,mean(ica,2));
corrRLNewDay1(n,2)=corr(tempRL,mean(iua,2));
corrRNewDay1(n,1)=corr(tempR,mean(ica,2));
corrRNewDay1(n,2)=corr(tempR,mean(iua,2));
corrLNewDay1(n,1)=corr(tempL,mean(ica,2));
corrLNewDay1(n,2)=corr(tempL,mean(iua,2));
end
cd(p)
figure,
A=[corrRLNewDay1 corrRNewDay1 corrLNewDay1];
bar([1:1:6],mean(A,1));
hold on
errorbar([1:1:6],mean(A,1),nansem(A,1),'.');

saveas(gcf,'corrTempOnlyCommonUniqueNewDay1.fig');

save('corrRLNewDay1.mat','corrRLNewDay1');
save('corrRNewDay1.mat','corrRNewDay1');
save('corrLNewDay1.mat','corrLNewDay1');

%% get an idea about the percetage of cells that only showed up once in each category

load('foldersAllN.mat');
perOnce=[];%first column is for uOld, second is for uNew, later are for uOld, uNewDay1-11 M2 method)
p=pwd;
for n=1:length(foldersAllN);
    load([foldersAllN{n},'\trackInfo.mat']);
     load([foldersAllN{n},'\trackLater.mat']);
     perOnce(n,1)=length(find(uOldTrack==1))/length(uOldTrack);
     perOnce(n,2)=length(find(uNewTrack==1))/length(uNewTrack);
     for m=1:length(trackLater);
         perOnce(n,m+2)=length(find(trackLater{m}==1))/length(trackLater{m});
     end
end
cd(p);

figure
M=nanmean(perOnce,1);
S=nansem(perOnce,1);
errorbar([1:1:length(M)],M,S);
title('percent of cells showing up once');
ylabel('percent');
save('perOnce.mat','perOnce');
saveas(gcf,'perOnce.fig');

%% change of cue cell patterns

load('foldersAllN.mat');
corrAllAllTempRLAllFolders={};
for n=1:15;
    corrAllAllTempRLAllFolders{n}=[]; %corr of all cells in all categories
end

idxOnlyOld=[];%cells only showed up in old env
idxOnly10DayNew=[];%cells showed up in 10 days in new

p=pwd;

for n=1:length(foldersAllN);
    cd(foldersAllN{n});
    load('corrAllAllTempRL.mat');
    for m=1:length(corrAllAllTempRL);
        a=corrAllAllTempRL{m};
    corrAllAllTempRLAllFolders{m}(end+1:end+size(a,1),:)=a;
    end
     load([foldersAllN{n},'\trackInfo.mat']);
     aa=zeros(length(uOldTrack),1);
     aa(uOldTrack==1)=1;
     idxOnlyOld(end+1:end+length(aa),1)=aa;
      bb=zeros(length(uNewTrack),1);
       bb(uNewTrack==10)=1;
     idxOnly10DayNew(end+1:end+length(bb),1)=bb;
end

cd(p);
save('corrAllAllTempRLAllFolders.mat','corrAllAllTempRLAllFolders');
save('idxOnlyOld.mat','idxOnlyOld');
save('idxOnly10DayNew.mat','idxOnly10DayNew');

%just plot unique

figure,
A=corrAllAllTempRLAllFolders{1};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:11],M,S,'r')
hold on
A=corrAllAllTempRLAllFolders{4};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:11],M,S,'g')

A=corrAllAllTempRLAllFolders{3};
M=nanmean(A,1);
S=nansem(A,1);
hold on
errorbar([1:1:11],M,S,'b')

A=corrAllAllTempRLAllFolders{1}(find(idxOnlyOld),1);
M=nanmean(A,1);
S=nansem(A,1);
hold on,
errorbar(1,M,S,'b.','MarkerSize',10)

hold on
A=corrAllAllTempRLAllFolders{4}(find(idxOnly10DayNew),:);
M=nanmean(A,1);
S=nansem(A,1);
hold on,
errorbar([1:1:11],M,S,'m')
xlim([0 12]);

saveas(gcf,'corrTemp_CommonUniqueNewOld.fig');
%% look at the activity correlation to template: RL, R, L
p=pwd;
nBinFieldsAllAll={};

load('foldersAllN.mat');
for n=1:15;
nBinFieldsAllAll{n}=[];
end

for n=1:length(foldersAllN);
    disp(n)
    cd(foldersAllN{n});
    load('nBinFieldsAll.mat');
    for m=1:length(nBinFieldsAll);
        nBinFieldsAllAll{m}(end+1:end+size(nBinFieldsAll{m},1),:)=nBinFieldsAll{m};
        
    end
end
        
cd(p)

save('nBinFieldsAllAll.mat','nBinFieldsAllAll');


c=gray(length(nBinFieldsAllAll)-1);

figure

subplot(141)
A=nBinFieldsAllAll{1};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M(1:end),S(1:end),'r');
hold on
A=nBinFieldsAllAll{2};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M(1:end),S(1:end),'m'); %new env temp only

hold on
A=nBinFieldsAllAll{3};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M(1:end),S(1:end),'g'); %new env temp only
hold on
A=nBinFieldsAllAll{4};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M(1:end),S(1:end),'b'); %new env temp only

for n=3:length(nBinFieldsAllAll)-4;
   
        A=nBinFieldsAllAll{n+4};
M=nanmean(A,1);
S=nansem(A,1);
hold on
errorbar([n:1:length(M)],M(n:end),S(n:end),'Color',c(n,:));
   
end

title('nBinFieldsAllAll');

subplot(142)
idx=[1 2 3 4 7 8 9 10 11 12 13 14 15];
A=nBinFieldsAllAll(idx);
M=[];
S=[];
for n=1:length(idx);
    B=A{n}; %new env temp only
    B=reshape(B,[1,size(B,1)*size(B,2)]);
    M(n)=nanmean(B);
    S(n)=nansem(B);
end

bar([1:1:length(idx)],M);
hold on
errorbar([1:1:length(idx)],M,S,'.');
title('nBinFieldsAllAll mean days');


subplot(143)
A=nBinFieldsAllAll{1};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M(1:end),S(1:end),'r'); %new env temp only
hold on
A=nBinFieldsAllAll{2};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M(1:end),S(1:end),'r'); %new env temp only

hold on
A=nBinFieldsAllAll{5};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M(1:end),S(1:end),'g'); %new env temp only
hold on
A=nBinFieldsAllAll{6};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M(1:end),S(1:end),'b'); %new env temp only

for n=3:length(nBinFieldsAllAll)-4;
   
        A=nBinFieldsAllAll{n+4};
M=nanmean(A,1);
S=nansem(A,1);
hold on
errorbar([n:1:length(M)],M(n:end),S(n:end),'Color',c(n,:));
   
end
title('nBinFieldsAllAll M2');

subplot(144)
idx=[1 2 5 6 7 8 9 10 11 12 13 14 15];
A=nBinFieldsAllAll(idx);
M=[];
S=[];
for n=1:length(idx);
    B=A{n}; %new env temp only
    B=reshape(B,[1,size(B,1)*size(B,2)]);
    M(n)=nanmean(B);
    S(n)=nansem(B);
end

bar([1:1:length(idx)],M);
hold on
errorbar([1:1:length(idx)],M,S,'.');

title('nBinFieldsAllAll mean days M2');
%% PERCENTAGE OF EACH CELL TYPE: using the M2 method

load('foldersAllN.mat');
pCellTypes=[];

%first number is the percentage of grid, cue r, cue l, speed p, speed n,
%others
p=pwd;
for n=1:length(foldersAllN);
    disp(n)
    cd(foldersAllN{n})
    load('cellTypeNEachDayInAll.mat');
    load('trackIDsAll.mat');
    idxTrackIDsAll=[1 5 6 7:1:15];%because this one include c0 and c1 so has one more element than cellTypeNEachDayAll
    
    idxCellType=[1 4 5 6:1:14];
    allCellN=[];
    
    for i=1:length(idxTrackIDsAll);
        allCellN(i)=length(trackIDsAll{idxTrackIDsAll(i)});
    end
    allCellN=sum(allCellN);
    cellType=cellTypeNEachDayInAll(:,idxCellType);
    cellType=sum(cellType,2);
    pCellTypes(1:5,n)=cellType./allCellN;
    pCellTypes(6,n)=1-sum(cellType(1:3)./allCellN);
end
cd(p)
save('pCellTypes.mat','pCellTypes');
figure,
bar([1:1:5],mean(pCellTypes(1:5,:),2));
hold on
errorbar([1:1:5],mean(pCellTypes(1:5,:),2),nansem(pCellTypes(1:5,:),2),'k.');
saveas(gcf,'pCellTypes.fig');