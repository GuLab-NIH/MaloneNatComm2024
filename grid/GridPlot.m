load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\indicesAllNewCueThresh.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\corrAllFOVs.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\dfofSigNormAllFOVs.mat');

i=indicesAllNewCueThresh.GridIdx;
corrAllGrid.dfof{1}=corrAllFOVs.dfof{1}(i,:);
corrAllGrid.dfof{2}=corrAllFOVs.dfof{2}(i,:);
corrAllGrid.dfof{3}=corrAllFOVs.dfof{3}(i,:);

corrAllGrid.dfofSig{1}=corrAllFOVs.dfofSig{1}(i,:);
corrAllGrid.dfofSig{2}=corrAllFOVs.dfofSig{2}(i,:);


corrAllGrid.dfofSig{3}=corrAllFOVs.dfofSig{3}(i,:);

corrAllGrid.dfofMeanToMean=mean(corrAllGrid.dfof{1},1);
corrAllGrid.dfofMeanToNext=mean(corrAllGrid.dfof{2},1);
corrAllGrid.dfofMeanToOthers=mean(corrAllGrid.dfof{3},1);

corrAllGrid.dfofSigMeanToMean=mean(corrAllGrid.dfofSig{1},1);
corrAllGrid.dfofSigMeanToNext=mean(corrAllGrid.dfofSig{2},1);
corrAllGrid.dfofSigMeanToOthers=mean(corrAllGrid.dfofSig{3},1);

save('corrAllGrid.mat','corrAllGrid');
        
figure
subplot(141)
title('dfof');
plot(corrAllGrid.dfofMeanToMean);
hold on
plot(corrAllGrid.dfofMeanToNext);
hold on
plot(corrAllGrid.dfofMeanToOthers);
legend('mean','next','others','Location','SouthEast');

subplot(142)
title('dfof sig');
plot(corrAllGrid.dfofSigMeanToMean);
hold on
plot(corrAllGrid.dfofSigMeanToNext);
hold on
plot(corrAllGrid.dfofSigMeanToOthers);
legend('mean','next','others','Location','SouthEast');

%only use dfofSig Mean To Mean and do shuffle
A=corrAllGrid.dfofSig{3};
MTempRL=corrAllGrid.dfofSigMeanToOthers;
S=std(A,1)/sqrt(size(A,1));

%shuffle
nShuffle=100;
MShuffle=[];
for n=1:nShuffle;  
    shuffleDfofSig=[];
    for m=1:size(A,1);
        shuffleDfofSig(m,:)=A(m,randperm(size(A,2)));
    end
    MShuffle(n,:)=mean(shuffleDfofSig,1);
end

isSig=[];
for n=1:size(MShuffle,2);
    if MTempRL(n)>=prctile(MShuffle(:,n),95);
        isSig(n)=1;
    elseif MTempRL(n)<=prctile(MShuffle(:,n),5);
        isSig(n)=1;
    else
        isSig(n)=0;
    end
end

subplot(143)
errorbar([1:1:length(MTempRL)],MTempRL,S);
hold on
errorbar([1:1:length(MTempRL)],mean(MShuffle,1),std(MShuffle,1))
title('corr run-By-Run NE only');
xlabel('Days');
ylabel('Mean correlation');
for n=1:length(MTempRL);
    if isSig(n)==1;
        plot(n,MTempRL(n)+0.01,'r*');
    end
end
xlim([0 length(MTempRL)+1]);

[r,p]=corr([1:1:length(MTempRL)-1]',MTempRL(2:end)');
%r=0.7334
%p=0.0158

title(['sig meanToOthers new only,p=',num2str(p)]);


sigToNDay1=[];%significance compared to day1 in new env
day1=2;%day1 is the first column of "A".
A=corrAllGrid.dfofSig{3};

for n=1:size(A,2);
    if n==day1;
        sigToNDay1(n,1)=nan;
    else
        [sigToNDay1(n,1),~]=ttest2(A(:,day1),A(:,n));
    end
end

subplot(144)
errorbar([1:1:length(MTempRL)],MTempRL,S);
% hold on
% errorbar([1:1:length(MTempRL)],mean(MShuffle,1),std(MShuffle,1))
title('corr run-By-Run NE only');
xlabel('Days');
ylabel('Mean correlation');
for n=1:length(MTempRL);
    if sigToNDay1(n)==1;
        hold on
        plot(n,MTempRL(n)+0.01,'r*');
    end
end
xlim([0 length(MTempRL)+1]);

saveas(gcf,'corrMeanGrid.fig');
%%
corrAllGrid.dfofNorm{1}=corrAllFOVs.dfofNorm{1}(i,:);
corrAllGrid.dfofNorm{2}=corrAllFOVs.dfofNorm{2}(i,:);
corrAllGrid.dfofNorm{3}=corrAllFOVs.dfofNorm{3}(i,:);

corrAllGrid.dfofSigNorm{1}=corrAllFOVs.dfofSigNorm{1}(i,:);
corrAllGrid.dfofSigNorm{2}=corrAllFOVs.dfofSigNorm{2}(i,:);


corrAllGrid.dfofSigNorm{3}=corrAllFOVs.dfofSigNorm{3}(i,:);

a=corrAllGrid.dfofNorm{1};
b=[];
for n=1:size(a,2);
c=a(:,n);
   d=c(~isnan(c));
   b(1,n)=mean(d);
end
corrAllGrid.dfofNormMeanToMean=b;
a=corrAllGrid.dfofNorm{2};
b=[];
for n=1:size(a,2);
c=a(:,n);
   d=c(~isnan(c));
   b(1,n)=mean(d);
end
corrAllGrid.dfofNormMeanToNext=b;
a=corrAllGrid.dfofNorm{3};
b=[];
for n=1:size(a,2);
c=a(:,n);
   d=c(~isnan(c));
   b(1,n)=mean(d);
end
corrAllGrid.dfofNormMeanToOthers=b;


a=corrAllGrid.dfofSigNorm{1};
b=[];
for n=1:size(a,2);
c=a(:,n);
   d=c(~isnan(c));
   b(1,n)=mean(d);
end
corrAllGrid.dfofSigNormMeanToMean=b;
a=corrAllGrid.dfofSigNorm{2};
b=[];
for n=1:size(a,2);
c=a(:,n);
   d=c(~isnan(c));
   b(1,n)=mean(d);
end
corrAllGrid.dfofSigNormMeanToNext=b;
a=corrAllGrid.dfofSigNorm{3};
b=[];
for n=1:size(a,2);
c=a(:,n);
   d=c(~isnan(c));
   b(1,n)=mean(d);
end
corrAllGrid.dfofSigNormMeanToOthers=b;


save('corrAllGrid.mat','corrAllGrid');

%%
load('corrAllGrid.mat');
figure
subplot(241)
plot(corrAllGrid.dfofNormMeanToMean);
hold on
plot(corrAllGrid.dfofNormMeanToNext);
hold on
plot(corrAllGrid.dfofNormMeanToOthers);
legend('mean','next','others','Location','SouthEast');
title('dfof');

subplot(242)
plot(corrAllGrid.dfofSigNormMeanToMean);
hold on
plot(corrAllGrid.dfofSigNormMeanToNext);
hold on
plot(corrAllGrid.dfofSigNormMeanToOthers);
legend('mean','next','others','Location','SouthEast');
title('dfof sig');

%only use dfof Mean To Mean and do shuffle
A=corrAllGrid.dfofNorm{3};
MTempRL=corrAllGrid.dfofNormMeanToOthers;
S=std(A,1)/sqrt(size(A,1));

%shuffle
nShuffle=100;
MShuffle=[];
for n=1:nShuffle;  
    shuffleDfofSig=[];
    for m=1:size(A,1);
        shuffleDfofSig(m,:)=A(m,randperm(size(A,2)));
    end
    MShuffle(n,:)=mean(shuffleDfofSig,1);
end

isSig=[];
for n=1:size(MShuffle,2);
    if MTempRL(n)>=prctile(MShuffle(:,n),95);
        isSig(n)=1;
    elseif MTempRL(n)<=prctile(MShuffle(:,n),5);
        isSig(n)=1;
    else
        isSig(n)=0;
    end
end

subplot(243)
errorbar([1:1:length(MTempRL)],MTempRL,S);
hold on
errorbar([1:1:length(MTempRL)],mean(MShuffle,1),std(MShuffle,1))
% title('dfof meanToNext');
xlabel('Days');
ylabel('Mean correlation');
for n=1:length(MTempRL);
    if isSig(n)==1;
        plot(n,MTempRL(n)+0.01,'r*');
    end
end
xlim([0 length(MTempRL)+1]);

[r,p]=corr([1:1:length(MTempRL)-1]',MTempRL(2:end)');
%r=0.7334
%p=0.0391

title(['dfof meanToOthers,p=',num2str(p)]);


sigToNDay1=[];%significance compared to day1 in new env
day1=2;%day1 is the first column of "A".
A=corrAllGrid.dfofNorm{3};

for n=1:size(A,2);
    if n==day1;
        sigToNDay1(n,1)=nan;
    else
        [sigToNDay1(n,1),~]=ttest2(A(:,day1),A(:,n));
    end
end

subplot(244)
errorbar([1:1:length(MTempRL)],MTempRL,S);
% hold on
% errorbar([1:1:length(MTempRL)],mean(MShuffle,1),std(MShuffle,1))
title('dfof mean to others');
xlabel('Days');
ylabel('Mean correlation');
for n=1:length(MTempRL);
    if sigToNDay1(n)==1;
        hold on
        plot(n,MTempRL(n)+0.01,'r*');
    end
end
xlim([0 length(MTempRL)+1]);


%only use dfofSig Mean To others and do shuffle
A=corrAllGrid.dfofSigNorm{1};
MTempRL=corrAllGrid.dfofSigNormMeanToMean;
S=std(A,1)/sqrt(size(A,1));

%shuffle
nShuffle=100;
MShuffle=[];
for n=1:nShuffle;  
    shuffleDfofSig=[];
    for m=1:size(A,1);
        shuffleDfofSig(m,:)=A(m,randperm(size(A,2)));
    end
    MShuffle(n,:)=mean(shuffleDfofSig,1);
end

isSig=[];
for n=1:size(MShuffle,2);
    if MTempRL(n)>=prctile(MShuffle(:,n),95);
        isSig(n)=1;
    elseif MTempRL(n)<=prctile(MShuffle(:,n),5);
        isSig(n)=1;
    else
        isSig(n)=0;
    end
end

subplot(245)
errorbar([1:1:length(MTempRL)],MTempRL,S);
hold on
errorbar([1:1:length(MTempRL)],mean(MShuffle,1),std(MShuffle,1))
title('dfof mean to next');
xlabel('Days');
ylabel('Mean correlation');
for n=1:length(MTempRL);
    if isSig(n)==1;
        plot(n,MTempRL(n)+0.01,'r*');
    end
end
xlim([0 length(MTempRL)+1]);

[r,p]=corr([1:1:length(MTempRL)-1]',MTempRL(2:end)');
%r=0.7334
%p=0.0158

title(['sig meanToMean,p=',num2str(p)]);


sigToNDay1=[];%significance compared to day1 in new env
day1=2;%day1 is the first column of "A".
A=corrAllGrid.dfofSigNorm{3};

for n=1:size(A,2);
    if n==day1;
        sigToNDay1(n,1)=nan;
    else
        [sigToNDay1(n,1),~]=ttest2(A(:,day1),A(:,n));
    end
end

subplot(246)
errorbar([1:1:length(MTempRL)],MTempRL,S);
% hold on
% errorbar([1:1:length(MTempRL)],mean(MShuffle,1),std(MShuffle,1))
title('sig mean to mean');
xlabel('Days');
ylabel('Mean correlation');
for n=1:length(MTempRL);
    if sigToNDay1(n)==1;
        hold on
        plot(n,MTempRL(n)+0.01,'r*');
    end
end
xlim([0 length(MTempRL)+1]);

%only use dfof Mean To Mean and do shuffle
A=corrAllGrid.dfofNorm{3};
MTempRL=corrAllGrid.dfofNormMeanToOthers;
S=std(A,1)/sqrt(size(A,1));

%shuffle
nShuffle=100;
MShuffle=[];
for n=1:nShuffle;  
    shuffleDfofSig=[];
    for m=1:size(A,1);
        shuffleDfofSig(m,:)=A(m,randperm(size(A,2)));
    end
    MShuffle(n,:)=mean(shuffleDfofSig,1);
end

isSig=[];
for n=1:size(MShuffle,2);
    if MTempRL(n)>=prctile(MShuffle(:,n),95);
        isSig(n)=1;
    elseif MTempRL(n)<=prctile(MShuffle(:,n),5);
        isSig(n)=1;
    else
        isSig(n)=0;
    end
end

subplot(247)
errorbar([1:1:length(MTempRL)],MTempRL,S);
hold on
errorbar([1:1:length(MTempRL)],mean(MShuffle,1),std(MShuffle,1))
% title('corr run-By-Run NE only');
xlabel('Days');
ylabel('Mean correlation');
for n=1:length(MTempRL);
    if isSig(n)==1;
        plot(n,MTempRL(n)+0.01,'r*');
    end
end
xlim([0 length(MTempRL)+1]);

[r,p]=corr([1:1:length(MTempRL)-1]',MTempRL(2:end)');
%r=0.7334
%p=0.0158

title(['sig meanToOthers,p=',num2str(p)]);


sigToNDay1=[];%significance compared to day1 in new env
day1=2;%day1 is the first column of "A".
A=corrAllGrid.dfofNorm{3};

for n=1:size(A,2);
    if n==day1;
        sigToNDay1(n,1)=nan;
    else
        [sigToNDay1(n,1),~]=ttest2(A(:,day1),A(:,n));
    end
end

subplot(248)
errorbar([1:1:length(MTempRL)],MTempRL,S);
% hold on
% errorbar([1:1:length(MTempRL)],mean(MShuffle,1),std(MShuffle,1))
title('sig mean to others');
xlabel('Days');
ylabel('Mean correlation');
for n=1:length(MTempRL);
    if sigToNDay1(n)==1;
        hold on
        plot(n,MTempRL(n)+0.01,'r*');
    end
end
xlim([0 length(MTempRL)+1]);

tightfig

saveas(gcf,'corrMeanAllCellNormToItself.fig');



% %% make a few shuffles of grid cells to make sure that this is not due to more cells than cue cells
% load('corrAllGrid.mat');
% N=10;
% NGrid=size(corrAllGrid.dfofSig{1},1);
% NGridRand=20;%this is the same number of cue cell R (cue cell R has more cells)
% figure
% for n=1:N;
%     i=randperm(NGrid);
%     i=i(1:NGridRand);
%     C=corrAllGrid.dfofSigNorm{1}(i,:);
%     M=mean(C,1);
%     S=std(C,1)/sqrt(NGridRand);
%    [r,p]=corr([1:1:length(M)-1]',M(2:end)');
%     subplot(2,5,n);
%     errorbar([1:1:length(M)-1],M(2:end),S(2:end));
%     title(['p=',num2str(p)]);
% end
% saveas(gcf,'randGridWithEqualNumberOfCueR.fig');
% 
% N=100;
% allP=[];
% for n=1:N;
%     i=randperm(NGrid);
%     i=i(1:NGridRand);
%     C=corrAllGrid.dfofSigNorm{1}(i,:);
%     M=mean(C,1);
%     S=std(C,1)/sqrt(NGridRand);
%    [~,allP(n)]=corr([1:1:length(M)-1]',M(2:end)');
%     
% end
% 
% figure,plot(allP,'r.');
% hold on
% plot([0 100],[0.05 0.05],'g')
% title(['p values of ',num2str(N),' random shuffles, >0.05 is ',num2str(100*length(find(allP>0.05))/N),'%']);
% ylabel('P value');
% xlabel('Number of shuffles')
% save('allP.mat','allP');
% saveas(gcf,'randGridWithEqualNumberOfCueRMoreShuffle.fig');
% 
% 
% %% make a few shuffles of grid cells to make sure that this is not due to more cells than cue cells
% load('corrAllGrid.mat');
% N=10;
% NGrid=size(corrAllGrid.dfofSig{1},1);
% NGridRand=87;%this is the same number of cue cell R (cue cell R has more cells)
% figure
% for n=1:N;
%     i=randperm(NGrid);
%     i=i(1:NGridRand);
%     C=corrAllGrid.dfofSigNorm{1}(i,:);
%     M=mean(C,1);
%     S=std(C,1)/sqrt(NGridRand);
%    [r,p]=corr([1:1:length(M)-1]',M(2:end)');
%     subplot(2,5,n);
%     errorbar([1:1:length(M)-1],M(2:end),S(2:end));
%     title(['p=',num2str(p)]);
% end
% saveas(gcf,'randGridWithEqualNumberOfCueR.fig');
% 
% N=100;
% allP=[];
% for n=1:N;
%     i=randperm(NGrid);
%     i=i(1:NGridRand);
%     C=corrAllGrid.dfofSigNorm{1}(i,:);
%     M=mean(C,1);
%     S=std(C,1)/sqrt(NGridRand);
%    [~,allP(n)]=corr([1:1:length(M)-1]',M(2:end)');
%     
% end
% 
% figure,plot(allP,'r.');
% hold on
% plot([0 100],[0.05 0.05],'g')
% title(['p values of ',num2str(N),' random shuffles, >0.05 is ',num2str(100*length(find(allP>0.05))/N),'%']);
% ylabel('P value');
% xlabel('Number of shuffles')
% save('allPMoreCells.mat','allP');
% saveas(gcf,'randGridWithEqualNumberOfCueRMoreShuffle.fig');
% 
% %this means even there is an even number of cue cells and grid cells, grid
% %cells still show this feature.
%%
i=indicesAllNewCueThresh.GridIdx;

dfofSigNormGrid=[];
for n=1:length(dfofSigNormAllFOVs);
    dfofSigNormGrid{n}=dfofSigNormAllFOVs{n}(i,:);
end
save('dfofSigNormGrid.mat','dfofSigNormGrid');
%%
%% sort activity according to tempRL

%sort the correlation of the activity with template (two sides adding up
%together)
%load a two side template
load('E:\ID20201118\20201224\loc1\pcaica\cueAnalysisNew_sig\tempRL.mat');
temp=tempRL;

%sort baesd on Dfof
lagsDfof=[];
%use the cue cell code to calculate lags
%sort cells based on the first day in new env
lagBins=60;%shift freely
day=1;
for n=1:size(dfofSigNormGrid{day},1)
[~,lagsDfof(n,1)]=findCorrLag(dfofSigNormGrid{day}(n,:),temp,lagBins);
end

[~,orderTempRL1]=sort(lagsDfof);
dfofSigNormGridSort_TempRL1={};
for n=1:length(dfofSigNormGrid);
    dfofSigNormGridSort_TempRL1{n}=dfofSigNormGrid{n}(orderTempRL1,:);
end

figure
for n=1:length(dfofSigNormGrid);
    subplot(1,length(dfofSigNormGrid),n);
    imagesc(dfofSigNormGridSort_TempRL1{n})
    title(['Day',num2str(n)]);
    axis off
end
tightfig;

save('orderTempRL1.mat','orderTempRL1');
save('dfofSigNormGridSort_TempRL1.mat','dfofSigNormGridSort_TempRL1')
saveas(gcf,'dfofSigNormGridSort_TempRL1.fig');

%% for the matrix above (sorted based on day1) we look at the matrix correlation

matrixCorr=[]; %matrix n to matrix n+1
load('dfofSigNormGridSort_TempRL1.mat');
a=dfofSigNormGridSort_TempRL1;
for n=1:length(a)-1;
    matrixCorr(n,1)=corr2(a{n},a{n+1});
end

figure,plot(matrixCorr,'r');
hold on
plot(matrixCorr,'r.','MarkerSize',10);
[r,p]=corr(matrixCorr(2:end),[1:1:length(matrixCorr)-1]');
title(['Matrix correlation One to Next p=',num2str(p)]);

ylabel('Matrix correlation');
xlabel('Matrix number');

saveas(gcf,'matrixCorrelation_TempRL1.fig');
save('matrixCorrTempRL1.mat','matrixCorr');

%% sort activity according to tempRL

%sort the correlation of the activity with template (two sides adding up
%together)
%load a two side template
load('E:\ID20201118\20201226\2ndloc1\pcaica\cueAnalysisNew_sig\tempRL.mat');
temp=tempRL;

%sort baesd on Dfof
lagsDfof=[];
%use the cue cell code to calculate lags
%sort cells based on the first day in new env
lagBins=60;%shift freely
day=2;
for n=1:size(dfofSigNormGrid{day},1)
[~,lagsDfof(n,1)]=findCorrLag(dfofSigNormGrid{day}(n,:),temp,lagBins);
end

[~,orderTempRL2]=sort(lagsDfof);
dfofSigNormGridSort_TempRL2={};
for n=1:length(dfofSigNormGrid);
    dfofSigNormGridSort_TempRL2{n}=dfofSigNormGrid{n}(orderTempRL2,:);
end

figure
for n=1:length(dfofSigNormGrid);
     subplot(1,length(dfofSigNormGrid),n);
    imagesc(dfofSigNormGridSort_TempRL2{n})
    title(['Day',num2str(n)]);
    axis off
end
tightfig;

save('orderTempRL2.mat','orderTempRL2');
save('dfofSigNormGridSort_TempRL2.mat','dfofSigNormGridSort_TempRL2')
saveas(gcf,'dfofSigNormGridSort_TempRL2.fig');

%% for the matrix above (sorted based on day2) we look at the matrix correlation

matrixCorr=[]; %matrix n to matrix n+1
load('dfofSigNormGridSort_TempRL2.mat');
a=dfofSigNormGridSort_TempRL2;
for n=1:length(a)-1;
    matrixCorr(n,1)=corr2(a{n},a{n+1});
end

figure,plot(matrixCorr,'r');
hold on
plot(matrixCorr,'r.','MarkerSize',10);
[r,p]=corr(matrixCorr(2:end),[1:1:length(matrixCorr)-1]');
title(['Matrix correlation One to Next p=',num2str(p)]);
ylabel('Matrix correlation');
xlabel('Matrix number');

saveas(gcf,'matrixCorrelation_TempRL2.fig');
save('matrixCorrTempRL2.mat','matrixCorr');


%% sort them accoriding to heighest peak
lagsDfof=[];
%use the cue cell code to calculate lags
%sort cells based on the first day in new env
day=1;
for n=1:size(dfofSigNormGrid{day},1)
[~,lagsDfof(n,1)]=max(dfofSigNormGrid{day}(n,:));
end

[~,orderMax1]=sort(lagsDfof);
dfofSigNormGridSort_Max1={};
for n=1:length(dfofSigNormGrid);
    dfofSigNormGridSort_Max1{n}=dfofSigNormGrid{n}(orderMax1,:);  
end

figure
for n=1:length(dfofSigNormGrid);
     subplot(1,length(dfofSigNormGrid),n);
    imagesc(dfofSigNormGridSort_Max1{n})
     title(['Day',num2str(n)]);
    axis off
end
tightfig;

save('orderMax1.mat','orderMax1');
saveas(gcf,'dfofSigNormGridSort_Max1.fig');
save('dfofSigNormGridSort_Max1.mat','dfofSigNormGridSort_Max1')

%% for the matrix above (sorted based on day1) we look at the matrix correlation

matrixCorr=[]; %matrix n to matrix n+1
load('dfofSigNormGridSort_Max1.mat');
a=dfofSigNormGridSort_Max1;
for n=1:length(a)-1;
    matrixCorr(n,1)=corr2(a{n},a{n+1});
end

figure,plot(matrixCorr,'r');
hold on
plot(matrixCorr,'r.','MarkerSize',10);
[r,p]=corr(matrixCorr(2:end),[1:1:length(matrixCorr)-1]');
title(['Matrix correlation One to Next p=',num2str(p)]);
ylabel('Matrix correlation');
xlabel('Matrix number');

saveas(gcf,'matrixCorrelation_Max1.fig');
save('matrixCorrMax1.mat','matrixCorr');

%% sort them accoriding to heighest peak
lagsDfof=[];
%use the cue cell code to calculate lags
%sort cells based on the first day in new env
day=2;
for n=1:size(dfofSigNormGrid{day},1)
[~,lagsDfof(n,1)]=max(dfofSigNormGrid{day}(n,:));
end

[~,orderMax2]=sort(lagsDfof);
dfofSigNormGridSort_Max2={};
for n=1:length(dfofSigNormGrid);
    dfofSigNormGridSort_Max2{n}=dfofSigNormGrid{n}(orderMax2,:);  
end

figure
for n=1:length(dfofSigNormGrid);
    subplot(1,length(dfofSigNormGrid),n);
    imagesc(dfofSigNormGridSort_Max2{n})
     title(['Day',num2str(n)]);
    axis off
end
tightfig;

save('orderMax2.mat','orderMax2');
saveas(gcf,'dfofSigNormGridSort_Max2.fig');
save('dfofSigNormGridSort_Max2.mat','dfofSigNormGridSort_Max2')

%% for the matrix above (sorted based on day2) we look at the matrix correlation

matrixCorr=[]; %matrix n to matrix n+1
load('dfofSigNormGridSort_Max2.mat');
a=dfofSigNormGridSort_Max2;
for n=1:length(a)-1;
    matrixCorr(n,1)=corr2(a{n},a{n+1});
end

figure,plot(matrixCorr,'r');
hold on
plot(matrixCorr,'r.','MarkerSize',10);
[r,p]=corr(matrixCorr(2:end),[1:1:length(matrixCorr)-1]');
title(['Matrix correlation One to Next p=',num2str(p)]);
ylabel('Matrix correlation');
xlabel('Matrix number');

saveas(gcf,'matrixCorrelation_Max2.fig');
save('matrixCorrMax2.mat','matrixCorr');

%% sort them accoriding to heighest peak
lagsDfof=[];
%use the cue cell code to calculate lags
%sort cells based on the first day in new env
day=11;
for n=1:size(dfofSigNormGrid{day},1)
[~,lagsDfof(n,1)]=max(dfofSigNormGrid{day}(n,:));
end

[~,orderMax11]=sort(lagsDfof);
dfofSigNormGridSort_Max11={};
for n=1:length(dfofSigNormGrid);
    dfofSigNormGridSort_Max11{n}=dfofSigNormGrid{n}(orderMax11,:);  
end

figure
for n=1:length(dfofSigNormGrid);
    subplot(1,length(dfofSigNormGrid),n);
    imagesc(dfofSigNormGridSort_Max11{n})
     title(['Day',num2str(n)]);
    axis off
end
tightfig;

save('orderMax11.mat','orderMax11');
saveas(gcf,'dfofSigNormGridSort_Max11.fig');
save('dfofSigNormGridSort_Max11.mat','dfofSigNormGridSort_Max11')

%% for the matrix above (sorted based on day2) we look at the matrix correlation

matrixCorr=[]; %matrix n to matrix n+1
load('dfofSigNormGridSort_Max11.mat');
a=dfofSigNormGridSort_Max11;
for n=1:length(a)-1;
    matrixCorr(n,1)=corr2(a{n},a{n+1});
end

figure,plot(matrixCorr,'r');
hold on
plot(matrixCorr,'r.','MarkerSize',10);
[r,p]=corr(matrixCorr(2:end),[1:1:length(matrixCorr)-1]');
title(['Matrix correlation One to Next p=',num2str(p)]);
ylabel('Matrix correlation');
xlabel('Matrix number');

saveas(gcf,'matrixCorrelation_Max11.fig');
save('matrixCorrMax11.mat','matrixCorr');
%% plot the order of cells
binWidth=7;
bins=[0:binWidth:ceil(size(dfofSigNormGrid{1},1)/binWidth)*binWidth];
[y1,n1]=histc(orderTempRL1,bins);
[y2,n2]=histc(orderTempRL2,bins);

MTempRL=zeros(length(bins)-1,length(bins)-1);
for xx=1:length(bins)-1;
    for yy=1:length(bins)-1;
        A=(n1==xx)&(n2==yy);
        MTempRL(yy,xx)=length(find(A));
    end
end

[y1,n1]=histc(orderMax1,bins);
[y2,n2]=histc(orderMax2,bins);
MMax=zeros(length(bins)-1,length(bins)-1);
for xx=1:length(bins)-1;
    for yy=1:length(bins)-1;
        A=(n1==xx)&(n2==yy);
        MMax(yy,xx)=length(find(A));
    end
end

figure,
subplot(121)
imagesc(MTempRL)
colormap(flipud(gray))
title('TempRL order');
axis equal
axis tight
xlabel('sort day1');
ylabel('sort day2');

subplot(122)
imagesc(MMax)
colormap(flipud(gray))
title('MaxPeak order');
axis equal
axis tight
xlabel('sort day1');
ylabel('sort day2');


saveas(gcf,'cellOrders.fig');

%% now determine the cell activitiy correlation to it own activity across days
load('dfofSigNormGrid.mat');
A=dfofSigNormGrid;
corrIndivi=[];%corr between this day and the next day. Each column is one day and each row is one cell
for n=1:length(A)-1;
    for m=1:size(A{1},1);
        corrIndivi(m,n)=corr(A{n}(m,:)',A{n+1}(m,:)');
    end
end

M=mean(corrIndivi,1);
S=std(corrIndivi,1);
figure,errorbar([1:1:length(M)],M,S,'r-')
[r,p]=corr([1:1:length(M)-1]', M(2:end)');
ylabel('Activity correlation n to n+1')
xlabel('Days');
title(['New Env increase corr p= ',num2str(p)]);
saveas(gcf,'individualCellCorrAcrossDays.fig');

save('corrIndivi.mat','corrIndivi');

%% field distribution
%%
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\indicesAllNewCueThresh.mat');
i=indicesAllNewCueThresh.GridIdx;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\fieldDistri.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\fieldDistriAmp.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\fieldDistriAmpNorm.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\foldersAllIndvFOV.mat');
fieldDistriGrid={};
fieldDistriAmpGrid={};
fieldDistriAmpNormGrid={};
fieldDistriCatGrid={};
fieldDistriAmpSumGrid=[];
fieldDistriAmpNormSumGrid=[];
fieldDistriGridEachBin=[];%in this data, the number of field discovered in each bin was calculated


for n=1:length(foldersAllIndvFOV{1});
    disp(n)
 fieldDistriGrid{n}={};
 fieldDistriCatGrid{n}=[];
fieldDistriAmpGrid{n}=[];
fieldDistriAmpNormGrid{n}=[];

fieldDistriGrid{n}=fieldDistri{n}(i);
fieldDistriCatGrid{n}=cell2mat(fieldDistriGrid{n});

fieldDistriAmpGrid{n}=fieldDistriAmp{n}(i,:);
fieldDistriAmpNormGrid{n}=fieldDistriAmpNorm{n}(i,:);

fieldDistriAmpGrid{n}(isnan(fieldDistriAmpGrid{n}))=0;
fieldDistriAmpNormGrid{n}(isnan(fieldDistriAmpNormGrid{n}))=0;

fieldDistriAmpSumGrid(n,:)=sum(fieldDistriAmpGrid{n},1);
fieldDistriAmpNormSumGrid(n,:)=sum(fieldDistriAmpNormGrid{n},1);

for m=1:max(fieldDistriCatGrid{1});%200bins for 10m track
    fieldDistriGridEachBin(n,m)=length(find(fieldDistriCatGrid{n}==m));
end
end


cd('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\gridCellsNewCueThresh');
save('fieldDistriGrid.mat','fieldDistriGrid');
save('fieldDistriCatGrid.mat','fieldDistriCatGrid');
save('fieldDistriAmpGrid.mat','fieldDistriAmpGrid');
save('fieldDistriAmpNorm.mat','fieldDistriAmpNormGrid');
save('fieldDistriAmpSumGrid.mat','fieldDistriAmpSumGrid');
save('fieldDistriAmpNormSumGrid.mat','fieldDistriAmpNormSumGrid');
save('fieldDistriGridEachBin.mat','fieldDistriGridEachBin');
%% plot them: distribution
load('E:\ID20201118\20201224\loc1\pcaica\cueAnalysisNew_sig\tempRL.mat');
oldTempRL=tempRL;
load('E:\ID20201118\20210103\1stloc1\pcaica\cueAnalysisNew_sig\tempRL');
figure,
for n=1:length(fieldDistriGrid);
    subplot(length(fieldDistriCatGrid),2,2*(n-1)+1);
    plot(fieldDistriGridEachBin(n,:))
    
    hold on
    if n==1;
        plot(oldTempRL*40);
    else
    plot(tempRL*40);
    hold on
       plot([890/5 890/5],[0 40],'r');
    end
    title(['Day',num2str(n)]);
    xlim([0 200]);
%     ylim([0 40]);
    axis off
    title(['Day',num2str(n)]);
end

r1=[];
p1=[];
rnext=[];
pnext=[];

for n=1:length(fieldDistriCatGrid)-1;
    [rnext(n,1),pnext(n,1)]=kstest2(fieldDistriCatGrid{n},fieldDistriCatGrid{n+1});
end

for n=1:length(fieldDistriCatGrid);
    subplot(length(fieldDistriCatGrid),2,2*n);
    [p,x]=ksdensity(fieldDistriCatGrid{n},'width',0.5);
    plot(x,p);
    hold on
    if n==1;
        plot(oldTempRL*max(p));
    else
    plot(tempRL*max(p));
     hold on
       plot([890/5 890/5],[0 max(p)],'r');
    end
    xlim([0 200]);
     ylim([0 0.015]);
     axis off
     [r1(n,1),p1(n,1)]=kstest2(fieldDistriCatGrid{1},fieldDistriCatGrid{n});
     if n>1;
     if rnext(n-1,1)==1;
         title('Sig Dif with previous day','Color','r');
     else
         title('Insig with previous day','Color','k');
     end
     end
     
end
  
tightfig;    
saveas(gcf,'fieldDistriGridbution.fig');
saveas(gcf,'fieldDistriGridbution.jpg');

%look at whether the cell activity is more consistent with the cue template

corrWithTemp=[];
% binnedOldTempRL=oldTempRL(1:2:length(oldTempRL)-1)+oldTempRL(1:2:length(oldTempRL));
% binnedTempRL=tempRL(1:2:length(tempRL)-1)+tempRL(1:2:length(tempRL));

for n=1:length(fieldDistriCatGrid);
    if n==1;
    corrWithTemp(n,1)=corr(fieldDistriGridEachBin(n,:)',oldTempRL);
    else
        corrWithTemp(n,1)=corr(fieldDistriGridEachBin(n,:)',tempRL);
    end
end
figure,plot(corrWithTemp);
xlabel('days')
ylabel('correlation to temp');
[r,p]=corr(corrWithTemp(2:end),[1:1:10]');
title(['p = ',num2str(p)]);
saveas(gcf,'distriCorrToTemp.fig');
saveas(gcf,'distriCorrToTemp.jpg');
%% plot them: amplitude
corrAmp=[];%correlation to the next activity norm adding up together

for n=1:length(fieldDistri)-1;
    corrAmp(n)=corr(fieldDistriAmpSumGrid(n,:)',fieldDistriAmpSumGrid(n+1,:)');
end

corrAmpNorm=[];%correlation to the next activity norm adding up together

for n=1:length(fieldDistri)-1;
      corrAmpNorm(n)=corr(fieldDistriAmpNormSumGrid(n,:)',fieldDistriAmpNormSumGrid(n+1,:)');
end

figure;

for n=1:length(fieldDistri);
      subplot(length(fieldDistri),2,2*(n-1)+1);
      plot(fieldDistriAmpSumGrid(n,:));
      hold on
      if n==1;
         plot(oldTempRL*max(fieldDistriAmpSumGrid(1,:)));  
      else
       plot(tempRL*max(fieldDistriAmpSumGrid(2,:)));
       hold on
       plot([890/5 890/5],[0 max(fieldDistriAmpSumGrid(2,:))],'r');
      end
     if n<length(fieldDistri);
               title(['corr to next = ',num2str(corrAmp(n))])
     end        
      axis off
    
      subplot(length(fieldDistri),2,2*n);
        plot(fieldDistriAmpNormSumGrid(n,:));
         hold on
      if n==1;
         plot(oldTempRL*max(fieldDistriAmpNormSumGrid(1,:)));  
      else
       plot(tempRL*max(fieldDistriAmpNormSumGrid(2,:)));
       hold on
       plot([890/5 890/5],[0 max(fieldDistriAmpNormSumGrid(2,:))],'r');
      end
        if n<length(fieldDistri);
               title(['corr to next = ',num2str(corrAmpNorm(n))])
            
        end
                
axis off
end
tightfig;

saveas(gcf,'fieldDistriAmpGrid.fig');
saveas(gcf,'fieldDistriAmpGrid.jpg');


%%
x=[];
p=[];

for n=1:length(fieldDistriGrid);   
[x(n,:),p(n,:)]=ksdensity(fieldDistriCatGrid{n},'width',5);
end


%another way is to count each bin

figure,

subplot(411)
plot(p(1,:)*5-2.5,mean(x(4:11,:),1)-mean(x(2:3,:),1))
hold on
plot([2.5:5:1000],tempRL*0.0009);
hold on
plot([890 890],[0 0.0009],'r');
xlim([-15*5 215*5]);
title('Distribution ksdensity Day (4-11)-(2-3)');


subplot(412)
plot([1:1:200],mean(fieldDistriGridEachBin(4:11,:),1)-mean(fieldDistriGridEachBin(2:3,:),1))
hold on
plot([1:1:200],tempRL*20);
hold on
plot([890/5 890/5],[0 20],'r');
xlim([0 200]);
title('Distribution each bin Day (4-11)-(2-3)');


subplot(413)
plot([1:1:size(fieldDistriAmpSumGrid,2)]*5-2.5,mean(fieldDistriAmpSumGrid(4:11,:),1)-mean(fieldDistriAmpSumGrid(2:3,:),1))
hold on
plot([2.5:5:1000],tempRL*2)
hold on
plot([890 890],[0 2],'r');
xlim([-15*5 215*5]);
title('Amplitude Day (4-11)-(2-3)');

subplot(414)
plot([1:1:size(fieldDistriAmpNormSumGrid,2)]*5-2.5,mean(fieldDistriAmpNormSumGrid(4:11,:),1)-mean(fieldDistriAmpNormSumGrid(2:3,:),1))
hold on
plot([2.5:5:1000],tempRL*4)
hold on
plot([890 890],[0 4],'r');
xlim([-15*5 215*5]);
title('Amplitude norm Day (4-11)-(2-3)');
saveas(gcf,'disriAmpBeforeAfterLearning.fig');
saveas(gcf,'disriAmpBeforeAfterLearning.jpg');

%% test whether any of these changes are significant
%working on field distribution in each bin
load('fieldDistriGridEachBin.mat');
a=fieldDistriGridEachBin;
aa=mean(a(4:11,:),1)-mean(a(2:3,:),1);
sig=[];
ss=[];
%randomly permute numbers in each row
NShuffle=100;
for n=1:NShuffle;
    s=[];
    for m=1:size(a,2);
        ac=a(:,m); 
        s(:,m)=ac(randperm(length(ac)));%permute across days
    end
    ss(n,:)=mean(s(4:11,:),1)-mean(s(2:3,:),1);

end
d=[];
for n=1:length(aa);
    d(n)=length(find(ss(:,n)<aa(n)))/NShuffle';
    if d(n)>0.95 || d(n)<0.05;
        sig(n)=1;
    else
        sig(n)=0;
    end
end
sigIdxFieldDistri=find(sig==1);
figure,plot(aa,'k');
hold on
plot(sigIdxFieldDistri,aa(sigIdxFieldDistri),'r.','MarkerSize',10);

saveas(gcf,'disriAmpBeforeAfterLearning_sig.fig');
save('sigIdxFieldDistri.mat','sigIdxFieldDistri');

%%
figure,

subplot(411)
plot(p(1,:)*5-2.5,mean(x(4:11,:),1)-mean(x(1,:),1))
hold on
plot([2.5:5:1000],tempRL*0.0009);
hold on
plot([890 890],[0 0.0009],'r');
xlim([-15*5 215*5]);
title('Distribution ksdensity Day (4-11)-(1)');


subplot(412)
plot([1:1:200],mean(fieldDistriGridEachBin(4:11,:),1)-mean(fieldDistriGridEachBin(1,:),1))
hold on
plot([1:1:200],tempRL*20);
hold on
plot([890/5 890/5],[0 20],'r');
xlim([0 200]);
title('Distribution each bin Day (4-11)-(1)');


subplot(413)
plot([1:1:size(fieldDistriAmpSumGrid,2)]*5-2.5,mean(fieldDistriAmpSumGrid(4:11,:),1)-mean(fieldDistriAmpSumGrid(1,:),1))
hold on
plot([2.5:5:1000],tempRL*2)
hold on
plot([890 890],[0 2],'r');
xlim([-15*5 215*5]);
title('Amplitude Day (4-11)-(1)');

subplot(414)
plot([1:1:size(fieldDistriAmpNormSumGrid,2)]*5-2.5,mean(fieldDistriAmpNormSumGrid(4:11,:),1)-mean(fieldDistriAmpNormSumGrid(1,:),1))
hold on
plot([2.5:5:1000],tempRL*4)
hold on
plot([890 890],[0 4],'r');
xlim([-15*5 215*5]);
title('Amplitude norm Day (4-11)-(1)');
saveas(gcf,'disriAmpOldNewEnv.fig');
saveas(gcf,'disriAmpOldNewEnv.jpg');

%% plot for the arrangement of fields: to the template all infields are 1, and nonfields are zero

load('E:\learningAnalysis\summaryManyMice_includingOldEnv\indicesAllNewCueThresh.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnv\fieldLocOnly.mat');

i=indicesAll.GridIdx;

fieldLocOnlyGrid={};

for n=1:length(fieldDistri);
    fieldLocOnlyGrid{n}=fieldLocOnly{n}(i,:);%each cell is has one row
    
end
save('fieldLocOnlyGrid.mat','fieldLocOnlyGrid');
%% order field loc only (0 and 1) activity based on different templates
%% sort activity according to tempRL

%sort the correlation of the activity with template (two sides adding up
%together)
%load a two side template
load('E:\ID20201118\20201224\loc1\pcaica\cueAnalysisNew_sig\tempRL.mat');
temp=tempRL;

%sort baesd on Dfof
lagsDfof=[];
%use the cue cell code to calculate lags
%sort cells based on the first day in new env
lagBins=60;
day=1;
for n=1:size(fieldLocOnlyGrid{day},1)
[~,lagsDfof(n,1)]=findCorrLag(fieldLocOnlyGrid{day}(n,:),temp,lagBins);
end

[~,orderLocTempRL1]=sort(lagsDfof);
fieldLocOnlySort_TempRL1={};
for n=1:length(fieldLocOnlyGrid);
    fieldLocOnlySort_TempRL1{n}=fieldLocOnlyGrid{n}(orderLocTempRL1,:);
end

figure
for n=1:length(fieldLocOnlyGrid);
    subplot(1,length(fieldLocOnlyGrid),n);
    imagesc(fieldLocOnlySort_TempRL1{n})
    colormap(flipud(gray));
    title(['Day',num2str(n)]);
    axis off
end
tightfig;

save('orderLocTempRL1.mat','orderLocTempRL1');
save('fieldLocOnlySort_TempRL1.mat','fieldLocOnlySort_TempRL1')
saveas(gcf,'fieldLocOnlySort_TempRL1.fig');

%% for the matrix above (sorted based on day1) we look at the matrix correlation

matrixCorrLocOnly=[]; %matrix n to matrix n+1
load('fieldLocOnlySort_TempRL1.mat');
a=fieldLocOnlySort_TempRL1;
for n=1:length(a)-1;
    matrixCorrLocOnly(n,1)=corr2(a{n},a{n+1});
end

figure,plot(matrixCorrLocOnly,'r');
hold on
plot(matrixCorrLocOnly,'r.','MarkerSize',10);
title('Matrix correlation One to Next');
ylabel('Matrix correlation');
xlabel('Matrix number');

% [r,p]=corr(matrixCorrLocOnly(2:end),[1:1:9]')
% r =
% 
%     0.7725
% 
% 
% p =
% 
%     0.0147
saveas(gcf,'matrixCorrLocOnlyelationLocOnly_TempRL1.fig');
save('matrixCorrLocOnlyTempRL1.mat','matrixCorrLocOnly');

%% sort activity according to tempRL

%sort the correlation of the activity with template (two sides adding up
%together)
%load a two side template
load('E:\ID20201118\20201226\2ndloc1\pcaica\cueAnalysisNew_sig\tempRL.mat');
temp=tempRL;

%sort baesd on Dfof
lagsDfof=[];
%use the cue cell code to calculate lags
%sort cells based on the first day in new env
lagBins=60;
day=2;
for n=1:size(fieldLocOnlyGrid{day},1)
[~,lagsDfof(n,1)]=findCorrLag(fieldLocOnlyGrid{day}(n,:),temp,lagBins);
end

[~,orderLocTempRL2]=sort(lagsDfof);
fieldLocOnlySort_TempRL2={};
for n=1:length(fieldLocOnlyGrid);
    fieldLocOnlySort_TempRL2{n}=fieldLocOnlyGrid{n}(orderLocTempRL2,:);
end

figure
for n=1:length(fieldLocOnlyGrid);
     subplot(1,length(fieldLocOnlyGrid),n);
    imagesc(fieldLocOnlySort_TempRL2{n})
     colormap(flipud(gray));
    title(['Day',num2str(n)]);
    axis off
end
tightfig;

save('orderLocTempRL2.mat','orderLocTempRL2');
save('fieldLocOnlySort_TempRL2.mat','fieldLocOnlySort_TempRL2')
saveas(gcf,'fieldLocOnlySort_TempRL2.fig');
%% for the matrix above (sorted based on day2) we look at the matrix correlation

matrixCorrLocOnly=[]; %matrix n to matrix n+1
load('fieldLocOnlySort_TempRL2.mat');
a=fieldLocOnlySort_TempRL2;
for n=1:length(a)-1;
    matrixCorrLocOnly(n,1)=corr2(a{n},a{n+1});
end

figure,plot(matrixCorrLocOnly,'r');
hold on
plot(matrixCorrLocOnly,'r.','MarkerSize',10);
title('Matrix correlation One to Next');
ylabel('Matrix correlation');
xlabel('Matrix number');

saveas(gcf,'matrixCorrLocOnlyelationLocOnly_TempRL2.fig');
save('matrixCorrLocOnlyTempRL2.mat','matrixCorrLocOnly');

%% sort activity according to tempRL

%sort the correlation of the activity with template (two sides adding up
%together)
%load a two side template
load('E:\ID20201118\20201226\2ndloc1\pcaica\cueAnalysisNew_sig\tempRL.mat');
temp=tempRL;

%sort baesd on Dfof
lagsDfof=[];
%use the cue cell code to calculate lags
%sort cells based on the first day in new env
lagBins=60;
day=11;
for n=1:size(fieldLocOnlyGrid{day},1)
[~,lagsDfof(n,1)]=findCorrLag(fieldLocOnlyGrid{day}(n,:),temp,lagBins);
end

[~,orderLocTempRL11]=sort(lagsDfof);
fieldLocOnlySort_TempRL11={};
for n=1:length(fieldLocOnlyGrid);
    fieldLocOnlySort_TempRL11{n}=fieldLocOnlyGrid{n}(orderLocTempRL11,:);
end

figure
for n=1:length(fieldLocOnlyGrid);
     subplot(1,length(fieldLocOnlyGrid),n);
    imagesc(fieldLocOnlySort_TempRL11{n})
     colormap(flipud(gray));
    title(['Day',num2str(n)]);
    axis off
end
tightfig;

save('orderLocTempRL11.mat','orderLocTempRL11');
save('fieldLocOnlySort_TempRL11.mat','fieldLocOnlySort_TempRL11')
saveas(gcf,'fieldLocOnlySort_TempRL11.fig');
%% for the matrix above (sorted based on day2) we look at the matrix correlation

matrixCorrLocOnly=[]; %matrix n to matrix n+1
load('fieldLocOnlySort_TempRL11.mat');
a=fieldLocOnlySort_TempRL11;
for n=1:length(a)-1;
    matrixCorrLocOnly(n,1)=corr2(a{n},a{n+1});
end

figure,plot(matrixCorrLocOnly,'r');
hold on
plot(matrixCorrLocOnly,'r.','MarkerSize',10);
title('Matrix correlation One to Next');
ylabel('Matrix correlation');
xlabel('Matrix number');

saveas(gcf,'matrixCorrLocOnlyelationLocOnly_TempRL11.fig');
save('matrixCorrLocOnlyTempRL11.mat','matrixCorrLocOnly');

%% plot the order of cells for old and new env
binWidth=7;
bins=[0:binWidth:ceil(size(fieldLocOnlyGrid{1},1)/binWidth)*binWidth];
%for old and new
[y1,n1]=histc(orderLocTempRL1,bins);
[y2,n2]=histc(orderLocTempRL2,bins);

MTempRL=zeros(length(bins)-1,length(bins)-1);
for xx=1:length(bins)-1;
    for yy=1:length(bins)-1;
        A=(n1==xx)&(n2==yy);
        MTempRL(yy,xx)=length(find(A));
    end
end

[y1,n1]=histc(orderLocTempRL1,bins);
[y2,n2]=histc(orderLocTempRL11,bins);
MTempRL2=zeros(length(bins)-1,length(bins)-1);
for xx=1:length(bins)-1;
    
    for yy=1:length(bins)-1;
        A=(n1==xx)&(n2==yy);
        MTempRL2(yy,xx)=length(find(A));
    end
end

figure,
subplot(121)
imagesc(MTempRL)
colormap(flipud(gray))
% colormap(parula)
title('TempRL order: old and new');
axis equal
axis tight
xlabel('sort old');
ylabel('sort day1new');

subplot(122)
imagesc(MTempRL2)
colormap(flipud(gray))
% colormap(cool)
title('TempRL order: new1 and new11');
axis equal
axis tight
xlabel('sort day1new');
ylabel('sort day10new');
saveas(gcf,'cellOrdersFieldLoc.fig');
saveas(gcf,'cellOrdersFieldLoc.jpg');

%% now determine the cell activitiy correlation to it own activity across days
A=fieldLocOnlyGrid;
corrLocIndivi=[];%corr between this day and the next day. Each column is one day and each row is one cell
for n=1:length(A)-1;
    for m=1:size(A{1},1);
        corrLocIndivi(m,n)=corr(A{n}(m,:)',A{n+1}(m,:)');
    end
    
end

M=[];
S=[];
for n=1:length(A)-1;
    B=corrLocIndivi(:,n);
    B=B(~isnan(B));
    M(n)=mean(B);
    S(n)=std(B);
end
figure,errorbar([1:1:length(M)],M,S,'r-')
[r,p]=corr([1:1:length(M)-1]', M(2:end)');
ylabel('Activity correlation n to n+1')
xlabel('Days');
title(['New Env increase corr p= ',num2str(p)]);
saveas(gcf,'individualCellCorrAcrossDaysLocOnly.fig');

%% look at the cell's max field width and spacing and in out field ratio

cd ..\
load('fieldInfoUsingAllCells.mat')
load('indicesAllNewCueThresh.mat');
cd('gridCellsNewCueThresh');
i=indicesAllNewCueThresh.GridIdx;

minSpacingGrid=minSpacing(i,:);
maxWidthGrid=maxWidth(i,:);
ratioInOtherGrid=ratioInOther(i,:);

%normalize all numbers based on the data on the first day in new env
minSpacingGridNorm=minSpacingNorm(i,:);
maxWidthGridNorm=maxWidthNorm(i,:);
ratioInOtherGridNorm=ratioInOtherNorm(i,:);

save('fieldInfoUsingAllCellsGrid.mat','minSpacingGrid','maxWidthGrid','ratioInOtherGrid','minSpacingGridNorm','maxWidthGridNorm','ratioInOtherGridNorm');
%%

figure,
subplot(231);
a=minSpacingGrid;
M=[];
S=[];
for n=1:size(a,2);
    b=a(:,n);
    b=b(~isnan(b));
     b=b(~isinf(b));
M(n)=mean(b);
S(n)=std(b)/sqrt(length(b));
end
errorbar([1:1:length(M)],M,S,'r');
hold on
plot([1:1:length(M)],M,'r.');
[r,p]=corr([1:1:length(M)-1]',M(2:end)');
title(['minSpacing p=',num2str(p)])
xlabel('spacing cm');
ylabel('cell number');

            
subplot(232);
a=maxWidthGrid;
M=[];
S=[];
for n=1:size(a,2);
    b=a(:,n);
    b=b(~isnan(b));
     b=b(~isinf(b));
M(n)=mean(b);
S(n)=std(b)/sqrt(length(b));
end
errorbar([1:1:length(M)],M,S,'r');
hold on
plot([1:1:length(M)],M,'r.');
[r,p]=corr([1:1:length(M)-1]',M(2:end)');
title(['maxWidth p=',num2str(p)])
xlabel('width cm');
ylabel('cell number');

subplot(233);
a=ratioInOtherGrid;
M=[];
S=[];
for n=1:size(a,2);
    b=a(:,n);
    b=b(~isnan(b));
      b=b(~isinf(b));
      b=b(b<500);
M(n)=mean(b);
S(n)=std(b)/sqrt(length(b));
end
errorbar([1:1:length(M)],M,S,'r');
hold on
plot([1:1:length(M)],M,'r.');
[r,p]=corr([1:1:length(M)-1]',M(2:end)');
title(['ratioInVSOther p=',num2str(p)])
xlabel('ratio');
ylabel('cell number');
saveas(gcf,'fieldInfoUsingAllCells.fig');

subplot(234);
a=minSpacingGridNorm;
M=[];
S=[];
for n=1:size(a,2);
    b=a(:,n);
    b=b(~isnan(b));
     b=b(~isinf(b));
M(n)=mean(b);
S(n)=std(b)/sqrt(length(b));
end
errorbar([1:1:length(M)],M,S,'r');
hold on
plot([1:1:length(M)],M,'r.');
[r,p]=corr([1:1:length(M)-1]',M(2:end)');
title(['minSpacing norm p=',num2str(p)])
xlabel('spacing cm');
ylabel('cell number');

            
subplot(235);
a=maxWidthGridNorm;
M=[];
S=[];
for n=1:size(a,2);
    b=a(:,n);
    b=b(~isnan(b));
     b=b(~isinf(b));
M(n)=mean(b);
S(n)=std(b)/sqrt(length(b));
end
errorbar([1:1:length(M)],M,S,'r');
hold on
plot([1:1:length(M)],M,'r.');
[r,p]=corr([1:1:length(M)-1]',M(2:end)');
title(['maxWidth norm p=',num2str(p)])
xlabel('width cm');
ylabel('cell number');

subplot(236);
a=ratioInOtherGridNorm;
M=[];
S=[];
for n=1:size(a,2);
    b=a(:,n);
    b=b(~isnan(b));
      b=b(~isinf(b));
      b=b(b<500);
M(n)=mean(b);
S(n)=std(b)/sqrt(length(b));
end
errorbar([1:1:length(M)],M,S,'r');
hold on
plot([1:1:length(M)],M,'r.');
[r,p]=corr([1:1:length(M)-1]',M(2:end)');
title(['ratioInVSOther norm p=',num2str(p)])
xlabel('ratio');
ylabel('cell number');

%
saveas(gcf,'fieldInfoUsingAllCellsGrid.fig');
saveas(gcf,'fieldInfoUsingAllCellsGrid.jpg');

%% look at the cell's max field width and spacing and in out field ratio
cd ..\
load('fieldInfoUsingAllCellsCorrected.mat')
load('indicesAllNewCueThresh.mat');
cd('gridCellsNewCueThresh');
i=indicesAllNewCueThresh.GridIdx;

minSpacingGrid=minSpacing(i,:);
maxWidthGrid=maxWidth(i,:);
ratioInOtherGrid=ratioInOther(i,:);

%normalize all numbers based on the data on the first day in new env
minSpacingGridNorm=minSpacingNorm(i,:);
maxWidthGridNorm=maxWidthNorm(i,:);
ratioInOtherGridNorm=ratioInOtherNorm(i,:);

save('fieldInfoUsingAllCellsCorrectedGrid.mat','minSpacingGrid','maxWidthGrid','ratioInOtherGrid','minSpacingGridNorm','maxWidthGridNorm','ratioInOtherGridNorm');
%%

figure,
subplot(231);
a=minSpacingGrid;
M=[];
S=[];
for n=1:size(a,2);
    b=a(:,n);
    b=b(~isnan(b));
     b=b(~isinf(b));
M(n)=mean(b);
S(n)=std(b)/sqrt(length(b));
end
errorbar([1:1:length(M)],M,S,'r');
hold on
plot([1:1:length(M)],M,'r.');
[r,p]=corr([1:1:length(M)-1]',M(2:end)');
title(['minSpacing p=',num2str(p)])
xlabel('spacing cm');
ylabel('cell number');

            
subplot(232);
a=maxWidthGrid;
M=[];
S=[];
for n=1:size(a,2);
    b=a(:,n);
    b=b(~isnan(b));
     b=b(~isinf(b));
M(n)=mean(b);
S(n)=std(b)/sqrt(length(b));
end
errorbar([1:1:length(M)],M,S,'r');
hold on
plot([1:1:length(M)],M,'r.');
[r,p]=corr([1:1:length(M)-1]',M(2:end)');
title(['maxWidth p=',num2str(p)])
xlabel('width cm');
ylabel('cell number');

subplot(233);
a=ratioInOtherGrid;
M=[];
S=[];
for n=1:size(a,2);
    b=a(:,n);
    b=b(~isnan(b));
      b=b(~isinf(b));
      b=b(b<500);
M(n)=mean(b);
S(n)=std(b)/sqrt(length(b));
end
errorbar([1:1:length(M)],M,S,'r');
hold on
plot([1:1:length(M)],M,'r.');
[r,p]=corr([1:1:length(M)-1]',M(2:end)');
title(['ratioInVSOther p=',num2str(p)])
xlabel('ratio');
ylabel('cell number');
saveas(gcf,'fieldInfoUsingAllCells.fig');

subplot(234);
a=minSpacingGridNorm;
M=[];
S=[];
for n=1:size(a,2);
    b=a(:,n);
    b=b(~isnan(b));
     b=b(~isinf(b));
M(n)=mean(b);
S(n)=std(b)/sqrt(length(b));
end
errorbar([1:1:length(M)],M,S,'r');
hold on
plot([1:1:length(M)],M,'r.');
[r,p]=corr([1:1:length(M)-1]',M(2:end)');
title(['minSpacing norm p=',num2str(p)])
xlabel('spacing cm');
ylabel('cell number');

            
subplot(235);
a=maxWidthGridNorm;
M=[];
S=[];
for n=1:size(a,2);
    b=a(:,n);
    b=b(~isnan(b));
     b=b(~isinf(b));
M(n)=mean(b);
S(n)=std(b)/sqrt(length(b));
end
errorbar([1:1:length(M)],M,S,'r');
hold on
plot([1:1:length(M)],M,'r.');
[r,p]=corr([1:1:length(M)-1]',M(2:end)');
title(['maxWidth norm p=',num2str(p)])
xlabel('width cm');
ylabel('cell number');

subplot(236);
a=ratioInOtherGridNorm;
M=[];
S=[];
for n=1:size(a,2);
    b=a(:,n);
    b=b(~isnan(b));
      b=b(~isinf(b));
      b=b(b<500);
M(n)=mean(b);
S(n)=std(b)/sqrt(length(b));
end
errorbar([1:1:length(M)],M,S,'r');
hold on
plot([1:1:length(M)],M,'r.');
[r,p]=corr([1:1:length(M)-1]',M(2:end)');
title(['ratioInVSOther norm p=',num2str(p)])
xlabel('ratio');
ylabel('cell number');

%
saveas(gcf,'fieldInfoUsingAllCellsCorrectedGrid.fig');
saveas(gcf,'fieldInfoUsingAllCellsCorrectedGrid.jpg');


%%
%% look at the cell's max field width and spacing and in out field ratio

cd ..\
load('fieldInfoUsingAllCellsMean.mat')
load('indicesAllNewCueThresh.mat');
cd('gridCellsNewCueThresh');
i=indicesAllNewCueThresh.GridIdx;

meanSpacingGrid=meanSpacing(i,:);
meanWidthGrid=meanWidth(i,:);
ratioOtherInGrid=ratioOtherIn(i,:);

%normalize all numbers based on the data on the first day in new env
meanSpacingGridNorm=meanSpacingNorm(i,:);
meanWidthGridNorm=meanWidthNorm(i,:);
ratioOtherInGridNorm=ratioOtherInNorm(i,:);

save('fieldInfoUsingAllCellsGridMean.mat','meanSpacingGrid','meanWidthGrid','ratioOtherInGrid','meanSpacingGridNorm','meanWidthGridNorm','ratioOtherInGridNorm');
%%

figure,
subplot(231);
a=meanSpacingGrid;
M=[];
S=[];
for n=1:size(a,2);
    b=a(:,n);
    b=b(~isnan(b));
     b=b(~isinf(b));
M(n)=mean(b);
S(n)=std(b)/sqrt(length(b));
end
errorbar([1:1:length(M)],M,S,'r');
hold on
plot([1:1:length(M)],M,'r.');
[r,p]=corr([1:1:length(M)-1]',M(2:end)');
title(['meanSpacing p=',num2str(p)])
xlabel('spacing cm');
ylabel('cell number');

            
subplot(232);
a=meanWidthGrid;
M=[];
S=[];
for n=1:size(a,2);
    b=a(:,n);
    b=b(~isnan(b));
     b=b(~isinf(b));
M(n)=mean(b);
S(n)=std(b)/sqrt(length(b));
end
errorbar([1:1:length(M)],M,S,'r');
hold on
plot([1:1:length(M)],M,'r.');
[r,p]=corr([1:1:length(M)-1]',M(2:end)');
title(['meanWidth p=',num2str(p)])
xlabel('width cm');
ylabel('cell number');

subplot(233);
a=ratioOtherInGrid;
M=[];
S=[];
for n=1:size(a,2);
    b=a(:,n);
    b=b(~isnan(b));
      b=b(~isinf(b));
      b=b(b<500);
M(n)=mean(b);
S(n)=std(b)/sqrt(length(b));
end
errorbar([1:1:length(M)],M,S,'r');
hold on
plot([1:1:length(M)],M,'r.');
[r,p]=corr([1:1:length(M)-1]',M(2:end)');
title(['ratioOtherVSIn p=',num2str(p)])
xlabel('ratio');
ylabel('cell number');
saveas(gcf,'fieldInfoUsingAllCells.fig');

subplot(234);
a=meanSpacingGridNorm;
M=[];
S=[];
for n=1:size(a,2);
    b=a(:,n);
    b=b(~isnan(b));
     b=b(~isinf(b));
M(n)=mean(b);
S(n)=std(b)/sqrt(length(b));
end
errorbar([1:1:length(M)],M,S,'r');
hold on
plot([1:1:length(M)],M,'r.');
[r,p]=corr([1:1:length(M)-1]',M(2:end)');
title(['meanSpacing norm p=',num2str(p)])
xlabel('spacing cm');
ylabel('cell number');

            
subplot(235);
a=meanWidthGridNorm;
M=[];
S=[];
for n=1:size(a,2);
    b=a(:,n);
    b=b(~isnan(b));
     b=b(~isinf(b));
M(n)=mean(b);
S(n)=std(b)/sqrt(length(b));
end
errorbar([1:1:length(M)],M,S,'r');
hold on
plot([1:1:length(M)],M,'r.');
[r,p]=corr([1:1:length(M)-1]',M(2:end)');
title(['meanWidth norm p=',num2str(p)])
xlabel('width cm');
ylabel('cell number');

subplot(236);
a=ratioOtherInGridNorm;
M=[];
S=[];
for n=1:size(a,2);
    b=a(:,n);
    b=b(~isnan(b));
      b=b(~isinf(b));
      b=b(b<500);
M(n)=mean(b);
S(n)=std(b)/sqrt(length(b));
end
errorbar([1:1:length(M)],M,S,'r');
hold on
plot([1:1:length(M)],M,'r.');
[r,p]=corr([1:1:length(M)-1]',M(2:end)');
title(['ratioOtherVSIn norm p=',num2str(p)])
xlabel('ratio');
ylabel('cell number');

%
saveas(gcf,'fieldInfoUsingAllCellsGridMean.fig');
saveas(gcf,'fieldInfoUsingAllCellsGridMean.jpg');

%% look at the cell's max field width and spacing and in out field ratio
cd ..\
load('fieldInfoUsingAllCellsCorrectedMean.mat')
load('indicesAllNewCueThresh.mat');
cd('gridCellsNewCueThresh');
i=indicesAllNewCueThresh.GridIdx;

meanSpacingGrid=meanSpacing(i,:);
meanWidthGrid=meanWidth(i,:);
ratioOtherInGrid=ratioOtherIn(i,:);

%normalize all numbers based on the data on the first day in new env
meanSpacingGridNorm=meanSpacingNorm(i,:);
meanWidthGridNorm=meanWidthNorm(i,:);
ratioOtherInGridNorm=ratioOtherInNorm(i,:);

save('fieldInfoUsingAllCellsCorrectedGridMean.mat','meanSpacingGrid','meanWidthGrid','ratioOtherInGrid','meanSpacingGridNorm','meanWidthGridNorm','ratioOtherInGridNorm');
%%

figure,
subplot(231);
a=meanSpacingGrid;
M=[];
S=[];
for n=1:size(a,2);
    b=a(:,n);
    b=b(~isnan(b));
     b=b(~isinf(b));
M(n)=mean(b);
S(n)=std(b)/sqrt(length(b));
end
errorbar([1:1:length(M)],M,S,'r');
hold on
plot([1:1:length(M)],M,'r.');
[r,p]=corr([1:1:length(M)-1]',M(2:end)');
title(['meanSpacing p=',num2str(p)])
xlabel('spacing cm');
ylabel('cell number');

            
subplot(232);
a=meanWidthGrid;
M=[];
S=[];
for n=1:size(a,2);
    b=a(:,n);
    b=b(~isnan(b));
     b=b(~isinf(b));
M(n)=mean(b);
S(n)=std(b)/sqrt(length(b));
end
errorbar([1:1:length(M)],M,S,'r');
hold on
plot([1:1:length(M)],M,'r.');
[r,p]=corr([1:1:length(M)-1]',M(2:end)');
title(['meanWidth p=',num2str(p)])
xlabel('width cm');
ylabel('cell number');

subplot(233);
a=ratioOtherInGrid;
M=[];
S=[];
for n=1:size(a,2);
    b=a(:,n);
    b=b(~isnan(b));
      b=b(~isinf(b));
      b=b(b<500);
M(n)=mean(b);
S(n)=std(b)/sqrt(length(b));
end
errorbar([1:1:length(M)],M,S,'r');
hold on
plot([1:1:length(M)],M,'r.');
[r,p]=corr([1:1:length(M)-1]',M(2:end)');
title(['ratioOtherVSIn p=',num2str(p)])
xlabel('ratio');
ylabel('cell number');
saveas(gcf,'fieldInfoUsingAllCells.fig');

subplot(234);
a=meanSpacingGridNorm;
M=[];
S=[];
for n=1:size(a,2);
    b=a(:,n);
    b=b(~isnan(b));
     b=b(~isinf(b));
M(n)=mean(b);
S(n)=std(b)/sqrt(length(b));
end
errorbar([1:1:length(M)],M,S,'r');
hold on
plot([1:1:length(M)],M,'r.');
[r,p]=corr([1:1:length(M)-1]',M(2:end)');
title(['meanSpacing norm p=',num2str(p)])
xlabel('spacing cm');
ylabel('cell number');

            
subplot(235);
a=meanWidthGridNorm;
M=[];
S=[];
for n=1:size(a,2);
    b=a(:,n);
    b=b(~isnan(b));
     b=b(~isinf(b));
M(n)=mean(b);
S(n)=std(b)/sqrt(length(b));
end
errorbar([1:1:length(M)],M,S,'r');
hold on
plot([1:1:length(M)],M,'r.');
[r,p]=corr([1:1:length(M)-1]',M(2:end)');
title(['meanWidth norm p=',num2str(p)])
xlabel('width cm');
ylabel('cell number');

subplot(236);
a=ratioOtherInGridNorm;
M=[];
S=[];
for n=1:size(a,2);
    b=a(:,n);
    b=b(~isnan(b));
      b=b(~isinf(b));
      b=b(b<500);
M(n)=mean(b);
S(n)=std(b)/sqrt(length(b));
end
errorbar([1:1:length(M)],M,S,'r');
hold on
plot([1:1:length(M)],M,'r.');
[r,p]=corr([1:1:length(M)-1]',M(2:end)');
title(['ratioOtherVSIn norm p=',num2str(p)])
xlabel('ratio');
ylabel('cell number');

%
saveas(gcf,'fieldInfoUsingAllCellsCorrectedGridMean.fig');
saveas(gcf,'fieldInfoUsingAllCellsCorrectedGridMean.jpg');

%% in out field ratio
cd ..\
load('fieldRatioInOut.mat')
load('indicesAllNewCueThresh.mat');
cd('gridCellsNewCueThresh');
i=indicesAllNewCueThresh.GridIdx;

ratioInOutGrid=ratioInOut(i,:);
ratioInOutCorrectedGrid=ratioInOutCorrected(i,:);

ratioInOutNormGrid=ratioInOutNorm(i,:);
ratioInOutCorrectedNormGrid=ratioInOutCorrectedNorm(i,:);


save('fieldRatioInOutGrid.mat','ratioInOutGrid','ratioInOutCorrectedGrid','ratioInOutNormGrid','ratioInOutCorrectedNormGrid');


%%
figure,


subplot(221);
a=ratioInOutGrid;
M=[];
S=[];
for n=1:size(a,2);
    b=a(:,n);
    b=b(~isnan(b));
     b=b(~isinf(b));
M(n)=mean(b);
S(n)=std(b)/sqrt(length(b));
end
errorbar([1:1:length(M)],M,S,'r');
hold on
plot([1:1:length(M)],M,'r.');
[r,p]=corr([1:1:length(M)-1]',M(2:end)');
ylabel('ratio in out');
xlabel('day');
title('ratio in out');


subplot(222);
a=ratioInOutCorrectedGrid;
M=[];
S=[];
for n=1:size(a,2);
    b=a(:,n);
    b=b(~isnan(b));
     b=b(~isinf(b));
M(n)=mean(b);
S(n)=std(b)/sqrt(length(b));
end
errorbar([1:1:length(M)],M,S,'r');
hold on
plot([1:1:length(M)],M,'r.');
[r,p]=corr([1:1:length(M)-1]',M(2:end)');
ylabel('ratio in out corrected');
xlabel('day');
title('ratio in out corrcted');



subplot(223);
a=ratioInOutNormGrid;
M=[];
S=[];
for n=1:size(a,2);
    b=a(:,n);
    b=b(~isnan(b));
     b=b(~isinf(b));
M(n)=mean(b);
S(n)=std(b)/sqrt(length(b));
end
errorbar([1:1:length(M)],M,S,'r');
hold on
plot([1:1:length(M)],M,'r.');
[r,p]=corr([1:1:length(M)-1]',M(2:end)');
ylabel('ratio in out Norm');
xlabel('day');
title('ratio in out Norm');


subplot(224);
a=ratioInOutCorrectedNormGrid;
M=[];
S=[];
for n=1:size(a,2);
    b=a(:,n);
    b=b(~isnan(b));
     b=b(~isinf(b));
M(n)=mean(b);
S(n)=std(b)/sqrt(length(b));
end
errorbar([1:1:length(M)],M,S,'r');
hold on
plot([1:1:length(M)],M,'r.');
[r,p]=corr([1:1:length(M)-1]',M(2:end)');
ylabel('ratio in out corrcted Norm');
xlabel('day');
title('ratio in out corrcted Norm');

saveas(gcf,'fieldRatioInOut.fig');
%% spatial information
p=pwd;
cd ..\
load('indicesAllNewCueThresh.mat');
load('SI.mat');
cd(p);
i=indicesAllNewCueThresh.GridIdx;
SIG1=SI1(i,:);
SIG2=SI2(i,:);
SIG3=SI3(i,:);

figure,
A=SIG1;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));

subplot(131)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('spatial information');
title(['SI original p=',num2str(p)]);

A=SIG2;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));

subplot(132)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('spatial information');
title(['SI 1cm/s p=',num2str(p)]);

A=SIG3;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));

subplot(133)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('spatial information');
title(['SI new 1cm/s p=',num2str(p)]);

saveas(gcf,'spatialInformation.fig');
save('SIG.mat','SIG1','SIG2','SIG3');
%% grid cell phase relationship during learning: use correlation to approximate
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\foldersAll.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\dfofSigNormAllFOVs.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\indicesAllNewCueThresh.mat');
p=pwd;
NCommonCellsEachFOV=[];

for n=1:length(foldersAll);
    cd(foldersAll{n})
    load('commonCells.mat');
    NCommonCellsEachFOV(n)=size(commonCells,1);
end

cd(p);
NCommonCellsEachFOV=[0 NCommonCellsEachFOV];
%cell number range
cellIdxRange=[];

for n=2:length(NCommonCellsEachFOV);
 cellIdxRange(n-1,2)=sum(NCommonCellsEachFOV(1:n));
 cellIdxRange(n-1,1)=sum(NCommonCellsEachFOV(1:n-1))+1;
end

gridFOVs=[];

for n=1:length(indicesAllNewCueThresh.GridIdx);
    a=[];
    i=indicesAllNewCueThresh.GridIdx(n);
    for m=1:size(cellIdxRange,1);
        if i>=cellIdxRange(m,1) && i<=cellIdxRange(m,2);
            a(m)=1;
        else
            a(m)=0;
        end
    end
    gridFOVs(n)=find(a==1);
    
end

gridIdxFOVs={}; %group grid cells that are in the same FOV
availableFOVs=unique(gridFOVs);
for n=1:length(availableFOVs);
    gridIdxFOVs{n}=indicesAllNewCueThresh.GridIdx(gridFOVs==availableFOVs(n));
end

%now get the dfof of these cells
gridDfofSigNormFOV={};

for n=1:length(gridIdxFOVs);
    gridDfofSigNormFOV{n}={};
    for m=1:length(dfofSigNormAllFOVs);
       gridDfofSigNormFOV{n}{m}= dfofSigNormAllFOVs{m}(gridIdxFOVs{n},:);
    end
end
        

% now compute the pairwise correlation

gridDfofSigNormFOVCorr={};

for n=1:length(gridDfofSigNormFOV);
    if length(gridIdxFOVs{n})>1;
    gridDfofSigNormFOVCorr{n}=[];%each column is one day, each row is one pair    
    A=gridDfofSigNormFOV{n};
    for m=1:length(A);
   C=[];     
AA=A{m};
for i=1:size(AA,1)-1;
    ai=AA(i,:);
    for ii=i+1:size(AA,1);
        aii=AA(ii,:);
        C(end+1,1)=corr(ai',aii');
    end
end

gridDfofSigNormFOVCorr{n}(:,m)=C;
    end
    end
end

gridDfofSigNormFOVCorrAll=[];

for n=1:length(gridDfofSigNormFOVCorr);
    gridDfofSigNormFOVCorrAll(end+1:end+size(gridDfofSigNormFOVCorr{n},1),:)=gridDfofSigNormFOVCorr{n};
end

gridDfofSigNormFOVCorrAllNorm=[];
for n=1:size(gridDfofSigNormFOVCorrAll,1);
    gridDfofSigNormFOVCorrAllNorm(n,:)=gridDfofSigNormFOVCorrAll(n,:)-gridDfofSigNormFOVCorrAll(n,2);
end

figure
subplot(211)
M=mean(gridDfofSigNormFOVCorrAll,1);
S=nansem(gridDfofSigNormFOVCorrAll,1);

errorbar([1:1:11],M,S)
title('pairwise mean activity corr, no normalization');

subplot(212)
M=mean(gridDfofSigNormFOVCorrAllNorm,1);
S=nansem(gridDfofSigNormFOVCorrAllNorm,1);

errorbar([1:1:11],M,S)
title('pairwise mean activity corr, normalization');

saveas(gcf,'pairwiseMeanActivityCorrelation.fig')

%% correlation using run by run and also dfof_sig
cd ..\
load('foldersAll.mat');
cd('gridCellsNewCueThresh');
CCAllFOVs=[];
CCAllFOVsNorm=[];
CCAllDfofSigFOVs=[];
CCAllDfofSigFOVsNorm=[];
p=pwd;
for n=1:length(foldersAll);
    cd(foldersAll{n});
    cd('reIdentifyCueCellsUniThreshAddingNewMice\gridActivityPairwiseCorr');
    clear CCAll
    clear CCAllNorm
    clear CCAllDfofSig
    clear CCAllDfofSigNorm
    if exist('CCAll.mat')>0;
        load('CCAll.mat');
          load('CCAllNorm.mat');
          load('CCAllDfofSig.mat');
          load('CCAllDfofSigNorm.mat');
   
    
    CCAllFOVs(end+1:end+size(CCAll,1),:)=CCAll;
    CCAllFOVsNorm(end+1:end+size(CCAllNorm,1),:)=CCAllNorm;
    CCAllDfofSigFOVs(end+1:end+size(CCAllDfofSig,1),:)=CCAllDfofSig;
    CCAllDfofSigFOVsNorm(end+1:end+size(CCAllDfofSigNorm,1),:)=CCAllDfofSigNorm;
    end
end

cd(p);

figure
subplot(221)
M=mean(CCAllFOVs,1);
S=nansem(CCAllFOVs,1);

errorbar([1:1:11],M,S)
title('pairwise RBR activity corr, no normalization');

subplot(223)
M=mean(CCAllFOVsNorm,1);
S=nansem(CCAllFOVsNorm,1);

errorbar([1:1:11],M,S)
title('pairwise RBR corr, normalization');

  
subplot(222)
M=mean(CCAllDfofSigFOVs,1);
S=nansem(CCAllDfofSigFOVs,1);

errorbar([1:1:11],M,S)
title('pairwise DFOF SIG corr, no normalization');

subplot(224)
M=mean(CCAllDfofSigFOVsNorm,1);
S=nansem(CCAllDfofSigFOVsNorm,1);

errorbar([1:1:11],M,S)
title('pairwise DFOF SIG corr, normalization');

  saveas(gcf,'pairwiseRBRAndDfofSigCorrelation.fig')
  
%% fluorescence

%3 calculations
pp=pwd;
cd ..\
load('foldersAll.mat');
cd(pp);
mF=[];% mean fluorescence. %each row is a cell and each column is a day
mDfof=[];% mean dfof, initial calculation,using 1cm/s speed thresh. %each row is a cell and each column is a day
mDfofSig=[];% mean dfof sig, new calculation,using 1cm/s speed thresh. %each row is a cell and each column is a day
mDfofMean=[]; %mean of mean dfof
mDfofSigMean=[]; %mean of mean dfof sig

for n=1:length(foldersAll);
  load([foldersAll{n} '\' 'commonCells.mat']);
  load([foldersAll{n} '\reIdentifyCueCellsUniThreshAddingNewMice\indices.mat']);
  load([foldersAll{n} '\' 'useFolders.mat']);
  disp(n)
  A1=[];
  A2=[];
  A3=[];
  A4=[];
  A5=[];
  commonGrid=commonCells(indices.iGrid,:);
  for m=1:size(commonGrid,2);
       disp(m)
      a=commonGrid(:,m);
      cd(useFolders{m});
      
      mm='dfof_*';
d=dir(mm);
for i=1:length(d);
    load(d(i).name);
end

    mm='F_*';
d=dir(mm);
for i=1:length(d);
    load(d(i).name);
end

    mm='dfofaveragesmooth_*';
d=dir(mm);
for i=1:length(d);
    load(d(i).name);
end
     F=F';
     dfof=dfof';
     dfof_sig=dfof_sig';
 
     dfofaveragesmooth=dfofaveragesmooth';
     dfofaveragesmooth_sig=dfofaveragesmooth_sig';
     
      A1(:,m)=mean(F(a,:),2);
       A2(:,m)=mean(dfof(a,:),2);
       A3(:,m)=mean(dfof_sig(a,:),2);
       A4(:,m)=mean(dfofaveragesmooth(a,:),2);
       A5(:,m)=mean( dfofaveragesmooth_sig(a,:),2);
  end
  mF(end+1:end+size(A1,1),:)=A1;
  mDfof(end+1:end+size(A2,1),:)=A2;
  mDfofSig(end+1:end+size(A3,1),:)=A3;
  mDfofMean(end+1:end+size(A2,1),:)=A4;
  mDfofSigMean(end+1:end+size(A3,1),:)=A5;
end


figure,
A=mF;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));

subplot(151)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean F');
title(['mean F p=',num2str(p)]);

A=mDfof;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));

subplot(152)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean dfof');
title(['mean dfof p=',num2str(p)]);

A=mDfofSig;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));

subplot(153)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean dfof sig');
title(['mean dfof sig p=',num2str(p)]);

A=mDfofMean;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));

subplot(154)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean dfof mean');
title(['mean dfof mean p=',num2str(p)]);

A=mDfofSigMean;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));

subplot(155)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean dfof sig mean');
title(['mean dfof sig mean p=',num2str(p)]);
tightfig;

cd(pp);
%
saveas(gcf,'allF.fig');
save('allF.mat','mF','mDfof','mDfofSig','mDfofMean','mDfofSigMean');

     %% sig trans info for grid cells
 cd ..\
 load('allSigTransSP.mat');
 load('indicesAllNewCueThresh.mat');
i=indicesAllNewCueThresh.GridIdx;

 cd('gridCellsNewCueThresh');
 
 mFreqSP=mFreqSP(i,:);
mFreqRunSP=mFreqRunSP(i,:);
mDisCoverSP=mDisCoverSP(i,:);
mDisCoverRunSP=mDisCoverRunSP(i,:);
mAmpSP=mAmpSP(i,:);
mSigIntegSP=mSigIntegSP(i,:);
mDurSP=mDurSP(i,:);
mSigIntegPerRunSP=mSigIntegPerRunSP(i,:);

 figure,
A=mFreqSP;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));

subplot(241)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean Freq');
title(['mean Freq p=',num2str(p)]);

A=mFreqRunSP;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));

subplot(242)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean Freq run');
title(['mean Freq run p=',num2str(p)]);

A=mDisCoverSP;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));

subplot(243)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean Dis cover');
title(['mean Dis cover p=',num2str(p)]);

A=mDisCoverRunSP;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));

subplot(244)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean Dis cover Run');
title(['mean Dis cover Run p=',num2str(p)]);

A=mAmpSP;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));

subplot(245)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean Amp');
title(['mean Amp p=',num2str(p)]);

A=mSigIntegSP;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));

subplot(246)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean SigInteg');
title(['mean SigInteg p=',num2str(p)]);

A=mDurSP;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));

subplot(247)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean Dur');
title(['mean Dur p=',num2str(p)]);

A=mSigIntegPerRunSP;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));

subplot(248)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean sig integ per run');
title(['mean sig integ per run p=',num2str(p)]);


saveas(gcf,'allSigTransSP.fig');
save('allSigTransSP.mat','mFreqSP','mFreqRunSP','mDisCoverSP','mDisCoverRunSP','mAmpSP','mSigIntegSP','mDurSP','mSigIntegPerRunSP');

    
 %%    
load('fieldDistriGrid.mat');
NField={};
MField=[];
SField=[];
for n=1:length(fieldDistriGrid);
    NField{n}=[];
    for m=1:length(fieldDistriGrid{n});
        
NField{n}(m)=length(fieldDistriGrid{n}{m});
    end
    MField(n,1)=nanmean(NField{n});
    SField(n,1)=nansem(NField{n},2);
end

figure,
errorbar([1:1:11],MField,SField)

% load('fieldDistriUC.mat');
% NFieldUC={};
% MFieldUC=[];
% SFieldUC=[];
% for n=1:length(fieldDistriUC);
%     NFieldUC{n}=[];
%     for m=1:length(fieldDistriUC{n});
%         
% NFieldUC{n}(m)=length(fieldDistriUC{n}{m});
%     end
%     MFieldUC(n,1)=nanmean(NFieldUC{n});
%     SFieldUC(n,1)=nansem(NFieldUC{n},2);
% end
% 
% hold on
% errorbar([1:1:11],MFieldUC,SFieldUC)
title('binsWithFieldsGrid');
% legend('common','uncommon','Location','Northeast');
saveas(gcf,'binsWithFieldsGrid.fig');
save('NField.mat','NField');
save('MField.mat','MField');
save('SField.mat','SField');
  
%% gaining fields
load('fieldDistriGrid.mat');
% load('fieldDistriUC.mat');

fieldDistri_zeroOne={};
% fieldDistriUC_zeroOne={};

for n=1:length(fieldDistriGrid);
    fieldDistri_zeroOne{n}=zeros(length(fieldDistriGrid{n}),200);
%     fieldDistriUC_zeroOne{n}=zeros(length(fieldDistriUC{n}),200);
    for m=1:length(fieldDistriGrid{n});
        i=fieldDistriGrid{n}{m};
        if ~isnan(i);
        fieldDistri_zeroOne{n}(m,i)=1;
        end
    end
%     for m=1:length(fieldDistriUC{n});
%         i=fieldDistriUC{n}{m};
%          if ~isnan(i);
%         fieldDistriUC_zeroOne{n}(m,i)=1;
%          end
%     end
end
save('fieldDistri_zeroOne.mat','fieldDistri_zeroOne');
% save('fieldDistriUC_zeroOne.mat','fieldDistriUC_zeroOne');

fieldDistriDiff=[];
b=[2 3];%day 1 and 2 in new env (2 and 3 in data)
a=[8 9 10 11];%AFter learning days
fieldDistriBeforeLearning=[];
fieldDistriAfterLearning=[];
for n=1:length(fieldDistriGrid{1});
    beforeLearning=[];
    for i=1:length(b);
        beforeLearning(end+1,:)=fieldDistri_zeroOne{b(i)}(n,:);
    end
    afterLearning=[];
    for i=1:length(a);
        afterLearning(end+1,:)=fieldDistri_zeroOne{a(i)}(n,:);
    end
    fieldDistriDiff(n,:)=mean(afterLearning,1)-mean(beforeLearning);
    fieldDistriBeforeLearning(n,:)=mean(beforeLearning,1);
     fieldDistriAfterLearning(n,:)=mean(afterLearning,1);
end

figure

subplot(211)
errorbar([1:1:200],mean(fieldDistriBeforeLearning,1),nansem(fieldDistriBeforeLearning,1));
hold on
errorbar([1:1:200],mean(fieldDistriAfterLearning,1),nansem(fieldDistriAfterLearning,1));
title('field distri before after individual cells');

subplot(212)
errorbar([1:1:200],mean(fieldDistriDiff,1),nansem(fieldDistriDiff,1));
title('field distri diff indiv cells')

saveas(gcf,'fieldDiffIndividualCells.fig');
save('fieldDistriDiff.mat','fieldDistriDiff')
save('fieldDistriBeforeLearning.mat','fieldDistriBeforeLearning')
save('fieldDistriAfterLearning.mat','fieldDistriAfterLearning')
    
   
%% calculate fractional changes of fieldDistriDiff 
fieldDistriDiffFraction=[];
load('fieldDistriDiff.mat');
load('fieldDistriBeforeLearning.mat');
%calculate fractional changes

%each case is different, do it one by one

for n=1:size(fieldDistriDiff,1);
   fieldDistriDiffFraction(n,:)=fieldDistriDiff(n,:)/mean(fieldDistriBeforeLearning(n,:));
end

save('fieldDistriDiffFraction.mat','fieldDistriDiffFraction');
    

