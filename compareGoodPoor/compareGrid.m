%compareing good and bad in different aspects

%% compare run by run consistency: use sig: run by run to mean

%compare grid and all in genera

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\gridCellsNewCueThresh\corrAllGrid.mat');
G=corrAllGrid;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\gridCellsNewCueThresh\corrAllGrid.mat');
B=corrAllGrid;

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\corrAllFOVs.mat');
GG=corrAllFOVs;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\corrAllFOVs.mat');
BB=corrAllFOVs;


%only use dfofSig Mean To Mean and do shuffle

figure
subplot(131)
A=G.dfofSig{1};
MTempRL=G.dfofSigMeanToMean;
S=std(A,1)/sqrt(size(A,1));
[~,p1]=corr([1:1:10]',MTempRL(2:end)');
errorbar([1:1:length(MTempRL)],MTempRL,S,'g','LineWidth',2);

A=B.dfofSig{1};
MTempRL=B.dfofSigMeanToMean;
S=std(A,1)/sqrt(size(A,1));
[~,p2]=corr([1:1:10]',MTempRL(2:end)');
hold on
errorbar([1:1:length(MTempRL)],MTempRL,S,'m','LineWidth',2);

A=GG.dfofSig{1};
MTempRL=GG.dfofSigMeanToMean;
S=std(A,1)/sqrt(size(A,1));
[~,p1]=corr([1:1:10]',MTempRL(2:end)');
hold on
errorbar([1:1:length(MTempRL)],MTempRL,S,'k');

A=BB.dfofSig{1};
MTempRL=BB.dfofSigMeanToMean;
S=std(A,1)/sqrt(size(A,1));
[~,p2]=corr([1:1:10]',MTempRL(2:end)');
hold on
errorbar([1:1:length(MTempRL)],MTempRL,S,'Color',[0.6 0.6 0.6]);


sig=[];
for n=1:length(MTempRL);
    g=G.dfofSig{1}(:,n);
    b=B.dfofSig{1}(:,n);
    [~,p]=ttest2(g,b);
    if p<0.05;
        sig(n)=1;
    else
        sig(n)=0;
    end
end

xlim([0.5 11.5]);
title(['corrAbsValues grid and all good p=',num2str(p1),'bad p=',num2str(p2)])

%% compare grid nongrid in good and bad learners

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\gridCellsNewCueThresh\corrAllGrid.mat');
GG=corrAllGrid;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\nonGridNewCueThresh\corrAllNGrid.mat');
GN=corrAllNGrid;

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\gridCellsNewCueThresh\corrAllGrid.mat');
BG=corrAllGrid;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\nonGridNewCueThresh\corrAllNGrid.mat');
BN=corrAllNGrid;

subplot(132)
A=GG.dfofSig{1};
MTempRL=GG.dfofSigMeanToMean;
S=std(A,1)/sqrt(size(A,1));

hold on
errorbar([1:1:length(MTempRL)],MTempRL,S,'g-');

A=GN.dfofSig{1};
MTempRL=GN.dfofSigMeanToMean;
S=std(A,1)/sqrt(size(A,1));
hold on
errorbar([1:1:length(MTempRL)],MTempRL,S,'g--');

A=BG.dfofSig{1};
MTempRL=BG.dfofSigMeanToMean;
S=std(A,1)/sqrt(size(A,1));

hold on
errorbar([1:1:length(MTempRL)],MTempRL,S,'m-');

A=BN.dfofSig{1};
MTempRL=BN.dfofSigMeanToMean;
S=std(A,1)/sqrt(size(A,1));
hold on
errorbar([1:1:length(MTempRL)],MTempRL,S,'m--');

title('commonNonCommonCorr');
xlabel('day');
ylabel('corr');
xlim([0.5 11.5]);

%% compare NORMALIZED run by run consistency: use sig: run by run to mean
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\gridCellsNewCueThresh\corrAllGrid.mat');
G=corrAllGrid;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\gridCellsNewCueThresh\corrAllGrid.mat');
B=corrAllGrid;

%only use dfofSig Mean To Mean and do shuffle

subplot(133)
A=G.dfofSigNorm{1};
MTempRL=G.dfofSigNormMeanToMean;
S=std(A,1)/sqrt(size(A,1));
[~,p1]=corr([1:1:10]',MTempRL(2:end)');
errorbar([1:1:length(MTempRL)],MTempRL,S,'g');
A=B.dfofSigNorm{1};
MTempRL=B.dfofSigNormMeanToMean;
S=std(A,1)/sqrt(size(A,1));
[~,p2]=corr([1:1:10]',MTempRL(2:end)');
hold on
errorbar([1:1:length(MTempRL)],MTempRL,S,'m');

sig=[];
for n=1:length(MTempRL);
    g=G.dfofSigNorm{1}(:,n);
    b=B.dfofSigNorm{1}(:,n);
      [~,p]=ttest2(g,b);
    if p<0.05;
        sig(n)=1;
    else
        sig(n)=0;
    end

end

title(['corrNormValue good p=',num2str(p1),'bad p=',num2str(p2)])
xlim([0.5 11.5]);
tightfig;

saveas(gcf,'corrAbsNormValue.fig');


% %%
% %look at the correlation of this with behavior, use 10 run data
% 
% %activity
% act=[G.dfofSigMeanToMean B.dfofSigMeanToMean];
% load('E:\learningAnalysis\behaviorNewBatch\predLickDataBad.mat');
% load('E:\learningAnalysis\behaviorNewBatch\predLickDataGood.mat');
% lick=[predLick10RunsMeanGood(1:11) predLick10RunsMeanBad(1:11)];
% load('E:\learningAnalysis\behaviorNewBatch\slowDownDataBad.mat');
% load('E:\learningAnalysis\behaviorNewBatch\slowDownDataGood.mat');
% slow=[slowDown10RunsMeanGood(1:11) slowDown10RunsMeanBad(1:11)];
% 
% figure,
% subplot(231)
% plot(act,lick,'r.','MarkerSize',10)
% [~,p]=corr(act',lick');
% title(['all lick p=',num2str(p)])
% 
% subplot(234)
% plot(act,slow,'r.','MarkerSize',10)
% [~,p]=corr(act',slow');
% title(['all slow p=',num2str(p)])
% 
% subplot(232)
% plot(act(1:11),lick(1:11),'r.','MarkerSize',10)
% [~,p]=corr(act(1:11)',lick(1:11)');
% title(['Good lick p=',num2str(p)])
% 
% subplot(235)
% plot(act(1:11),slow(1:11),'r.','MarkerSize',10)
% [~,p]=corr(act(1:11)',slow(1:11)');
% title(['Good slow p=',num2str(p)])
% 
% subplot(233)
% plot(act(12:end),lick(12:end),'r.','MarkerSize',10)
% [~,p]=corr(act(12:end)',lick(12:end)');
% title(['Good lick p=',num2str(p)])
% 
% subplot(236)
% plot(act(12:end),slow(12:end),'r.','MarkerSize',10)
% [~,p]=corr(act(12:end)',slow(12:end)');
% title(['Good slow p=',num2str(p)])
% 
% saveas(gcf,'absCorrAndBehavior.fig');

% %% look at matrix corrlations
% 
% load('E:\learningAnalysis\summaryManyMice_includingOldEnvGood\gridCells\dfofSigNormGridSort_Max1.mat');
% G=dfofSigNormGridSort_Max1;
% load('E:\learningAnalysis\summaryManyMice_includingOldEnvBad\gridCells\dfofSigNormGridSort_Max1.mat');
% B=dfofSigNormGridSort_Max1;
% figure
% for n=1:length(G);
%     subplot(2,length(G),n);
%     imagesc(G{n})
%     
%     title(['Day',num2str(n)]);
%     axis off
% end
% for n=1:length(B);
%     subplot(2,length(G),length(G)+n);
%     imagesc(B{n})
%     
% %     title(['Day',num2str(n)]);
%     axis off
% end
% 
% tightfig;
% saveas(gcf,'sortedByMaxOld_fig');
% %% compare matrix correlation
% 
% %grid, non grid, in good and bad learnears
% 
% load('E:\learningAnalysis\summaryManyMice_includingOldEnvGood\gridCells\dfofSigNormGridSort_Max2.mat');
% G=dfofSigNormGridSort_Max2;
% 
% load('E:\learningAnalysis\summaryManyMice_includingOldEnvBad\gridCells\dfofSigNormGridSort_Max2.mat');
% B=dfofSigNormGridSort_Max2;
% figure
% for n=1:length(G);
%     subplot(2,length(G),n);
%     imagesc(G{n})
%     
%     title(['Day',num2str(n)]);
%     axis off
% end
% for n=1:length(B);
%     subplot(2,length(G),length(G)+n);
%     imagesc(B{n})
%     
% %     title(['Day',num2str(n)]);
%     axis off
% end
% 
% tightfig;
% saveas(gcf,'sortedByMaxNew1_fig');
% 
% %%
% 
% load('E:\learningAnalysis\summaryManyMice_includingOldEnvGood\gridCells\dfofSigNormGridSort_Max11.mat');
% G=dfofSigNormGridSort_Max11;
% load('E:\learningAnalysis\summaryManyMice_includingOldEnvBad\gridCells\dfofSigNormGridSort_Max11.mat');
% B=dfofSigNormGridSort_Max11;
% figure
% for n=1:length(G);
%     subplot(2,length(G),n);
%     imagesc(G{n})
%     
%     title(['Day',num2str(n)]);
%     axis off
% end
% for n=1:length(B);
%     subplot(2,length(G),length(G)+n);
%     imagesc(B{n})
%     
% %     title(['Day',num2str(n)]);
%     axis off
% end
% 
% tightfig;
% saveas(gcf,'sortedByMaxNew11_fig');

%% matrix correlation: left whole matrix, right: individual cells
%compare grid non grid, in good and bad learners
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\gridCellsNewCueThresh\matrixCorrMax1.mat');
G=matrixCorr;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\nonGridNewCueThresh\matrixCorrMax1.mat');
GN=matrixCorr;

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\gridCellsNewCueThresh\matrixCorrMax1.mat');
B=matrixCorr;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\nonGridNewCueThresh\matrixCorrMax1.mat');
BN=matrixCorr;

figure
subplot(121);
plot(G,'g');
[r,p1]=corr(G(2:end),[1:1:length(G)-1]');
hold on
plot(GN,'g--');
[r,p11]=corr(GN(2:end),[1:1:length(GN)-1]');
hold on
plot(B,'m');
[r,p2]=corr(B(2:end),[1:1:length(B)-1]');
plot(BN,'m--');
[r,p12]=corr(BN(2:end),[1:1:length(BN)-1]');
xlim([0 11])
ylim([-0.05 0.7])
title(['matrix corr good p=',num2str(p1),'bad p=',num2str(p2)]);

breakyaxis([0.1 0.35]);


%%

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\gridCellsNewCueThresh\corrIndivi.mat');
G=corrIndivi;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\nonGridNewCueThresh\corrIndivi.mat');
GN=corrIndivi;

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\gridCellsNewCueThresh\corrIndivi.mat');
B=corrIndivi;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\nonGridNewCueThresh\corrIndivi.mat');
BN=corrIndivi;

subplot(122)
M=mean(G,1);
S=std(G,1)/sqrt(size(G,1));
errorbar([1:1:length(M)],M,S,'g-')
[r,p1]=corr([1:1:length(M)-1]', M(2:end)');

M=mean(GN,1);
S=std(GN,1)/sqrt(size(GN,1));
hold on
errorbar([1:1:length(M)],M,S,'g--')
[r,p1]=corr([1:1:length(M)-1]', M(2:end)');

M=mean(B,1);
S=std(B,1)/sqrt(size(B,1));
hold on
errorbar([1:1:length(M)],M,S,'m-')
[r,p2]=corr([1:1:length(M)-1]', M(2:end)');
title(['indiv corr good p=',num2str(p1),'bad p=',num2str(p2)]);

M=mean(BN,1);
S=std(BN,1)/sqrt(size(BN,1));
hold on
errorbar([1:1:length(M)],M,S,'m--')
[r,p2]=corr([1:1:length(M)-1]', M(2:end)');
title(['indiv corr good p=',num2str(p1),'bad p=',num2str(p2)]);
xlim([0 11])
ylim([-0.1 0.75])
breakyaxis([0.05 0.35]);

sig=[];
for n=1:size(G,2);
    g=G(:,n);
    b=B(:,n);
    [sig(n,1),~]=ttest2(g,b);
end

saveas(gcf,'matrixAndIndiviCorr.fig');
%% field distribution
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\gridCellsNewCueThresh\fieldDistriGridEachBin.mat');
G=fieldDistriGridEachBin;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\gridCellsNewCueThresh\fieldDistriGridEachBin.mat');
B=fieldDistriGridEachBin;
load('E:\ID20201118\20210103\1stloc1\pcaica\cueAnalysisNew_sig\tempRL.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\gridCellsNewCueThresh\dfofSigNormGrid.mat');
NG=size(dfofSigNormGrid{1},1);
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\gridCellsNewCueThresh\dfofSigNormGrid.mat');
NB=size(dfofSigNormGrid{1},1);


figure,
subplot(211);
plot([1:1:200],(mean(G(4:11,:),1)-mean(G(2:3,:),1))/NG,'g')
hold on
plot([1:1:200],(mean(B(4:11,:),1)-mean(B(2:3,:),1))/NB,'m')

plot([1:1:200],tempRL*0.1);
hold on
plot([890/5 890/5],[0 0.1],'r');
xlim([0 200]);
title('fieldDistri');

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\gridCellsNewCueThresh\fieldDistriAmpSumGrid.mat');
G=fieldDistriAmpSumGrid;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\gridCellsNewCueThresh\fieldDistriAmpSumGrid.mat');
B=fieldDistriAmpSumGrid;

subplot(212);
plot([1:1:200],(mean(G(4:11,:),1)-mean(G(2:3,:),1))/NG,'g')
hold on
plot([1:1:200],(mean(B(4:11,:),1)-mean(B(2:3,:),1))/NB,'m')

plot([1:1:200],tempRL*0.01);
hold on
plot([890/5 890/5],[0 0.01],'r');
xlim([0 200]);
title('fieldAmp');
saveas(gcf,'fieldDistriAmp_beforeAfterLearn.fig');

%% test field distri significance

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\gridCellsNewCueThresh\fieldDistriGridEachBin.mat');
G=fieldDistriGridEachBin;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\gridCellsNewCueThresh\sigIdxFieldDistri.mat');
GSig=sigIdxFieldDistri;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\gridCellsNewCueThresh\fieldDistriGridEachBin.mat');
B=fieldDistriGridEachBin;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\gridCellsNewCueThresh\sigIdxFieldDistri.mat');
BSig=sigIdxFieldDistri;
load('E:\ID20201118\20210103\1stloc1\pcaica\cueAnalysisNew_sig\tempRL.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\gridCellsNewCueThresh\dfofSigNormGrid.mat');
NG=size(dfofSigNormGrid{1},1);
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\gridCellsNewCueThresh\dfofSigNormGrid.mat');
NB=size(dfofSigNormGrid{1},1);

figure,
a=(mean(G(4:11,:),1)-mean(G(2:3,:),1))/NG;
plot([1:1:200],a,'g')
hold on
plot(GSig,a(GSig),'r.','MarkerSize',10);

a=(mean(B(4:11,:),1)-mean(B(2:3,:),1))/NG;
plot([1:1:200],a,'m')
hold on
plot(BSig,a(BSig),'r.','MarkerSize',10);

plot([1:1:200],tempRL*0.1);
hold on
plot([890/5 890/5],[0 0.1],'r');
xlim([0 200]);
title('fieldDistri');
saveas(gcf,'fieldDistri__beforeAfterLearn_sig.fig');

%% plotting all activity together
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\gridCellsNewCueThresh\dfofSigNormGridSort_Max1.mat');
G1=dfofSigNormGridSort_Max1;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\gridCellsNewCueThresh\dfofSigNormGridSort_Max2.mat');
G2=dfofSigNormGridSort_Max2;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\gridCellsNewCueThresh\dfofSigNormGridSort_Max11.mat');
G11=dfofSigNormGridSort_Max11;

figure
for n=1:length(G1);
    subplot(3,length(G1),n);
    imagesc(G1{n})
    
    title(['Day',num2str(n)]);
    axis off
end

for n=1:length(G2);
    subplot(3,length(G2),length(G2)+n);
    imagesc(G2{n})
    
%     title(['Day',num2str(n)]);
    axis off
end

for n=1:length(G11);
    subplot(3,length(G11),length(G11)*2+n);
    imagesc(G11{n})
    
%     title(['Day',num2str(n)]);
    axis off
end
tightfig;

saveas(gcf,'sortedActivity_Good.fig');

%% plotting all activity together
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\gridCellsNewCueThresh\dfofSigNormGridSort_Max1.mat');
B1=dfofSigNormGridSort_Max1;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\gridCellsNewCueThresh\dfofSigNormGridSort_Max2.mat');
B2=dfofSigNormGridSort_Max2;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\gridCellsNewCueThresh\dfofSigNormGridSort_Max11.mat');
B11=dfofSigNormGridSort_Max11;

figure
for n=1:length(B1);
    subplot(3,length(B1),n);
    imagesc(B1{n})
    
    title(['Day',num2str(n)]);
    axis off
end

for n=1:length(B2);
    subplot(3,length(B2),length(B2)+n);
    imagesc(B2{n})
    
%     title(['Day',num2str(n)]);
    axis off
end

for n=1:length(B11);
    subplot(3,length(B11),length(B11)*2+n);
    imagesc(B11{n})
    
%     title(['Day',num2str(n)]);
    axis off
end
tightfig;

saveas(gcf,'sortedActivity_Bad.fig');
%% plot dfof together
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\gridCellsNewCueThresh\allF.mat');

GF=mF;
GDfof=mDfof;
GDfofSig=mDfofSig;
GDfofMean=mDfofMean;
GDfofSigMean=mDfofSigMean;

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\gridCellsNewCueThresh\allF.mat');

BF=mF;
BDfof=mDfof;
BDfofSig=mDfofSig;
BDfofMean=mDfofMean;
BDfofSigMean=mDfofSigMean;

figure,

subplot(151)
A=GF;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));
errorbar([1:1:11],M,S,'g')

A=BF;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));
hold on
errorbar([1:1:11],M,S,'m')

xlabel('days')
ylabel('mean F');
title(['mean F']);


subplot(152)
A=GDfof;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));
errorbar([1:1:11],M,S,'g')

A=BDfof;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));
hold on
errorbar([1:1:11],M,S,'m')

xlabel('days')
ylabel('mean dfof');
title(['mean dfof']);

subplot(153)
A=GDfofSig;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));
errorbar([1:1:11],M,S,'g')

A=BDfofSig;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));
hold on
errorbar([1:1:11],M,S,'m')
xlabel('days')
ylabel('mean dfof sig');
title(['mean dfof sig']);

subplot(154)
A=GDfofMean;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));
errorbar([1:1:11],M,S,'g')

A=BDfofMean;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));
hold on
errorbar([1:1:11],M,S,'m')
xlabel('days')
ylabel('mean dfof mean');
title(['mean dfof mean']);

subplot(155)
A=GDfofSigMean;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));
errorbar([1:1:11],M,S,'g')

A=BDfofSigMean;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));
hold on
errorbar([1:1:11],M,S,'m')

xlabel('days')
ylabel('mean dfof sig mean');
title(['mean dfof sig mean']);
tightfig;
saveas(gcf,'allF.fig');

%%
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\gridCellsNewCueThresh\allSigTransSP.mat');

GmFreq=(mFreqSP);
GmFreqRun=(mFreqRunSP);
GmDisCover=(mDisCoverSP);
GmDisCoverRun=(mDisCoverRunSP);
GmAmp=(mAmpSP);
GmSigInteg=(mSigIntegSP);
GmDur=(mDurSP);
GmSigIntegPerRun=(mSigIntegPerRunSP);

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\gridCellsNewCueThresh\allSigTransSP.mat');

BmFreq=(mFreqSP);
BmFreqRun=(mFreqRunSP);
BmDisCover=(mDisCoverSP);
BmDisCoverRun=(mDisCoverRunSP);
BmAmp=(mAmpSP);
BmSigInteg=(mSigIntegSP);
BmDur=(mDurSP);
BmSigIntegPerRun=(mSigIntegPerRunSP);

sig=[];

figure,
A=GmFreq;
B=BmFreq;
M1=nanmean(A,1);
S1=nansem(A,1);

M2=nanmean(B,1);
S2=nansem(B,1);
N=1;
for n=1:size(A,2);
    [sig(N,n),~]=ttest2(A(:,n),B(:,n));
end

subplot(241)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')
xlabel('days')
ylabel('mean Freq n/cm');
title(['mean Freq n/cm']);

A=GmFreqRun;
B=BmFreqRun;
M1=nanmean(A,1);
S1=nansem(A,1);

M2=nanmean(B,1);
S2=nansem(B,1);
N=2;
for n=1:size(A,2);
    [sig(N,n),~]=ttest2(A(:,n),B(:,n));
end

subplot(242)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')
xlabel('days')
ylabel('mean Freq Run n/cm');
title(['mean Freq Run n/cm']);

A=GmDisCover;
B=BmDisCover;
M1=nanmean(A,1);
S1=nansem(A,1);

M2=nanmean(B,1);
S2=nansem(B,1);
N=3;
for n=1:size(A,2);
    [sig(N,n),~]=ttest2(A(:,n),B(:,n));
end

subplot(243)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')
xlabel('days')
ylabel('mean DisCover %');
title(['mean DisCover %W']);


A=GmDisCoverRun;
B=BmDisCoverRun;
M1=nanmean(A,1);
S1=nansem(A,1);

M2=nanmean(B,1);
S2=nansem(B,1);
N=4;
for n=1:size(A,2);
    [sig(N,n),~]=ttest2(A(:,n),B(:,n));
end

subplot(244)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')
xlabel('days')
ylabel('mean DisCover run %');
title(['mean DisCover run %']);


A=GmAmp;
B=BmAmp;
M1=nanmean(A,1);
S1=nansem(A,1);

M2=nanmean(B,1);
S2=nansem(B,1);
N=5;
for n=1:size(A,2);
    [sig(N,n),~]=ttest2(A(:,n),B(:,n));
end
subplot(245)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')
xlabel('days');
ylabel('mean Amp');
title(['mean Amp']);

A=GmSigInteg;
B=BmSigInteg;
M1=nanmean(A,1);
S1=nansem(A,1);

M2=nanmean(B,1);
S2=nansem(B,1);
N=6;
for n=1:size(A,2);
    [sig(N,n),~]=ttest2(A(:,n),B(:,n));
end
subplot(246)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')
xlabel('days');
ylabel('mean SigInteg');
title(['mean SigInteg']);

A=GmDur;
B=BmDur;
M1=nanmean(A,1);
S1=nansem(A,1);

M2=nanmean(B,1);
S2=nansem(B,1);
N=7;
for n=1:size(A,2);
    [sig(N,n),~]=ttest2(A(:,n),B(:,n));
end
subplot(247)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')

xlabel('days')
ylabel('mean Duration (sec)');
title(['mean Duration (sec)']);


A=GmSigIntegPerRun;
B=BmSigIntegPerRun;
M1=nanmean(A,1);
S1=nansem(A,1);

M2=nanmean(B,1);
S2=nansem(B,1);
N=8;
for n=1:size(A,2);
    [sig(N,n),~]=ttest2(A(:,n),B(:,n));
end
subplot(248)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')

xlabel('days')
ylabel('mean sigIntegPerRun');
title(['mean sigIntegPerRun']);

tightfig;
saveas(gcf,'allSigTransSP.fig');

%% grid scales
clear all
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\gridCellsNewCueThresh\fieldInfoUsingAllCellsCorrectedGrid.mat');
G=minSpacingGrid;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\gridCellsNewCueThresh\fieldInfoUsingAllCellsCorrectedGrid.mat');
B=minSpacingGrid;

figure

a=G;
% a(a>200)=nan;
M=[];
S=[];
for n=1:size(a,2);
    b=a(:,n);
    b=b(~isnan(b));
      b=b(~isinf(b));
%       b=b(b<500);
M(n)=mean(b);
S(n)=std(b)/sqrt(length(b));
end
errorbar([1:1:length(M)],M,S,'g');

[r1,p1]=corr([1:1:10]',M(2:end)','tail','left');


a=B;
% a(a>200)=nan;
M=[];
S=[];
for n=1:size(a,2);
    b=a(:,n);
    b=b(~isnan(b));
      b=b(~isinf(b));
%       b=b(b<500);
M(n)=mean(b);
S(n)=std(b)/sqrt(length(b));
end
hold on
errorbar([1:1:length(M)],M,S,'m');

[r2,p2]=corr([1:1:10]',M(2:end)');

title(['gridScale, p1=', num2str(p1),' p2=',num2str(p2)]);

saveas(gcf,'grid scale.fig');

%% 

%% replot the above: include good and bad learners before after learning and now the differences were calculated according to the onw group change
%LOAD CUE TEMPLATE
load('E:\ID20201118\20201226\2ndloc1\pcaica\cueAnalysisNew_sig\tempRL.mat');
tempNRL=tempRL;


load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\gridCellsNewCueThresh\fieldDistriAfterLearning.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\gridCellsNewCueThresh\fieldDistriBeforeLearning.mat');

GA=fieldDistriAfterLearning;
GB=fieldDistriBeforeLearning;

figure,

subplot(211);

%plot cue template
tempNRL(tempNRL==1)= 0.3606;
tempNRL(tempNRL==0)= 0.1848;

plot(tempNRL,'color',[0.2 0.2 0.2])
hold on
plot([890/5 890/5],[min(tempNRL) max(tempNRL)],'b');
hold on
[lineOut, fillOut] = semshade(GB,0.2,'k',[1:1:200],1);
hold on
  [lineOut, fillOut] = semshade(GA,0.2,'g',[1:1:200],1);
      
sigG=[];
mG=[];

for n=1:200;
    [sigG(n),~]=ttest2(GA(:,n),GB(:,n));
    mG(n)=mean(GA(:,n))-mean(GB(:,n));
end

GAM=mean(GA,1);
GBM=mean(GB,1);

for n=1:length(sigG);
    if sigG(n)==1;
        if mG(n)>0;
            hold on
            plot(n,GAM(n),'r.','MarkerSize',10);
        else
            hold on
            plot(n,GBM(n),'r.','MarkerSize',10);
        end
    end
end
axis tight

title('good leaners before after learning');

load('E:\ID20201118\20201226\2ndloc1\pcaica\cueAnalysisNew_sig\tempRL.mat');
tempNRL=tempRL;


load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\gridCellsNewCueThresh\fieldDistriAfterLearning.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\gridCellsNewCueThresh\fieldDistriBeforeLearning.mat');

BA=fieldDistriAfterLearning;
BB=fieldDistriBeforeLearning;

subplot(212);
%plot cue template
tempNRL(tempNRL==1)= 0.3747;
tempNRL(tempNRL==0)= 0.1369;

plot(tempNRL,'color',[0.2 0.2 0.2])
hold on
plot([890/5 890/5],[min(tempNRL) max(tempNRL)],'b');

hold on
[lineOut, fillOut] = semshade(BB,0.2,'k',[1:1:200],1);
hold on
  [lineOut, fillOut] = semshade(BA,0.2,'m',[1:1:200],1);
      
sigB=[];
mB=[];

for n=1:200;
    [sigB(n),~]=ttest2(BA(:,n),BB(:,n));
    mB(n)=mean(BA(:,n))-mean(BB(:,n));
end

BAM=mean(BA,1);
BBM=mean(BB,1);

for n=1:length(sigB);
    if sigB(n)==1;
        if mB(n)<0;
            hold on
            plot(n,BAM(n),'r.','MarkerSize',10);
        else
            hold on
            plot(n,BBM(n),'r.','MarkerSize',10);
        end
    end
end
axis tight
title('bad leaners before after learning');

saveas(gcf,'fieldDistriDiffIndividualCellsGoodBad.fig');

%% frid fields numbers

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\gridCellsNewCueThresh\nFieldRealGrid.mat');
G=nFieldRealGrid;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\gridCellsNewCueThresh\nFieldRealGrid.mat');
B=nFieldRealGrid;

figure,
errorbar([1:1:11],nanmean(G,1),nansem(G,1),'g');
hold on
errorbar([1:1:11],nanmean(B,1),nansem(B,1),'m');
[pAnova,pMC] = anovaRM2W(G,B);
[r1,p1]=ttest2(G(:,1),G(:,2));
[r2,p2]=ttest2(B(:,1),B(:,2));
title(['nField p',num2str(pAnova),'p1=',num2str(p1),'p2=',num2str(p2)]);
saveas(gcf,'nField.fig');

%% 
%% dfof_sig of cue cells
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\allFSPNew.mat');
GmF=(mFSP);
GmDfof=(mDfofSP);
GmDfofSig=(mDfofSigSP);
GmDfofMean=(mDfofMean);
GmDfofSigMean=(mDfofSigMean);
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\indicesAllNewCueThresh.mat');
iG=indicesAllNewCueThresh.GridIdx;

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\allFSPNew.mat');
BmF=(mFSP);
BmDfof=(mDfofSP);
BmDfofSig=(mDfofSigSP);
BmDfofMean=(mDfofMean);
BmDfofSigMean=(mDfofSigMean);
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\indicesAllNewCueThresh.mat');
iB=indicesAllNewCueThresh.GridIdx;

GmDfofSigCue=GmDfofSig(iG,:);
BmDfofSigCue=BmDfofSig(iB,:);

A=GmDfofSigCue;
B=BmDfofSigCue;
[pAnova,pMC] = anovaRM2W(A,B);%p=0.0027
[r1,p1]=ttest(A(:,1),A(:,2));%p1=0.0028
[r2,p2]=ttest(B(:,1),B(:,2));%p2=4.8511e04

M1=mean(A,1);
S1=std(A,1)/sqrt(size(A,1));

M2=mean(B,1);
S2=std(B,1)/sqrt(size(B,1));
N=3;
for n=1:size(A,2);
    [sig(N,n),~]=ttest2(A(:,n),B(:,n));
end
figure
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')

xlabel('days')
ylabel('mean DfofSig');
title(['mean DfofSig grid cells']);

saveas(gcf,'dfofSigGridCells.fig');

%% plotting grid cell fluorescence change during learning

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\allFSPNew.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\indicesAllNewCueThresh.mat');
i=indicesAllNewCueThresh.GridIdx;

GmF=(mFSP(i,:));
GmDfof=(mDfofSP(i,:));
GmDfofSig=(mDfofSigSP(i,:));
GmDfofMean=(mDfofMean(i,:));
GmDfofSigMean=(mDfofSigMean(i,:));

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\allFSPNew.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\indicesAllNewCueThresh.mat');
i=indicesAllNewCueThresh.GridIdx;

BmF=(mFSP(i,:));
BmDfof=(mDfofSP(i,:));
BmDfofSig=(mDfofSigSP(i,:));
BmDfofMean=(mDfofMean(i,:));
BmDfofSigMean=(mDfofSigMean(i,:));
sig=[];

figure,
A=GmF;
B=BmF;
M1=mean(A,1);
S1=std(A,1)/sqrt(size(A,1));

M2=mean(B,1);
S2=std(B,1)/sqrt(size(B,1));
N=1;
for n=1:size(A,2);
    [sig(N,n),~]=ttest2(A(:,n),B(:,n));
end

subplot(151)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')
xlabel('days')
ylabel('mean F');
title(['mean F']);

A=GmDfof;
B=BmDfof;
M1=mean(A,1);
S1=std(A,1)/sqrt(size(A,1));

M2=mean(B,1);
S2=std(B,1)/sqrt(size(B,1));
N=2;
for n=1:size(A,2);
    [sig(N,n),~]=ttest2(A(:,n),B(:,n));
end
subplot(152)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')
xlabel('days');
ylabel('mean Dfof');
title(['mean Dfof']);

A=GmDfofSig;
B=BmDfofSig;
[pAnova,pMC] = anovaRM2W(A,B);%p=0.0027
[r1,p1]=ttest(A(:,1),A(:,2));%p1=0.0028
[r2,p2]=ttest(B(:,1),B(:,2));%p2=4.8511e04

M1=mean(A,1);
S1=std(A,1)/sqrt(size(A,1));

M2=mean(B,1);
S2=std(B,1)/sqrt(size(B,1));
N=3;
for n=1:size(A,2);
    [sig(N,n),~]=ttest2(A(:,n),B(:,n));
end
subplot(153)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')

xlabel('days')
ylabel('mean DfofSig');
title(['mean DfofSig']);

A=GmDfofMean;
B=BmDfofMean;
M1=mean(A,1);
S1=std(A,1)/sqrt(size(A,1));

M2=mean(B,1);
S2=std(B,1)/sqrt(size(B,1));
N=4;
for n=1:size(A,2);
    [sig(N,n),~]=ttest2(A(:,n),B(:,n));
end
subplot(154)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')
xlabel('days');
ylabel('mean DfofMean');
title(['mean DfofMean']);

A=GmDfofSigMean;
B=BmDfofSigMean;
M1=mean(A,1);
S1=std(A,1)/sqrt(size(A,1));

M2=mean(B,1);
S2=std(B,1)/sqrt(size(B,1));
N=5;
for n=1:size(A,2);
    [sig(N,n),~]=ttest2(A(:,n),B(:,n));
end
subplot(155)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')

xlabel('days')
ylabel('mean DfofSigMean');
title(['mean DfofSigMean']);

tightfig;

%
saveas(gcf,'allFCommonCellsSPNew_Grid.fig');
