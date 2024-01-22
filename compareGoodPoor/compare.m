%compareing good and bad in different aspects

%% compare run by run consistency: use sig: run by run to mean
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\corrAllFOVs.mat');
G=corrAllFOVs;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\corrAllFOVs.mat');
B=corrAllFOVs;

%only use dfofSig Mean To next and do shuffle

figure

A=G.dfofSig{3};
MTempRL=G.dfofSigMeanToOthers;
S1=std(A,1)/sqrt(size(A,1));
[~,p1]=corr([1:1:10]',MTempRL(2:end)','tail','right');
% %two tail
% [~,p1]=corr([1:1:10]',MTempRL(2:end)');%p1= 0.0396


errorbar([1:1:length(MTempRL)],MTempRL,S1,'g');
g=A;

A=B.dfofSig{3};
MTempRL=B.dfofSigMeanToOthers;
S1=std(A,1)/sqrt(size(A,1));
[~,p2]=corr([1:1:10]',MTempRL(2:end)','tail','left');

% %two tail
% [~,p2]=corr([1:1:10]',MTempRL(2:end)');%p2= 0.0396

hold on
errorbar([1:1:length(MTempRL)],MTempRL,S1,'m');
b=A;

[pAnova,pMC] = anovaRM2W(g,b)

sig=[];
for n=1:length(MTempRL);
    g=G.dfofSig{3}(:,n);
    b=B.dfofSig{3}(:,n);
    [sig(n,1),~]=ttest2(g,b);
end

title(['good p=',num2str(p1),'bad p=',num2str(p2)])

saveas(gcf,'corrAbsValues.fig');

%% compare with non tracked cells

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\corrAllFOVs.mat');
GC=corrAllFOVs;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\corrAllFOVUncommon.mat');
GU=corrAllFOVUncommon;

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\corrAllFOVs.mat');
BC=corrAllFOVs;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\corrAllFOVUncommon.mat');
BU=corrAllFOVUncommon

figure
A=GC.dfofSig{3};
MTempRL=GC.dfofSigMeanToOthers;
S1=std(A,1)/sqrt(size(A,1));

hold on
errorbar([1:1:length(MTempRL)],MTempRL,S1,'g-');
hold on
errorbar([1:1:length(MTempRL)],GU.dfofSigMeanToNext,GU.dfofSigMeanToNext_sem,'g--')
[r,p]=corr([1:1:length(MTempRL)-1]',MTempRL(2:end)');
[r,p]=corr([1:1:length(MTempRL)-1]',GU.dfofSigMeanToNext(2:end)');

A=BC.dfofSig{3};
MTempRL=BC.dfofSigMeanToOthers;
S1=std(A,1)/sqrt(size(A,1));

hold on
errorbar([1:1:length(MTempRL)],MTempRL,S1,'m-');
hold on
errorbar([1:1:length(MTempRL)],BU.dfofSigMeanToNext,BU.dfofSigMeanToNext_sem,'m--')
[r,p]=corr([1:1:length(MTempRL)-1]',MTempRL(2:end)');
[r,p]=corr([1:1:length(MTempRL)-1]',BU.dfofSigMeanToNext(2:end)');

title('commonNonCommonCorr');
xlabel('day');
ylabel('corr');
xlim([0.5 11.5]);
saveas(gcf,'commonNonCommonCorr.fig');

%% compare NORMALIZED run by run consistency: use sig: run by run to mean
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\corrAllFOVs.mat');
G=corrAllFOVs;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\corrAllFOVs.mat');
B=corrAllFOVs;

%only use dfofSig Mean To Mean and do shuffle

figure
A=G.dfofSigNorm{3};
MTempRL=G.dfofSigNormMeanToOthers;
S1=std(A,1)/sqrt(size(A,1));
[~,p1]=corr([1:1:10]',MTempRL(2:end)');
errorbar([1:1:length(MTempRL)],MTempRL,S1,'g');

A=B.dfofSigNorm{3};
MTempRL=B.dfofSigNormMeanToOthers;
S1=std(A,1)/sqrt(size(A,1));
[~,p2]=corr([1:1:10]',MTempRL(2:end)');
hold on
errorbar([1:1:length(MTempRL)],MTempRL,S1,'m');
sig=[];
for n=1:length(MTempRL);
    g=G.dfofSigNorm{1}(:,n);
    b=B.dfofSigNorm{1}(:,n);
    [sig(n,1),~]=ttest2(g,b);
end

title(['good p=',num2str(p1),'bad p=',num2str(p2)])
xlim([0.5 11.5]);
saveas(gcf,'corrNormValue.fig');


%%
% %look at the correlation of this with behavior, use 10 run data
% 
% %activity
% act=[G.dfofSigMeanToOthers B.dfofSigMeanToOthers];
% load('E:\learningAnalysis\behaviorNewBatch\predLickDataBad.mat');
% load('E:\learningAnalysis\behaviorNewBatch\predLickDataGood.mat');
% lick=[predLick10RunsMeanGood(1:11) predLick10RunsMeanBad(1:11)];
% load('E:\learningAnalysis\behaviorNewBatch\slowDownDataBad.mat');
% load('E:\learningAnalysis\behaviorNewBatch\slowDownDataGood.mat');
% slow=[slowDown10RunsMeanGood(1:11) slowDown10RunsMeanBad(1:11)];
% 
% figure,
% subplot(231)
% plot(act,lick,'r.')
% [~,p]=corr(act',lick');
% title(['all lick p=',num2str(p)])
% 
% subplot(234)
% plot(act,slow,'r.')
% [~,p]=corr(act',slow');
% title(['all slow p=',num2str(p)])
% 
% subplot(232)
% plot(act(1:11),lick(1:11),'r.')
% [~,p]=corr(act(1:11)',lick(1:11)');
% title(['Good lick p=',num2str(p)])
% 
% subplot(235)
% plot(act(1:11),slow(1:11),'r.')
% [~,p]=corr(act(1:11)',slow(1:11)');
% title(['Good slow p=',num2str(p)])
% 
% subplot(233)
% plot(act(12:end),lick(12:end),'r.')
% [~,p]=corr(act(12:end)',lick(12:end)');
% title(['Good lick p=',num2str(p)])
% 
% subplot(236)
% plot(act(12:end),slow(12:end),'r.')
% [~,p]=corr(act(12:end)',slow(12:end)');
% title(['Good slow p=',num2str(p)])
% 
% saveas(gcf,'absCorrAndBehavior.fig');

%% look at matrix corrlations

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\dfofSigNormAllFOVsSort_Max1.mat');
G=dfofSigNormAllFOVsSort_Max1;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\dfofSigNormAllFOVsSort_Max1.mat');
B=dfofSigNormAllFOVsSort_Max1;
figure
for n=1:length(G);
    subplot(2,length(G),n);
    imagesc(G{n})
    
    title(['Day',num2str(n)]);
    axis off
end
for n=1:length(B);
    subplot(2,length(G),length(G)+n);
    imagesc(B{n})
    
%     title(['Day',num2str(n)]);
    axis off
end

tightfig;
saveas(gcf,'sortedByMaxOld_fig');
%%

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\dfofSigNormAllFOVsSort_Max2.mat');
G=dfofSigNormAllFOVsSort_Max2;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\dfofSigNormAllFOVsSort_Max2.mat');
B=dfofSigNormAllFOVsSort_Max2;
figure
for n=1:length(G);
    subplot(2,length(G),n);
    imagesc(G{n})
    
    title(['Day',num2str(n)]);
    axis off
end
for n=1:length(B);
    subplot(2,length(G),length(G)+n);
    imagesc(B{n})
    
%     title(['Day',num2str(n)]);
    axis off
end

tightfig;
saveas(gcf,'sortedByMaxNew1_fig');

%%

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\dfofSigNormAllFOVsSort_Max11.mat');
G=dfofSigNormAllFOVsSort_Max11;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\dfofSigNormAllFOVsSort_Max11.mat');
B=dfofSigNormAllFOVsSort_Max11;
figure
for n=1:length(G);
    subplot(2,length(G),n);
    imagesc(G{n})
    
    title(['Day',num2str(n)]);
    axis off
end
for n=1:length(B);
    subplot(2,length(G),length(G)+n);
    imagesc(B{n})
    
%     title(['Day',num2str(n)]);
    axis off
end

tightfig;
saveas(gcf,'sortedByMaxNew11_fig');

%%
%plotting all of them on the same figure
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\dfofSigNormAllFOVsSort_Max1.mat');
G1=dfofSigNormAllFOVsSort_Max1;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\dfofSigNormAllFOVsSort_Max1.mat');
B1=dfofSigNormAllFOVsSort_Max1;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\dfofSigNormAllFOVsSort_Max2.mat');
G2=dfofSigNormAllFOVsSort_Max2;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\dfofSigNormAllFOVsSort_Max2.mat');
B2=dfofSigNormAllFOVsSort_Max2;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\dfofSigNormAllFOVsSort_Max11.mat');
G11=dfofSigNormAllFOVsSort_Max11;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\dfofSigNormAllFOVsSort_Max11.mat');
B11=dfofSigNormAllFOVsSort_Max11;

N=length(G1);
figure
for n=1:N;
    subplot(6,N,n);
    imagesc(G1{n});
    axis off
    subplot(6,N,n+N);
    imagesc(G2{n});
     axis off
     subplot(6,N,n+2*N);
    imagesc(G11{n});
     axis off
     
     subplot(6,N,n+3*N);
    imagesc(B1{n});
    axis off
    subplot(6,N,n+4*N);
    imagesc(B2{n});
     axis off
     subplot(6,N,n+5*N);
    imagesc(B11{n});
     axis off
end
tightfig;
saveas(gcf,'sortedActivityAll.fig');



%% matrix correlation: left whole matrix, right: individual cells
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\matrixCorrMax1.mat');
G=matrixCorr;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\matrixCorrMax1.mat');
B=matrixCorr;

figure
subplot(121);
plot(G,'g.','MarkerSize',10);
hold on
plot([2:1:length(G)],G(2:end),'g');
[r,p1]=corr(G(2:end),[1:1:length(G)-1]','tail','right');%p1=1.8383e-05
% do two tail
[r,p1]=corr(G(2:end),[1:1:length(G)-1]');%p1= 3.6766e-05

hold on
plot(B,'m.','MarkerSize',10);
hold on
plot([2:1:length(B)],B(2:end),'m');
xlim([0 11])
ylim([-0.05 0.7])
breakyaxis([0.05 0.45]);
[r,p2]=corr(B(2:end),[1:1:length(B)-1]');
title(['matrix corr good p=',num2str(p1),'bad p=',num2str(p2)]);

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\corrIndivi.mat');
G=corrIndivi;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\corrIndivi.mat');
B=corrIndivi;

subplot(122)
M1=mean(G,1);
S1=std(G,1)/sqrt(size(G,1));
hold on
plot([1:1:length(M1)],M1,'g.','MarkerSize',10);
hold on
errorbar([1:1:length(M1)],M1(1:end),S1(1:end),'g.')
hold on
plot([2:1:length(M1)],M1(2:end),'g-');

[r,p1]=corr([1:1:length(M1)-1]', M1(2:end)','tail','right');%p1=9.4462e-06
% do two tail
[r,p1]=corr([1:1:length(M1)-1]', M1(2:end)');%p1=1.8892e-05

M1=mean(B,1);
S1=std(B,1)/sqrt(size(B,1));
hold on
plot([1:1:length(M1)],M1,'m.','MarkerSize',10);
hold on
errorbar([1:1:length(M1)],M1(1:end),S1(1:end),'m.')
hold on
plot([2:1:length(M1)],M1(2:end),'m-');

[r,p2]=corr([1:1:length(M1)-1]', M1(2:end)');
title(['indiv corr good p=',num2str(p1),'bad p=',num2str(p2)]);
xlim([0 11])
ylim([-0.1 0.7])
breakyaxis([0 0.45]);

sig=[];
for n=1:size(G,2);
    g=G(:,n);
    b=B(:,n);
    [sig(n,1),~]=ttest2(g,b);
end

saveas(gcf,'matrixAndIndiviCorr.fig');
%% field distribution
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\fieldDistriEachBin.mat');
G=fieldDistriEachBin;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\fieldDistriEachBin.mat');
B=fieldDistriEachBin;
load('E:\ID20201118\20210103\1stloc1\pcaica\cueAnalysisNew_sig\tempRL.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\dfofSigNormAllFOVs.mat');
NG=size(dfofSigNormAllFOVs{1},1);
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\dfofSigNormAllFOVs.mat');
NB=size(dfofSigNormAllFOVs{1},1);


figure,
subplot(211);
plot([1:1:200],(mean(G(4:11,:),1)-mean(G(2:3,:),1))/NG,'g')
hold on
plot([1:1:200],(mean(B(4:11,:),1)-mean(B(2:3,:),1))/NB,'m')

plot([1:1:200],tempRL*0.06);
hold on
plot([890/5 890/5],[0 0.06],'r');
xlim([0 200]);
title('fieldDistri');

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\fieldDistriAmpSum.mat');
G=fieldDistriAmpSum;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\fieldDistriAmpSum.mat');
B=fieldDistriAmpSum;

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

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\fieldDistriEachBin.mat');
G=fieldDistriEachBin;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\sigIdxFieldDistri.mat');
GSig=sigIdxFieldDistri;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\fieldDistriEachBin.mat');
B=fieldDistriEachBin;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\sigIdxFieldDistri.mat');
BSig=sigIdxFieldDistri;
load('E:\ID20201118\20210103\1stloc1\pcaica\cueAnalysisNew_sig\tempRL.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\dfofSigNormAllFOVs.mat');
NG=size(dfofSigNormAllFOVs{1},1);
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\dfofSigNormAllFOVs.mat');
NB=size(dfofSigNormAllFOVs{1},1);

figure,
a=(mean(G(4:11,:),1)-mean(G(2:3,:),1))/NG;
plot([1:1:200],a,'g')
hold on
plot(GSig,a(GSig),'r.','MarkerSize',10);

a=(mean(B(4:11,:),1)-mean(B(2:3,:),1))/NG;
plot([1:1:200],a,'m')
hold on
plot(BSig,a(BSig),'r.','MarkerSize',10);

plot([1:1:200],tempRL*0.06);
hold on
plot([890/5 890/5],[0 0.06],'r');
xlim([0 200]);
title('fieldDistri');
saveas(gcf,'fieldDistri__beforeAfterLearn_sig.fig');

%% plotting grid cells multiple sessions: REMOVE THE LAST TWO: BECAUSE ONLY ONE CELL
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\indicesAllNewCueThresh.mat');
G=indicesAllNewCueThresh;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\indicesAllNewCueThresh.mat');
B=indicesAllNewCueThresh;
NG=size(G.Grid,1);
BG=size(B.Grid,1);

%GRID
GG=G.GridIdxLowerThreshs;
BB=B.GridIdxLowerThreshs;
GGL=[];
BBL=[];
for n=1:length(GG);
    GGL(n)=length(GG{n})/NG;
   BBL(n)=length(BB{n})/BG;
end


figure
subplot(251)
plot([1:1:length(GGL)],GGL,'g-','LineWidth',1);
hold on
plot([1:1:length(GGL)],GGL,'g.','MarkerSize',10);
hold on
plot([1:1:length(BBL)],BBL,'m-','LineWidth',1);
plot([1:1:length(BBL)],BBL,'m.','MarkerSize',10);
[r,p]=ttest(GGL,BBL,'tail','right');
title(['grid paired t, p=', num2str(p)]);
xlim([0.5 11.5]);



%CUE L  
GG=G.CueLIdxLowerThreshs;
BB=B.CueLIdxLowerThreshs;
GGL=[];
BBL=[];
for n=1:length(GG);
    GGL(n)=length(GG{n})/NG;
   BBL(n)=length(BB{n})/BG;
end


subplot(253)
plot([1:1:length(GGL)],GGL,'g-','LineWidth',1);
hold on
plot([1:1:length(GGL)],GGL,'g.','MarkerSize',10);
hold on
plot([1:1:length(BBL)],BBL,'m-','LineWidth',1);
plot([1:1:length(BBL)],BBL,'m.','MarkerSize',10);
[r,p]=ttest(GGL,BBL,'tail','right');
title(['cueL paired t, p=', num2str(p)]);
xlim([0.5 11.5]);

%CUE R
GG=G.CueRIdxLowerThreshs;
BB=B.CueRIdxLowerThreshs;
GGL=[];
BBL=[];
for n=1:length(GG);
    GGL(n)=length(GG{n})/NG;
   BBL(n)=length(BB{n})/BG;
end

subplot(252)
plot([1:1:length(GGL)],GGL,'g-','LineWidth',1);
hold on
plot([1:1:length(GGL)],GGL,'g.','MarkerSize',10);
hold on
plot([1:1:length(BBL)],BBL,'m-','LineWidth',1);
plot([1:1:length(BBL)],BBL,'m.','MarkerSize',10);
[r,p]=ttest(GGL,BBL,'tail','right');
title(['cueR paired t, p=', num2str(p)]);
xlim([0.5 11.5]);


%speed p
GG=G.SpeedPIdxLowerThreshs;
BB=B.SpeedPIdxLowerThreshs;
GGL=[];
BBL=[];
for n=1:length(GG);
    GGL(n)=length(GG{n})/NG;
   BBL(n)=length(BB{n})/BG;
end

subplot(254)
plot([1:1:length(GGL)],GGL,'g-','LineWidth',1);
hold on
plot([1:1:length(GGL)],GGL,'g.','MarkerSize',10);
hold on
plot([1:1:length(BBL)],BBL,'m-','LineWidth',1);
plot([1:1:length(BBL)],BBL,'m.','MarkerSize',10);
[r,p]=ttest(GGL,BBL,'tail','right');
title(['speedP paired t, p=', num2str(p)]);
xlim([0.5 11.5]);

%speed N
GG=G.SpeedNIdxLowerThreshs;
BB=B.SpeedNIdxLowerThreshs;
GGL=[];
BBL=[];
for n=1:length(GG);
    GGL(n)=length(GG{n})/NG;
   BBL(n)=length(BB{n})/BG;
end

subplot(255)
plot([1:1:length(GGL)],GGL,'g-','LineWidth',1);
hold on
plot([1:1:length(GGL)],GGL,'g.','MarkerSize',10);
hold on
plot([1:1:length(BBL)],BBL,'m-','LineWidth',1);
plot([1:1:length(BBL)],BBL,'m.','MarkerSize',10);
[r,p]=ttest(GGL,BBL,'tail','right');
title(['speedN paired t, p=', num2str(p)]);
xlim([0.5 11.5]);


%% go through individual FOVs
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\indicesAllNewCueThresh.mat');
G=indicesAllNewCueThresh;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\indicesAllNewCueThresh.mat');
B=indicesAllNewCueThresh;


GG=G.perGridAllThreshIndivFOVs(1:23,:);
BB=B.perGridAllThreshIndivFOVs;
subplot(256)
errorbar([1:1:11],mean(GG,1),std(GG,1)/sqrt(size(GG,1)),'g');
hold on
errorbar([1:1:11],mean(BB,1),std(BB,1)/sqrt(size(BB,1)),'m');
xlim([0.5 11.5]);
sigG=[];
for n=1:size(GG,2);
    [~,p]=ttest2(GG(:,n),BB(:,n),'tail','right');
    if p<0.05;
        sigG(n,1)=1;
    else
        sigG(n,1)=0;
    end
end
[r,p]=ttest(mean(GG,1),mean(BB,1),'tail','right');
title(['grid indiv mean paired p=',num2str(p)]);


GG=G.perCueRAllThreshIndivFOVs(1:23,:);
BB=B.perCueRAllThreshIndivFOVs;
subplot(257)
errorbar([1:1:11],mean(GG,1),std(GG,1)/sqrt(size(GG,1)),'g');
hold on
errorbar([1:1:11],mean(BB,1),std(BB,1)/sqrt(size(BB,1)),'m');
xlim([0.5 11.5]);
sigCR=[];
for n=1:size(GG,2);
    [~,p]=ttest2(GG(:,n),BB(:,n),'tail','left');
    if p<0.05;
        sigCR(n,1)=1;
    else
        sigCR(n,1)=0;
    end
end

[r,p]=ttest(mean(GG,1),mean(BB,1),'tail','right');
title(['cueR indiv mean paired p=',num2str(p)]);



GG=G.perCueLAllThreshIndivFOVs(1:23,:);
BB=B.perCueLAllThreshIndivFOVs;
[pAnova,pMC] = anovaRM2W(GG,BB)%P= 0.3344
subplot(258)
errorbar([1:1:11],mean(GG,1),std(GG,1)/sqrt(size(GG,1)),'g');
hold on
errorbar([1:1:11],mean(BB,1),std(BB,1)/sqrt(size(BB,1)),'m');
xlim([0.5 11.5]);
sigCL=[];
for n=1:size(GG,2);
    [~,p]=ttest2(GG(:,n),BB(:,n));
    if p<0.05;
        sigCL(n,1)=1;
    else
        sigCL(n,1)=0;
    end
end
[r,p]=ttest(mean(GG,1),mean(BB,1),'tail','right');
title(['cueL indiv mean paired p=',num2str(p)]);

GG=G.perSpeedPAllThreshIndivFOVs(1:23,:);
BB=B.perSpeedPAllThreshIndivFOVs;
subplot(259)
errorbar([1:1:11],mean(GG,1),std(GG,1)/sqrt(size(GG,1)),'g');
hold on
errorbar([1:1:11],mean(BB,1),std(BB,1)/sqrt(size(BB,1)),'m');
xlim([0.5 11.5]);
sigCL=[];
for n=1:size(GG,2);
    [~,p]=ttest2(GG(:,n),BB(:,n));
    if p<0.05;
        sigCL(n,1)=1;
    else
        sigCL(n,1)=0;
    end
end
[r,p]=ttest(mean(GG,1),mean(BB,1),'tail','right');
title(['SpeedP indiv mean paired p=',num2str(p)]);

GG=G.perSpeedNAllThreshIndivFOVs(1:23,:);
BB=B.perSpeedNAllThreshIndivFOVs;
subplot(2,5,10)
errorbar([1:1:11],mean(GG,1),std(GG,1)/sqrt(size(GG,1)),'g');
hold on
errorbar([1:1:11],mean(BB,1),std(BB,1)/sqrt(size(BB,1)),'m');
xlim([0.5 11.5]);
sigCL=[];
for n=1:size(GG,2);
    [~,p]=ttest2(GG(:,n),BB(:,n));
    if p<0.05;
        sigCL(n,1)=1;
    else
        sigCL(n,1)=0;
    end
end
[r,p]=ttest(mean(GG,1),mean(BB,1),'tail','right');
title(['SpeedN indiv mean paired p=',num2str(p)]);

saveas(gcf,'stableCells.fig');

%using the above data to compare the cells identified in different
%sessions: at least one session, or at least 6 sessions
%%
figure,

A1=G.perGridAllThreshIndivFOVs(1:23,1);
A2=B.perGridAllThreshIndivFOVs(:,1);
A3=G.perGridAllThreshIndivFOVs(1:23,6);
A4=B.perGridAllThreshIndivFOVs(:,6);

M(1)=mean(A1);
M(2)=mean(A2);
M(3)=mean(A3);
M(4)=mean(A4);
S(1)=nansem(A1,1);
S(2)=nansem(A2,1);
S(3)=nansem(A3,1);
S(4)=nansem(A4,1);

[r,p1]=ttest2(A1,A2);

[r,p2]=ttest2(A3,A4);
subplot(121)
bar([1:1:4],M);
hold on
errorbar([1:1:4],M,S,'.');
title(['grid p1=',num2str(p1),' p2=',num2str(p2)])

A1=G.perCueRAllThreshIndivFOVs(1:23,1);
A2=B.perCueRAllThreshIndivFOVs(:,1);
A3=G.perCueRAllThreshIndivFOVs(1:23,6);
A4=B.perCueRAllThreshIndivFOVs(:,6);

M(1)=mean(A1);
M(2)=mean(A2);
M(3)=mean(A3);
M(4)=mean(A4);
S(1)=nansem(A1,1);
S(2)=nansem(A2,1);
S(3)=nansem(A3,1);
S(4)=nansem(A4,1);

[r,p1]=ttest2(A1,A2);

[r,p2]=ttest2(A3,A4);
subplot(122)
bar([1:1:4],M);
hold on
errorbar([1:1:4],M,S,'.');
title(['cueR p1=',num2str(p1),' p2=',num2str(p2)])
saveas(gcf,'stableCells_1timeAnd6Times.fig');

%% just plot the percentage of cells identified six times for grid, cueR, cueL
A1=G.perGridAllThreshIndivFOVs(1:23,6);
A2=B.perGridAllThreshIndivFOVs(:,6);
A3=G.perCueRAllThreshIndivFOVs(1:23,6);
A4=B.perCueRAllThreshIndivFOVs(:,6);
A5=G.perCueLAllThreshIndivFOVs(1:23,6);
A6=B.perCueLAllThreshIndivFOVs(:,6);

M(1)=mean(A1);
M(2)=mean(A2);
M(3)=mean(A3);
M(4)=mean(A4);
M(5)=mean(A5);
M(6)=mean(A6);
S(1)=nansem(A1,1);
S(2)=nansem(A2,1);
S(3)=nansem(A3,1);
S(4)=nansem(A4,1);
S(5)=nansem(A5,1);
S(6)=nansem(A6,1);
figure
bar([1:1:6],M);
hold on
errorbar([1:1:6],M,S,'k.');
[r,p1]=ttest2(A1,A2);
[r,p2]=ttest2(A3,A4);
[r,p3]=ttest2(A5,A6);
[r,p4]=ttest2(A3,A5);
[r,p5]=ttest2(A4,A6);

title(['Grid cueR cueL p1=',num2str(p1),' p2=',num2str(p2),' p3=',num2str(p3)])

saveas(gcf,'stableCells_6Times.fig');
%% plotting spatial information
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\SI.mat');
GSI1=(SI1);
GSI2=(SI2);
% GSI3=(SI3);
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\SI.mat');
BSI1=(SI1);
BSI2=(SI2);
% BSI3=(SI3);

figure,
A=GSI1;
B=BSI1;
M1=mean(A,1);
S1=std(A,1)/sqrt(size(A,1));

M2=mean(B,1);
S2=std(B,1)/sqrt(size(B,1));

subplot(121)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')
xlabel('days')
ylabel('spatial information');
title(['SI original thresh']);

A=GSI2;
B=BSI2;
M1=mean(A,1);
S1=std(A,1)/sqrt(size(A,1));

M2=mean(B,1);
S2=std(B,1)/sqrt(size(B,1));

subplot(122)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')
xlabel('days');
ylabel('spatial information');
title(['SI 1cm/s']);

% A=GSI3;
% B=BSI3;
% M1=mean(A,1);
% S1=std(A,1)/sqrt(size(A,1));
% 
% M2=mean(B,1);
% S2=std(B,1)/sqrt(size(B,1));
% 
% subplot(133)
% errorbar([1:1:11],M1,S1,'g')
% hold on
% errorbar([1:1:11],M2,S2,'m')
% 
% xlabel('days')
% ylabel('spatial information');
% title(['SI new 1cm/s']);

saveas(gcf,'spatialInformation.fig');

%% plotting fluorscence
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\SI.mat');
GSI1=(SI1);
GSI2=(SI2);
GSI3=(SI3);
GSI4=(SI4);

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\SI.mat');
BSI1=(SI1);
BSI2=(SI2);
BSI3=(SI3);
BSI4=(SI4);

% GSI2N=[];
% for n=1:size(GSI2,1);
%     a=GSI2(n,:);
%     GSI2N(n,:)=a/a(2);
% end
% 
% BSI2N=[];
% for n=1:size(BSI2,1);
%     a=BSI2(n,:);
%     BSI2N(n,:)=a/a(2);
% end


figure,
A=GSI1;
B=BSI1;
M1=mean(A,1);
S1=std(A,1)/sqrt(size(A,1));

M2=mean(B,1);
S2=std(B,1)/sqrt(size(B,1));

subplot(221)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')
xlabel('days')
ylabel('spatial information');
title(['SI method1']);

A=GSI2;
B=BSI2;
M1=mean(A,1);
S1=std(A,1)/sqrt(size(A,1));

M2=mean(B,1);
S2=std(B,1)/sqrt(size(B,1));

subplot(222)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')
xlabel('days');
ylabel('spatial information');
title(['SI method1 1cms']);


A=GSI3;
B=BSI3;
M1=mean(A,1);
S1=std(A,1)/sqrt(size(A,1));

M2=mean(B,1);
S2=std(B,1)/sqrt(size(B,1));

subplot(223)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')
xlabel('days')
ylabel('spatial information');
title(['SI method2']);

A=GSI4;
B=BSI4;
M1=mean(A,1);
S1=std(A,1)/sqrt(size(A,1));

M2=mean(B,1);
S2=std(B,1)/sqrt(size(B,1));

subplot(224)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')
xlabel('days');
ylabel('spatial information');
title(['SI method2 1cms']);

saveas(gcf,'spatialInformation.fig');

%% plotting fluorscence


load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\allF.mat');
GmF=(mF);
GmDfof=(mDfof);
GmDfofSig=(mDfofSig);
GmDfofMean=(mDfofMean);
GmDfofSigMean=(mDfofSigMean);

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\allF.mat');
BmF=(mF);
BmDfof=(mDfof);
BmDfofSig=(mDfofSig);
BmDfofMean=(mDfofMean);
BmDfofSigMean=(mDfofSigMean);
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
M1=mean(A,1);
S1=std(A,1)/sqrt(size(A,1));

M2=mean(B,1);
S2=std(B,1)/sqrt(size(B,1));
N=3;
for n=1:size(A,2);
    [sig(N,n),~]=ttest2(A(:,n),B(:,n),'tail','right');
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
saveas(gcf,'allFCommonCells.fig');
%% normlize dfof using tbe last day in env

A=[];
B=GmDfofSig;
for n=1:size(B,1);
    A(n,:)=B(n,:)/B(n,1);
end

M=nanmean(A,1);
S=nansem(A,1);

figure, errorbar([1:1:11],M,S,'g');

A=[];
B=BmDfofSig;
for n=1:size(B,1);
    A(n,:)=B(n,:)/B(n,1);
end

M=nanmean(A,1);
S=nansem(A,1);
hold on

errorbar([1:1:11],M,S,'m');
title('normDfofSig common speedThresh');

saveas(gcf,'normDfofSig_common.fig')

%% plotting fluorscence


load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\allFSP.mat');
GmF=(mFSP);
GmDfof=(mDfofSP);
GmDfofSig=(mDfofSigSP);
GmDfofMean=(mDfofMean);
GmDfofSigMean=(mDfofSigMean);

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\allFSP.mat');
BmF=(mFSP);
BmDfof=(mDfofSP);
BmDfofSig=(mDfofSigSP);
BmDfofMean=(mDfofMean);
BmDfofSigMean=(mDfofSigMean);
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
saveas(gcf,'allFCommonCellsSP.fig');
%% normlize dfof using tbe last day in env

A=[];
B=GmDfofSig;
for n=1:size(B,1);
    A(n,:)=B(n,:)/B(n,1);
end

M=nanmean(A,1);
S=nansem(A,1);

figure, errorbar([1:1:11],M,S,'g');

A=[];
B=BmDfofSig;
for n=1:size(B,1);
    A(n,:)=B(n,:)/B(n,1);
end

M=nanmean(A,1);
S=nansem(A,1);
hold on

errorbar([1:1:11],M,S,'m');
title('normDfofSig common speedThresh');

saveas(gcf,'normDfofSig_common_speedThresh.fig')


%% plotting fluorscence


load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\allFSPNew.mat');
GmF=(mFSP);
GmDfof=(mDfofSP);
GmDfofSig=(mDfofSigSP);
GmDfofMean=(mDfofMean);
GmDfofSigMean=(mDfofSigMean);

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\allFSPNew.mat');
BmF=(mFSP);
BmDfof=(mDfofSP);
BmDfofSig=(mDfofSigSP);
BmDfofMean=(mDfofMean);
BmDfofSigMean=(mDfofSigMean);
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
saveas(gcf,'allFCommonCellsSPNew.fig');


%% normlize dfof using tbe last day in env

A=[];
B=GmDfofSig;
for n=1:size(B,1);
    A(n,:)=B(n,:)/B(n,1);
end

sigG(1)=nan;

for n=2:11;
    [~,sigG(n)]=ttest(A(:,1),A(:,n));
end
M=nanmean(A,1);
S=nansem(A,1);

figure, errorbar([1:1:11],M,S,'g');

A=[];
B=BmDfofSig;
for n=1:size(B,1);
    A(n,:)=B(n,:)/B(n,1);
end
sigB(1)=nan;

for n=2:11;
    [~,sigB(n)]=ttest(A(:,1),A(:,n));
end
M=nanmean(A,1);
S=nansem(A,1);
hold on

errorbar([1:1:11],M,S,'m');
title('normDfofSig common speedThresh');

saveas(gcf,'normDfofSig_common_speedThreshNew.fig')


%% plotting fluorscence


load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\allFAllCells.mat');
GmF=(mFAllCells);
GmDfof=(mDfofAllCells);
GmDfofSig=(mDfofSigAllCells);
GmDfofMean=(mDfofMeanAllCells);
GmDfofSigMean=(mDfofSigMeanAllCells);

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\allFAllCells.mat');
BmF=(mFAllCells);
BmDfof=(mDfofAllCells);
BmDfofSig=(mDfofSigAllCells);
BmDfofMean=(mDfofMeanAllCells);
BmDfofSigMean=(mDfofSigMeanAllCells);
sig=[];

figure,
A=GmF;
B=BmF;

M1=[];
S1=[];
M2=[];
S2=[];
for n=1:length(A);
M1(n)=mean(A{n});
S1(n)=nansem(A{n},1);
M2(n)=mean(B{n});
S2(n)=nansem(B{n},1);
end

N=1;
for n=1:length(A);
    [sig(N,n),~]=ttest2(A{n},B{n});
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
M1=[];
S1=[];
M2=[];
S2=[];
for n=1:length(A);
M1(n)=mean(A{n});
S1(n)=nansem(A{n},1);
M2(n)=mean(B{n});
S2(n)=nansem(B{n},1);
end
N=2;
for n=1:length(A);
    [sig(N,n),~]=ttest2(A{n},B{n});
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
M1=[];
S1=[];
M2=[];
S2=[];
for n=1:length(A);
M1(n)=mean(A{n});
S1(n)=nansem(A{n},1);
M2(n)=mean(B{n});
S2(n)=nansem(B{n},1);
end
N=3;
for n=1:length(A);
    [sig(N,n),~]=ttest2(A{n},B{n});
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
M1=[];
S1=[];
M2=[];
S2=[];
for n=1:length(A);
M1(n)=mean(A{n});
S1(n)=nansem(A{n},1);
M2(n)=mean(B{n});
S2(n)=nansem(B{n},1);
end
N=4;
for n=1:length(A);
    [sig(N,n),~]=ttest2(A{n},B{n});
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
M1=[];
S1=[];
M2=[];
S2=[];
for n=1:length(A);
M1(n)=mean(A{n});
S1(n)=nansem(A{n},1);
M2(n)=mean(B{n});
S2(n)=nansem(B{n},1);
end
N=5;
for n=1:length(A);
    [sig(N,n),~]=ttest2(A{n},B{n});
end

subplot(155)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')

xlabel('days')
ylabel('mean DfofSigMean');
title(['mean DfofSigMean']);

tightfig;
saveas(gcf,'allFAllCells.fig');


%% plotting fluorscence


load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\allFAllCellsSP.mat');
GmF=(mFAllCellsSP);
GmDfof=(mDfofAllCellsSP);
GmDfofSig=(mDfofSigAllCellsSP);
GmDfofMean=(mDfofMeanAllCells);
GmDfofSigMean=(mDfofSigMeanAllCells);

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\allFAllCellsSP.mat');
BmF=(mFAllCellsSP);
BmDfof=(mDfofAllCellsSP);
BmDfofSig=(mDfofSigAllCellsSP);
BmDfofMean=(mDfofMeanAllCells);
BmDfofSigMean=(mDfofSigMeanAllCells);
sig=[];

figure,
A=GmF;
B=BmF;

M1=[];
S1=[];
M2=[];
S2=[];
for n=1:length(A);
M1(n)=mean(A{n});
S1(n)=nansem(A{n},1);
M2(n)=mean(B{n});
S2(n)=nansem(B{n},1);
end

N=1;
for n=1:length(A);
    [sig(N,n),~]=ttest2(A{n},B{n});
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
M1=[];
S1=[];
M2=[];
S2=[];
for n=1:length(A);
M1(n)=mean(A{n});
S1(n)=nansem(A{n},1);
M2(n)=mean(B{n});
S2(n)=nansem(B{n},1);
end
N=2;
for n=1:length(A);
    [sig(N,n),~]=ttest2(A{n},B{n});
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
M1=[];
S1=[];
M2=[];
S2=[];
for n=1:length(A);
M1(n)=mean(A{n});
S1(n)=nansem(A{n},1);
M2(n)=mean(B{n});
S2(n)=nansem(B{n},1);
end
N=3;
for n=1:length(A);
    [sig(N,n),~]=ttest2(A{n},B{n});
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
M1=[];
S1=[];
M2=[];
S2=[];
for n=1:length(A);
M1(n)=mean(A{n});
S1(n)=nansem(A{n},1);
M2(n)=mean(B{n});
S2(n)=nansem(B{n},1);
end
N=4;
for n=1:length(A);
    [sig(N,n),~]=ttest2(A{n},B{n});
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
M1=[];
S1=[];
M2=[];
S2=[];
for n=1:length(A);
M1(n)=mean(A{n});
S1(n)=nansem(A{n},1);
M2(n)=mean(B{n});
S2(n)=nansem(B{n},1);
end
N=5;
for n=1:length(A);
    [sig(N,n),~]=ttest2(A{n},B{n});
end

subplot(155)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')

xlabel('days')
ylabel('mean DfofSigMean');
title(['mean DfofSigMean']);

tightfig;
saveas(gcf,'allFAllCellsSP.fig');

%% plotting fluorscence


load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\allFAllCellsSPNew.mat');
GmF=(mFAllCellsSP);
GmDfof=(mDfofAllCellsSP);
GmDfofSig=(mDfofSigAllCellsSP);
GmDfofMean=(mDfofMeanAllCells);
GmDfofSigMean=(mDfofSigMeanAllCells);

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\allFAllCellsSPNew.mat');
BmF=(mFAllCellsSP);
BmDfof=(mDfofAllCellsSP);
BmDfofSig=(mDfofSigAllCellsSP);
BmDfofMean=(mDfofMeanAllCells);
BmDfofSigMean=(mDfofSigMeanAllCells);
sig=[];

figure,
A=GmF;
B=BmF;

M1=[];
S1=[];
M2=[];
S2=[];
for n=1:length(A);
M1(n)=mean(A{n});
S1(n)=nansem(A{n},1);
M2(n)=mean(B{n});
S2(n)=nansem(B{n},1);
end

N=1;
for n=1:length(A);
    [sig(N,n),~]=ttest2(A{n},B{n});
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
M1=[];
S1=[];
M2=[];
S2=[];
for n=1:length(A);
M1(n)=mean(A{n});
S1(n)=nansem(A{n},1);
M2(n)=mean(B{n});
S2(n)=nansem(B{n},1);
end
N=2;
for n=1:length(A);
    [sig(N,n),~]=ttest2(A{n},B{n});
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
M1=[];
S1=[];
M2=[];
S2=[];
for n=1:length(A);
M1(n)=mean(A{n});
S1(n)=nansem(A{n},1);
M2(n)=mean(B{n});
S2(n)=nansem(B{n},1);
end
N=3;
for n=1:length(A);
    [sig(N,n),~]=ttest2(A{n},B{n});
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
M1=[];
S1=[];
M2=[];
S2=[];
for n=1:length(A);
M1(n)=mean(A{n});
S1(n)=nansem(A{n},1);
M2(n)=mean(B{n});
S2(n)=nansem(B{n},1);
end
N=4;
for n=1:length(A);
    [sig(N,n),~]=ttest2(A{n},B{n});
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
M1=[];
S1=[];
M2=[];
S2=[];
for n=1:length(A);
M1(n)=mean(A{n});
S1(n)=nansem(A{n},1);
M2(n)=mean(B{n});
S2(n)=nansem(B{n},1);
end
N=5;
for n=1:length(A);
    [sig(N,n),~]=ttest2(A{n},B{n});
end

subplot(155)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')

xlabel('days')
ylabel('mean DfofSigMean');
title(['mean DfofSigMean']);

tightfig;

saveas(gcf,'allFAllCellsSPNew.fig');
%% mean FOV: florecemce of all cells addded together: 
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\allFAllCellsSPFOV.mat');
GmF=(mFAllCellsSPFOV);
GmDfof=(mDfofAllCellsSPFOV);
GmDfofSig=(mDfofSigAllCellsSPFOV);
GmDfofMean=(mDfofMeanAllCellsFOV);
GmDfofSigMean=(mDfofSigMeanAllCellsFOV);

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\allFAllCellsSPFOV.mat');
BmF=(mFAllCellsSPFOV);
BmDfof=(mDfofAllCellsSPFOV);
BmDfofSig=(mDfofSigAllCellsSPFOV);
BmDfofMean=(mDfofMeanAllCellsFOV);
BmDfofSigMean=(mDfofSigMeanAllCellsFOV);
sig=[];

figure,
A=GmF(2:end,:);
B=BmF;

M1=nanmean(A,1);
S1=nansem(A,1);
M2=nanmean(B,1);
S2=nansem(A,1);


% N=1;
% for n=1:length(A);
%     [sig(N,n),~]=ttest2(A{n},B{n});
% end

subplot(151)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')
xlabel('days')
ylabel('mean F');
title(['mean F']);

A=GmDfof(2:end,:);
B=BmDfof;
M1=nanmean(A,1);
S1=nansem(A,1);
M2=nanmean(B,1);
S2=nansem(A,1);
% N=2;
% for n=1:length(A);
%     [sig(N,n),~]=ttest2(A{n},B{n});
% end
subplot(152)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')
xlabel('days');
ylabel('mean Dfof');
title(['mean Dfof']);

A=GmDfofSig(2:end,:);
B=BmDfofSig;
M1=nanmean(A,1);
S1=nansem(A,1);
M2=nanmean(B,1);
S2=nansem(A,1);
% N=3;
% for n=1:length(A);
%     [sig(N,n),~]=ttest2(A{n},B{n});
% end
subplot(153)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')

xlabel('days')
ylabel('mean DfofSig');
title(['mean DfofSig']);

A=GmDfofMean(2:end,:);
B=BmDfofMean;
M1=nanmean(A,1);
S1=nansem(A,1);
M2=nanmean(B,1);
S2=nansem(A,1);
% N=4;
% for n=1:length(A);
%     [sig(N,n),~]=ttest2(A{n},B{n});
% end
% 
subplot(154)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')
xlabel('days');
ylabel('mean DfofMean');
title(['mean DfofMean']);

A=GmDfofSigMean(2:end,:);
B=BmDfofSigMean;
M1=nanmean(A,1);
S1=nansem(A,1);
M2=nanmean(B,1);
S2=nansem(A,1);
% N=5;
% for n=1:length(A);
%     [sig(N,n),~]=ttest2(A{n},B{n});
% end

subplot(155)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')

xlabel('days')
ylabel('mean DfofSigMean');
title(['mean DfofSigMean']);

tightfig;
saveas(gcf,'allFAllCells_meanFov.fig');

%% normlize dfof using tbe last day in env

A=[];
B=GmDfofSig;
for n=1:size(B,1);
    A(n,:)=B(n,:)/B(n,1);
end

M=nanmean(A,1);
S=nansem(A,1);

figure, errorbar([1:1:11],M,S,'g');

A=[];
B=BmDfofSig;
for n=1:size(B,1);
    A(n,:)=B(n,:)/B(n,1);
end

M=nanmean(A,1);
S=nansem(A,1);
hold on

errorbar([1:1:11],M,S,'m');
title('DfofSig FOVs allCells');

saveas(gcf,'normDfofSigFOVsAllCells.fig')
%% MEAN FOV: florecemce of all cells addded together: all uncommon cells
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\allFAllCellsSPUCFOV.mat');
GmF=(mFAllCellsSPUCFOV);
GmDfof=(mDfofAllCellsSPUCFOV);
GmDfofSig=(mDfofSigAllCellsSPUCFOV);
GmDfofMean=(mDfofMeanAllCellsUCFOV);
GmDfofSigMean=(mDfofSigMeanAllCellsUCFOV);

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\allFAllCellsSPUCFOV.mat');
BmF=(mFAllCellsSPUCFOV);
BmDfof=(mDfofAllCellsSPUCFOV);
BmDfofSig=(mDfofSigAllCellsSPUCFOV);
BmDfofMean=(mDfofMeanAllCellsUCFOV);
BmDfofSigMean=(mDfofSigMeanAllCellsUCFOV);
sig=[];

figure,
A=GmF(2:end,:);
B=BmF;

M1=nanmean(A,1);
S1=nansem(A,1);
M2=nanmean(B,1);
S2=nansem(A,1);


% N=1;
% for n=1:length(A);
%     [sig(N,n),~]=ttest2(A{n},B{n});
% end

subplot(151)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')
xlabel('days')
ylabel('mean F');
title(['mean F']);

A=GmDfof(2:end,:);
B=BmDfof;
M1=nanmean(A,1);
S1=nansem(A,1);
M2=nanmean(B,1);
S2=nansem(A,1);
% N=2;
% for n=1:length(A);
%     [sig(N,n),~]=ttest2(A{n},B{n});
% end
subplot(152)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')
xlabel('days');
ylabel('mean Dfof');
title(['mean Dfof']);

A=GmDfofSig(2:end,:);
B=BmDfofSig;
M1=nanmean(A,1);
S1=nansem(A,1);
M2=nanmean(B,1);
S2=nansem(A,1);
% N=3;
% for n=1:length(A);
%     [sig(N,n),~]=ttest2(A{n},B{n});
% end
subplot(153)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')

xlabel('days')
ylabel('mean DfofSig');
title(['mean DfofSig']);

A=GmDfofMean(2:end,:);
B=BmDfofMean;
M1=nanmean(A,1);
S1=nansem(A,1);
M2=nanmean(B,1);
S2=nansem(A,1);
% N=4;
% for n=1:length(A);
%     [sig(N,n),~]=ttest2(A{n},B{n});
% end
% 
subplot(154)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')
xlabel('days');
ylabel('mean DfofMean');
title(['mean DfofMean']);

A=GmDfofSigMean(2:end,:);
B=BmDfofSigMean;
M1=nanmean(A,1);
S1=nansem(A,1);
M2=nanmean(B,1);
S2=nansem(A,1);
% N=5;
% for n=1:length(A);
%     [sig(N,n),~]=ttest2(A{n},B{n});
% end

subplot(155)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')

xlabel('days')
ylabel('mean DfofSigMean');
title(['mean DfofSigMean']);

tightfig;
saveas(gcf,'allFAllCells_meanFovUncommon.fig');

%% normlize dfof using tbe last day in env

A=[];
B=GmDfofSig;
for n=1:size(B,1);
    A(n,:)=B(n,:)/B(n,1);
end

M=nanmean(A,1);
S=nansem(A,1);

figure, errorbar([1:1:11],M,S,'g');

A=[];
B=BmDfofSig;
for n=1:size(B,1);
    A(n,:)=B(n,:)/B(n,1);
end

M=nanmean(A,1);
S=nansem(A,1);
hold on

errorbar([1:1:11],M,S,'m');
title('sumDfofSig FOVs Uncommon');

saveas(gcf,'normDfofSigFOVsUncommon.fig')

%% sum: florecemce of all cells addded together: all uncommon cells
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\allFAllCellsSumSPFOV.mat');
GmF=(mFAllCellsSumSPFOV);
GmDfof=(mDfofAllCellsSumSPFOV);
GmDfofSig=(mDfofSigAllCellsSumSPFOV);
GmDfofMean=(mDfofMeanAllCellsSumFOV);
GmDfofSigMean=(mDfofSigMeanAllCellsSumFOV);

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\allFAllCellsSumSPFOV.mat');
BmF=(mFAllCellsSumSPFOV);
BmDfof=(mDfofAllCellsSumSPFOV);
BmDfofSig=(mDfofSigAllCellsSumSPFOV);
BmDfofMean=(mDfofMeanAllCellsSumFOV);
BmDfofSigMean=(mDfofSigMeanAllCellsSumFOV);
sig=[];

figure,
A=GmF(2:end,:);
B=BmF;

M1=nanmean(A,1);
S1=nansem(A,1);
M2=nanmean(B,1);
S2=nansem(A,1);


% N=1;
% for n=1:length(A);
%     [sig(N,n),~]=ttest2(A{n},B{n});
% end

subplot(151)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')
xlabel('days')
ylabel('mean F');
title(['mean F']);

A=GmDfof(2:end,:);
B=BmDfof;
M1=nanmean(A,1);
S1=nansem(A,1);
M2=nanmean(B,1);
S2=nansem(A,1);
% N=2;
% for n=1:length(A);
%     [sig(N,n),~]=ttest2(A{n},B{n});
% end
subplot(152)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')
xlabel('days');
ylabel('mean Dfof');
title(['mean Dfof']);

A=GmDfofSig(2:end,:);
B=BmDfofSig;
M1=nanmean(A,1);
S1=nansem(A,1);
M2=nanmean(B,1);
S2=nansem(A,1);
% N=3;
% for n=1:length(A);
%     [sig(N,n),~]=ttest2(A{n},B{n});
% end
subplot(153)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')

xlabel('days')
ylabel('mean DfofSig');
title(['mean DfofSig']);

A=GmDfofMean(2:end,:);
B=BmDfofMean;
M1=nanmean(A,1);
S1=nansem(A,1);
M2=nanmean(B,1);
S2=nansem(A,1);
% N=4;
% for n=1:length(A);
%     [sig(N,n),~]=ttest2(A{n},B{n});
% end
% 
subplot(154)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')
xlabel('days');
ylabel('mean DfofMean');
title(['mean DfofMean']);

A=GmDfofSigMean(2:end,:);
B=BmDfofSigMean;
M1=nanmean(A,1);
S1=nansem(A,1);
M2=nanmean(B,1);
S2=nansem(A,1);
% N=5;
% for n=1:length(A);
%     [sig(N,n),~]=ttest2(A{n},B{n});
% end

subplot(155)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')

xlabel('days')
ylabel('mean DfofSigMean');
title(['mean DfofSigMean']);

tightfig;
saveas(gcf,'allFAllCells_meanSumFov.fig');

%% normlize dfof using tbe last day in env

A=[];
B=GmDfofSig;
for n=1:size(B,1);
    A(n,:)=B(n,:)/B(n,1);
end

M=nanmean(A,1);
S=nansem(A,1);

figure, errorbar([1:1:11],M,S,'g');

A=[];
B=BmDfofSig;
for n=1:size(B,1);
    A(n,:)=B(n,:)/B(n,1);
end

M=nanmean(A,1);
S=nansem(A,1);
hold on

errorbar([1:1:11],M,S,'m');
title('DfofSig FOVs allCells');

saveas(gcf,'normDfofSigFOVsSumAllCells.fig')


%% mean FOV: florecemce of all cells addded together: new speed thresh method: keep sig trans that the start was not below threshold
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\allFAllCellsSPFOVNew.mat');
GmF=(mFAllCellsSPFOV);
GmDfof=(mDfofAllCellsSPFOV);
GmDfofSig=(mDfofSigAllCellsSPFOV);
GmDfofMean=(mDfofMeanAllCellsFOV);
GmDfofSigMean=(mDfofSigMeanAllCellsFOV);

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\allFAllCellsSPFOVNew.mat');
BmF=(mFAllCellsSPFOV);
BmDfof=(mDfofAllCellsSPFOV);
BmDfofSig=(mDfofSigAllCellsSPFOV);
BmDfofMean=(mDfofMeanAllCellsFOV);
BmDfofSigMean=(mDfofSigMeanAllCellsFOV);
sig=[];

figure,
A=GmF(2:end,:);
B=BmF;

M1=nanmean(A,1);
S1=nansem(A,1);
M2=nanmean(B,1);
S2=nansem(A,1);


% N=1;
% for n=1:length(A);
%     [sig(N,n),~]=ttest2(A{n},B{n});
% end

subplot(151)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')
xlabel('days')
ylabel('mean F');
title(['mean F']);

A=GmDfof(2:end,:);
B=BmDfof;
M1=nanmean(A,1);
S1=nansem(A,1);
M2=nanmean(B,1);
S2=nansem(A,1);
% N=2;
% for n=1:length(A);
%     [sig(N,n),~]=ttest2(A{n},B{n});
% end
subplot(152)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')
xlabel('days');
ylabel('mean Dfof');
title(['mean Dfof']);

A=GmDfofSig(2:end,:);
B=BmDfofSig;
M1=nanmean(A,1);
S1=nansem(A,1);
M2=nanmean(B,1);
S2=nansem(A,1);
% N=3;
% for n=1:length(A);
%     [sig(N,n),~]=ttest2(A{n},B{n});
% end
subplot(153)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')

xlabel('days')
ylabel('mean DfofSig');
title(['mean DfofSig']);

A=GmDfofMean(2:end,:);
B=BmDfofMean;
M1=nanmean(A,1);
S1=nansem(A,1);
M2=nanmean(B,1);
S2=nansem(A,1);
% N=4;
% for n=1:length(A);
%     [sig(N,n),~]=ttest2(A{n},B{n});
% end
% 
subplot(154)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')
xlabel('days');
ylabel('mean DfofMean');
title(['mean DfofMean']);

A=GmDfofSigMean(2:end,:);
B=BmDfofSigMean;
M1=nanmean(A,1);
S1=nansem(A,1);
M2=nanmean(B,1);
S2=nansem(A,1);
% N=5;
% for n=1:length(A);
%     [sig(N,n),~]=ttest2(A{n},B{n});
% end

subplot(155)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')

xlabel('days')
ylabel('mean DfofSigMean');
title(['mean DfofSigMean']);

tightfig;
saveas(gcf,'allFAllCells_meanFovNew.fig');

%% normlize dfof using tbe last day in env

A=[];
B=GmDfofSig;
for n=1:size(B,1);
    A(n,:)=B(n,:)/B(n,1);
end

M=nanmean(A,1);
S=nansem(A,1);

figure, errorbar([1:1:11],M,S,'g');

A=[];
B=BmDfofSig;
for n=1:size(B,1);
    A(n,:)=B(n,:)/B(n,1);
end

M=nanmean(A,1);
S=nansem(A,1);
hold on

errorbar([1:1:11],M,S,'m');
title('DfofSig FOVs allCells');

saveas(gcf,'normDfofSigFOVsAllCellsNew.fig')


%% sig trans features

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\allSigTrans.mat');
GmFreq=(mFreq);
GmFreqRun=(mFreqRun);
GmDisCover=(mDisCover);
GmDisCoverRun=(mDisCoverRun);
GmAmp=(mAmp);
GmSigInteg=(mSigInteg);
GmDur=(mDur);
GmSigIntegPerRun=(mSigIntegPerRun);

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\allSigTrans.mat');
BmFreq=(mFreq);
BmFreqRun=(mFreqRun);
BmDisCover=(mDisCover);
BmDisCoverRun=(mDisCoverRun);
BmAmp=(mAmp);
BmSigInteg=(mSigInteg);
BmDur=(mDur);
BmSigIntegPerRun=(mSigIntegPerRun);

sig=[];

figure,
A=GmFreq;
B=BmFreq;
M1=mean(A,1);
S1=std(A,1)/sqrt(size(A,1));

M2=mean(B,1);
S2=std(B,1)/sqrt(size(B,1));
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
M1=mean(A,1);
S1=std(A,1)/sqrt(size(A,1));

M2=mean(B,1);
S2=std(B,1)/sqrt(size(B,1));
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
M1=mean(A,1);
S1=std(A,1)/sqrt(size(A,1));

M2=mean(B,1);
S2=std(B,1)/sqrt(size(B,1));
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
M1=mean(A,1);
S1=std(A,1)/sqrt(size(A,1));

M2=mean(B,1);
S2=std(B,1)/sqrt(size(B,1));
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
M1=mean(A,1);
S1=std(A,1)/sqrt(size(A,1));

M2=mean(B,1);
S2=std(B,1)/sqrt(size(B,1));
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
M1=mean(A,1);
S1=std(A,1)/sqrt(size(A,1));

M2=mean(B,1);
S2=std(B,1)/sqrt(size(B,1));
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
M1=mean(A,1);
S1=std(A,1)/sqrt(size(A,1));

M2=mean(B,1);
S2=std(B,1)/sqrt(size(B,1));
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
M1=mean(A,1);
S1=std(A,1)/sqrt(size(A,1));

M2=mean(B,1);
S2=std(B,1)/sqrt(size(B,1));
N=7;
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
saveas(gcf,'allSigTrans.fig');


%% sig trans features all cells

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\allSigTransAllCells.mat');
GmFreq=(mFreqAllCells);
GmFreqRun=(mFreqRunAllCells);
GmDisCover=(mDisCoverAllCells);
GmDisCoverRun=(mDisCoverRunAllCells);
GmAmp=(mAmpAllCells);
GmSigInteg=(mSigIntegAllCells);
GmDur=(mDurAllCells);
GmSigIntegPerRun=(mSigIntegPerRunAllCells);

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\allSigTransAllCells.mat');
BmFreq=(mFreqAllCells);
BmFreqRun=(mFreqRunAllCells);
BmDisCover=(mDisCoverAllCells);
BmDisCoverRun=(mDisCoverRunAllCells);
BmAmp=(mAmpAllCells);
BmSigInteg=(mSigIntegAllCells);
BmDur=(mDurAllCells);
BmSigIntegPerRun=(mSigIntegPerRunAllCells);

sig=[];

figure,
A=GmFreq;
B=BmFreq;

for n=1:length(A);
M1(n)=nanmean(A{n});
S1(n)=nansem(A{n},1);
M2(n)=nanmean(B{n});
S2(n)=nansem(B{n},1);
end

N=1;
for n=1:length(A);
    [sig(N,n),~]=ttest2(A{n},B{n});
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

for n=1:length(A);
M1(n)=nanmean(A{n});
S1(n)=nansem(A{n},1);
M2(n)=nanmean(B{n});
S2(n)=nansem(B{n},1);
end

N=2;
for n=1:length(A);
    [sig(N,n),~]=ttest2(A{n},B{n});
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

for n=1:length(A);
M1(n)=nanmean(A{n});
S1(n)=nansem(A{n},1);
M2(n)=nanmean(B{n});
S2(n)=nansem(B{n},1);
end

N=3;
for n=1:length(A);
    [sig(N,n),~]=ttest2(A{n},B{n});
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

for n=1:length(A);
M1(n)=nanmean(A{n});
S1(n)=nansem(A{n},1);
M2(n)=nanmean(B{n});
S2(n)=nansem(B{n},1);
end

N=4;
for n=1:length(A);
    [sig(N,n),~]=ttest2(A{n},B{n});
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

for n=1:length(A);
M1(n)=nanmean(A{n});
S1(n)=nansem(A{n},1);
M2(n)=nanmean(B{n});
S2(n)=nansem(B{n},1);
end

N=5;
for n=1:length(A);
    [sig(N,n),~]=ttest2(A{n},B{n});
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

for n=1:length(A);
M1(n)=nanmean(A{n});
S1(n)=nansem(A{n},1);
M2(n)=nanmean(B{n});
S2(n)=nansem(B{n},1);
end

N=6;
for n=1:length(A);
    [sig(N,n),~]=ttest2(A{n},B{n});
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

for n=1:length(A);
M1(n)=nanmean(A{n});
S1(n)=nansem(A{n},1);
M2(n)=nanmean(B{n});
S2(n)=nansem(B{n},1);
end

N=7;
for n=1:length(A);
    [sig(N,n),~]=ttest2(A{n},B{n});
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

for n=1:length(A);
M1(n)=nanmean(A{n});
S1(n)=nansem(A{n},1);
M2(n)=nanmean(B{n});
S2(n)=nansem(B{n},1);
end

N=8;
for n=1:length(A);
    [sig(N,n),~]=ttest2(A{n},B{n});
end

subplot(248)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')

xlabel('days')
ylabel('mean sigIntegPerRun');
title(['mean sigIntegPerRun']);


tightfig;
saveas(gcf,'allSigTransAllCells.fig');

%% sig trans features

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\allSigTransSP.mat');
GmFreq=(mFreqSP);
GmFreqRun=(mFreqRunSP);
GmDisCover=(mDisCoverSP);
GmDisCoverRun=(mDisCoverRunSP);
GmAmp=(mAmpSP);
GmSigInteg=(mSigIntegSP);
GmDur=(mDurSP);
GmSigIntegPerRun=(mSigIntegPerRunSP);

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\allSigTransSP.mat');
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


%% sig trans features all cells

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\allSigTransAllCellsSP.mat');
GmFreq=(mFreqAllCellsSP);
GmFreqRun=(mFreqRunAllCellsSP);
GmDisCover=(mDisCoverAllCellsSP);
GmDisCoverRun=(mDisCoverRunAllCellsSP);
GmAmp=(mAmpAllCellsSP);
GmSigInteg=(mSigIntegAllCellsSP);
GmDur=(mDurAllCellsSP);
GmSigIntegPerRun=(mSigIntegPerRunAllCellsSP);

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\allSigTransAllCellsSP.mat');
BmFreq=(mFreqAllCellsSP);
BmFreqRun=(mFreqRunAllCellsSP);
BmDisCover=(mDisCoverAllCellsSP);
BmDisCoverRun=(mDisCoverRunAllCellsSP);
BmAmp=(mAmpAllCellsSP);
BmSigInteg=(mSigIntegAllCellsSP);
BmDur=(mDurAllCellsSP);
BmSigIntegPerRun=(mSigIntegPerRunAllCellsSP);

sig=[];

figure,
A=GmFreq;
B=BmFreq;

for n=1:length(A);
M1(n)=nanmean(A{n});
S1(n)=nansem(A{n},1);
M2(n)=nanmean(B{n});
S2(n)=nansem(B{n},1);
end

N=1;
for n=1:length(A);
    [sig(N,n),~]=ttest2(A{n},B{n});
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

for n=1:length(A);
M1(n)=nanmean(A{n});
S1(n)=nansem(A{n},1);
M2(n)=nanmean(B{n});
S2(n)=nansem(B{n},1);
end

N=2;
for n=1:length(A);
    [sig(N,n),~]=ttest2(A{n},B{n});
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

for n=1:length(A);
M1(n)=nanmean(A{n});
S1(n)=nansem(A{n},1);
M2(n)=nanmean(B{n});
S2(n)=nansem(B{n},1);
end

N=3;
for n=1:length(A);
    [sig(N,n),~]=ttest2(A{n},B{n});
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

for n=1:length(A);
M1(n)=nanmean(A{n});
S1(n)=nansem(A{n},1);
M2(n)=nanmean(B{n});
S2(n)=nansem(B{n},1);
end

N=4;
for n=1:length(A);
    [sig(N,n),~]=ttest2(A{n},B{n});
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

for n=1:length(A);
M1(n)=nanmean(A{n});
S1(n)=nansem(A{n},1);
M2(n)=nanmean(B{n});
S2(n)=nansem(B{n},1);
end

N=5;
for n=1:length(A);
    [sig(N,n),~]=ttest2(A{n},B{n});
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

for n=1:length(A);
M1(n)=nanmean(A{n});
S1(n)=nansem(A{n},1);
M2(n)=nanmean(B{n});
S2(n)=nansem(B{n},1);
end

N=6;
for n=1:length(A);
    [sig(N,n),~]=ttest2(A{n},B{n});
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

for n=1:length(A);
M1(n)=nanmean(A{n});
S1(n)=nansem(A{n},1);
M2(n)=nanmean(B{n});
S2(n)=nansem(B{n},1);
end

N=7;
for n=1:length(A);
    [sig(N,n),~]=ttest2(A{n},B{n});
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

for n=1:length(A);
M1(n)=nanmean(A{n});
S1(n)=nansem(A{n},1);
M2(n)=nanmean(B{n});
S2(n)=nansem(B{n},1);
end

N=8;
for n=1:length(A);
    [sig(N,n),~]=ttest2(A{n},B{n});
end

subplot(248)
errorbar([1:1:11],M1,S1,'g')
hold on
errorbar([1:1:11],M2,S2,'m')

xlabel('days')
ylabel('mean sigIntegPerRun');
title(['mean sigIntegPerRun']);


tightfig;
saveas(gcf,'allSigTransAllCellsSP.fig');


%% plotting field distributions

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\fieldDistriEachBin.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\fieldDistriEachBinUC.mat');
G=fieldDistriEachBin;
GUC=fieldDistriEachBinUC;%this number has already been normalized by number of cells
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\NCommon.mat');
% load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\NUncommon.mat');
NG=NCommon;
% NUG=NUncommon;

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\fieldDistriEachBin.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\fieldDistriEachBinUC.mat');
B=fieldDistriEachBin;
BUC=fieldDistriEachBinUC;%this number has already been normalized by number of cells
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\NCommon.mat');
% load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\NUncommon.mat');
NB=NCommon;
% NUB=NUncommon;


figure
for n=1:size(G,1);
    subplot(2,11,n);
    hold on
    plot(G(n,:)/NG,'g');
    hold on
    plot(B(n,:)/NB,'m');
    title(['day C',num2str(n-1)]);
    
      subplot(2,11,size(G,1)+n);
    hold on
    plot(GUC(n,:),'g');
    hold on
    plot(BUC(n,:),'m');
    title(['day UC',num2str(n-1)]);
end
tightfig

saveas(gcf,'fieldDistributionIndividualDays.fig');

figure
for n=1:size(G,1);
    subplot(2,11,n);
    hold on
    plot(G(n,:)/NG,'g');
     hold on
    plot(GUC(n,:),'k');
    title(['day Good',num2str(n-1)]);
    
      subplot(2,11,size(G,1)+n);
      hold on
    plot(B(n,:)/NB,'m');
 
    hold on
    plot(BUC(n,:),'k');
    title(['day Bad',num2str(n-1)]);
end
tightfig
saveas(gcf,'fieldDistributionIndividualDaysCommonUncommon.fig');

%% correlation with tempalte
load('E:\ID20201118\20201226\2ndloc1\pcaica\cueAnalysisNew_sig\tempRL.mat');
tempNRL=tempRL;
load('E:\ID20201118\20201226\2ndloc1\pcaica\cueAnalysisNew_sig\tempR.mat');
tempNR=tempR;
load('E:\ID20201118\20201226\2ndloc1\pcaica\cueAnalysisNew_sig\tempL.mat');
tempNL=tempL;


load('E:\ID20201118\20201224\loc1\pcaica\cueAnalysisNew_sig\tempRL.mat');
tempORL=tempRL;
load('E:\ID20201118\20201224\loc1\pcaica\cueAnalysisNew_sig\tempR.mat');
tempOR=tempR;
load('E:\ID20201118\20201224\loc1\pcaica\cueAnalysisNew_sig\tempL.mat');
tempOL=tempL;


tempR={};%containing: first: old env, 2-end: new env.
tempRL={};
tempL={};
tempR{1}=tempOR;
tempRL{1}=tempORL;
tempL{1}=tempOL;
for n=2:11;
    tempR{n}=tempNR;
    tempRL{n}=tempNRL;
    tempL{n}=tempNL;
end

corrTempRLG=[];
corrTempRG=[];
corrTempLG=[];
corrTempRLB=[];
corrTempRB=[];
corrTempLB=[];

for n=1:length(tempR);
    corrTempRLG(n)=corr(tempRL{n},G(n,:)'/NG);
    corrTempRG(n)=corr(tempR{n},G(n,:)'/NG);
    corrTempLG(n)=corr(tempL{n},G(n,:)'/NG);
    corrTempRLB(n)=corr(tempRL{n},B(n,:)'/NB);
    corrTempRB(n)=corr(tempR{n},B(n,:)'/NB);
    corrTempLB(n)=corr(tempL{n},B(n,:)'/NB);
end

figure,
subplot(231)
title('corr to temp RL');
plot(corrTempRLG,'g');
hold on
plot(corrTempRLB,'m');
ylim([-0.5 0.4]);
title('C L and R');

subplot(232)
title('corr to temp R');
plot(corrTempRG,'g');
hold on
plot(corrTempRB,'m');
ylim([-0.5 0.4]);
title('C R');

subplot(233)
title('corr to temp R');
plot(corrTempLG,'g');
hold on
plot(corrTempLB,'m');
ylim([-0.5 0.4]);
title('C L');

corrTempRLGUC=[];
corrTempRGUC=[];
corrTempLGUC=[];
corrTempRLBUC=[];
corrTempRBUC=[];
corrTempLBUC=[];

for n=1:length(tempR);
    corrTempRLGUC(n)=corr(tempRL{n},GUC(n,:)');
    corrTempRGUC(n)=corr(tempR{n},GUC(n,:)');
    corrTempLGUC(n)=corr(tempL{n},GUC(n,:)');
    corrTempRLBUC(n)=corr(tempRL{n},BUC(n,:)');
    corrTempRBUC(n)=corr(tempR{n},BUC(n,:)');
    corrTempLBUC(n)=corr(tempL{n},BUC(n,:)');
end

subplot(234)
title('corr to temp RL');
plot(corrTempRLGUC,'g');
hold on
plot(corrTempRLBUC,'m');
ylim([-0.5 0.4]);
title('UC L and R');

subplot(235)
title('corr to temp R');
plot(corrTempRGUC,'g');
hold on
plot(corrTempRBUC,'m');
ylim([-0.5 0.4]);
title('UC R');

subplot(236)
title('corr to temp R');
plot(corrTempLGUC,'g');
hold on
plot(corrTempLBUC,'m');
ylim([-0.5 0.4]);
title('UC L');
saveas(gcf,'fieldCorrToEnv.fig');

%% correlation with tempalte: individual cells
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\corrToEnvR.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\corrToEnvL.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\corrToEnvRL.mat');

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\corrToEnvRUC.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\corrToEnvLUC.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\corrToEnvRLUC.mat');

GR=corrToEnvR;
GL=corrToEnvL;
GRL=corrToEnvRL;
GUCR=corrToEnvRUC;
GUCL=corrToEnvLUC;
GUCRL=corrToEnvRLUC;

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\corrToEnvR.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\corrToEnvL.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\corrToEnvRL.mat');

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\corrToEnvRUC.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\corrToEnvLUC.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\corrToEnvRLUC.mat');

BR=corrToEnvR;
BL=corrToEnvL;
BRL=corrToEnvRL;
BUCR=corrToEnvRUC;
BUCL=corrToEnvLUC;
BUCRL=corrToEnvRLUC;
%
figure,
subplot(231)
A=GRL;
M=nanmean(A,2);
S=nansem(A,2);
A1=BRL;
M1=nanmean(A1,2);
S1=nansem(A1,2);
errorbar([1:1:11],M,S,'g');
hold on
errorbar([1:1:11],M1,S1,'m');
title('C L and R');
ylim([-0.07 0.07])


subplot(232)
A=GR;
M=nanmean(A,2);
S=nansem(A,2);
A1=BR;
M1=nanmean(A1,2);
S1=nansem(A1,2);
errorbar([1:1:11],M,S,'g');
hold on
errorbar([1:1:11],M1,S1,'m');

title('C R');
ylim([-0.07 0.07])

subplot(233)
A=GL;
M=nanmean(A,2);
S=nansem(A,2);
A1=BL;
M1=nanmean(A1,2);
S1=nansem(A1,2);
errorbar([1:1:11],M,S,'g');
hold on
errorbar([1:1:11],M1,S1,'m');
title('C L');
ylim([-0.07 0.07])

subplot(234)
S=[];
M=[];
A=GUCRL;
for n=1:length(A);
M(n)=nanmean(A{n},2);
S(n)=nansem(A{n},2);
end
hold on
errorbar([1:1:11],M,S,'g');
S1=[];
M1=[];
A1=BUCRL;
for n=1:length(A);
M1(n)=nanmean(A1{n},2);
S1(n)=nansem(A1{n},2);
end
hold on
errorbar([1:1:11],M1,S1,'m');
title('UC L and R');
ylim([-0.07 0.07])

subplot(235)
S=[];
M=[];
A=GUCR;
for n=1:length(A);
M(n)=nanmean(A{n},2);
S(n)=nansem(A{n},2);
end
hold on
errorbar([1:1:11],M,S,'g');
S1=[];
M1=[];
A1=BUCR;
for n=1:length(A);
M1(n)=nanmean(A1{n},2);
S1(n)=nansem(A1{n},2);
end
hold on
errorbar([1:1:11],M1,S1,'m');
title('UC R');
ylim([-0.07 0.07])

subplot(236)
S=[];
M=[];
A=GUCL;
for n=1:length(A);
M(n)=nanmean(A{n},2);
S(n)=nansem(A{n},2);
end
hold on
errorbar([1:1:11],M,S,'g');
S1=[];
M1=[];
A1=BUCL;
for n=1:length(A);
M1(n)=nanmean(A1{n},2);
S1(n)=nansem(A1{n},2);
end
hold on
errorbar([1:1:11],M1,S1,'m');
title('UC L');
ylim([-0.07 0.07])

saveas(gcf,'fieldCorrToEnv_indivCells.fig');

%% correlation with tempalte: individual FOVs
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\corrTempRCAllFOVs.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\corrTempLCAllFOVs.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\corrTempRLCAllFOVs.mat');

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\corrTempRUCAllFOVs.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\corrTempLUCAllFOVs.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\corrTempRLUCAllFOVs.mat');

GR=corrTempRCAllFOVs;
GL=corrTempLCAllFOVs;
GRL=corrTempRLCAllFOVs;
GUCR=corrTempRUCAllFOVs;
GUCL=corrTempLUCAllFOVs;
GUCRL=corrTempRLUCAllFOVs;

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\corrTempRCAllFOVs.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\corrTempLCAllFOVs.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\corrTempRLCAllFOVs.mat');

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\corrTempRUCAllFOVs.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\corrTempLUCAllFOVs.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\corrTempRLUCAllFOVs.mat');

BR=corrTempRCAllFOVs;
BL=corrTempLCAllFOVs;
BRL=corrTempRLCAllFOVs;
BUCR=corrTempRUCAllFOVs;
BUCL=corrTempLUCAllFOVs;
BUCRL=corrTempRLUCAllFOVs;
%
figure,
subplot(231)
A=GRL;
M=nanmean(A,1);
S=nansem(A,1);
A1=BRL;
M1=nanmean(A1,1);
S1=nansem(A1,1);
errorbar([1:1:11],M,S,'g');
hold on
errorbar([1:1:11],M1,S1,'m');
title('C L and R');


subplot(232)
A=GR;
M=nanmean(A,1);
S=nansem(A,1);
A1=BR;
M1=nanmean(A1,1);
S1=nansem(A1,1);
errorbar([1:1:11],M,S,'g');
hold on
errorbar([1:1:11],M1,S1,'m');

title('C R');

subplot(233)
A=GL;
M=nanmean(A,1);
S=nansem(A,1);
A1=BL;
M1=nanmean(A1,1);
S1=nansem(A1,1);
errorbar([1:1:11],M,S,'g');
hold on
errorbar([1:1:11],M1,S1,'m');
title('C L');

subplot(234)
A=GUCRL;
M=nanmean(A,1);
S=nansem(A,1);
A1=BUCRL;
M1=nanmean(A1,1);
S1=nansem(A1,1);
errorbar([1:1:11],M,S,'g');
hold on
errorbar([1:1:11],M1,S1,'m');
title('UC L and R');

subplot(235)
A=GUCR;
M=nanmean(A,1);
S=nansem(A,1);
A1=BUCR;
M1=nanmean(A1,1);
S1=nansem(A1,1);
errorbar([1:1:11],M,S,'g');
hold on
errorbar([1:1:11],M1,S1,'m');
title('UC R');

subplot(236)
A=GUCL;
M=nanmean(A,1);
S=nansem(A,1);
A1=BUCL;
M1=nanmean(A1,1);
S1=nansem(A1,1);
errorbar([1:1:11],M,S,'g');
hold on
errorbar([1:1:11],M1,S1,'m');
title('UC L');

saveas(gcf,'fieldCorrTemp_indivFOVs.fig');

%% correlation with tempalte: individual Mice
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\corrTempRCAllMice.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\corrTempLCAllMice.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\corrTempRLCAllMice.mat');

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\corrTempRUCAllMice.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\corrTempLUCAllMice.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\corrTempRLUCAllMice.mat');

GR=corrTempRCAllMice;
GL=corrTempLCAllMice;
GRL=corrTempRLCAllMice;
GUCR=corrTempRUCAllMice;
GUCL=corrTempLUCAllMice;
GUCRL=corrTempRLUCAllMice;

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\corrTempRCAllMice.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\corrTempLCAllMice.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\corrTempRLCAllMice.mat');

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\corrTempRUCAllMice.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\corrTempLUCAllMice.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\corrTempRLUCAllMice.mat');

BR=corrTempRCAllMice;
BL=corrTempLCAllMice;
BRL=corrTempRLCAllMice;
BUCR=corrTempRUCAllMice;
BUCL=corrTempLUCAllMice;
BUCRL=corrTempRLUCAllMice;
%
figure,
subplot(231)
A=GRL;
M=nanmean(A,1);
S=nansem(A,1);
A1=BRL;
M1=nanmean(A1,1);
S1=nansem(A1,1);
errorbar([1:1:11],M,S,'g');
hold on
errorbar([1:1:11],M1,S1,'m');
title('C L and R');


subplot(232)
A=GR;
M=nanmean(A,1);
S=nansem(A,1);
A1=BR;
M1=nanmean(A1,1);
S1=nansem(A1,1);
errorbar([1:1:11],M,S,'g');
hold on
errorbar([1:1:11],M1,S1,'m');

title('C R');

subplot(233)
A=GL;
M=nanmean(A,1);
S=nansem(A,1);
A1=BL;
M1=nanmean(A1,1);
S1=nansem(A1,1);
errorbar([1:1:11],M,S,'g');
hold on
errorbar([1:1:11],M1,S1,'m');
title('C L');

subplot(234)
A=GUCRL;
M=nanmean(A,1);
S=nansem(A,1);
A1=BUCRL;
M1=nanmean(A1,1);
S1=nansem(A1,1);
errorbar([1:1:11],M,S,'g');
hold on
errorbar([1:1:11],M1,S1,'m');
title('UC L and R');

subplot(235)
A=GUCR;
M=nanmean(A,1);
S=nansem(A,1);
A1=BUCR;
M1=nanmean(A1,1);
S1=nansem(A1,1);
errorbar([1:1:11],M,S,'g');
hold on
errorbar([1:1:11],M1,S1,'m');
title('UC R');

subplot(236)
A=GUCL;
M=nanmean(A,1);
S=nansem(A,1);
A1=BUCL;
M1=nanmean(A1,1);
S1=nansem(A1,1);
errorbar([1:1:11],M,S,'g');
hold on
errorbar([1:1:11],M1,S1,'m');
title('UC L');

saveas(gcf,'fieldCorrTemp_indivMice.fig');



%% plot mean number of fields
%COMMON CELLS
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\MField.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\SField.mat');
GM=MField;
GS=SField;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\MField.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\SField.mat');
BM=MField;
BS=SField;

figure,
subplot(121)
errorbar([1:1:11],GM,GS,'g');
hold on
errorbar([1:1:11],BM,BS,'m');
title('commonCellNFields');

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\MFieldUC.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\SFieldUC.mat');
GM=MFieldUC;
GS=SFieldUC;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\MFieldUC.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\SFieldUC.mat');
BM=MFieldUC;
BS=SFieldUC;
ylim([34 65])

subplot(122)
errorbar([1:1:11],GM,GS,'g');
hold on
errorbar([1:1:11],BM,BS,'m');
title('unCommonCellNFields');
ylim([34 65])

saveas(gcf,'NBinsWithFields.fig');

%% compare field chagnes
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\fieldDistriDiff.mat');
G=fieldDistriDiff;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\fieldDistriDiff.mat');
B=fieldDistriDiff;

figure,
[lineOut, fillOut] = semshade(G,0.2,'g',[1:1:200],1);
hold on
[lineOut, fillOut] = semshade(B,0.2,'m',[1:1:200],1);
hold on
plot([1 200],[0 0],'k--');
    
sig=[];
m=[];

for n=1:200;
    [sig(n),~]=ttest2(G(:,n),B(:,n));
    m(n)=mean(G(:,n))-mean(B(:,n));
end

GM=mean(G,1);
BM=mean(B,1);

for n=1:length(sig);
    if sig(n)==1;
        if m(n)>0;
            hold on
            plot(n,GM(n),'r.','MarkerSize',10);
        else
            hold on
            plot(n,BM(n),'r.','MarkerSize',10);
        end
    end
end

saveas(gcf,'fieldDistriDiffIndividualCells.fig');

%% replot the above: include good and bad learners before after learning and now the differences were calculated according to the onw group change
%LOAD CUE TEMPLATE
load('E:\ID20201118\20201226\2ndloc1\pcaica\cueAnalysisNew_sig\tempRL.mat');
tempNRL=tempRL;


load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\fieldDistriAfterLearning.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\fieldDistriBeforeLearning.mat');

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


load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\fieldDistriAfterLearning.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\fieldDistriBeforeLearning.mat');

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

%% plot "fieldDistriDiffIndividualCells" with significance above

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\fieldDistriDiff.mat');
G=fieldDistriDiff;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\fieldDistriDiff.mat');
B=fieldDistriDiff;

figure,
[lineOut, fillOut] = semshade(G,0.2,'g',[1:1:200],1);
hold on
[lineOut, fillOut] = semshade(B,0.2,'m',[1:1:200],1);
hold on
plot([1 200],[0 0],'k--');

GM=mean(G,1);
BM=mean(B,1);

for n=1:length(sigG);
    if sigG(n)==1;
   
            plot(n,GM(n),'r.','MarkerSize',10);

        end
    end

for n=1:length(sigB);
    if sigB(n)==1;
   
            plot(n,BM(n),'r.','MarkerSize',10);

        end
end

    axis tight
    
    saveas(gcf,'fieldDistriDiffIndividualCells_diffSig.fig');
    
    
    %% plot "fieldDistriDiffIndividualCells" with significance above: this fraction is 

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\fieldDistriDiffFraction.mat');
G=fieldDistriDiffFraction;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\fieldDistriDiffFraction.mat');
B=fieldDistriDiffFraction;

% %normalized by mean
% load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\fieldDistriBeforeLearning.mat');
% GB=fieldDistriBeforeLearning;
% load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\fieldDistriBeforeLearning.mat');
% BB=fieldDistriBeforeLearning;
% 
% GB=mean(GB,1);
% BB=mean(BB,1);
% 
% G=G/sum(GB);
% B=B/sum(BB);

figure,
[lineOut, fillOut] = semshade(G,0.2,'g',[1:1:200],1);
hold on
[lineOut, fillOut] = semshade(B,0.2,'m',[1:1:200],1);
hold on
plot([1 200],[0 0],'k--');

GM=mean(G,1);
BM=mean(B,1);

for n=1:length(sigG);
    if sigG(n)==1;
   
            plot(n,GM(n),'r.','MarkerSize',10);

        end
    end

for n=1:length(sigB);
    if sigB(n)==1;
   
            plot(n,BM(n),'r.','MarkerSize',10);

        end
end

    axis tight
    
    saveas(gcf,'fieldDistriDiffFractionIndividualCells_diffSig.fig');

    
    
%% correlation with tempalte: individual cells
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\corrToEnvRDfof.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\corrToEnvLDfof.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\corrToEnvRLDfof.mat');


GR=corrToEnvRDfof;
GL=corrToEnvLDfof;
GRL=corrToEnvRLDfof;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\corrToEnvRDfof.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\corrToEnvLDfof.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\corrToEnvRLDfof.mat');

BR=corrToEnvRDfof;
BL=corrToEnvLDfof;
BRL=corrToEnvRLDfof;
%
figure,
subplot(131)
A=GRL;
M=nanmean(A,2);
S=nansem(A,2);
A1=BRL;
M1=nanmean(A1,2);
S1=nansem(A1,2);
errorbar([1:1:11],M,S,'g');
hold on
errorbar([1:1:11],M1,S1,'m');
title('C L and R');
ylim([-0.07 0.07])


subplot(132)
A=GR;
M=nanmean(A,2);
S=nansem(A,2);
A1=BR;
M1=nanmean(A1,2);
S1=nansem(A1,2);
errorbar([1:1:11],M,S,'g');
hold on
errorbar([1:1:11],M1,S1,'m');

title('C R');
ylim([-0.07 0.07])

subplot(133)
A=GL;
M=nanmean(A,2);
S=nansem(A,2);
A1=BL;
M1=nanmean(A1,2);
S1=nansem(A1,2);
errorbar([1:1:11],M,S,'g');
hold on
errorbar([1:1:11],M1,S1,'m');
title('C L');
ylim([-0.07 0.07])

% subplot(234)
% S=[];
% M=[];
% A=GUCRL;
% for n=1:length(A);
% M(n)=nanmean(A{n},2);
% S(n)=nansem(A{n},2);
% end
% hold on
% errorbar([1:1:11],M,S,'g');
% S1=[];
% M1=[];
% A1=BUCRL;
% for n=1:length(A);
% M1(n)=nanmean(A1{n},2);
% S1(n)=nansem(A1{n},2);
% end
% hold on
% errorbar([1:1:11],M1,S1,'m');
% title('UC L and R');
% ylim([-0.07 0.07])
% 
% subplot(235)
% S=[];
% M=[];
% A=GUCR;
% for n=1:length(A);
% M(n)=nanmean(A{n},2);
% S(n)=nansem(A{n},2);
% end
% hold on
% errorbar([1:1:11],M,S,'g');
% S1=[];
% M1=[];
% A1=BUCR;
% for n=1:length(A);
% M1(n)=nanmean(A1{n},2);
% S1(n)=nansem(A1{n},2);
% end
% hold on
% errorbar([1:1:11],M1,S1,'m');
% title('UC R');
% ylim([-0.07 0.07])
% 
% subplot(236)
% S=[];
% M=[];
% A=GUCL;
% for n=1:length(A);
% M(n)=nanmean(A{n},2);
% S(n)=nansem(A{n},2);
% end
% hold on
% errorbar([1:1:11],M,S,'g');
% S1=[];
% M1=[];
% A1=BUCL;
% for n=1:length(A);
% M1(n)=nanmean(A1{n},2);
% S1(n)=nansem(A1{n},2);
% end
% hold on
% errorbar([1:1:11],M1,S1,'m');
% title('UC L');
% ylim([-0.07 0.07])

saveas(gcf,'dfofCorrToEnv_indivCells.fig');

%% comparing good and bad learners by randomly selecting 300 cells several times in each learner group
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\fieldDistriUC_zeroOne.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\fieldDistri_zeroOne.mat');
G=fieldDistri_zeroOne;
GUC=fieldDistriUC_zeroOne;

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\fieldDistriUC_zeroOne.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\fieldDistri_zeroOne.mat');
B=fieldDistri_zeroOne;
BUC=fieldDistriUC_zeroOne;

nG=size(G{1},1);
nB=size(B{1},1);

NCell=200;
idxG=[];%indices of good learner cells, each column is one shuffle
idxB=[];%indices of bad learner cells, each column is one shuffle

GS={};%bin distribution of G SHUFFLE, each cell is one shuffle
BS={};%bin distribution of B SHUFFLE, each cell is one shuffle

NShuffle=100;
for n=1:NShuffle;
    iG=randperm(nG);
    iB=randperm(nB);
    idxG(:,n)=iG(1:NCell);
    idxB(:,n)=iB(1:NCell);
end

for m=1:length(G);
    GS{m}=[];
    BS{m}=[];
    
    GG=[];
    BB=[];
for n=1:NShuffle;   
      
        GG(n,:)=sum(G{m}(idxG(:,n),:));
        BB(n,:)=sum(B{m}(idxB(:,n),:));        
    end
    GS{m}=GG;
    BS{m}=BB;
end
    

load('E:\ID20201118\20201226\2ndloc1\pcaica\cueAnalysisNew_sig\tempRL.mat');
tempNRL=tempRL;
load('E:\ID20201118\20201226\2ndloc1\pcaica\cueAnalysisNew_sig\tempR.mat');
tempNR=tempR;
load('E:\ID20201118\20201226\2ndloc1\pcaica\cueAnalysisNew_sig\tempL.mat');
tempNL=tempL;


load('E:\ID20201118\20201224\loc1\pcaica\cueAnalysisNew_sig\tempRL.mat');
tempORL=tempRL;
load('E:\ID20201118\20201224\loc1\pcaica\cueAnalysisNew_sig\tempR.mat');
tempOR=tempR;
load('E:\ID20201118\20201224\loc1\pcaica\cueAnalysisNew_sig\tempL.mat');
tempOL=tempL;


tempR={};%containing: first: old env, 2-end: new env.
tempRL={};
tempL={};
tempR{1}=tempOR;
tempRL{1}=tempORL;
tempL{1}=tempOL;
for n=2:11;
    tempR{n}=tempNR;
    tempRL{n}=tempNRL;
    tempL{n}=tempNL;
end

corrRLG=[];
corrRG=[];
corrLG=[];

corrRLB=[];
corrRB=[];
corrLB=[];

for n=1:NShuffle;
    for m=1:length(tempR);
        R=tempR{m};
        L=tempL{m};
        RL=tempRL{m};
        g=GS{m}(n,:)';
        b=BS{m}(n,:)';
        
       corrRLG(n,m)=corr(RL,g);
       corrRG(n,m)=corr(R,g);
       corrLG(n,m)=corr(L,g); 
       
       corrRLB(n,m)=corr(RL,b);
       corrRB(n,m)=corr(R,b);
       corrLB(n,m)=corr(L,b);
    end
end

figure
subplot(131);
A=corrRLG;
M=mean(A,1);
S=nansem(A,1);
errorbar([1:1:11],M,S,'g');

A=corrRLB;
M=mean(A,1);
S=nansem(A,1);
hold on
errorbar([1:1:11],M,S,'m');
title('RL');

subplot(132);
A=corrRG;
M=mean(A,1);
S=nansem(A,1);
errorbar([1:1:11],M,S,'g');

A=corrRB;
M=mean(A,1);
S=nansem(A,1);
hold on
errorbar([1:1:11],M,S,'m');
title('R');

subplot(133);
A=corrLG;
M=mean(A,1);
S=nansem(A,1);
errorbar([1:1:11],M,S,'g');

A=corrLB;
M=mean(A,1);
S=nansem(A,1);
hold on
errorbar([1:1:11],M,S,'m');
title('L');

sigRL=[];
sigR=[];
sigL=[];

for n=1:11;
    [sigRL(n),~]=ttest2(corrRLG(:,n),corrRLB(:,n));
    [sigR(n),~]=ttest2(corrRG(:,n),corrRB(:,n));
    [sigL(n),~]=ttest2(corrLG(:,n),corrLB(:,n));
end


saveas(gcf,'correlationToTemp_randomCells.fig');

%plot field distribution


%% in field and other ratio

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\fieldInfoUsingAllCells.mat');
G=ratioInOther;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\fieldInfoUsingAllCells.mat');
B=ratioInOther;
[pAnova,pMC] = anovaRM2W(G,B)
figure
subplot(121)
a=G;
a(a>259.452)=nan;
M=[];
S=[];
for n=1:size(a,2);
    b=a(:,n);
    b=b(~isnan(b));
      b=b(~isinf(b));
      b=b(b<259.452);
M(n)=mean(b);
S(n)=std(b)/sqrt(length(b));
end
errorbar([1:1:length(M)],M,S,'g');


a=B;
a(a>234.824)=nan;
M1=[];
S=[];
for n=1:size(a,2);
    b=a(:,n);
    b=b(~isnan(b));
      b=b(~isinf(b));
      b=b(b<234.824);
M1(n)=mean(b);
S(n)=std(b)/sqrt(length(b));
end
hold on
errorbar([1:1:length(M1)],M1,S,'m');
title('in field and other');

%NORM
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\fieldInfoUsingAllCells.mat');
G=ratioInOtherNorm;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\fieldInfoUsingAllCells.mat');
B=ratioInOtherNorm;;

subplot(122)

a=G;
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
errorbar([1:1:length(M)],M,S,'g');


a=B;
a(a>100)=nan;
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
hold on
errorbar([1:1:length(M)],M,S,'m');

title('norm');

saveas(gcf,'inFieldVSOthersRatio.fig');


%% compare in field and other F diferences
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\meanFFieldInOthers.mat');
GI=dfofFieldIn;
GO=dfofFieldOther;
GIN=dfofFieldInNorm;
GON=dfofFieldOtherNorm;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\meanFFieldInOthers.mat');
BI=dfofFieldIn;
BO=dfofFieldOther;
BIN=dfofFieldInNorm;
BON=dfofFieldOtherNorm;

figure,
subplot(221);
a=GI;
M=nanmean(a,1);
S=nansem(a,1);

errorbar([1:1:length(M)],M,S,'g');
a=BI;
M=nanmean(a,1);
S=nansem(a,1);

hold on
errorbar([1:1:length(M)],M,S,'m');

% [r,p]=corr([1:1:length(M)-1]',M(2:end)');
title(['in field dfof mean'])
ylabel('dfof mean');
xlabel('day');


subplot(222);
a=GO;
M=nanmean(a,1);
S=nansem(a,1);

errorbar([1:1:length(M)],M,S,'g');
a=BO;
M=nanmean(a,1);
S=nansem(a,1);

hold on
errorbar([1:1:length(M)],M,S,'m');

% [r,p]=corr([1:1:length(M)-1]',M(2:end)');
title(['OTHER NON field dfof mean'])
ylabel('dfof mean');
xlabel('day');

subplot(223);
a=GIN;
M=nanmean(a,1);
S=nansem(a,1);

errorbar([1:1:length(M)],M,S,'g');
a=BIN;
M=nanmean(a,1);
S=nansem(a,1);

hold on
errorbar([1:1:length(M)],M,S,'m');

% [r,p]=corr([1:1:length(M)-1]',M(2:end)');
title(['in field dfof mean NORM'])
ylabel('dfof mean');
xlabel('day');


subplot(224);
a=GON;
M=nanmean(a,1);
S=nansem(a,1);

errorbar([1:1:length(M)],M,S,'g');
a=BON;
M=nanmean(a,1);
S=nansem(a,1);

hold on
errorbar([1:1:length(M)],M,S,'m');

% [r,p]=corr([1:1:length(M)-1]',M(2:end)');
title(['OTHER NON field dfof mean NORM'])
ylabel('dfof mean');
xlabel('day');

saveas(gcf,'inOutFieldMeanF.fig')

%% other to in

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\fieldInfoUsingAllCellsMean.mat');
G=ratioOtherIn;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\fieldInfoUsingAllCellsMean.mat');
B=ratioOtherIn;

figure
subplot(121)
a=G;
% a(a>200)=nan;
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
errorbar([1:1:length(M)],M,S,'g');


a=B;
% a(a>200)=nan;
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
hold on
errorbar([1:1:length(M)],M,S,'m');
title('other AND in field');

%NORM
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\fieldInfoUsingAllCellsMean.mat');
G=ratioOtherInNorm;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\fieldInfoUsingAllCellsMean.mat');
B=ratioOtherInNorm;

subplot(122)

a=G;
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
errorbar([1:1:length(M)],M,S,'g');


a=B;
a(a>100)=nan;
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
hold on
errorbar([1:1:length(M)],M,S,'m');

title('norm');

saveas(gcf,'OthersVSinFieldRatio.fig');



%% IN AND OUT
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\fieldRatioInOut.mat');
G=ratioInOut;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\fieldRatioInOut.mat');
B=ratioInOut;

figure
subplot(121)
a=G;
% a(a>200)=nan;
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
errorbar([1:1:length(M)],M,S,'g');


a=B;
% a(a>200)=nan;
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
hold on
errorbar([1:1:length(M)],M,S,'m');
title('other AND in field');

%NORM
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\fieldRatioInOut.mat');
G=ratioInOutNorm;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\fieldRatioInOut.mat');
B=ratioInOutNorm;

subplot(122)

a=G;
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
errorbar([1:1:length(M)],M,S,'g');


a=B;
a(a>100)=nan;
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
hold on
errorbar([1:1:length(M)],M,S,'m');

title('norm');

saveas(gcf,'FieldInVSOutRatio.fig');

%% correlation to behavior

% %% correlate to behavior: run by run
% load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\corrAllFOVs.mat');
% G=corrAllFOVs;
% load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\corrAllFOVs.mat');
% B=corrAllFOVs;
% 
% %activity
% act=[G.dfofSigMeanToOthers B.dfofSigMeanToOthers];
% load('E:\learningAnalysis\behaviorNewBatch\predLickDataBad.mat');
% load('E:\learningAnalysis\behaviorNewBatch\predLickDataGood.mat');
% lick=[predLick10RunsMeanGood(1:11) predLick10RunsMeanBad(1:11)];
% load('E:\learningAnalysis\behaviorNewBatch\slowDownDataBad.mat');
% load('E:\learningAnalysis\behaviorNewBatch\slowDownDataGood.mat');
% slow=[slowDown10RunsMeanGood(1:11) slowDown10RunsMeanBad(1:11)];
% 
% figure,
% % subplot(231)
% % plot(act,lick,'r.','MarkerSize',10)
% % [~,p]=corr(act',lick');
% % title(['all lick p=',num2str(p)])
% % 
% % subplot(234)
% % plot(act,slow,'r.','MarkerSize',10)
% % [~,p]=corr(act',slow');
% % title(['all slow p=',num2str(p)])
% 
% subplot(221)
% plot(act(1:11),lick(1:11),'r.','MarkerSize',10)
% [~,p]=corr(act(1:11)',lick(1:11)');
% title(['Good lick p=',num2str(p)])
% 
% subplot(223)
% plot(act(1:11),slow(1:11),'r.','MarkerSize',10)
% [~,p]=corr(act(1:11)',slow(1:11)');
% title(['Good slow p=',num2str(p)])
% 
% subplot(222)
% plot(act(12:end),lick(12:end),'r.','MarkerSize',10)
% [~,p]=corr(act(12:end)',lick(12:end)');
% title(['Bad lick p=',num2str(p)])
% 
% subplot(224)
% plot(act(12:end),slow(12:end),'r.','MarkerSize',10)
% [~,p]=corr(act(12:end)',slow(12:end)');
% title(['Bad slow p=',num2str(p)])
% 
% saveas(gcf,'absCorrRunByRunAndBehavior.fig');

% %% correlate to behavior: matrix correlation
% load('E:\learningAnalysis\summaryManyMice_includingOldEnvGood\matrixCorrMax1.mat');
% G=matrixCorr;
% load('E:\learningAnalysis\summaryManyMice_includingOldEnvBad\matrixCorrMax1.mat');
% B=matrixCorr;
% 
% %activity
% act=[G(2:end)' B(2:end)'];
% load('E:\learningAnalysis\behaviorNewBatch\predLickDataBad.mat');
% load('E:\learningAnalysis\behaviorNewBatch\predLickDataGood.mat');
% lick=[predLick10RunsMeanGood(3:11) predLick10RunsMeanBad(3:11)];
% load('E:\learningAnalysis\behaviorNewBatch\slowDownDataBad.mat');
% load('E:\learningAnalysis\behaviorNewBatch\slowDownDataGood.mat');
% slow=[slowDown10RunsMeanGood(3:11) slowDown10RunsMeanBad(3:11)];
% 
% figure,
% % subplot(231)
% % plot(act,lick,'r.','MarkerSize',10)
% % [~,p]=corr(act',lick');
% % title(['all lick p=',num2str(p)])
% % 
% % subplot(234)
% % plot(act,slow,'r.','MarkerSize',10)
% % [~,p]=corr(act',slow');
% % title(['all slow p=',num2str(p)])
% 
% subplot(221)
% plot(act(1:9),lick(1:9),'r.','MarkerSize',10)
% [~,p]=corr(act(1:9)',lick(1:9)');
% title(['Good lick p=',num2str(p)])
% 
% subplot(223)
% plot(act(1:9),slow(1:9),'r.','MarkerSize',10)
% [~,p]=corr(act(1:9)',slow(1:9)');
% title(['Good slow p=',num2str(p)])
% 
% subplot(222)
% plot(act(10:end),lick(10:end),'r.','MarkerSize',10)
% [~,p]=corr(act(10:end)',lick(10:end)');
% title(['bad lick p=',num2str(p)])
% 
% subplot(224)
% plot(act(10:end),slow(10:end),'r.','MarkerSize',10)
% [~,p]=corr(act(10:end)',slow(10:end)');
% title(['bad slow p=',num2str(p)])
% 
% saveas(gcf,'matrixCorrAndBehavior.fig');

% %% correlate to behavior: individual correlation
% load('E:\learningAnalysis\summaryManyMice_includingOldEnvGood\corrIndivi.mat');
% G=corrIndivi;
% load('E:\learningAnalysis\summaryManyMice_includingOldEnvBad\corrIndivi.mat');
% B=corrIndivi;
% 
% %activity
% act=[mean(G(:,2:end),1) mean(B(:,2:end),1)];
% load('E:\learningAnalysis\behaviorNewBatch\predLickDataBad.mat');
% load('E:\learningAnalysis\behaviorNewBatch\predLickDataGood.mat');
% lick=[predLick10RunsMeanGood(3:11) predLick10RunsMeanBad(3:11)];
% load('E:\learningAnalysis\behaviorNewBatch\slowDownDataBad.mat');
% load('E:\learningAnalysis\behaviorNewBatch\slowDownDataGood.mat');
% slow=[slowDown10RunsMeanGood(3:11) slowDown10RunsMeanBad(3:11)];
% 
% figure,
% % subplot(231)
% % plot(act,lick,'r.','MarkerSize',10)
% % [~,p]=corr(act',lick');
% % title(['all lick p=',num2str(p)])
% % 
% % subplot(234)
% % plot(act,slow,'r.','MarkerSize',10)
% % [~,p]=corr(act',slow');
% % title(['all slow p=',num2str(p)])
% 
% subplot(221)
% plot(act(1:9),lick(1:9),'r.','MarkerSize',10)
% [~,p]=corr(act(1:9)',lick(1:9)');
% title(['Good lick p=',num2str(p)])
% 
% subplot(223)
% plot(act(1:9),slow(1:9),'r.','MarkerSize',10)
% [~,p]=corr(act(1:9)',slow(1:9)');
% title(['Good slow p=',num2str(p)])
% 
% subplot(222)
% plot(act(10:end),lick(10:end),'r.','MarkerSize',10)
% [~,p]=corr(act(10:end)',lick(10:end)');
% title(['bad lick p=',num2str(p)])
% 
% subplot(224)
% plot(act(10:end),slow(10:end),'r.','MarkerSize',10)
% [~,p]=corr(act(10:end)',slow(10:end)');
% title(['bad slow p=',num2str(p)])
% 
% saveas(gcf,'indiviCorrAndBehavior.fig');


%% plotting consistency run by run with behavior in good and bad learners
% load('E:\learningAnalysis\summaryManyMice_includingOldEnvGood\allRunAllCellLearn0.mat');
% load('E:\learningAnalysis\summaryManyMice_includingOldEnvGood\allRunAllCellLearn1.mat');
% load('E:\learningAnalysis\summaryManyMice_includingOldEnvGood\allRunAllCellLearn2.mat');
% A=allRunAllCellLearn2;
% B=allRunAllCellLearn1;
% C=allRunAllCellLearn0;
% 
% load('E:\learningAnalysis\summaryManyMice_includingOldEnvBad\allRunAllCellLearn0.mat');
% load('E:\learningAnalysis\summaryManyMice_includingOldEnvBad\allRunAllCellLearn1.mat');
% load('E:\learningAnalysis\summaryManyMice_includingOldEnvBad\allRunAllCellLearn2.mat');
% D=allRunAllCellLearn2;
% E=allRunAllCellLearn1;
% F=allRunAllCellLearn0;
% 
%  semA=std(A)/sqrt(length(A));
%      meanA=mean(A);
%      semB=std(B)/sqrt(length(B));
%      meanB=mean(B);
%       semC=std(C)/sqrt(length(C));
%  meanC=mean(C);
%   semD=std(D)/sqrt(length(D));
%      meanD=mean(D);
%         semE=std(E)/sqrt(length(E));
%      meanE=mean(E);
%        semF=std(F)/sqrt(length(F));
%      meanF=mean(F);
%      
%         [~,p1]=ttest2(A,B)
%      [~,p2]=ttest2(A,C)
%       [~,p3]=ttest2(D,E)
%       [~,p4]=ttest2(D,F)
%       
%        figure  
%     name={'GL2';'GL1';'GL0';'BL2';'BL1';'BL0'};
%      bar([1 2 3 4 5 6],[meanA meanB meanC meanD meanE meanF]);
%      set(gca,'xticklabel',name);
%      hold on
%      errorbar([1 2 3 4 5 6],[meanA meanB meanC meanD meanE meanF],[semA semB semC semD semE semF],'.');
%      title(['p1=',num2str(p1),' p2=',num2str(p2),' p3=',num2str(p3),' p4=',num2str(p4)]);
% 
%       saveas(gcf,'corr_allRunAllCellsLearn.fig');
%        saveas(gcf,'corr_allRunAllCellsLearn.eps');
       
      
       
%      %% plotting run by run consistency correlating with behavior at run by run basis across days
%      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      %this analysis turns out not good because if we randomly pick a few
%      %numbers of runs in learned runs, we also got lower consistency, which
%      %is similar to no learn runs. So the difference between learn and
%      %unlearn runs may not be due to the ral differences.
%      
%      %%%%%%%%%%%%%%%%%%%%%%%%%% 
%        %only look at new env
%        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        %good learners
%        load('E:\learningAnalysis\summaryManyMice_includingOldEnvGood\allCorrCellLearn12DaysToOwn.mat');
%        load('E:\learningAnalysis\summaryManyMice_includingOldEnvGood\allCorrCellLearn0DaysToOwn.mat');
%        
%         figure
% 
%     M=[];
%   S=[];
%   D=allCorrCellLearn12DaysToOwn;
%   
%   for n=1:length(D);
%       M(n)=mean(D{n},'omitnan');
%       S(n)=nanstd(D{n})/sqrt(length(D{n}));
%   end
%   hold on
%   errorbar([1:1:length(D)-1],M(2:end),S(2:end),'g');
%   [~,p1]=corr([1:1:length(M)-1]',M(2:end)');
% 
%       M=[];
%   S=[];
%   D=allCorrCellLearn0DaysToOwn;
%   
%   for n=1:length(D);
% %       if n==9;
% %           D{n}=D{n}(1:end~=95);
% %       end
%       M(n)=mean(D{n},'omitnan');
%       S(n)=nanstd(D{n})/sqrt(length(D{n}));
%       
%   end
%   hold on
%   errorbar([1:1:length(D)-1],M(2:end),S(2:end),'g--');
%     [~,p2]=corr([1:1:length(M)-1]',M(2:end)');
% 
%   sig=[];
%   p=[];
%   for n=2:11;
%       [~,p(n)]=ttest2(allCorrCellLearn12DaysToOwn{n},allCorrCellLearn0DaysToOwn{n},'tail','right');
%       
%       if p(n)<0.05;
%           sig(n)=1;
%           
%           hold on
%           plot(n-1,mean(allCorrCellLearn12DaysToOwn{n},'omitnan')*1.05,'r*');
%       else 
%           sig(n)=0;
%       end
%   end
%   
% 
% %   legend('learn12','learn0','Location','Southeast');
% %     title(['learn vs nolearn indiv corr p1=',num2str(p1),' p2=',num2str(p2)])
%     
%      load('E:\learningAnalysis\summaryManyMice_includingOldEnvBad\allCorrCellLearn12DaysToOwn.mat');
%        load('E:\learningAnalysis\summaryManyMice_includingOldEnvBad\allCorrCellLearn0DaysToOwn.mat');
%        
%     
%     M=[];
%   S=[];
%   D=allCorrCellLearn12DaysToOwn;
%   
%   for n=1:length(D);
%       M(n)=mean(D{n},'omitnan');
%       S(n)=nanstd(D{n})/sqrt(length(D{n}));
%   end
%   hold on
%   errorbar([1:1:length(D)-1],M(2:end),S(2:end),'m');
%   [~,p3]=corr([1:1:length(M)-1]',M(2:end)');
% 
%       M=[];
%   S=[];
%   D=allCorrCellLearn0DaysToOwn;
%   
%   for n=1:length(D);
%       M(n)=mean(D{n},'omitnan');
%       S(n)=nanstd(D{n})/sqrt(length(D{n}));
%   end
%   hold on
%   errorbar([1:1:length(D)-1],M(2:end),S(2:end),'m--');
%     [~,p4]=corr([1:1:length(M)-1]',M(2:end)');
% 
%   sig=[];
%   p=[];
%   for n=2:11;
%       [~,p(n)]=ttest2(allCorrCellLearn12DaysToOwn{n},allCorrCellLearn0DaysToOwn{n},'tail','right');
%       
%       if p(n)<0.05;
%           sig(n)=1;
%           
%           hold on
%           plot(n-1,mean(allCorrCellLearn12DaysToOwn{n},'omitnan')*1.05,'r*');
%       else 
%           sig(n)=0;
%       end
%   end
%   
%  xlim([0 11]);
% %   legend('learn12','learn0','Location','Southeast');
%     title(['p1=',num2str(p1),' p2=',num2str(p2),' p3=',num2str(p3),' p4=',num2str(p4)])
% saveas(gcf,'runByRunCorrLearnNoLearnRuns_notGood.fig');

