%plotting days tracked togethr
p=pwd;
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

% cd('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells');
% load('daysAllNorm.mat');
% 
% GN=daysAllNorm;
% daysMeanNormG=[];
% daysSEMNormG=[];
% daysSTDNormG=[];
% for n=1:length(daysAllNorm);
%     daysMeanNormG(n)=mean(daysAllNorm{n});
%     daysSEMNormG(n)=std(daysAllNorm{n})/sqrt(length(daysAllNorm{n}));
%     daysSTDNormG(n)=std(daysAllNorm{n});
% end

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

% cd('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\commonUniqueCells');
% load('daysAllNorm.mat');
% 
% BN=daysAllNorm;
% daysMeanNormB=[];
% daysSEMNormB=[];
% daysSTDNormB=[];
% for n=1:length(daysAllNorm);
%     daysMeanNormB(n)=mean(daysAllNorm{n});
%     daysSEMNormB(n)=std(daysAllNorm{n})/sqrt(length(daysAllNorm{n}));
%     daysSTDNormB(n)=std(daysAllNorm{n});
% end

cd(p);

sigCellDay=[];
for n=1:length(G);
    [sigCellDay(n),~]=ttest2(G{n},B{n});
end
% 
% sigCellDayNorm=[];
% for n=1:length(GN)-1;
%     [sigCellDayNorm(n),~]=ttest2(GN{n},BN{n});
% end



%% cell percentage
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\nCellsAllPer.mat');

G=nCellsAllPer;

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\commonUniqueCells\nCellsAllPer.mat');

B=nCellsAllPer;

cd(p);

sigCell=[];
for n=1:size(G,2);
    [sigCell(n),~]=ttest2(G(:,n),B(:,n));
end

[sigCellCommon,p]=ttest2(G(:,1),B(:,1),'Tail','left');%p=0.0461
[sigCellUNew,p]=ttest2(G(:,3),B(:,3),'Tail','right');%sig!!
%% cell percentage method 2
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\nCellsAllPerM2.mat');

GM=nCellsAllPerM2;

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\commonUniqueCells\nCellsAllPerM2.mat');
BM=nCellsAllPerM2;


sigCellM=[];
for n=1:size(GM,2);
    [sigCellM(n),~]=ttest2(GM(:,n),BM(:,n),'Tail','left');
end

%%
figure
subplot(221)
idx=[1 2 3 6:1:14];
errorbar([1:1:length(idx)],daysMeanG(idx),daysSEMG(idx),'g.')
hold on
errorbar([1:1:length(idx)],daysMeanB(idx),daysSEMB(idx),'m.')

legend('good','bad','Location','Northeast')
title('days tracked');

subplot(222)
idx=[1 4 5 6:1:14];
errorbar([1:1:length(idx)],daysMeanG(idx),daysSEMG(idx),'g.')
hold on
errorbar([1:1:length(idx)],daysMeanB(idx),daysSEMB(idx),'m.')

legend('good','bad','Location','Northeast')
title('days tracked M2');

% subplot(222)
% errorbar([1:1:length(daysSEMNormG)],daysMeanNormG,daysSEMNormG,'g.')
% hold on
% errorbar([1:1:length(daysSEMNormB)],daysMeanNormB,daysSEMNormB,'m.')
% legend('good','bad','Location','Northeast')
% title('days tracked norm');

subplot(223)
errorbar([1:1:size(G,2)],mean(G,1),std(G,1)/sqrt(size(G,1)),'g.');
hold on
errorbar([1:1:size(B,2)],mean(B,1),std(B,1)/sqrt(size(B,1)),'m.');
legend('good','bad','Location','Northeast')
title('cell percentage');

subplot(224)
errorbar([1:1:size(GM,2)],mean(GM,1),std(GM,1)/sqrt(size(GM,1)),'g.');
hold on
errorbar([1:1:size(BM,2)],mean(BM,1),std(BM,1)/sqrt(size(BM,1)),'m.');
legend('good','bad','Location','Northeast')
title('cell percentage M2');

saveas(gcf,'compare.fig');



%% correlation of cells appreared on different days
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\corrAllAllS.mat');
G=corrAllAllS;
idx=[1 2 3 6 7 8 9 10 11 12 13 14];
%cell c01, a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10

AG=G(idx);
M=[];
S=[];
for n=1:length(idx);
    B=AG{n};
    B=reshape(B,[1,size(B,1)*size(B,2)]);
    M(n)=nanmean(B);
    S(n)=nansem(B);
end

figure,

subplot(132)
errorbar([1:1:length(idx)],M,S,'g');

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\commonUniqueCells\corrAllAllS.mat');
B=corrAllAllS;
idx=[1 2 3 6 7 8 9 10 11 12 13 14];

AB=B(idx);
M=[];
S=[];
for n=1:length(idx);
    B=AB{n};
    B=reshape(B,[1,size(B,1)*size(B,2)]);
    M(n)=nanmean(B);
    S(n)=nansem(B);
end

hold on
errorbar([1:1:length(idx)],M,S,'m');
legend('good','bad','Location','Northeast')

title('correlation of all types');
sig=[];
for n=1:length(idx);
    sig(n)=ttest2(nanmean(AG{n},2),nanmean(AB{n},2));
end


load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\corrAllAllS.mat');
G=corrAllAllS;
idx=[1 4 5 6 7 8 9 10 11 12 13 14];

AG=G(idx);
M=[];
S=[];
for n=1:length(idx);
    B=AG{n};
    B=reshape(B,[1,size(B,1)*size(B,2)]);
    M(n)=nanmean(B);
    S(n)=nansem(B);
end


subplot(133)
errorbar([1:1:length(idx)],M,S,'g');

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\commonUniqueCells\corrAllAllS.mat');
B=corrAllAllS;
idx=[1 4 5 6 7 8 9 10 11 12 13 14];

AB=B(idx);
M=[];
S=[];
for n=1:length(idx);
    B=AB{n};
    B=reshape(B,[1,size(B,1)*size(B,2)]);
    M(n)=nanmean(B);
    S(n)=nansem(B);
end

hold on
errorbar([1:1:length(idx)],M,S,'m');
legend('good','bad','Location','Northeast')

title('correlation of all types M2');

sigM2=[];
for n=1:length(idx);
    [~,sigM2(n)]=ttest2(nanmean(AG{n},2),nanmean(AB{n},2));
end

subplot(131)
errorbar([1:1:11],nanmean(AG{1},1),nansem(AG{1},1),'g')
hold on
errorbar([1:1:11],nanmean(AB{1},1),nansem(AB{1},1),'m')
title('commonColdNewOnly');

saveas(gcf,'compareCorrelation.fig');

%% correlation of mean activity cells appreared on different days
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\corrMeanActAllAllS.mat');
G=corrMeanActAllAllS;
idx=[1 2 3 6 7 8 9 10 11 12 13 14];

AG=G(idx);
M=[];
S=[];
for n=1:length(idx);
    B=AG{n}(:,2:end);%remove the env switch
    B=reshape(B,[1,size(B,1)*size(B,2)]);
    M(n)=nanmean(B);
    S(n)=nansem(B);
end

figure,

subplot(132)
errorbar([1:1:length(idx)],M,S,'g.');

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\commonUniqueCells\corrMeanActAllAllS.mat');
B=corrMeanActAllAllS;
idx=[1 2 3 6 7 8 9 10 11 12 13 14];

AB=B(idx);
M=[];
S=[];
for n=1:length(idx);
    B=AB{n}(:,2:end);%remove the env switch
    B=reshape(B,[1,size(B,1)*size(B,2)]);
    M(n)=nanmean(B);
    S(n)=nansem(B);
end

hold on
errorbar([1:1:length(idx)],M,S,'m.');
legend('good','bad','Location','Northeast')

title('correlation of MeanAct all types');
sig=[];
for n=1:length(idx);
    sig(n)=ttest2(nanmean(AG{n},2),nanmean(AB{n},2));
end


load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\corrMeanActAllAllS.mat');
G=corrMeanActAllAllS;
idx=[1 4 5 6 7 8 9 10 11 12 13 14];

AG=G(idx);
M=[];
S=[];
for n=1:length(idx);
    B=AG{n}(:,2:end);%remove the env switch
    B=reshape(B,[1,size(B,1)*size(B,2)]);
    M(n)=nanmean(B);
    S(n)=nansem(B);
end


subplot(133)
errorbar([1:1:length(idx)],M,S,'g.');

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\commonUniqueCells\corrMeanActAllAllS.mat');
B=corrMeanActAllAllS;
idx=[1 4 5 6 7 8 9 10 11 12 13 14];

AB=B(idx);
M=[];
S=[];
for n=1:length(idx);
    B=AB{n}(:,2:end);%remove the env switch
    B=reshape(B,[1,size(B,1)*size(B,2)]);
    M(n)=nanmean(B);
    S(n)=nansem(B);
end

hold on
errorbar([1:1:length(idx)],M,S,'m.');
legend('good','bad','Location','Northeast')

title('correlation of MeanAct all types M2');

sigM2=[];
for n=1:length(idx);
    sigM2(n)=ttest2(nanmean(AG{n},2),nanmean(AB{n},2));
end

subplot(131)
errorbar([1:1:10],nanmean(AG{1},1),nansem(AG{1},1),'g')
hold on
errorbar([1:1:10],nanmean(AB{1},1),nansem(AB{1},1),'m')
title('commonColdNewOnly');

saveas(gcf,'compareCorrelationMeanAct.fig');

%% CUE SCore and speed

%max cue score on the day they appear
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\maxCueScoresAppearDayAll.mat');
G=maxCueScoresAppearDayAll;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\commonUniqueCells\maxCueScoresAppearDayAll.mat');
B=maxCueScoresAppearDayAll;

GM=[];
GS=[];
for n=1:length(G);
    A=G{n};
    if ~isempty(A);
    GM(n)=nanmean(A);
    GS(n)=nansem(A);
    else
        GM(n)=nan;
        GS(n)=nan;
    end
end


BM=[];
BS=[];
for n=1:length(B);
    A=B{n};
    if ~isempty(A);
    BM(n)=nanmean(A);
    BS(n)=nansem(A);
    else
        BM(n)=nan;
        BS(n)=nan;
    end
end


figure,
subplot(441)
errorbar([1:1:length(GM)],GM,GS,'g');
hold on
errorbar([1:1:length(BM)],BM,BS,'m');
title('max cue score on appear day');

%all cue R scores
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\cueScoresRAllDaysAll.mat');
G=cueScoresRAllDaysAll;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\commonUniqueCells\cueScoresRAllDaysAll.mat');
B=cueScoresRAllDaysAll;

subplot(442);
A=G{2};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M,S,'g');

hold on
A=B{2};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M,S,'m');
title('common old new cue R')

subplot(446)
A=G{4};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M,S,'g');

hold on
A=B{4};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M,S,'m');
title('unique new cue R')

subplot(4,4,10)
A=G{6};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M,S,'g');

hold on
A=B{6};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M,S,'m');
title('unique new M2 cue R')

subplot(4,4,14)
M=[];
S=[];
for n=1:length(G);
    BB=G{n}(:,2:end);%remove the env switch
    BB=reshape(BB,[1,size(BB,1)*size(BB,2)]);
    M(n)=nanmean(BB);
    S(n)=nansem(BB);
end
errorbar([1:1:length(M)],M,S,'g.');

hold on
M=[];
S=[];
for n=1:length(B);
    BB=B{n}(:,2:end);%remove the env switch
    BB=reshape(BB,[1,size(BB,1)*size(BB,2)]);
    M(n)=nanmean(BB);
    S(n)=nansem(BB);
end
errorbar([1:1:length(M)],M,S,'m.');


%all cue L scores
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\cueScoresLAllDaysAll.mat');
G=cueScoresLAllDaysAll;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\commonUniqueCells\cueScoresLAllDaysAll.mat');
B=cueScoresLAllDaysAll;

subplot(443);
A=G{2};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M,S,'g');

hold on
A=B{2};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M,S,'m');
title('common old new cue L')

subplot(447)
A=G{4};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M,S,'g');

hold on
A=B{4};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M,S,'m');
title('unique new cue L')

subplot(4,4,11)
A=G{6};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M,S,'g');

hold on
A=B{6};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M,S,'m');
title('unique new M2 cue L')

subplot(4,4,15)
M=[];
S=[];
for n=1:length(G);
    BB=G{n}(:,2:end);%remove the env switch
    BB=reshape(BB,[1,size(BB,1)*size(BB,2)]);
    M(n)=nanmean(BB);
    S(n)=nansem(BB);
end
errorbar([1:1:length(M)],M,S,'g.');

hold on
M=[];
S=[];
for n=1:length(B);
    BB=B{n}(:,2:end);%remove the env switch
    BB=reshape(BB,[1,size(BB,1)*size(BB,2)]);
    M(n)=nanmean(BB);
    S(n)=nansem(BB);
end
errorbar([1:1:length(M)],M,S,'m.');


%all speed scores
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\speedScoresAllDaysAll.mat');
G=speedScoresAllDaysAll;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\commonUniqueCells\speedScoresAllDaysAll.mat');
B=speedScoresAllDaysAll;

subplot(444);
A=G{2};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M,S,'g');

hold on
A=B{2};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M,S,'m');
title('common old new speed')

subplot(448)
A=G{4};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M,S,'g');

hold on
A=B{4};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M,S,'m');
title('unique new speed')

subplot(4,4,12)
A=G{6};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M,S,'g');

hold on
A=B{6};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:length(M)],M,S,'m');
title('unique new M2 speed')

subplot(4,4,16)
M=[];
S=[];
for n=1:length(G);
    BB=G{n}(:,2:end);%remove the env switch
    BB=reshape(BB,[1,size(BB,1)*size(BB,2)]);
    M(n)=nanmean(BB);
    S(n)=nansem(BB);
end
errorbar([1:1:length(M)],M,S,'g.');

hold on
M=[];
S=[];
for n=1:length(B);
    BB=B{n}(:,2:end);%remove the env switch
    BB=reshape(BB,[1,size(BB,1)*size(BB,2)]);
    M(n)=nanmean(BB);
    S(n)=nansem(BB);
end
errorbar([1:1:length(M)],M,S,'m.');

saveas(gcf,'compareCueScore.fig');

%% PLOTTING CELL TYPES IN EACH CATEGORY

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\cellTypesPerAll.mat');
G=cellTypesPerAll;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvbADMoreMice\commonUniqueCells\cellTypesPerAll.mat');
B=cellTypesPerAll;

figure
for n=1:length(G);
    subplot(1,5,n)
   A=G{n};
   M=nanmean(A,1);
   S=nansem(A,1);
   
   errorbar([1:1:length(M)],M,S,'g');
   
   A=B{n};
   M=nanmean(A,1);
   S=nansem(A,1);
   hold on
   errorbar([1:1:length(M)],M,S,'m');
   
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

saveas(gcf,'cellTypesPer.fig');

%% new cell types: cells in at least half sessions are grid cells
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\cellTypesNewAllFinalPercentAllFOVs.mat');
G=cellTypesNewAllFinalPercentAllFOVs;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvbADMoreMice\commonUniqueCells\cellTypesNewAllFinalPercentAllFOVs.mat');
B=cellTypesNewAllFinalPercentAllFOVs;

figure
for n=1:length(B);
     A1=G{n};
     A2=B{n};
     subplot(1,5,n);
 errorbar([1:1:15],nanmean(A1,1),nansem(A1,1),'g-');
 hold on
  errorbar([1:1:15],nanmean(A2,1),nansem(A2,1),'m-');
  
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

saveas(gcf,'cellTypesPerMoreThanHalfTrackedDays.fig');

%%
%only look at shared and unique cells in new env

%cell types on their appread day
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\cellTypesPerAll.mat');
G=cellTypesPerAll;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvbADMoreMice\commonUniqueCells\cellTypesPerAll.mat');
B=cellTypesPerAll;

%use the shared cells in new env
i=[2 4];%shared, unique old and unique new
figure
sigCellTypePer=[];
for n=1:length(G);
    subplot(5,2,2*(n-1)+1);
    A1=G{n}(:,i);
    A2=B{n}(:,i);
    
    names={'cON';'uO';'uN'};
    errorbar([1 2],mean(A1,1),nansem(A1,1),'g');
    hold on
     errorbar([1 2],mean(A2,1),nansem(A2,1),'m');
     set(gca,'xtick',[1:3],'xticklabel',names);
      [~,sigCellTypePer(n,1)]=ttest(A1(:,1),A1(:,2),'tail','right');
          [~,sigCellTypePer(n,2)]=ttest(A2(:,1),A2(:,2),'tail','right');
     if n==1;
        title(['g ',num2str(sigCellTypePer(n,1)),',',num2str(sigCellTypePer(n,2))]);
    elseif n==2;
        title(['cR ',num2str(sigCellTypePer(n,1)),',',num2str(sigCellTypePer(n,2))]);
    elseif n==3;
        title(['cL ',num2str(sigCellTypePer(n,1)),',',num2str(sigCellTypePer(n,2))]);
         elseif n==4;
        title(['sP ',num2str(sigCellTypePer(n,1)),',',num2str(sigCellTypePer(n,2))]);
    else
        title(['sN ',num2str(sigCellTypePer(n,1)),',',num2str(sigCellTypePer(n,2))]);
     end 
  
        
     
end

%stable cells in at least half sessions
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\cellTypesNewAllFinalPercentAllFOVs.mat');
G=cellTypesNewAllFinalPercentAllFOVs;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvbADMoreMice\commonUniqueCells\cellTypesNewAllFinalPercentAllFOVs.mat');
B=cellTypesNewAllFinalPercentAllFOVs;

sigCellTypeStablePer=[];

for n=1:length(G);
   subplot(5,2,2*n);
    A1=G{n}(:,i);
    A2=B{n}(:,i);
    
    names={'cON';'uO'};
    errorbar([1:2],mean(A1,1),nansem(A1,1),'g');
    hold on
     errorbar([1:2],mean(A2,1),nansem(A2,1),'m');
     set(gca,'xtick',[1:2],'xticklabel',names)
  [~,sigCellTypeStablePer(n,1)]=ttest(A1(:,1),A1(:,2),'tail','right');
          [~,sigCellTypeStablePer(n,2)]=ttest(A2(:,1),A2(:,2),'tail','right'); 
          
          if n==1;
        title(['g s ',num2str(sigCellTypeStablePer(n,1)),',',num2str(sigCellTypeStablePer(n,2))]);
    elseif n==2;
        title(['cR s ',num2str(sigCellTypeStablePer(n,1)),',',num2str(sigCellTypeStablePer(n,2))]);
    elseif n==3;
        title(['cL s ',num2str(sigCellTypeStablePer(n,1)),',',num2str(sigCellTypeStablePer(n,2))]);
         elseif n==4;
        title(['sP s ',num2str(sigCellTypeStablePer(n,1)),',',num2str(sigCellTypeStablePer(n,2))]);
    else
        title(['sN s ',num2str(sigCellTypeStablePer(n,1)),',',num2str(sigCellTypeStablePer(n,2))]);
     end 
  
      
end

saveas(gcf,'compareCOldNewAndUniqueNew_cellTypes.fig');
    
%% new cell types: cells fractions on each day
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\cellTypesNewAllFractinEachDayAllFOVs.mat');
G=cellTypesNewAllFractinEachDayAllFOVs;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvbADMoreMice\commonUniqueCells\cellTypesNewAllFractinEachDayAllFOVs.mat');
B=cellTypesNewAllFractinEachDayAllFOVs;

figure
for n=1:length(B);
     A1=G{n};
     A2=B{n};
     subplot(1,5,n);
 errorbar([1:1:14],nanmean(A1,1),nansem(A1,1),'g.');
 hold on
  errorbar([1:1:14],nanmean(A2,1),nansem(A2,1),'m.');
  
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

saveas(gcf,'cellTypesFractionsEachDay.fig');  

%% plot the days of these above cells stay as that cell type in good and bad learners
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\cellTypeDayM.mat');
GM=cellTypeDayM;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvbADMoreMice\commonUniqueCells\cellTypeDayM.mat');
BM=cellTypeDayM;

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\cellTypeDaySEM.mat');
GSEM=cellTypeDaySEM;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvbADMoreMice\commonUniqueCells\cellTypeDaySEM.mat');
BSEM=cellTypeDaySEM;

figure
for n=1:size(GM,1);
    
    subplot(1,5,n);
    errorbar([1:1:15],GM(n,:),GSEM(n,:),'g');
    hold on
      errorbar([1:1:15],BM(n,:),BSEM(n,:),'m');
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

saveas(gcf,'cellTypesDayStayingThatCell.fig');  

%% plot the days of these above cells tracked among the 11 days
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\cellTypeDayTrackM.mat');
GM=cellTypeDayTrackM;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvbADMoreMice\commonUniqueCells\cellTypeDayTrackM.mat');
BM=cellTypeDayTrackM;

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\cellTypeDayTrackSEM.mat');
GSEM=cellTypeDayTrackSEM;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvbADMoreMice\commonUniqueCells\cellTypeDayTrackSEM.mat');
BSEM=cellTypeDayTrackSEM;

figure
for n=1:size(GM,1);
    
    subplot(1,5,n);
    errorbar([1:1:15],GM(n,:),GSEM(n,:),'g.');
    hold on
      errorbar([1:1:15],BM(n,:),BSEM(n,:),'m.');
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

saveas(gcf,'cellTypesDayTracks.fig'); 

%% plotting pie graph for the fraction of grid and cue cells in each category
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\cellTypesNewAllFractinEachDayAllFOVs.mat');
G=cellTypesNewAllFractinEachDayAllFOVs;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvbADMoreMice\commonUniqueCells\cellTypesNewAllFractinEachDayAllFOVs.mat');
B=cellTypesNewAllFractinEachDayAllFOVs;

figure,
subplot(121);
pie(mean(G{1}));
colormap(summer);
title('good grid');
subplot(122);
pie(nanmean(G{2}));
colormap(summer);
title('good cueR');
saveas(gcf,'gridCueRGood.fig');

figure
subplot(121)
pie(mean(B{1}));
colormap(spring);
title('bad grid');
subplot(122);
pie(nanmean(B{2}));
colormap(spring);
title('bad cueR');
saveas(gcf,'gridCueRBad.fig');

%% correlation of cells appreared on different days
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\corrAllAllSTempRL.mat');
G=corrAllAllSTempRL;
idx=[1 2 3 4 7 8 9 10 11 12 13 14 15];

AG=G(idx);
M=[];
S=[];
for n=1:length(idx);
    B=AG{n}(:,2:end);
    B=reshape(B,[1,size(B,1)*size(B,2)]);
    M(n)=nanmean(B);
    S(n)=nansem(B);
end

figure,

subplot(321)
errorbar([1:1:length(idx)],M,S,'g.');

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\commonUniqueCells\corrAllAllSTempRL.mat');
B=corrAllAllSTempRL;
idx=[1 2 3 4 7 8 9 10 11 12 13 14 15];

AB=B(idx);
M=[];
S=[];
for n=1:length(idx);
    B=AB{n}(:,2:end);
    B=reshape(B,[1,size(B,1)*size(B,2)]);
    M(n)=nanmean(B);
    S(n)=nansem(B);
end

hold on
errorbar([1:1:length(idx)],M,S,'m.');
legend('good','bad','Location','Northeast')

title('correlationNewTempRL of all types');
sig=[];
for n=1:length(idx);
    sig(n)=ttest2(nanmean(AG{n},2),nanmean(AB{n},2));
end


load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\corrAllAllSTempRL.mat');
G=corrAllAllSTempRL;
idx=[1 2 4 5 6 7 8 9 10 11 12 13 14 15];

AG=G(idx);
M=[];
S=[];
for n=1:length(idx);
    B=AG{n}(:,2:end);
    B=reshape(B,[1,size(B,1)*size(B,2)]);
    M(n)=nanmean(B);
    S(n)=nansem(B);
end


subplot(322)
errorbar([1:1:length(idx)],M,S,'g.');

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\commonUniqueCells\corrAllAllSTempRL.mat');
B=corrAllAllSTempRL;
idx=[1 2 4 5 6 7 8 9 10 11 12 13 14 15];

AB=B(idx);
M=[];
S=[];
for n=1:length(idx);
    B=AB{n}(:,2:end);
    B=reshape(B,[1,size(B,1)*size(B,2)]);
    M(n)=nanmean(B);
    S(n)=nansem(B);
end

hold on
errorbar([1:1:length(idx)],M,S,'m.');
legend('good','bad','Location','Northeast')

title('correlationNewTempRL of all types M2');

sigM2=[];
for n=1:length(idx);
    sigM2(n)=ttest2(nanmean(AG{n},2),nanmean(AB{n},2));
end

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\corrAllAllSTempR.mat');
G=corrAllAllSTempR;
idx=[1 2 3 4 7 8 9 10 11 12 13 14 15];

AG=G(idx);
M=[];
S=[];
for n=1:length(idx);
    B=AG{n}(:,2:end);
    B=reshape(B,[1,size(B,1)*size(B,2)]);
    M(n)=nanmean(B);
    S(n)=nansem(B);
end



subplot(323)
errorbar([1:1:length(idx)],M,S,'g.');

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\commonUniqueCells\corrAllAllSTempR.mat');
B=corrAllAllSTempR;
idx=[1 2 3 4 7 8 9 10 11 12 13 14 15];

AB=B(idx);
M=[];
S=[];
for n=1:length(idx);
    B=AB{n}(:,2:end);
    B=reshape(B,[1,size(B,1)*size(B,2)]);
    M(n)=nanmean(B);
    S(n)=nansem(B);
end

hold on
errorbar([1:1:length(idx)],M,S,'m.');
legend('good','bad','Location','Northeast')

title('correlationNewTempR of all types');
sig=[];
for n=1:length(idx);
    sig(n)=ttest2(nanmean(AG{n},2),nanmean(AB{n},2));
end


load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\corrAllAllSTempR.mat');
G=corrAllAllSTempR;
idx=[1 2 4 5 6 7 8 9 10 11 12 13 14 15];

AG=G(idx);
M=[];
S=[];
for n=1:length(idx);
    B=AG{n}(:,2:end);
    B=reshape(B,[1,size(B,1)*size(B,2)]);
    M(n)=nanmean(B);
    S(n)=nansem(B);
end


subplot(324)
errorbar([1:1:length(idx)],M,S,'g.');

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\commonUniqueCells\corrAllAllSTempR.mat');
B=corrAllAllSTempR;
idx=[1 2 4 5 6 7 8 9 10 11 12 13 14 15];

AB=B(idx);
M=[];
S=[];
for n=1:length(idx);
    B=AB{n}(:,2:end);
    B=reshape(B,[1,size(B,1)*size(B,2)]);
    M(n)=nanmean(B);
    S(n)=nansem(B);
end

hold on
errorbar([1:1:length(idx)],M,S,'m.');
legend('good','bad','Location','Northeast')

title('correlationNewTempR of all types M2');

sigM2=[];
for n=1:length(idx);
    sigM2(n)=ttest2(nanmean(AG{n},2),nanmean(AB{n},2));
end

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\corrAllAllSTempL.mat');
G=corrAllAllSTempL;
idx=[1 2 3 4 7 8 9 10 11 12 13 14 15];

AG=G(idx);
M=[];
S=[];
for n=1:length(idx);
    B=AG{n}(:,2:end);
    B=reshape(B,[1,size(B,1)*size(B,2)]);
    M(n)=nanmean(B);
    S(n)=nansem(B);
end



subplot(325)
errorbar([1:1:length(idx)],M,S,'g.');

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\commonUniqueCells\corrAllAllSTempL.mat');
B=corrAllAllSTempL;
idx=[1 2 3 4 7 8 9 10 11 12 13 14 15];

AB=B(idx);
M=[];
S=[];
for n=1:length(idx);
    B=AB{n}(:,2:end);
    B=reshape(B,[1,size(B,1)*size(B,2)]);
    M(n)=nanmean(B);
    S(n)=nansem(B);
end

hold on
errorbar([1:1:length(idx)],M,S,'m.');
legend('good','bad','Location','Northeast')

title('correlationNewTempL of all types');
sig=[];
for n=1:length(idx);
    sig(n)=ttest2(nanmean(AG{n},2),nanmean(AB{n},2));
end


load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\corrAllAllSTempL.mat');
G=corrAllAllSTempL;
idx=[1 2 4 5 6 7 8 9 10 11 12 13 14 15];

AG=G(idx);
M=[];
S=[];
for n=1:length(idx);
    B=AG{n}(:,2:end);
    B=reshape(B,[1,size(B,1)*size(B,2)]);
    M(n)=nanmean(B);
    S(n)=nansem(B);
end


subplot(326)
errorbar([1:1:length(idx)],M,S,'g.');

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\commonUniqueCells\corrAllAllSTempL.mat');
B=corrAllAllSTempL;
idx=[1 2 4 5 6 7 8 9 10 11 12 13 14 15];

AB=B(idx);
M=[];
S=[];
for n=1:length(idx);
    B=AB{n}(:,2:end);
    B=reshape(B,[1,size(B,1)*size(B,2)]);
    M(n)=nanmean(B);
    S(n)=nansem(B);
end

hold on
errorbar([1:1:length(idx)],M,S,'m.');
legend('good','bad','Location','Northeast')

title('correlationNewTempL of all types M2');

sigM2=[];
for n=1:length(idx);
    sigM2(n)=ttest2(nanmean(AG{n},2),nanmean(AB{n},2));
end


saveas(gcf,'compareCorrelationToNewTemp.fig');

%% 
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\corrRLNewDay1.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\corrRNewDay1.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\corrLNewDay1.mat');
G=[corrRLNewDay1 corrRNewDay1 corrLNewDay1];
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\commonUniqueCells\corrRLNewDay1.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\commonUniqueCells\corrRNewDay1.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\commonUniqueCells\corrLNewDay1.mat');
B=[corrRLNewDay1 corrRNewDay1 corrLNewDay1];

figure,
bar([1:1:6],mean(G,1));
hold on
errorbar([1:1:6],mean(G,1),nansem(G,1),'g.');
hold on
bar([7:1:12],mean(B,1));
hold on
errorbar([7:1:12],mean(B,1),nansem(B,1),'m.');
saveas(gcf,'compareCorrelationToNewTempDay1Only.fig');

%% get an idea about the percetage of cells that only showed up once in each category
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\perOnce.mat');
G=perOnce;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\commonUniqueCells\perOnce.mat');
B=perOnce;

figure

M=nanmean(G,1);
S=nansem(G,1);
errorbar([1:1:length(M)],M,S,'g.');

M=nanmean(B,1);
S=nansem(B,1);
hold on
errorbar([1:1:length(M)],M,S,'m.');
title('percent of cells showing up once');
ylabel('percent');
saveas(gcf,'perOnce.fig');

%% different cell types: correlation with template
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\corrAllAllTempRLAllFolders.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\idxOnlyOld.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\idxOnly10DayNew.mat');

G=corrAllAllTempRLAllFolders;
GO=idxOnlyOld;
GN=idxOnly10DayNew;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\commonUniqueCells\corrAllAllTempRLAllFolders.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\commonUniqueCells\idxOnlyOld.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\commonUniqueCells\idxOnly10DayNew.mat');
B=corrAllAllTempRLAllFolders;
BO=idxOnlyOld;
BN=idxOnly10DayNew;

figure
subplot(131);
A=G{1};
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:11],M,S,'g')
A=B{1};
M=nanmean(A,1);
S=nansem(A,1);
hold on
errorbar([1:1:11],M,S,'m')

subplot(132);
A1=G{1}(:,2);
A2=G{4}(:,2);
A3=A2(find(GN));
A4=B{1}(:,2);
A5=B{4}(:,2);
A6=A5(find(BN));
M=[];
S=[];
M(1)=nanmean(A1);
M(2)=nanmean(A2);
M(3)=nanmean(A3);
M(4)=nanmean(A4);
M(5)=nanmean(A5);
M(6)=nanmean(A6);
S(1)=nansem(A1,1);
S(2)=nansem(A2,1);
S(3)=nansem(A3,1);
S(4)=nansem(A4,1);
S(5)=nansem(A5,1);
S(6)=nansem(A6,1);
bar([1:1:6],M);
hold on
errorbar([1:1:6],M,S,'.');
title('new env day1');

subplot(133);
A1=G{1}(:,1);
A2=G{3}(:,1);
A3=A2(find(GO));
A4=B{1}(:,1);
A5=B{3}(:,1);
A6=A5(find(BO));
M=[];
S=[];
M(1)=nanmean(A1);
M(2)=nanmean(A2);
M(3)=nanmean(A3);
M(4)=nanmean(A4);
M(5)=nanmean(A5);
M(6)=nanmean(A6);
S(1)=nansem(A1,1);
S(2)=nansem(A2,1);
S(3)=nansem(A3,1);
S(4)=nansem(A4,1);
S(5)=nansem(A5,1);
S(6)=nansem(A6,1);
bar([1:1:6],M);
hold on
errorbar([1:1:6],M,S,'.');
title('new env day1');

%% compare number of bins with fields
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\nBinFieldsAllAll.mat');
G=nBinFieldsAllAll;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\commonUniqueCells\nBinFieldsAllAll.mat');
B=nBinFieldsAllAll;

figure

subplot(121)
idx=[1 2 3 4 7 8 9 10 11 12 13 14 15];
A=G(idx);
M=[];
S=[];
for n=1:length(idx);
    BB=A{n}; %new env temp only
    BB=reshape(BB,[1,size(BB,1)*size(BB,2)]);
    M(n)=nanmean(BB);
    S(n)=nansem(BB);
end

hold on
errorbar([1:1:length(idx)],M,S,'g.');

A=B(idx);
M=[];
S=[];
for n=1:length(idx);
    BB=A{n}; %new env temp only
    BB=reshape(BB,[1,size(BB,1)*size(BB,2)]);
    M(n)=nanmean(BB);
    S(n)=nansem(BB);
end

hold on
errorbar([1:1:length(idx)],M,S,'m.');

title('nBinFieldsAllAll mean days');


subplot(122)
idx=[1 2 5 6 7 8 9 10 11 12 13 14 15];
A=G(idx);
M=[];
S=[];
for n=1:length(idx);
    BB=A{n}; %new env temp only
    BB=reshape(BB,[1,size(BB,1)*size(BB,2)]);
    M(n)=nanmean(BB);
    S(n)=nansem(BB);
end

hold on
errorbar([1:1:length(idx)],M,S,'g.');

A=B(idx);
M=[];
S=[];
for n=1:length(idx);
    BB=A{n}; %new env temp only
    BB=reshape(BB,[1,size(BB,1)*size(BB,2)]);
    M(n)=nanmean(BB);
    S(n)=nansem(BB);
end

hold on
errorbar([1:1:length(idx)],M,S,'m.');

saveas(gcf,'nBinFields.fig');
title('nBinFieldsAllAll mean days M2');

%% percentage of grid cue in good bad

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells\pCellTypes.mat');
G=pCellTypes;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\commonUniqueCells\pCellTypes.mat');
B=pCellTypes;

idxG=[1 4 7];
idxB=[2 5 8];
figure,
M=[];
M([1 3 5])=mean(G(1:3,:),2);
M([2 4 6])=mean(B(1:3,:),2);
S=[];
S([1 3 5])=nansem(G(1:3,:),2);
S([2 4 6])=nansem(B(1:3,:),2);
bar([1:1:6],M);
hold on
errorbar([1:1:6],M,S,'k.')
sig=[];
[~,sig(1)]=ttest2(G(1,:),B(1,:));
[~,sig(2)]=ttest2(G(2,:),B(2,:));
[~,sig(3)]=ttest2(G(3,:),B(3,:));

[~,sig(4)]=ttest2(G(2,:),G(3,:));
[~,sig(5)]=ttest2(B(2,:),B(3,:));

saveas(gcf,'pCellTypes.fig');