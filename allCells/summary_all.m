%put all corr together
foldersAll={};
% foldersAll{1}='E:\ID20201118\findCells\allCellsWith1210\analyzeCells';
% foldersAll{2}='E:\ID20210206\findCells\loc1\subsetOfFOVsWithAnotherOldNo11\analyzeCells';
% foldersAll{3}='E:\ID20210206\findCells\loc2\subsetFOVsWithAnotherOldNo11\analyzeCells';
% foldersAll{4}='E:\ID20210206\findCells\loc3\subsetFOVsWithOldNo11\analyzeCells';
foldersAll{1}='E:\ID20210207\findCells\loc1\subsetFOVOldNO11\analyzeCells';
foldersAll{2}='E:\ID20210207\findCells\loc2\subsetFOVWithOldNo11\analyzeCells';

% foldersAll{6}='E:\ID20210208\findCells\loc1\subsetFOVsWithOldNo11\analyzeCells';
% foldersAll{7}='E:\ID20210209\findCells\loc1\subsetFOVsWithOldNo11\analyzeCells';
% foldersAll{8}='E:\ID20210209\findCells\loc3RedidPCAon0311\subsetFOVsWithOldNo11\analyzeCells';
% foldersAll{9}='E:\ID20210413\findCells\loc1\subsetOld1_10New\analyzeCells';
% foldersAll{10}='E:\ID20210413\findCells\loc2\subsetOldAnd1_10\analyzeCells';
foldersAll{3}='E:\ID20210519_1\findCells\loc1\subsetFOVOldNew1_10\analyzeCells';
foldersAll{4}='E:\ID20210519_1\findCells\loc2\subsetFOVOld1-10\analyzeCells';
foldersAll{5}='E:\ID20210519_2\findCells\loc1\subsetFOVsOld1_10\analyzeCells';
foldersAll{6}='E:\ID20210519_2\findCells\loc2\subsetFOVOld1-10\analyzeCells';
foldersAll{7}='E:\ID20210519_2\findCells\loc3\subsetFOVOld1_10\analyzeCells';
foldersAll{8}='E:\ID20210811A\findCells\loc1\oldNew1_10\analyzeCells';
foldersAll{9}='E:\ID20210811A\findCells\loc3\oldNew1_10\analyzeCells';

save('foldersAll.mat','foldersAll');

%%

corrAllFOVs=[];
corrAllFOVs.dfof{1}=[];
corrAllFOVs.dfof{2}=[];
corrAllFOVs.dfof{3}=[];
corrAllFOVs.dfofSig{1}=[];
corrAllFOVs.dfofSig{2}=[];
corrAllFOVs.dfofSig{3}=[];

for n=1:length(foldersAll);
    clear corrAll
  load([foldersAll{n} '\' 'corrAll.mat']);
  disp(n)
  
      corrAllFOVs.dfof{1}(end+1:end+size(corrAll.dfof{1},1),:)=corrAll.dfof{1};
  corrAllFOVs.dfof{2}(end+1:end+size(corrAll.dfof{2},1),:)=corrAll.dfof{2};
  corrAllFOVs.dfof{3}(end+1:end+size(corrAll.dfof{3},1),:)=corrAll.dfof{3};
  corrAllFOVs.dfofSig{1}(end+1:end+size(corrAll.dfofSig{1},1),:)=corrAll.dfofSig{1};
  corrAllFOVs.dfofSig{2}(end+1:end+size(corrAll.dfofSig{2},1),:)=corrAll.dfofSig{2};
  corrAllFOVs.dfofSig{3}(end+1:end+size(corrAll.dfofSig{3},1),:)=corrAll.dfofSig{3}; 
  end




corrAllFOVs.dfofMeanToMean=nanmean(corrAllFOVs.dfof{1},1);
corrAllFOVs.dfofMeanToNext=nanmean(corrAllFOVs.dfof{2},1);
corrAllFOVs.dfofMeanToOthers=nanmean(corrAllFOVs.dfof{3},1);

corrAllFOVs.dfofSigMeanToMean=nanmean(corrAllFOVs.dfofSig{1},1);
corrAllFOVs.dfofSigMeanToNext=nanmean(corrAllFOVs.dfofSig{2},1);
corrAllFOVs.dfofSigMeanToOthers=nanmean(corrAllFOVs.dfofSig{3},1);

save('corrAllFOVs.mat','corrAllFOVs');
        
figure
subplot(141)
title('dfof');
plot(corrAllFOVs.dfofMeanToMean);
hold on
plot(corrAllFOVs.dfofMeanToNext);
hold on
plot(corrAllFOVs.dfofMeanToOthers);
legend('mean','next','others','Location','SouthEast');

subplot(142)
title('dfof sig');
plot(corrAllFOVs.dfofSigMeanToMean);
hold on
plot(corrAllFOVs.dfofSigMeanToNext);
hold on
plot(corrAllFOVs.dfofSigMeanToOthers);
legend('mean','next','others','Location','SouthEast');



%only use dfofSig Mean To Mean and do shuffle
A=corrAllFOVs.dfofSig{3};
MTempRL=corrAllFOVs.dfofSigMeanToOthers;
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

title(['sig meanToMean New only,p=',num2str(p)]);


sigToNDay1=[];%significance compared to day1 in new env
day1=2;%day1 is the first column of "A".
% A=corrAllFOVs.dfofSig{1};
A=corrAllFOVs.dfofSig{3};%use one run to the others

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


saveas(gcf,'corrMeanAllCell.fig');



%%

load('corrAllFOVs.mat');
corrAllFOVs.dfofNorm{1}=[];
corrAllFOVs.dfofNorm{2}=[];
corrAllFOVs.dfofNorm{3}=[];
corrAllFOVs.dfofSigNorm{1}=[];
corrAllFOVs.dfofSigNorm{2}=[];
corrAllFOVs.dfofSigNorm{3}=[];

for n=1:length(foldersAll);
    clear corrAll
  load([foldersAll{n} '\' 'corrAll.mat']);
  disp(n)
  
      corrAllFOVs.dfofNorm{1}(end+1:end+size(corrAll.dfofNorm{1},1),:)=corrAll.dfofNorm{1};
  corrAllFOVs.dfofNorm{2}(end+1:end+size(corrAll.dfofNorm{2},1),:)=corrAll.dfofNorm{2};
  corrAllFOVs.dfofNorm{3}(end+1:end+size(corrAll.dfofNorm{3},1),:)=corrAll.dfofNorm{3};
  corrAllFOVs.dfofSigNorm{1}(end+1:end+size(corrAll.dfofSigNorm{1},1),:)=corrAll.dfofSigNorm{1};
  corrAllFOVs.dfofSigNorm{2}(end+1:end+size(corrAll.dfofSigNorm{2},1),:)=corrAll.dfofSigNorm{2};
  corrAllFOVs.dfofSigNorm{3}(end+1:end+size(corrAll.dfofSigNorm{3},1),:)=corrAll.dfofSigNorm{3}; 
  end

corrAllFOVs.dfofNormMeanToMean=mean(corrAllFOVs.dfofNorm{1},1);
corrAllFOVs.dfofNormMeanToNext=mean(corrAllFOVs.dfofNorm{2},1);
corrAllFOVs.dfofNormMeanToOthers=mean(corrAllFOVs.dfofNorm{3},1);

corrAllFOVs.dfofSigNormMeanToMean=mean(corrAllFOVs.dfofSigNorm{1},1);
corrAllFOVs.dfofSigNormMeanToNext=mean(corrAllFOVs.dfofSigNorm{2},1);
corrAllFOVs.dfofSigNormMeanToOthers=mean(corrAllFOVs.dfofSigNorm{3},1);

save('corrAllFOVs.mat','corrAllFOVs');

%%

a=corrAllFOVs.dfofNorm{1};
b=[];
for n=1:size(a,2);
c=a(:,n);
   d=c(~isnan(c));
   b(1,n)=mean(d);
end
corrAllFOVs.dfofNormMeanToMean=b;
a=corrAllFOVs.dfofNorm{2};
b=[];
for n=1:size(a,2);
c=a(:,n);
   d=c(~isnan(c));
   b(1,n)=mean(d);
end
corrAllFOVs.dfofNormMeanToNext=b;
a=corrAllFOVs.dfofNorm{3};
b=[];
for n=1:size(a,2);
c=a(:,n);
   d=c(~isnan(c));
   b(1,n)=mean(d);
end
corrAllFOVs.dfofNormMeanToOthers=b;


a=corrAllFOVs.dfofSigNorm{1};
b=[];
for n=1:size(a,2);
c=a(:,n);
   d=c(~isnan(c));
   b(1,n)=mean(d);
end
corrAllFOVs.dfofSigNormMeanToMean=b;
a=corrAllFOVs.dfofSigNorm{2};
b=[];
for n=1:size(a,2);
c=a(:,n);
   d=c(~isnan(c));
   b(1,n)=mean(d);
end
corrAllFOVs.dfofSigNormMeanToNext=b;
a=corrAllFOVs.dfofSigNorm{3};
b=[];
for n=1:size(a,2);
c=a(:,n);
   d=c(~isnan(c));
   b(1,n)=mean(d);
end
corrAllFOVs.dfofSigNormMeanToOthers=b;


save('corrAllFOVs.mat','corrAllFOVs');

%%
load('corrAllFOVs.mat');
figure
subplot(241)
plot(corrAllFOVs.dfofNormMeanToMean);
hold on
plot(corrAllFOVs.dfofNormMeanToNext);
hold on
plot(corrAllFOVs.dfofNormMeanToOthers);
legend('mean','next','others','Location','SouthEast');
title('dfof');

subplot(242)
plot(corrAllFOVs.dfofSigNormMeanToMean);
hold on
plot(corrAllFOVs.dfofSigNormMeanToNext);
hold on
plot(corrAllFOVs.dfofSigNormMeanToOthers);
legend('mean','next','others','Location','SouthEast');
title('dfof sig');

%only use dfof Mean To Mean and do shuffle
A=corrAllFOVs.dfofNorm{3};
MTempRL=corrAllFOVs.dfofNormMeanToOthers;
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
A=corrAllFOVs.dfofNorm{3};

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
A=corrAllFOVs.dfofSigNorm{1};
MTempRL=corrAllFOVs.dfofSigNormMeanToMean;
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
A=corrAllFOVs.dfofSigNorm{3};

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
A=corrAllFOVs.dfofNorm{3};
MTempRL=corrAllFOVs.dfofNormMeanToOthers;
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
A=corrAllFOVs.dfofNorm{3};

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

%% plotting tracked and untracked together
load('foldersAll.mat');
corrAllFOVUncommon=[];
corrAllFOVUncommon.dfof{1}=[];
corrAllFOVUncommon.dfof{2}=[];
corrAllFOVUncommon.dfof{3}=[];
corrAllFOVUncommon.dfofSig{1}=[];
corrAllFOVUncommon.dfofSig{2}=[];
corrAllFOVUncommon.dfofSig{3}=[];

for m=1:11;
              corrAllFOVUncommon.dfof{1}{m}=[];
            corrAllFOVUncommon.dfof{2}{m}=[];
              corrAllFOVUncommon.dfof{3}{m}=[];
              corrAllFOVUncommon.dfofSig{1}{m}=[];
            corrAllFOVUncommon.dfofSig{2}{m}=[];
              corrAllFOVUncommon.dfofSig{3}{m}=[]; 
end
p=pwd;
for n=1:length(foldersAll);
    disp(n)
    cd(foldersAll{n})
    load('commonUncommon\uncommonCorr.mat');
    for m=1:length(uncommonCorr.dfof{1});

    corrAllFOVUncommon.dfof{1}{m}(end+1:end+length(uncommonCorr.dfof{1}{m}),1)=uncommonCorr.dfof{1}{m};
    corrAllFOVUncommon.dfof{2}{m}(end+1:end+length(uncommonCorr.dfof{2}{m}),1)=uncommonCorr.dfof{2}{m};
    corrAllFOVUncommon.dfof{3}{m}(end+1:end+length(uncommonCorr.dfof{3}{m}),1)=uncommonCorr.dfof{3}{m};
    corrAllFOVUncommon.dfofSig{1}{m}(end+1:end+length(uncommonCorr.dfofSig{1}{m}),1)=uncommonCorr.dfofSig{1}{m};
    corrAllFOVUncommon.dfofSig{2}{m}(end+1:end+length(uncommonCorr.dfofSig{2}{m}),1)=uncommonCorr.dfofSig{2}{m};
    corrAllFOVUncommon.dfofSig{3}{m}(end+1:end+length(uncommonCorr.dfofSig{3}{m}),1)=uncommonCorr.dfofSig{3}{m};
    end
end
cd(p);

for n=1:length(corrAllFOVUncommon.dfof{1});
    for m=1:length(corrAllFOVUncommon.dfof{1});
          a=corrAllFOVUncommon.dfof{1}{m};
          a=a(~isnan(a));
        corrAllFOVUncommon.dfofMeanToMean(1,m)=mean(a);
        corrAllFOVUncommon.dfofMeanToMean_sem(1,m)=std(a)/sqrt(length(a));
     a=corrAllFOVUncommon.dfof{2}{m};
          a=a(~isnan(a));
        corrAllFOVUncommon.dfofMeanToNext(1,m)=mean(a);
        corrAllFOVUncommon.dfofMeanToNext_sem(1,m)=std(a)/sqrt(length(a));
         a=corrAllFOVUncommon.dfof{3}{m};
          a=a(~isnan(a));
        corrAllFOVUncommon.dfofMeanToOthers(1,m)=mean(a);
         corrAllFOVUncommon.dfofMeanToOthers_sem(1,m)=std(a)/sqrt(length(a));
        
         a=corrAllFOVUncommon.dfofSig{1}{m};
          a=a(~isnan(a));
        corrAllFOVUncommon.dfofSigMeanToMean(1,m)=mean(a);
        corrAllFOVUncommon.dfofSigMeanToMean_sem(1,m)=std(a)/sqrt(length(a));
     a=corrAllFOVUncommon.dfofSig{2}{m};
          a=a(~isnan(a));
        corrAllFOVUncommon.dfofSigMeanToNext(1,m)=mean(a);
        corrAllFOVUncommon.dfofSigMeanToNext_sem(1,m)=std(a)/sqrt(length(a));
         a=corrAllFOVUncommon.dfofSig{3}{m};
          a=a(~isnan(a));
        corrAllFOVUncommon.dfofSigMeanToOthers(1,m)=mean(a);
        corrAllFOVUncommon.dfofSigMeanToOthers_sem(1,m)=std(a)/sqrt(length(a));
    end
end

save('corrAllFOVUncommon.mat','corrAllFOVUncommon')

%% plotting common uncommon together

load('corrAllFOVs.mat');
load('corrAllFOVUncommon.mat');

%use sig mean to mean

figure
%only use dfofSig Mean To Mean and do shuffle
A=corrAllFOVs.dfofSig{3};
MTempRL=corrAllFOVs.dfofSigMeanToOthers;
S=std(A,1)/sqrt(size(A,1));

% %shuffle
% nShuffle=100;
% MShuffle=[];
% for n=1:nShuffle;  
%     shuffleDfofSig=[];
%     for m=1:size(A,1);
%         shuffleDfofSig(m,:)=A(m,randperm(size(A,2)));
%     end
%     MShuffle(n,:)=mean(shuffleDfofSig,1);
% end
% 
% isSig=[];
% for n=1:size(MShuffle,2);
%     if MTempRL(n)>=prctile(MShuffle(:,n),95);
%         isSig(n)=1;
%     elseif MTempRL(n)<=prctile(MShuffle(:,n),5);
%         isSig(n)=1;
%     else
%         isSig(n)=0;
%     end
% end

errorbar([1:1:length(MTempRL)],MTempRL,S);
hold on
errorbar([1:1:length(MTempRL)],corrAllFOVUncommon.dfofSigMeanToOthers,corrAllFOVUncommon.dfofSigMeanToOthers_sem)
title('corr run-By-Run NE only');
xlabel('Days');
ylabel('Mean correlation');
% for n=1:length(MTempRL);
%     if isSig(n)==1;
%         plot(n,MTempRL(n)+0.01,'r*');
%     end
% end
% xlim([0 length(MTempRL)+1]);

[r,p]=corr([1:1:length(MTempRL)-1]',MTempRL(2:end)');
%r=0.7334
%p=0.0158
[r,p]=corr([1:1:length(MTempRL)-1]',corrAllFOVUncommon.dfofSigMeanToOthers(2:end)');

title(['sig meanToOthers New only,p=',num2str(p)]);

sig=[];

A=corrAllFOVs.dfofSig{3};
B=corrAllFOVUncommon.dfofSig{3};

for n=1:size(A,2);
    [~,p]=ttest2(A(:,n),B{n});
    if p<0.05;
        sig(n,1)=1;
    else
        sig(n,1)=0;
    end
end

saveas(gcf,'corr_commonUncommon.fig');

%% look at all cell activity
dfofSigNormAllFOVs{1}=[];
dfofSigNormAllFOVs{2}=[];
dfofSigNormAllFOVs{3}=[];
dfofSigNormAllFOVs{4}=[];
dfofSigNormAllFOVs{5}=[];
dfofSigNormAllFOVs{6}=[];
dfofSigNormAllFOVs{7}=[];
dfofSigNormAllFOVs{8}=[];
dfofSigNormAllFOVs{9}=[];
dfofSigNormAllFOVs{10}=[];
dfofSigNormAllFOVs{11}=[];

p=pwd;

for n=1:length(foldersAll);
  cd(foldersAll{n});
  disp(n)
  
cd ..\
load('allDfof_sigNorm.mat');
   for m=1:length(dfofSigNormAllFOVs);
      dfofSigNormAllFOVs{m}(end+1:end+size(allDfof_sigNorm{m},1),:)=allDfof_sigNorm{m};
   end
  end


cd(p);

save('dfofSigNormAllFOVs.mat','dfofSigNormAllFOVs');



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
for n=1:size(dfofSigNormAllFOVs{day},1)
[~,lagsDfof(n,1)]=findCorrLag(dfofSigNormAllFOVs{day}(n,:),temp,lagBins);
end

[~,orderTempRL1]=sort(lagsDfof);
dfofSigNormAllFOVsSort_TempRL1={};
for n=1:length(dfofSigNormAllFOVs);
    dfofSigNormAllFOVsSort_TempRL1{n}=dfofSigNormAllFOVs{n}(orderTempRL1,:);
end

figure
for n=1:length(dfofSigNormAllFOVs);
    subplot(1,length(dfofSigNormAllFOVs),n);
    imagesc(dfofSigNormAllFOVsSort_TempRL1{n})
    
    title(['Day',num2str(n)]);
    axis off
end
tightfig;

save('orderTempRL1.mat','orderTempRL1');
save('dfofSigNormAllFOVsSort_TempRL1.mat','dfofSigNormAllFOVsSort_TempRL1')
saveas(gcf,'dfofSigNormAllFOVsSort_TempRL1.fig');

%% for the matrix above (sorted based on day1) we look at the matrix correlation

matrixCorr=[]; %matrix n to matrix n+1
load('dfofSigNormAllFOVsSort_TempRL1.mat');
a=dfofSigNormAllFOVsSort_TempRL1;
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
lagBins=60;
day=2;
for n=1:size(dfofSigNormAllFOVs{day},1)
[~,lagsDfof(n,1)]=findCorrLag(dfofSigNormAllFOVs{day}(n,:),temp,lagBins);
end

[~,orderTempRL2]=sort(lagsDfof);
dfofSigNormAllFOVsSort_TempRL2={};
for n=1:length(dfofSigNormAllFOVs);
    dfofSigNormAllFOVsSort_TempRL2{n}=dfofSigNormAllFOVs{n}(orderTempRL2,:);
end

figure
for n=1:length(dfofSigNormAllFOVs);
     subplot(1,length(dfofSigNormAllFOVs),n);
    imagesc(dfofSigNormAllFOVsSort_TempRL2{n})
    title(['Day',num2str(n)]);
    axis off
end
tightfig;

save('orderTempRL2.mat','orderTempRL2');
save('dfofSigNormAllFOVsSort_TempRL2.mat','dfofSigNormAllFOVsSort_TempRL2')
saveas(gcf,'dfofSigNormAllFOVsSort_TempRL2.fig');
%% for the matrix above (sorted based on day2) we look at the matrix correlation

matrixCorr=[]; %matrix n to matrix n+1
load('dfofSigNormAllFOVsSort_TempRL2.mat');
a=dfofSigNormAllFOVsSort_TempRL2;
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
for n=1:size(dfofSigNormAllFOVs{day},1)
[~,lagsDfof(n,1)]=findCorrLag(dfofSigNormAllFOVs{day}(n,:),temp,lagBins);
end

[~,orderTempRL11]=sort(lagsDfof);
dfofSigNormAllFOVsSort_TempRL11={};
for n=1:length(dfofSigNormAllFOVs);
    dfofSigNormAllFOVsSort_TempRL11{n}=dfofSigNormAllFOVs{n}(orderTempRL11,:);
end

figure
for n=1:length(dfofSigNormAllFOVs);
     subplot(1,length(dfofSigNormAllFOVs),n);
    imagesc(dfofSigNormAllFOVsSort_TempRL11{n})
    title(['Day',num2str(n)]);
    axis off
end
tightfig;

save('orderTempRL11.mat','orderTempRL11');
save('dfofSigNormAllFOVsSort_TempRL11.mat','dfofSigNormAllFOVsSort_TempRL11')
saveas(gcf,'dfofSigNormAllFOVsSort_TempRL11.fig');
%% for the matrix above (sorted based on day2) we look at the matrix correlation

matrixCorr=[]; %matrix n to matrix n+1
load('dfofSigNormAllFOVsSort_TempRL11.mat');
a=dfofSigNormAllFOVsSort_TempRL11;
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

saveas(gcf,'matrixCorrelation_TempRL11.fig');
save('matrixCorrTempRL11.mat','matrixCorr');

%% sort them accoriding to heighest peak
lagsDfof=[];
%use the cue cell code to calculate lags
%sort cells based on the first day in new env
day=1;
for n=1:size(dfofSigNormAllFOVs{day},1)
[~,lagsDfof(n,1)]=max(dfofSigNormAllFOVs{day}(n,:));
end

[~,orderMax1]=sort(lagsDfof);
dfofSigNormAllFOVsSort_Max1={};
for n=1:length(dfofSigNormAllFOVs);
    dfofSigNormAllFOVsSort_Max1{n}=dfofSigNormAllFOVs{n}(orderMax1,:);  
end

figure
for n=1:length(dfofSigNormAllFOVs);
     subplot(1,length(dfofSigNormAllFOVs),n);
    imagesc(dfofSigNormAllFOVsSort_Max1{n})
     title(['Day',num2str(n)]);
    axis off
end
tightfig;

save('orderMax1.mat','orderMax1');
saveas(gcf,'dfofSigNormAllFOVsSort_Max1.fig');
save('dfofSigNormAllFOVsSort_Max1.mat','dfofSigNormAllFOVsSort_Max1')

%% for the matrix above (sorted based on day1) we look at the matrix correlation

matrixCorr=[]; %matrix n to matrix n+1
load('dfofSigNormAllFOVsSort_Max1.mat');
a=dfofSigNormAllFOVsSort_Max1;
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
for n=1:size(dfofSigNormAllFOVs{day},1)
[~,lagsDfof(n,1)]=max(dfofSigNormAllFOVs{day}(n,:));
end

[~,orderMax2]=sort(lagsDfof);
dfofSigNormAllFOVsSort_Max2={};
for n=1:length(dfofSigNormAllFOVs);
    dfofSigNormAllFOVsSort_Max2{n}=dfofSigNormAllFOVs{n}(orderMax2,:);  
end

figure
for n=1:length(dfofSigNormAllFOVs);
    subplot(1,length(dfofSigNormAllFOVs),n);
    imagesc(dfofSigNormAllFOVsSort_Max2{n})
     title(['Day',num2str(n)]);
    axis off
end
tightfig;

save('orderMax2.mat','orderMax2');
saveas(gcf,'dfofSigNormAllFOVsSort_Max2.fig');
save('dfofSigNormAllFOVsSort_Max2.mat','dfofSigNormAllFOVsSort_Max2')

%% for the matrix above (sorted based on day2) we look at the matrix correlation

matrixCorr=[]; %matrix n to matrix n+1
load('dfofSigNormAllFOVsSort_Max2.mat');
a=dfofSigNormAllFOVsSort_Max2;
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
for n=1:size(dfofSigNormAllFOVs{day},1)
[~,lagsDfof(n,1)]=max(dfofSigNormAllFOVs{day}(n,:));
end

[~,orderMax11]=sort(lagsDfof);
dfofSigNormAllFOVsSort_Max11={};

for n=1:length(dfofSigNormAllFOVs);
    dfofSigNormAllFOVsSort_Max11{n}=dfofSigNormAllFOVs{n}(orderMax11,:);  
end

figure
for n=1:length(dfofSigNormAllFOVs);
    subplot(1,length(dfofSigNormAllFOVs),n);
    imagesc(dfofSigNormAllFOVsSort_Max11{n})
     title(['Day',num2str(n)]);
    axis off
end
tightfig;

save('orderMax11.mat','orderMax11');
saveas(gcf,'dfofSigNormAllFOVsSort_Max11.fig');
save('dfofSigNormAllFOVsSort_Max11.mat','dfofSigNormAllFOVsSort_Max11')

%% for the matrix above (sorted based on day2) we look at the matrix correlation

matrixCorr=[]; %matrix n to matrix n+1
load('dfofSigNormAllFOVsSort_Max11.mat');
a=dfofSigNormAllFOVsSort_Max11;
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

%% plot the order of cells for old and new env
binWidth=25;
bins=[0:binWidth:ceil(size(dfofSigNormAllFOVs{1},1)/binWidth)*binWidth];
%for old and new
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
subplot(221)
imagesc(MTempRL)
colormap(flipud(gray))
% colormap(parula)
title('TempRL order: old and new');
axis equal
axis tight
xlabel('sort old');
ylabel('sort day1new');

subplot(222)
imagesc(MMax)
colormap(flipud(gray))
% colormap(parula)
title('MaxPeak order: old and new');
axis equal
axis tight
xlabel('sort old');
ylabel('sort day1new');

%for day 1 and 10 in new env
[y1,n1]=histc(orderTempRL2,bins);
[y2,n2]=histc(orderTempRL11,bins);

MTempRL=zeros(length(bins)-1,length(bins)-1);
for xx=1:length(bins)-1;
    for yy=1:length(bins)-1;
        A=(n1==xx)&(n2==yy);
        MTempRL(yy,xx)=length(find(A));
    end
end

[y1,n1]=histc(orderMax2,bins);
[y2,n2]=histc(orderMax11,bins);
MMax=zeros(length(bins)-1,length(bins)-1);
for xx=1:length(bins)-1;
    
    for yy=1:length(bins)-1;
        A=(n1==xx)&(n2==yy);
        MMax(yy,xx)=length(find(A));
    end
end

subplot(223)
imagesc(MTempRL)
colormap(flipud(gray))
% colormap(parula)
title('TempRL order: new 1 and 10');
axis equal
axis tight
xlabel('sort day1new');
ylabel('sort day10new');

subplot(224)
imagesc(MMax)
colormap(flipud(gray))
% colormap(parula)
title('MaxPeak order: new 1 and 10');
axis equal
axis tight
xlabel('sort day1new');
ylabel('sort day10new');

saveas(gcf,'cellOrders.fig');
saveas(gcf,'cellOrders.jpg');

%% now determine the cell activitiy correlation to it own activity across days
load('dfofSigNormAllFOVs.mat');
A=dfofSigNormAllFOVs;
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

% %% determine cell indices
% 
% indicesAll=[];
% indicesAll.Grid=[];
% indicesAll.CueL=[];
% indicesAll.CueR=[];
% indicesAll.SpeedP=[];
% indicesAll.SpeedN=[];
% indicesAll.GridIncludeCue=[];
% p=pwd;
% for n=1:length(foldersAll);
%     cd(foldersAll{n});
% 
%         load('indices.mat');
%         indicesAll.Grid(end+1:end+size(indices.commonCellsGrid,1),:)=indices.commonCellsGrid;
%         indicesAll.CueL(end+1:end+size(indices.commonCellsCueCellL,1),:)=indices.commonCellsCueCellL;
%         indicesAll.CueR(end+1:end+size(indices.commonCellsCueCellR,1),:)=indices.commonCellsCueCellR;
%         indicesAll.SpeedP(end+1:end+size(indices.commonCellsSpeedP,1),:)=indices.commonCellsSpeedP;
%         indicesAll.SpeedN(end+1:end+size(indices.commonCellsSpeedN,1),:)=indices.commonCellsSpeedN;
%         indicesAll.GridIncludeCue(end+1:end+size(indices.commonCellsGridIncludeCue,1),:)=indices.commonCellsGridIncludeCue;
% end
% 
%     
% cd(p);
% 
% %select the ones that are more than 50% is one cell type
% indicesAll.thresh=0.5;
% A=indicesAll.Grid;
% A=sum(A,2);
% indicesAll.GridIdx=find(A>size(indicesAll.Grid,2)*indicesAll.thresh);
% indicesAll.nonGridIdx=find(A<=size(indicesAll.Grid,2)*indicesAll.thresh);
% 
% A=indicesAll.GridIncludeCue;
% A=sum(A,2);
% indicesAll.GridIncludeCueIdx=find(A>size(indicesAll.GridIncludeCue,2)*indicesAll.thresh);
% 
% A=indicesAll.CueL;
% A=sum(A,2);
% indicesAll.CueLIdx=find(A>size(indicesAll.CueL,2)*indicesAll.thresh);
% 
% A=indicesAll.CueR;
% A=sum(A,2);
% indicesAll.CueRIdx=find(A>size(indicesAll.CueR,2)*indicesAll.thresh);
% 
% A=indicesAll.SpeedP;
% A=sum(A,2);
% indicesAll.SpeedPIdx=find(A>size(indicesAll.SpeedP,2)*indicesAll.thresh);
% 
% A=indicesAll.SpeedN;
% A=sum(A,2);
% indicesAll.SpeedNIdx=find(A>size(indicesAll.SpeedN,2)*indicesAll.thresh);
% 
% B=unique([indicesAll.GridIdx;indicesAll.CueLIdx;indicesAll.CueRIdx]);
% indicesAll.nonGridnonCueIdx=setdiff([1:1:size(indicesAll.Grid,1)],B);
% save('indicesAll.mat','indicesAll');
% 
% 
% figure,
% subplot(161);
% imagesc(indicesAll.Grid);
% title('grid cells');
% subplot(162);
% imagesc(indicesAll.GridIncludeCue);
% title('right cue cells');
% subplot(163);
% imagesc(indicesAll.CueR);
% title('right cue cells');
% subplot(164);
% imagesc(indicesAll.CueL);
% title('left cue cells');
% subplot(165);
% imagesc(indicesAll.SpeedP);
% title('positive speed cells');
% subplot(166);
% imagesc(indicesAll.SpeedN);
% title('negative speed cells');
% 
% saveas(gcf,'cellTypes.fig');


%% determine cell indices using new cue cell threshold
load('foldersAll.mat');
indicesAllNewCueThresh=[];
indicesAllNewCueThresh.Grid=[];
indicesAllNewCueThresh.CueL=[];
indicesAllNewCueThresh.CueR=[];
indicesAllNewCueThresh.SpeedP=[];
indicesAllNewCueThresh.SpeedN=[];
indicesAllNewCueThresh.GridIncludeCue=[];
p=pwd;
for n=1:length(foldersAll);
    cd(foldersAll{n});

        load('reIdentifyCueCellsUniThreshAddingNewMice\indices.mat');
        indicesAllNewCueThresh.Grid(end+1:end+size(indices.commonCellsGrid,1),:)=indices.commonCellsGrid;
        indicesAllNewCueThresh.CueL(end+1:end+size(indices.commonCellsCueCellL,1),:)=indices.commonCellsCueCellL;
        indicesAllNewCueThresh.CueR(end+1:end+size(indices.commonCellsCueCellR,1),:)=indices.commonCellsCueCellR;
        indicesAllNewCueThresh.SpeedP(end+1:end+size(indices.commonCellsSpeedP,1),:)=indices.commonCellsSpeedP;
        indicesAllNewCueThresh.SpeedN(end+1:end+size(indices.commonCellsSpeedN,1),:)=indices.commonCellsSpeedN;
        indicesAllNewCueThresh.GridIncludeCue(end+1:end+size(indices.commonCellsGridIncludeCue,1),:)=indices.commonCellsGridIncludeCue;
        indicesAllNewCueThresh.perGridAllThreshIndivFOVs(n,:)=indices.perGridAllThresh;
        indicesAllNewCueThresh.perCueLAllThreshIndivFOVs(n,:)=indices.perCueLAllThresh;
        indicesAllNewCueThresh.perCueRAllThreshIndivFOVs(n,:)=indices.perCueRAllThresh;
        indicesAllNewCueThresh.perSpeedPAllThreshIndivFOVs(n,:)=indices.perSpeedPAllThresh;
        indicesAllNewCueThresh.perSpeedNAllThreshIndivFOVs(n,:)=indices.perSpeedNAllThresh;

end

    
cd(p);

%select the ones that are more than 50% is one cell type
indicesAllNewCueThresh.thresh=size(indicesAllNewCueThresh.Grid,2)*0.5;
indicesAllNewCueThresh.allThreshs=[1:1:size(indicesAllNewCueThresh.Grid,2)];
A=indicesAllNewCueThresh.Grid;
A=sum(A,2);
indicesAllNewCueThresh.GridIdx=find(A>indicesAllNewCueThresh.thresh);
indicesAllNewCueThresh.nonGridIdx=find(A<=indicesAllNewCueThresh.thresh);

A=indicesAllNewCueThresh.Grid;
A=sum(A,2);
for n=1:length(indicesAllNewCueThresh.allThreshs);
indicesAllNewCueThresh.GridIdxLowerThreshs{n}=find(A>=indicesAllNewCueThresh.allThreshs(n));
end


A=indicesAllNewCueThresh.GridIncludeCue;
A=sum(A,2);
indicesAllNewCueThresh.GridIncludeCueIdx=find(A>indicesAllNewCueThresh.thresh);

A=indicesAllNewCueThresh.CueL;
A=sum(A,2);
indicesAllNewCueThresh.CueLIdx=find(A>indicesAllNewCueThresh.thresh);
A=indicesAllNewCueThresh.CueL;
A=sum(A,2);
for n=1:length(indicesAllNewCueThresh.allThreshs);
indicesAllNewCueThresh.CueLIdxLowerThreshs{n}=find(A>=indicesAllNewCueThresh.allThreshs(n));
end


A=indicesAllNewCueThresh.CueR;
A=sum(A,2);
indicesAllNewCueThresh.CueRIdx=find(A>indicesAllNewCueThresh.thresh);

A=indicesAllNewCueThresh.CueR;
A=sum(A,2);
for n=1:length(indicesAllNewCueThresh.allThreshs);
indicesAllNewCueThresh.CueRIdxLowerThreshs{n}=find(A>=indicesAllNewCueThresh.allThreshs(n));
end

A=indicesAllNewCueThresh.SpeedP;
A=sum(A,2);
indicesAllNewCueThresh.SpeedPIdx=find(A>indicesAllNewCueThresh.thresh);

A=indicesAllNewCueThresh.SpeedP;
A=sum(A,2);
for n=1:length(indicesAllNewCueThresh.allThreshs);
indicesAllNewCueThresh.SpeedPIdxLowerThreshs{n}=find(A>=indicesAllNewCueThresh.allThreshs(n));
end


A=indicesAllNewCueThresh.SpeedN;
A=sum(A,2);
indicesAllNewCueThresh.SpeedNIdx=find(A>indicesAllNewCueThresh.thresh);

A=indicesAllNewCueThresh.SpeedN;
A=sum(A,2);
for n=1:length(indicesAllNewCueThresh.allThreshs);
indicesAllNewCueThresh.SpeedNIdxLowerThreshs{n}=find(A>=indicesAllNewCueThresh.allThreshs(n));
end


B=unique([indicesAllNewCueThresh.GridIdx;indicesAllNewCueThresh.CueLIdx;indicesAllNewCueThresh.CueRIdx]);
indicesAllNewCueThresh.nonGridnonCueIdx=setdiff([1:1:size(indicesAllNewCueThresh.Grid,1)],B);
save('indicesAllNewCueThresh.mat','indicesAllNewCueThresh');


figure,
subplot(161);
imagesc(indicesAllNewCueThresh.Grid);
title('grid cells');
subplot(162);
imagesc(indicesAllNewCueThresh.GridIncludeCue);
title('right cue cells');
subplot(163);
imagesc(indicesAllNewCueThresh.CueR);

title('right cue cells');
subplot(164);
imagesc(indicesAllNewCueThresh.CueL);
title('left cue cells');
subplot(165);
imagesc(indicesAllNewCueThresh.SpeedP);
title('positive speed cells');
subplot(166);
imagesc(indicesAllNewCueThresh.SpeedN);
title('negative speed cells');

saveas(gcf,'cellTypesNewCueThresh.fig');



%% shift of all activity: according to cue score calculation
p=pwd;
load('foldersAll');
sigLagL=[];
sigLagR=[];

for n=1:length(foldersAll);
    disp(n)

    cd(foldersAll{n});
    load('useFolders.mat');
    load('commonCells.mat');   
  
    allLagLCommonCells=[];
     allLagRCommonCells=[];
    for m=1:length(useFolders);
        disp(m)
        cd(useFolders{m});
        
        load('allROIs.mat');
        allLagL=nan(size(roi,3),1);
        allLagR=nan(size(roi,3),1);
        load('cueAnalysisNew_sig\useIdx.mat');
        load('cueAnalysisNew_sig\newScoreShuffleTemplate\Left\cueCells.mat');
        allLagL(useIdx,1)=cueCells.realLags;
        allLagLCommonCells(:,m)=allLagL(commonCells(:,m));
        load('cueAnalysisNew_sig\newScoreShuffleTemplate\Right\cueCells.mat');
        allLagR(useIdx,1)=cueCells.realLags;
        allLagRCommonCells(:,m)=allLagR(commonCells(:,m));
    end
    
    sigLagL(end+1:end+size(allLagLCommonCells,1),:)=allLagLCommonCells;
   sigLagR(end+1:end+size(allLagRCommonCells,1),:)=allLagRCommonCells; 
end
cd(p);
save('sigLagL.mat','sigLagL');
save('sigLagR.mat','sigLagR');
        
%% distribution of fields

load('foldersAll.mat');
p=pwd;
foldersAllIndvFOV={};
commonCellsAllIdx={};

for n=1:length(foldersAll);
    cd(foldersAll{n});
    load('useFolders.mat');
    load('commonCells.mat');   
   foldersAllIndvFOV{n}=useFolders;
commonCellsAllIdx{n}=commonCells; 
    
end

cd(p);
save('foldersAllIndvFOV.mat','foldersAllIndvFOV');
save('commonCellsAllIdx.mat','commonCellsAllIdx');

%% 
%for all follwoing data, each cell is one day
fieldDistri={};%each cell is the bin number of all cells on that day, each neuron has one cell in this structure 
fieldDistriCat={};%the individual cells for each neuron in "fieldDistri" were contatinated together
fieldDistriAmp={};%field only activity on individual days
fieldDistriAmpNorm={};%each activity is normalized so that it's ranged from 0 to 1
fieldDistriAmpSum=[];%add eah cell of fieldDistriAmp together
fieldDistriAmpNormSum=[];%add eah cell of fieldDistriAmpNorm together
fieldDistriEachBin=[];%in this data, the number of field discovered in each bin was calculated

for n=1:length(foldersAllIndvFOV{1});
    disp(n)
 fieldDistri{n}={};
 fieldDistriCat{n}=[];
fieldDistriAmp{n}=[];
fieldDistriAmpNorm{n}=[];

for m=1:length(foldersAllIndvFOV);
     disp(m)
 cd(foldersAllIndvFOV{m}{n});
  B=commonCellsAllIdx{m}(:,n);
   load('PValueClassifier_KY2_6_sig\allCellsCorrected.mat');
for i=1:length(B);
        D=allCellsCorrected.inFieldBins{B(i)};
    fieldDistri{n}{end+1}= D';
   C=allCellsCorrected.dfofaveragesmoothFields(:,B(i));
   fieldDistriAmp{n}(:,end+1)=C;
   E=C/max(C);
  fieldDistriAmpNorm{n}(:,end+1)=E;
end
end
fieldDistriCat{n}=cell2mat(fieldDistri{n});
fieldDistriAmp{n}=fieldDistriAmp{n}';
fieldDistriAmpNorm{n}=fieldDistriAmpNorm{n}';

fieldDistriAmp{n}(isnan(fieldDistriAmp{n}))=0;
fieldDistriAmpNorm{n}(isnan(fieldDistriAmpNorm{n}))=0;

fieldDistriAmpSum(n,:)=sum(fieldDistriAmp{n},1);
fieldDistriAmpNormSum(n,:)=sum(fieldDistriAmpNorm{n},1);

for m=1:max(fieldDistriCat{1});%200bins for 10m track
    fieldDistriEachBin(n,m)=length(find(fieldDistriCat{n}==m));
end
end
cd(p);
save('fieldDistri.mat','fieldDistri');
save('fieldDistriCat.mat','fieldDistriCat');
save('fieldDistriAmp.mat','fieldDistriAmp');
save('fieldDistriAmpNorm.mat','fieldDistriAmpNorm');
save('fieldDistriAmpSum.mat','fieldDistriAmpSum');
save('fieldDistriAmpNormSum.mat','fieldDistriAmpNormSum');
save('fieldDistriEachBin.mat','fieldDistriEachBin');

%% plot them: distribution

load('E:\ID20201118\20201224\loc1\pcaica\cueAnalysisNew_sig\tempRL.mat');
oldTempRL=tempRL;
load('E:\ID20201118\20210103\1stloc1\pcaica\cueAnalysisNew_sig\tempRL');
figure,
for n=1:length(fieldDistri);
    subplot(length(fieldDistriCat),2,2*(n-1)+1);
    plot(fieldDistriEachBin(n,:))
    hold on
    if n==1;
        plot(oldTempRL*140);
    else
    plot(tempRL*140);
    hold on
       plot([890/5 890/5],[0 140],'r');
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

for n=1:length(fieldDistriCat)-1;
    [rnext(n,1),pnext(n,1)]=kstest2(fieldDistriCat{n},fieldDistriCat{n+1});
end

for n=1:length(fieldDistriCat);
    subplot(length(fieldDistriCat),2,2*n);
    [p,x]=ksdensity(fieldDistriCat{n},'width',0.5);
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
     [r1(n,1),p1(n,1)]=kstest2(fieldDistriCat{1},fieldDistriCat{n});
     if n>1;
     if rnext(n-1,1)==1;
         title('Sig Dif with previous day','Color','r');
     else
         title('Insig with previous day','Color','k');
     end
     end
     
end
  
tightfig;    
saveas(gcf,'fieldDistribution.fig');
saveas(gcf,'fieldDistribution.jpg');

%look at whether the cell activity is more consistent with the cue template

corrWithTemp=[];
% binnedOldTempRL=oldTempRL(1:2:length(oldTempRL)-1)+oldTempRL(1:2:length(oldTempRL));
% binnedTempRL=tempRL(1:2:length(tempRL)-1)+tempRL(1:2:length(tempRL));

for n=1:length(fieldDistriCat);
    [p,x]=ksdensity(fieldDistriCat{n},'width',1);
    if n==1;
    corrWithTemp(n,1)=corr(fieldDistriEachBin(n,:)',oldTempRL);
    else
        corrWithTemp(n,1)=corr(fieldDistriEachBin(n,:)',tempRL);
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
    corrAmp(n)=corr(fieldDistriAmpSum(n,:)',fieldDistriAmpSum(n+1,:)');
end

corrAmpNorm=[];%correlation to the next activity norm adding up together

for n=1:length(fieldDistri)-1;
      corrAmpNorm(n)=corr(fieldDistriAmpNormSum(n,:)',fieldDistriAmpNormSum(n+1,:)');
end

figure;

for n=1:length(fieldDistri);
      subplot(length(fieldDistri),2,2*(n-1)+1);
      plot(fieldDistriAmpSum(n,:));
      hold on
      if n==1;
         plot(oldTempRL*max(fieldDistriAmpSum(1,:)));  
      else
       plot(tempRL*max(fieldDistriAmpSum(2,:)));
       hold on
       plot([890/5 890/5],[0 max(fieldDistriAmpSum(2,:))],'r');
      end
     if n<length(useFolders);
               title(['corr to next = ',num2str(corrAmp(n))])
     end        
      axis off
    
      subplot(length(fieldDistri),2,2*n);
        plot(fieldDistriAmpNormSum(n,:));
         hold on
      if n==1;
         plot(oldTempRL*max(fieldDistriAmpNormSum(1,:)));  
      else
       plot(tempRL*max(fieldDistriAmpNormSum(2,:)));
       hold on
       plot([890/5 890/5],[0 max(fieldDistriAmpNormSum(2,:))],'r');
      end
        if n<length(useFolders);
               title(['corr to next = ',num2str(corrAmpNorm(n))])
            
        end
                
axis off
end
tightfig;

saveas(gcf,'fieldDistriAmp.fig');
saveas(gcf,'fieldDistriAmp.jpg');

%%
x=[];
p=[];

for n=1:length(fieldDistri);   
[x(n,:),p(n,:)]=ksdensity(fieldDistriCat{n},'width',5);
end

figure,

subplot(411)
plot(p(1,:)*5-2.5,mean(x(4:11,:),1)-mean(x(2:3,:),1))
hold on
plot([2.5:5:1000],tempRL*0.001)
plot([890 890],[0 0.001],'r');

xlim([-15*5 215*5]);
title('Distribution Day (4-11)-(2-3)');

subplot(412)
plot([1:1:200],mean(fieldDistriEachBin(4:11,:),1)-mean(fieldDistriEachBin(2:3,:),1))
hold on
plot([1:1:200],tempRL*20);
hold on
plot([890/5 890/5],[0 20],'r');
xlim([0 200]);
title('Distribution each bin Day (4-11)-(2-3)');


subplot(413)
plot([1:1:size(fieldDistriAmpSum,2)]*5-2.5,mean(fieldDistriAmpSum(4:11,:),1)-mean(fieldDistriAmpSum(2:3,:),1))
hold on
plot([2.5:5:1000],tempRL*5)
hold on
plot([890 890],[0 5],'r');

xlim([-15*5 215*5]);
title('Amplitude Day (4-11)-(2-3)');

subplot(414)
plot([1:1:size(fieldDistriAmpNormSum,2)]*5-2.5,mean(fieldDistriAmpNormSum(4:11,:),1)-mean(fieldDistriAmpNormSum(2:3,:),1))
hold on
plot([2.5:5:1000],tempRL*10)
hold on
plot([890 890],[0 10],'r');

xlim([-15*5 215*5]);
title('Amplitude norm Day (4-11)-(2-3)');
saveas(gcf,'disriAmpBeforeAfterLearning.fig');
saveas(gcf,'disriAmpBeforeAfterLearning.jpg');

%% test whether any of these changes are significant
%working on field distribution in each bin
load('fieldDistriEachBin.mat');
a=fieldDistriEachBin;
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
x=[];
p=[];

for n=1:length(fieldDistri);   
[x(n,:),p(n,:)]=ksdensity(fieldDistriCat{n},'width',5);
end

figure,

subplot(411)
plot(p(1,:)*5-2.5,mean(x(7:11,:),1)-mean(x(2:6,:),1))
hold on
plot([2.5:5:1000],tempRL*0.001)
plot([890 890],[0 0.001],'r');

xlim([-15*5 215*5]);
title('Distribution Day (7-11)-(2-6)');

subplot(412)
plot([1:1:200],mean(fieldDistriEachBin(7:11,:),1)-mean(fieldDistriEachBin(2:6,:),1))
hold on
plot([1:1:200],tempRL*20);
hold on
plot([890/5 890/5],[0 20],'r');
xlim([0 200]);
title('Distribution each bin Day (7-11)-(2-6)');


subplot(413)
plot([1:1:size(fieldDistriAmpSum,2)]*5-2.5,mean(fieldDistriAmpSum(7:11,:),1)-mean(fieldDistriAmpSum(2:6,:),1))
hold on
plot([2.5:5:1000],tempRL*5)
hold on
plot([890 890],[0 5],'r');

xlim([-15*5 215*5]);
title('Amplitude Day (7-11)-(2-6)');

subplot(414)
plot([1:1:size(fieldDistriAmpNormSum,2)]*5-2.5,mean(fieldDistriAmpNormSum(7:11,:),1)-mean(fieldDistriAmpNormSum(2:6,:),1))
hold on
plot([2.5:5:1000],tempRL*10)
hold on
plot([890 890],[0 10],'r');

xlim([-15*5 215*5]);
title('Amplitude norm Day (7-11)-(2-6)');
saveas(gcf,'disriAmpBeforeAfterLearning_diffDays.fig');
saveas(gcf,'disriAmpBeforeAfterLearning_diffDays.jpg');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%look at field distribution of tracked and untracked cells

p=pwd;
load('foldersAll.mat');
load('foldersAllIndvFOV.mat');

%for all follwoing data, each cell is one day
fieldDistriUC={};%each cell is the bin number of all cells on that day, each neuron has one cell in this structure 
fieldDistriCatUC={};%the individual cells for each neuron in "fieldDistri" were contatinated together
fieldDistriAmpUC={};%field only activity on individual days
fieldDistriAmpNormUC={};%each activity is normalized so that it's ranged from 0 to 1
fieldDistriAmpSumUC=[];%add eah cell of fieldDistriAmp together
fieldDistriAmpNormSumUC=[];%add eah cell of fieldDistriAmpNorm together
fieldDistriEachBinUC=[];%in this data, the number of field discovered in each bin was calculated

p=pwd;
uncommonCellsAll=[];
  
for n=1:length(foldersAll);
    cd(foldersAll{n});
   uncommonCellsAll{n}=[];
load('commonUncommon\uncommonCells.mat');
uncommonCellsAll{n}=uncommonCells; %n is each animal, within that, each cell is one day
end

cd(p)

NUncommonCellsAll=[];
for n=1:length(uncommonCellsAll{1});
    N=0;
    for m=1:length(foldersAll);
        N=N+length(uncommonCellsAll{m}{n});
    end
    NUncommonCellsAll(n,1)=N;
end


for n=1:length(foldersAllIndvFOV{1});
    disp(n)
 fieldDistriUC{n}={};
 fieldDistriCatUC{n}=[];
fieldDistriAmpUC{n}=[];
fieldDistriAmpNormUC{n}=[];

for m=1:length(foldersAllIndvFOV);
     disp(m)
 cd(foldersAllIndvFOV{m}{n});
  B=uncommonCellsAll{m}{n};%going through each animal across taht day
   load('PValueClassifier_KY2_6_sig\allCellsCorrected.mat');
for i=1:length(B);
        D=allCellsCorrected.inFieldBins{B(i)};
    fieldDistriUC{n}{end+1}=D';
   C=allCellsCorrected.dfofaveragesmoothFields(:,B(i));
   fieldDistriAmpUC{n}(:,end+1)=C;
   E=C/max(C);
  fieldDistriAmpNormUC{n}(:,end+1)=E;
end
end
fieldDistriCatUC{n}=cell2mat(fieldDistriUC{n});
fieldDistriAmpUC{n}=fieldDistriAmpUC{n}';
fieldDistriAmpNormUC{n}=fieldDistriAmpNormUC{n}';

fieldDistriAmpUC{n}(isnan(fieldDistriAmpUC{n}))=0;
fieldDistriAmpNormUC{n}(isnan(fieldDistriAmpNormUC{n}))=0;

fieldDistriAmpSumUC(n,:)=sum(fieldDistriAmpUC{n},1)/NUncommonCellsAll(n);%normalized by the number of cells on each day
fieldDistriAmpNormSumUC(n,:)=sum(fieldDistriAmpNormUC{n},1)/NUncommonCellsAll(n);%normalized by the number of cells on each day

for m=1:max(fieldDistriCatUC{1});%200bins for 10m track
    fieldDistriEachBinUC(n,m)=length(find(fieldDistriCatUC{n}==m));
end
end

%normalize field distribution in each bin
for n=1:size(fieldDistriEachBinUC,1);
    fieldDistriEachBinUC(n,:)= fieldDistriEachBinUC(n,:)/NUncommonCellsAll(n);
end


cd(p);
save('fieldDistriUC.mat','fieldDistriUC');
save('fieldDistriCatUC.mat','fieldDistriCatUC');
save('fieldDistriAmpUC.mat','fieldDistriAmpUC');
save('fieldDistriAmpNormUC.mat','fieldDistriAmpNormUC');
save('fieldDistriAmpSumUC.mat','fieldDistriAmpSumUC');
save('fieldDistriAmpNormSumUC.mat','fieldDistriAmpNormSumUC');
save('fieldDistriEachBinUC.mat','fieldDistriEachBinUC');

%% plot them: distribution

load('E:\ID20201118\20201224\loc1\pcaica\cueAnalysisNew_sig\tempRL.mat');
oldTempRL=tempRL;
load('E:\ID20201118\20210103\1stloc1\pcaica\cueAnalysisNew_sig\tempRL');
figure,
for n=1:length(fieldDistriUC);
    subplot(length(fieldDistriCatUC),2,2*(n-1)+1);
    plot(fieldDistriEachBinUC(n,:))
    hold on
    if n==1;
        plot(oldTempRL*140);
    else
    plot(tempRL*140);
    hold on
       plot([890/5 890/5],[0 140],'r');
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

for n=1:length(fieldDistriCatUC)-1;
    [rnext(n,1),pnext(n,1)]=kstest2(fieldDistriCatUC{n},fieldDistriCatUC{n+1});
end

for n=1:length(fieldDistriCatUC);
    subplot(length(fieldDistriCatUC),2,2*n);
    [p,x]=ksdensity(fieldDistriCatUC{n},'width',0.5);
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
     [r1(n,1),p1(n,1)]=kstest2(fieldDistriCatUC{1},fieldDistriCatUC{n});
     if n>1;
     if rnext(n-1,1)==1;
         title('Sig Dif with previous day','Color','r');
     else
         title('Insig with previous day','Color','k');
     end
     end
     
end
  
tightfig;    
saveas(gcf,'fieldDistributionUC.fig');
saveas(gcf,'fieldDistributionUC.jpg');

%look at whether the cell activity is more consistent with the cue template

corrWithTemp=[];
% binnedOldTempRL=oldTempRL(1:2:length(oldTempRL)-1)+oldTempRL(1:2:length(oldTempRL));
% binnedTempRL=tempRL(1:2:length(tempRL)-1)+tempRL(1:2:length(tempRL));

for n=1:length(fieldDistriCatUC);
    [p,x]=ksdensity(fieldDistriCatUC{n},'width',1);
    if n==1;
    corrWithTemp(n,1)=corr(fieldDistriEachBinUC(n,:)',oldTempRL);
    else
        corrWithTemp(n,1)=corr(fieldDistriEachBinUC(n,:)',tempRL);
    end
end
figure,plot(corrWithTemp);
xlabel('days')
ylabel('correlation to temp');
[r,p]=corr(corrWithTemp(2:end),[1:1:10]');
title(['p = ',num2str(p)]);
saveas(gcf,'distriCorrToTempUC.fig');
saveas(gcf,'distriCorrToTempUC.jpg');

%% plot them: amplitude
corrAmp=[];%correlation to the next activity norm adding up together

for n=1:length(fieldDistriUC)-1;
    corrAmp(n)=corr(fieldDistriAmpSumUC(n,:)',fieldDistriAmpSumUC(n+1,:)');
end

corrAmpNorm=[];%correlation to the next activity norm adding up together

for n=1:length(fieldDistriUC)-1;
      corrAmpNorm(n)=corr(fieldDistriAmpNormSumUC(n,:)',fieldDistriAmpNormSumUC(n+1,:)');
end

figure;

for n=1:length(fieldDistriUC);
      subplot(length(fieldDistriUC),2,2*(n-1)+1);
      plot(fieldDistriAmpSumUC(n,:));
      hold on
      if n==1;
         plot(oldTempRL*max(fieldDistriAmpSumUC(1,:)));  
      else
       plot(tempRL*max(fieldDistriAmpSumUC(2,:)));
       hold on
       plot([890/5 890/5],[0 max(fieldDistriAmpSumUC(2,:))],'r');
      end
     if n<length(useFolders);
               title(['corr to next = ',num2str(corrAmp(n))])
     end        
      axis off
    
      subplot(length(fieldDistriUC),2,2*n);
        plot(fieldDistriAmpNormSumUC(n,:));
         hold on
      if n==1;
         plot(oldTempRL*max(fieldDistriAmpNormSumUC(1,:)));  
      else
       plot(tempRL*max(fieldDistriAmpNormSumUC(2,:)));
       hold on
       plot([890/5 890/5],[0 max(fieldDistriAmpNormSumUC(2,:))],'r');
      end
        if n<length(useFolders);
               title(['corr to next = ',num2str(corrAmpNorm(n))])
            
        end
                
axis off
end
tightfig;

saveas(gcf,'fieldDistriAmpUC.fig');
saveas(gcf,'fieldDistriAmpUC.jpg');

%%
x=[];
p=[];

for n=1:length(fieldDistriUC);   
[x(n,:),p(n,:)]=ksdensity(fieldDistriCatUC{n},'width',5);
end

figure,

subplot(411)
plot(p(1,:)*5-2.5,mean(x(4:11,:),1)-mean(x(2:3,:),1))
hold on
plot([2.5:5:1000],tempRL*0.001)
plot([890 890],[0 0.001],'r');

xlim([-15*5 215*5]);
title('Distribution Day (4-11)-(2-3)');

subplot(412)
plot([1:1:200],mean(fieldDistriEachBinUC(4:11,:),1)-mean(fieldDistriEachBinUC(2:3,:),1))
hold on
plot([1:1:200],tempRL*0.03);
hold on
plot([890/5 890/5],[0 0.03],'r');
xlim([0 200]);
title('Distribution each bin Day (4-11)-(2-3)');


subplot(413)
plot([1:1:size(fieldDistriAmpSumUC,2)]*5-2.5,mean(fieldDistriAmpSumUC(4:11,:),1)-mean(fieldDistriAmpSumUC(2:3,:),1))
hold on
plot([2.5:5:1000],tempRL*0.005)
hold on
plot([890 890],[0 0.005],'r');

xlim([-15*5 215*5]);
title('Amplitude Day (4-11)-(2-3)');

subplot(414)
plot([1:1:size(fieldDistriAmpNormSumUC,2)]*5-2.5,mean(fieldDistriAmpNormSumUC(4:11,:),1)-mean(fieldDistriAmpNormSumUC(2:3,:),1))
hold on
plot([2.5:5:1000],tempRL*0.01)
hold on
plot([890 890],[0 0.01],'r');

xlim([-15*5 215*5]);
title('Amplitude norm Day (4-11)-(2-3)');
saveas(gcf,'disriAmpBeforeAfterLearningUC.fig');
saveas(gcf,'disriAmpBeforeAfterLearningUC.jpg');

%% test whether any of these changes are significant
%working on field distribution in each bin
load('fieldDistriEachBinUC.mat');
a=fieldDistriEachBinUC;
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
sigIdxFieldDistriUC=find(sig==1);
figure,plot(aa,'k');
hold on
plot(sigIdxFieldDistriUC,aa(sigIdxFieldDistriUC),'r.','MarkerSize',10);

saveas(gcf,'disriAmpBeforeAfterLearning_sigUC.fig');
save('sigIdxFieldDistriUC.mat','sigIdxFieldDistriUC');
        
%% 
load('fieldDistriEachBinUC.mat');
load('sigIdxFieldDistriUC.mat');
a=fieldDistriEachBinUC;
aa=mean(a(4:11,:),1)-mean(a(2:3,:),1);

allSigUC=zeros(200,1);
allSigUC(sigIdxFieldDistriUC)=1;
S=contiguous(allSigUC,[1]);
CS=S{1,2};

figure,plot(aa,'k');
for n=1:size(CS,1);
    idx=[CS(n,1):1:CS(n,2)];
    hold on
plot(idx,aa(idx),'r.','MarkerSize',5);
hold on 
plot(idx,aa(idx),'r','LineWidth',1);
end

load('fieldDistriEachBin.mat');
load('sigIdxFieldDistri.mat');
load('dfofSigNormAllFOVs.mat');
nCommonCell=size(dfofSigNormAllFOVs{1},1);
a=fieldDistriEachBin/nCommonCell;
aa=mean(a(4:11,:),1)-mean(a(2:3,:),1);

allSig=zeros(200,1);
allSig(sigIdxFieldDistri)=1;
S=contiguous(allSig,[1]);
CS=S{1,2};

hold on,plot(aa,'g');
for n=1:size(CS,1);
    idx=[CS(n,1):1:CS(n,2)];
    hold on
plot(idx,aa(idx),'r.','MarkerSize',5);
hold on 
plot(idx,aa(idx),'r','LineWidth',1);
end

saveas(gcf,'disriAmpBeforeAfterLearning_sigCommonUncommon.fig');
title('disriAmpBeforeAfterLearning_sigCommonUncommon');
% hold on,plot(aa,'g');
% hold on
% plot(sigIdxFieldDistri,aa(sigIdxFieldDistri),'r.','MarkerSize',10);





%% a better method is to look at each bin across days and ask which bin shows significant increase of filed scross days
% sig=[];
% load('fieldDistriAmpNormSum.mat');
% A=fieldDistriAmpNormSum(2:end,:);
% M=[];
% for n=1:5;
%  M(n,:)=mean(A(n*2-1:n*2,:),1);
% end
% day=[1:1:5]';
% for n=1:size(M,2);
%     [~,p]=corr(M(:,n),day);
%     if p<0.05
%         sig(n)=1;
%     else
%         sig(n)=0;
%     end
% end

%this is too strangent: not many bins are significant

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
plot([1:1:200],mean(fieldDistriEachBin(4:11,:),1)-mean(fieldDistriEachBin(1,:),1))
hold on
plot([1:1:200],tempRL*20);
hold on
plot([890/5 890/5],[0 20],'r');
xlim([0 200]);
title('Distribution each bin Day (4-11)-(1)');


subplot(413)
plot([1:1:size(fieldDistriAmpSum,2)]*5-2.5,mean(fieldDistriAmpSum(4:11,:),1)-mean(fieldDistriAmpSum(1,:),1))
hold on
plot([2.5:5:1000],tempRL*2)
hold on
plot([890 890],[0 2],'r');
xlim([-15*5 215*5]);
title('Amplitude Day (4-11)-(1)');

subplot(414)
plot([1:1:size(fieldDistriAmpNormSum,2)]*5-2.5,mean(fieldDistriAmpNormSum(4:11,:),1)-mean(fieldDistriAmpNormSum(1,:),1))
hold on
plot([2.5:5:1000],tempRL*4)
hold on
plot([890 890],[0 4],'r');
xlim([-15*5 215*5]);
title('Amplitude norm Day (4-11)-(1)');
saveas(gcf,'disriAmpOldNewEnv.fig');
saveas(gcf,'disriAmpOldNewEnv.jpg');

%% plot for the arrangement of fields: to the template all infields are 1, and nonfields are zero
fieldLocOnly={};

for n=1:length(fieldDistri);
    fieldLocOnly{n}=zeros(length(fieldDistri{1}),200);%each cell is has one row
    for m=1:length(fieldDistri{n});
        f=fieldDistri{n}{m};
        if ~isnan(f)&~isempty(f);
        fieldLocOnly{n}(m,f)=1;
        end
    end
end
save('fieldLocOnly.mat','fieldLocOnly');
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
for n=1:size(fieldLocOnly{day},1)
[~,lagsDfof(n,1)]=findCorrLag(fieldLocOnly{day}(n,:),temp,lagBins);
end

[~,orderLocTempRL1]=sort(lagsDfof);
fieldLocOnlySort_TempRL1={};
for n=1:length(fieldLocOnly);
    fieldLocOnlySort_TempRL1{n}=fieldLocOnly{n}(orderLocTempRL1,:);
end

figure
for n=1:length(fieldLocOnly);
    subplot(1,length(fieldLocOnly),n);
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

[r,p]=corr([1:1:9]',matrixCorrLocOnly(2:end));
figure,plot(matrixCorrLocOnly,'r');
hold on
plot(matrixCorrLocOnly,'r.','MarkerSize',10);
title(['Matrix correlation One to Next p=',num2str(p)]);
ylabel('Matrix correlation');
xlabel('Matrix number');

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
for n=1:size(fieldLocOnly{day},1)
[~,lagsDfof(n,1)]=findCorrLag(fieldLocOnly{day}(n,:),temp,lagBins);
end

[~,orderLocTempRL2]=sort(lagsDfof);
fieldLocOnlySort_TempRL2={};
for n=1:length(fieldLocOnly);
    fieldLocOnlySort_TempRL2{n}=fieldLocOnly{n}(orderLocTempRL2,:);
end

figure
for n=1:length(fieldLocOnly);
     subplot(1,length(fieldLocOnly),n);
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
[r,p]=corr([1:1:9]',matrixCorrLocOnly(2:end));

figure,plot(matrixCorrLocOnly,'r');
hold on
plot(matrixCorrLocOnly,'r.','MarkerSize',10);
title(['Matrix correlation One to Next p=',num2str(p)]);
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
for n=1:size(fieldLocOnly{day},1)
[~,lagsDfof(n,1)]=findCorrLag(fieldLocOnly{day}(n,:),temp,lagBins);
end

[~,orderLocTempRL11]=sort(lagsDfof);
fieldLocOnlySort_TempRL11={};
for n=1:length(fieldLocOnly);
    fieldLocOnlySort_TempRL11{n}=fieldLocOnly{n}(orderLocTempRL11,:);
end

figure
for n=1:length(fieldLocOnly);
     subplot(1,length(fieldLocOnly),n);
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
title(['Matrix correlation One to Next p=',num2str(p)]);
ylabel('Matrix correlation');
xlabel('Matrix number');

saveas(gcf,'matrixCorrLocOnlyelationLocOnly_TempRL11.fig');
save('matrixCorrLocOnlyTempRL11.mat','matrixCorrLocOnly');

%% plot the order of cells for old and new env
binWidth=22;
bins=[0:binWidth:ceil(size(dfofSigNormAllFOVs{1},1)/binWidth)*binWidth];
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
load('dfofSigNormAllFOVs.mat');
A=fieldLocOnly;
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
load('foldersAll.mat');
load('foldersAllIndvFOV.mat');
minSpacing=[];
maxWidth=[];
ratioInOther=[];
p=pwd;
for n=1:length(foldersAll);
    cd(foldersAll{n});
    load('commonCells.mat');
    load('useFolders.mat');
    s=[];
    w=[];
    r=[];
    for m=1:length(useFolders);
        cells=commonCells(:,m);
        cd(useFolders{m});
        load('PValueClassifier_KY2_6_sig\allCells.mat');
        s(1:length(cells),m)=allCells.minSpacing(cells);
        w(1:length(cells),m)=allCells.maxWidth(cells);
        %the ratioInOut in the data is in field and out fields, here we
        %want in field vs all others (out fields and non identified bins)
        for i=1:length(cells);
            inBins=allCells.inFieldBins{cells(i)};
            activity=allCells.dfofaveragesmooth(:,cells(i));
            otherBins=setdiff([1:1:length(activity)],inBins);
            if isempty(find(isnan(inBins)))&isempty(find(isnan(otherBins)));
            r(i,m)=mean(activity(inBins))/mean(activity(otherBins));
            else
                r(i,m)=nan;
            end
        end
    end
    minSpacing(end+1:end+size(s,1),:)=s;
  maxWidth(end+1:end+size(w,1),:)=w;
    ratioInOther(end+1:end+size(r,1),:)=r;
end
cd(p);

%normalize all numbers based on the data on the first day in new env
minSpacingNorm=[];
maxWidthNorm=[];
ratioInOtherNorm=[];
normDay=2;

for n=1:size(minSpacing,1);
    minSpacingNorm(n,:)=minSpacing(n,:)/minSpacing(n,normDay);
    maxWidthNorm(n,:)=maxWidth(n,:)/maxWidth(n,normDay);
ratioInOtherNorm(n,:)=ratioInOther(n,:)/ratioInOther(n,normDay);
end


save('fieldInfoUsingAllCells.mat','minSpacing','maxWidth','ratioInOther','minSpacingNorm','maxWidthNorm','ratioInOtherNorm');
%%

figure,
subplot(231);
a=minSpacing;
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
ylabel('spacing cm');
xlabel('day');

            
subplot(232);
a=maxWidth;
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
ylabel('width cm');
xlabel('day');

subplot(233);
a=ratioInOther;
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
ylabel('ratio');
xlabel('day');
% saveas(gcf,'fieldInfoUsingAllCells.fig');

subplot(234);
a=minSpacingNorm;
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
ylabel('spacing cm');
xlabel('day');

            
subplot(235);
a=maxWidthNorm;
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
ylabel('width cm');
xlabel('day');

subplot(236);
a=ratioInOtherNorm;
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
ylabel('ratio');
xlabel('day');

%
saveas(gcf,'fieldInfoUsingAllCells.fig');
% saveas(gcf,'fieldInfoUsingAllCells.jpg');

%% look at the cell's max field width and spacing and in out field ratio
load('foldersAll.mat');
load('foldersAllIndvFOV.mat');
minSpacing=[];
maxWidth=[];
ratioInOther=[];
p=pwd;
for n=1:length(foldersAll);
    cd(foldersAll{n});
    load('commonCells.mat');
    load('useFolders.mat');
    s=[];
    w=[];
    r=[];
    for m=1:length(useFolders);
        cells=commonCells(:,m);
        cd(useFolders{m});
        load('PValueClassifier_KY2_6_sig\allCellsCorrected.mat');
        s(1:length(cells),m)=allCellsCorrected.minSpacing(cells);
        w(1:length(cells),m)=allCellsCorrected.maxWidth(cells);
        %the ratioInOut in the data is in field and out fields, here we
        %want in field vs all others (out fields and non identified bins)
        for i=1:length(cells);
            inBins=allCellsCorrected.inFieldBins{cells(i)};
            activity=allCellsCorrected.dfofaveragesmooth(:,cells(i));
            otherBins=setdiff([1:1:length(activity)],inBins);
            if isempty(find(isnan(inBins)))&isempty(find(isnan(otherBins)));
            r(i,m)=mean(activity(inBins))/mean(activity(otherBins));
            else
                r(i,m)=nan;
            end
        end
    end
    minSpacing(end+1:end+size(s,1),:)=s;
  maxWidth(end+1:end+size(w,1),:)=w;
    ratioInOther(end+1:end+size(r,1),:)=r;
end
cd(p);
%normalize all numbers based on the data on the first day in new env
minSpacingNorm=[];
maxWidthNorm=[];
ratioInOtherNorm=[];
normDay=2;

for n=1:size(minSpacing,1);
    minSpacingNorm(n,:)=minSpacing(n,:)/minSpacing(n,normDay);
    maxWidthNorm(n,:)=maxWidth(n,:)/maxWidth(n,normDay);
ratioInOtherNorm(n,:)=ratioInOther(n,:)/ratioInOther(n,normDay);
end


save('fieldInfoUsingAllCellsCorrected.mat','minSpacing','maxWidth','ratioInOther','minSpacingNorm','maxWidthNorm','ratioInOtherNorm');
%%

figure,
subplot(231);
a=minSpacing;
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
ylabel('spacing cm');
xlabel('day');

            
subplot(232);
a=maxWidth;
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
ylabel('width cm');
xlabel('day');

subplot(233);
a=ratioInOther;
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
ylabel('ratio');
xlabel('day');
% saveas(gcf,'fieldInfoUsingAllCells.fig');

subplot(234);
a=minSpacingNorm;
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
ylabel('spacing cm');
xlabel('day');


            
subplot(235);
a=maxWidthNorm;
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
ylabel('width cm');
xlabel('day');

subplot(236);
a=ratioInOtherNorm;
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
ylabel('ratio');
xlabel('day');

%
saveas(gcf,'fieldInfoUsingAllCellsCorrected.fig');
saveas(gcf,'fieldInfoUsingAllCellsCorrected.jpg');


%% look at the cell's max field width and spacing and in out field ratio
load('foldersAll.mat');
load('foldersAllIndvFOV.mat');
meanSpacing=[];
meanWidth=[];
ratioOtherIn=[];
p=pwd;
for n=1:length(foldersAll);
    cd(foldersAll{n});
    load('commonCells.mat');
    load('useFolders.mat');
    s=[];
    w=[];
    r=[];
    for m=1:length(useFolders);
        cells=commonCells(:,m);
        cd(useFolders{m});
        load('PValueClassifier_KY2_6_sig\allCells.mat');
        
        mS=[];%mean spacing
        mW=[];%mean width
        for i=1:length(allCells.fieldSpacings);
            mS(i)=mean(allCells.fieldSpacings{i});
            mW(i)=mean(allCells.fieldWidths{i});
        end
        
        
        s(1:length(cells),m)=mS(cells);
        w(1:length(cells),m)=mW(cells);
        %the ratioInOut in the data is in field and out fields, here we
        %want in field vs all others (out fields and non identified bins)
        for i=1:length(cells);
            inBins=allCells.inFieldBins{cells(i)};
            activity=allCells.dfofaveragesmooth(:,cells(i));
            otherBins=setdiff([1:1:length(activity)],inBins);
            if isempty(find(isnan(inBins)))&isempty(find(isnan(otherBins)));
            r(i,m)=mean(activity(otherBins))/mean(activity(inBins));
            else
                r(i,m)=nan;
            end
        end
    end
    meanSpacing(end+1:end+size(s,1),:)=s;
  meanWidth(end+1:end+size(w,1),:)=w;
    ratioOtherIn(end+1:end+size(r,1),:)=r;
end
cd(p);

%normalize all numbers based on the data on the first day in new env
meanSpacingNorm=[];
meanWidthNorm=[];
ratioOtherInNorm=[];
normDay=2;

for n=1:size(meanSpacing,1);
    meanSpacingNorm(n,:)=meanSpacing(n,:)/meanSpacing(n,normDay);
    meanWidthNorm(n,:)=meanWidth(n,:)/meanWidth(n,normDay);
ratioOtherInNorm(n,:)=ratioOtherIn(n,:)/ratioOtherIn(n,normDay);
end


save('fieldInfoUsingAllCellsMean.mat','meanSpacing','meanWidth','ratioOtherIn','meanSpacingNorm','meanWidthNorm','ratioOtherInNorm');
%%

figure,
subplot(231);
a=meanSpacing;
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
ylabel('mean spacing cm');
xlabel('day');

            
subplot(232);
a=meanWidth;
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
ylabel('width cm');
xlabel('day');

subplot(233);
a=ratioOtherIn;
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
ylabel('ratio');
xlabel('day');
% saveas(gcf,'fieldInfoUsingAllCells.fig');

subplot(234);
a=meanSpacingNorm;
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
ylabel('spacing cm');
xlabel('day');

            
subplot(235);
a=meanWidthNorm;
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
ylabel('width cm');
xlabel('day');

subplot(236);
a=ratioOtherInNorm;
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
ylabel('ratio');
xlabel('day');

%
saveas(gcf,'fieldInfoUsingAllCellsMean.fig');
% saveas(gcf,'fieldInfoUsingAllCells.jpg');

%% look at the cell's max field width and spacing and in out field ratio
load('foldersAll.mat');
load('foldersAllIndvFOV.mat');
meanSpacing=[];
meanWidth=[];
ratioOtherIn=[];
p=pwd;
for n=1:length(foldersAll);
    cd(foldersAll{n});
    load('commonCells.mat');
    load('useFolders.mat');
    s=[];
    w=[];
    r=[];
    for m=1:length(useFolders);
        cells=commonCells(:,m);
        cd(useFolders{m});
        load('PValueClassifier_KY2_6_sig\allCellsCorrected.mat');
        
  mS=[];%mean spacing
        mW=[];%mean width
        for i=1:length(allCellsCorrected.fieldSpacings);
            mS(i)=mean(allCellsCorrected.fieldSpacings{i});
            mW(i)=mean(allCellsCorrected.fieldWidths{i});
        end
               
        s(1:length(cells),m)=mS(cells);
        w(1:length(cells),m)=mW(cells);
        %the ratioInOut in the data is in field and out fields, here we
        %want in field vs all others (out fields and non identified bins)
        for i=1:length(cells);
            inBins=allCellsCorrected.inFieldBins{cells(i)};
            activity=allCellsCorrected.dfofaveragesmooth(:,cells(i));
            otherBins=setdiff([1:1:length(activity)],inBins);
            if isempty(find(isnan(inBins)))&isempty(find(isnan(otherBins)));
            r(i,m)=mean(activity(otherBins))/mean(activity(inBins));
            else
                r(i,m)=nan;
            end
        end
    end
    meanSpacing(end+1:end+size(s,1),:)=s;
  meanWidth(end+1:end+size(w,1),:)=w;
    ratioOtherIn(end+1:end+size(r,1),:)=r;
end
cd(p);
%normalize all numbers based on the data on the first day in new env
meanSpacingNorm=[];
meanWidthNorm=[];
ratioOtherInNorm=[];
normDay=2;

for n=1:size(meanSpacing,1);
    meanSpacingNorm(n,:)=meanSpacing(n,:)/meanSpacing(n,normDay);
    meanWidthNorm(n,:)=meanWidth(n,:)/meanWidth(n,normDay);
ratioOtherInNorm(n,:)=ratioOtherIn(n,:)/ratioOtherIn(n,normDay);
end


save('fieldInfoUsingAllCellsCorrectedMean.mat','meanSpacing','meanWidth','ratioOtherIn','meanSpacingNorm','meanWidthNorm','ratioOtherInNorm');
%%

figure,
subplot(231);
a=meanSpacing;
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
ylabel('spacing cm');
xlabel('day');

            
subplot(232);
a=meanWidth;
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
ylabel('width cm');
xlabel('day');

subplot(233);
a=ratioOtherIn;
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
ylabel('ratio');
xlabel('day');
% saveas(gcf,'fieldInfoUsingAllCells.fig');

subplot(234);
a=meanSpacingNorm;
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
ylabel('spacing cm');
xlabel('day');


            
subplot(235);
a=meanWidthNorm;
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
ylabel('width cm');
xlabel('day');

subplot(236);
a=ratioOtherInNorm;
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
title(['ratioOtherVsIn norm p=',num2str(p)])
ylabel('ratio');
xlabel('day');

%
saveas(gcf,'fieldInfoUsingAllCellsCorrectedMean.fig');
saveas(gcf,'fieldInfoUsingAllCellsCorrectedMean.jpg');

%% use real in and out field ratio

load('foldersAll.mat');
load('foldersAllIndvFOV.mat');
ratioInOut=[];
ratioInOutCorrected=[];
p=pwd;

for n=1:length(foldersAll);
    cd(foldersAll{n});
    load('commonCells.mat');
    load('useFolders.mat');
  
    r=[];
    rc=[];
    for m=1:length(useFolders);
        cells=commonCells(:,m);
        cd(useFolders{m});
        load('PValueClassifier_KY2_6_sig\allCells.mat');
        load('PValueClassifier_KY2_6_sig\allCellsCorrected.mat');
  r(1:length(cells),m)=allCells.ratioInOut(cells);
  rc(1:length(cells),m)=allCellsCorrected.ratioInOut(cells);      
    end

    ratioInOut(end+1:end+size(r,1),:)=r;
    ratioInOutCorrected(end+1:end+size(r,1),:)=rc;
end
cd(p);

ratioInOutNorm=[];
ratioInOutCorrectedNorm=[];
normDay=2;
for n=1:size(ratioInOut,1);
   
ratioInOutNorm(n,:)=ratioInOut(n,:)/ratioInOut(n,normDay);
ratioInOutCorrectedNorm(n,:)=ratioInOutCorrected(n,:)/ratioInOutCorrected(n,normDay);
end

save('fieldRatioInOut.mat','ratioInOut','ratioInOutCorrected','ratioInOutNorm','ratioInOutCorrectedNorm');
%%
figure,


subplot(221);
a=ratioInOut;
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
a=ratioInOutCorrected;
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
a=ratioInOutNorm;
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
a=ratioInOutCorrectedNorm;
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

%% mean activity in and out
load('foldersAll.mat');
load('foldersAllIndvFOV.mat');

dfofFieldIn=[];
dfofFieldOther=[];
p=pwd;


for n=1:length(foldersAll);
    cd(foldersAll{n});
    load('commonCells.mat');
    load('useFolders.mat');
  
    ifield=[];
    nfield=[];
    for m=1:length(useFolders);
        cells=commonCells(:,m);
        cd(useFolders{m});
        load('PValueClassifier_KY2_6_sig\allCells.mat');
        
     for i=1:length(cells);
         inBins=allCells.inFieldBins{cells(i)};
            activity=allCells.dfofaveragesmooth(:,cells(i));
            otherBins=setdiff([1:1:length(activity)],inBins);
             if isempty(find(isnan(inBins)))&isempty(find(isnan(otherBins)));
            nfield(i,m)=nanmean(activity(otherBins));
            ifield(i,m)=nanmean(activity(inBins));
             else
               nfield(i,m)=nan;  
                ifield(i,m)=nan;  
             end
     end
            
    end

    dfofFieldIn(end+1:end+size(ifield,1),:)=ifield;
    dfofFieldOther(end+1:end+size(nfield,1),:)=nfield;
end
cd(p);

dfofFieldInNorm=[];
dfofFieldOtherNorm=[];
normDay=2;

for n=1:size(dfofFieldIn,1);
    dfofFieldInNorm(n,:)=dfofFieldIn(n,:)/dfofFieldIn(n,normDay);
    dfofFieldOtherNorm(n,:)=dfofFieldOther(n,:)/dfofFieldOther(n,normDay);
end

save('meanFFieldInOthers.mat','dfofFieldIn','dfofFieldOther','dfofFieldInNorm','dfofFieldOtherNorm');

figure,
subplot(121);
a=dfofFieldIn;
M=nanmean(a,1);
S=nansem(a,1);

errorbar([1:1:length(M)],M,S,'r');
a=dfofFieldOther;
M=nanmean(a,1);
S=nansem(a,1);

hold on
errorbar([1:1:length(M)],M,S,'k');

% [r,p]=corr([1:1:length(M)-1]',M(2:end)');
title(['in out field dfof mean',num2str(p)])
ylabel('dfof mean');
xlabel('day');

            
subplot(122);
a=dfofFieldInNorm;
M=nanmean(a,1);
S=nansem(a,1);

errorbar([1:1:length(M)],M,S,'r');
a=dfofFieldOtherNorm;
M=nanmean(a,1);
S=nansem(a,1);

hold on
errorbar([1:1:length(M)],M,S,'k');

% [r,p]=corr([1:1:length(M)-1]',M(2:end)');
title(['in out field dfof mean norm',num2str(p)])
ylabel('dfof mean norm');
xlabel('day');

saveas(gcf,'meanFFieldInOthers.fig');


%% run by run consistency with run by run behavior
p=pwd;
load('foldersAll.mat');

 allLickRunAllCellNoOldAll=[];
 allNoLickRunAllCellNoOldAll=[];
  allWrongLickRunAllCellNoOldAll=[];
  allSlowRunAllCellNoOldAll=[];
 allNoSlowRunAllCellNoOldAll=[];

for n=1:length(foldersAll);
    cd(foldersAll{n});
    cd('corrBehaviorRunByRun')
     load('allLickRunAllCellNoOld.mat');
 load('allNoLickRunAllCellNoOld.mat');
   load('allWrongLickRunAllCellNoOld.mat');
  load('allSlowRunAllCellNoOld.mat');
 load('allNoSlowRunAllCellNoOld.mat');
 
  A=allLickRunAllCellNoOld;
B=allNoLickRunAllCellNoOld;
C=allWrongLickRunAllCellNoOld;
D=allSlowRunAllCellNoOld;
E=allNoSlowRunAllCellNoOld;

   allLickRunAllCellNoOldAll(end+1:end+length(A),1)=A;
     allNoLickRunAllCellNoOldAll(end+1:end+length(B),1)=B;
      allWrongLickRunAllCellNoOldAll(end+1:end+length(C),1)=C;
           allSlowRunAllCellNoOldAll(end+1:end+length(D),1)=D;
           allNoSlowRunAllCellNoOldAll(end+1:end+length(E),1)=E;
end

cd(p)

 A=allLickRunAllCellNoOldAll;
B=allNoLickRunAllCellNoOldAll;
C=allWrongLickRunAllCellNoOldAll;
D=allSlowRunAllCellNoOldAll;
E=allNoSlowRunAllCellNoOldAll;

 semA=std(A)/sqrt(length(A));
     meanA=mean(A);
     semB=std(B)/sqrt(length(B));
     meanB=mean(B);
      semC=std(C)/sqrt(length(C));
 meanC=mean(C);
  semD=std(D)/sqrt(length(D));
     meanD=mean(D);
        semE=std(E)/sqrt(length(E));
     meanE=mean(E);
    [~,p1]=ttest2(A,B);
     [~,p2]=ttest2(A,C);
      [~,p3]=ttest2(D,E);
     
     figure  
 name={'L';'NL';'WL';'S';'NS'};
     bar([1 2 3 4 5],[meanA meanB meanC meanD meanE]);
     set(gca,'xticklabel',name);
     hold on
     errorbar([1 2 3 4 5],[meanA meanB meanC meanD meanE],[semA semB semC semD semE],'.');
     title(['p1=',num2str(p1),' p2=',num2str(p2),' p3=',num2str(p3)]);
     
     saveas(gcf,'corr_allRunAllCellsLickSlow.fig');
save('allLickRunAllCellNoOldAll.mat','allLickRunAllCellNoOldAll');
save('allNoLickRunAllCellNoOldAll.mat','allNoLickRunAllCellNoOldAll');
save('allWrongLickRunAllCellNoOldAll.mat','allWrongLickRunAllCellNoOldAll');
save('allSlowRunAllCellNoOldAll.mat','allSlowRunAllCellNoOldAll');
save('allNoSlowRunAllCellNoOldAll.mat','allNoSlowRunAllCellNoOldAll');

     %% focus on learn2 (slow lick both) learn1 (slow or lick) and learn0 (on slow no lick): no old env
     
p=pwd;
load('foldersAll.mat');

 allRunAllCellLearn2=[];
 allRunAllCellLearn1=[];
  allRunAllCellLearn0=[];


for n=1:length(foldersAll);

    cd(foldersAll{n});
    cd('corrBehaviorRunByRun')
     load('allRunAllCellLearn2NoOld.mat');
 load('allRunAllCellLearn1NoOld.mat');
   load('allRunAllCellLearn0NoOld.mat');

  A=allRunAllCellLearn2NoOld;
B=allRunAllCellLearn1NoOld;
C=allRunAllCellLearn0NoOld;


   allRunAllCellLearn2(end+1:end+length(A),1)=A;
     allRunAllCellLearn1(end+1:end+length(B),1)=B;
      allRunAllCellLearn0(end+1:end+length(C),1)=C;
   
end

cd(p)

 A=allRunAllCellLearn2;
B=allRunAllCellLearn1;
C=allRunAllCellLearn0;


 semA=std(A)/sqrt(length(A));
     meanA=mean(A);
     semB=std(B)/sqrt(length(B));
     meanB=mean(B);
      semC=std(C)/sqrt(length(C));
 meanC=mean(C);

    [~,p1]=ttest2(A,B);
     [~,p2]=ttest2(A,C);
    
     
    figure  
      name={'L2';'L1';'L0'};
  bar([1 2 3],[meanA meanB meanC]);
     set(gca,'xticklabel',name);
     hold on
     errorbar([1 2 3],[meanA meanB meanC],[semA semB semC],'.');
     title([' p1=',num2str(p1),' p2=',num2str(p2)]);
 saveas(gcf,'corr_allRunAllCellsLearn.fig');
     
save('allRunAllCellLearn2.mat','allRunAllCellLearn2');
save('allRunAllCellLearn1.mat','allRunAllCellLearn1');
save('allRunAllCellLearn0.mat','allRunAllCellLearn0');


%% plot them as a funciton of days
p=pwd;
load('foldersAll.mat');

 allCorrCellLearn2Days={};
 allCorrCellLearn1Days={};
  allCorrCellLearn12Days={};
  allCorrCellLearn0Days={};
   for n=1:11
     allCorrCellLearn2Days{n}=[];
     allCorrCellLearn1Days{n}=[];
     allCorrCellLearn12Days{n}=[];
     allCorrCellLearn0Days{n}=[];
   end
     
  for n=1:length(foldersAll);
    cd(foldersAll{n});
    cd('corrBehaviorRunByRun')
     load('corrToMeanRunLearn0.mat');
 load('corrToMeanRunLearn1.mat');
   load('corrToMeanRunLearn2.mat');
   
for m=1:11;
    A=mean(corrToMeanRunLearn2{m},1,'omitnan');
  B=mean(corrToMeanRunLearn1{m},1,'omitnan');
  C=mean([corrToMeanRunLearn2{m};corrToMeanRunLearn1{m}],1,'omitnan');
    D=mean(corrToMeanRunLearn0{m},1,'omitnan');
    if ~isnan(A);
    allCorrCellLearn2Days{m}(end+1:end+length(A))=A;
    end
        if ~isnan(B);
    allCorrCellLearn1Days{m}(end+1:end+length(B))=B;
        end
         if ~isnan(C);
    allCorrCellLearn12Days{m}(end+1:end+length(C))=C;
         end
         if ~isnan(D)
allCorrCellLearn0Days{m}(end+1:end+length(D))=D;
         end
end
  end
  
     
  
  cd(p);
  
  %%
  figure,
  
  subplot(121);
  M=[];
  S=[];
  D=allCorrCellLearn2Days;
  
  for n=1:length(D);
      M(n)=mean(D{n});
      S(n)=std(D{n})/sqrt(length(D{n}));
  end
  hold on
  errorbar([1:1:length(D)],M,S);
  
    M=[];
  S=[];
  D=allCorrCellLearn1Days;
  
  for n=1:length(D);
      M(n)=mean(D{n});
      S(n)=std(D{n})/sqrt(length(D{n}));
  end
  hold on
  errorbar([1:1:length(D)],M,S);
  
  M=[];
  S=[];
  D=allCorrCellLearn12Days;
  
  for n=1:length(D);
      M(n)=mean(D{n});
      S(n)=std(D{n})/sqrt(length(D{n}));
  end
  hold on
  errorbar([1:1:length(D)],M,S);
  

      M=[];
  S=[];
  D=allCorrCellLearn0Days;
  
  for n=1:length(D);
      M(n)=mean(D{n});
      S(n)=std(D{n})/sqrt(length(D{n}));
  end
  hold on
  errorbar([1:1:length(D)],M,S);
  
  legend('learn2','learn1','learn12','learn0','Location','Southwest');
  title('all learning conditions')
  
    subplot(122);

  M=[];
  S=[];
  D=allCorrCellLearn12Days;
  
  for n=1:length(D);
      M(n)=mean(D{n});
      S(n)=std(D{n})/sqrt(length(D{n}));
  end
  hold on
  errorbar([1:1:length(D)],M,S);
  [~,p1]=corr([1:1:length(M)-1]',M(2:end)');

      M=[];
  S=[];
  D=allCorrCellLearn0Days;
  
  for n=1:length(D);
      M(n)=mean(D{n});
      S(n)=std(D{n})/sqrt(length(D{n}));
  end
  hold on
  errorbar([1:1:length(D)],M,S);
    [~,p2]=corr([1:1:length(M)-1]',M(2:end)');

  sig=[];
  p=[];
  for n=1:11;
      [~,p(n)]=ttest2(allCorrCellLearn12Days{n},allCorrCellLearn0Days{n});
      
      if p(n)<0.05;
          sig(n)=1;
          
          hold on
          plot(n,mean(allCorrCellLearn12Days{n})*1.03,'r*');
      else 
          sig(n)=0;
      end
  end
  
  
  legend('learn12','learn0','Location','Southwest');
    title(['learn vs nolearn p1=',num2str(p1),' p2=',num2str(p2)])
    saveas(gcf,'corr_allRunAllCellsLearnWithDays.fig');
    
%% look at the correlation or learn and no learn runs within themselves
p=pwd;
load('foldersAll.mat');

 allCorrCellLearn12DaysToOwn={};%this is learn (including 1 and 2) correlation on their own (one run to the rest)
  allCorrCellLearn0DaysToOwn={};%this is no learn (including 0) correlation on their own (one run to the rest)
   allCorrCellLearnRand12DaysToOwn={};%random number of learn days (same number as no learn runs)
   
   allCorrCellLearn12DaysToOwnMean={};%this is learn (including 1 and 2) correlation on their own (one run to the rest)
  allCorrCellLearn0DaysToOwnMean={};%this is no learn (including 0) correlation on their own (one run to the rest)
allCorrCellLearnRand12DaysToOwnMean={};%random number of learn days (same number as no learn runs)
   
  
   for n=1:11
     allCorrCellLearn12DaysToOwn{n}=[];
     allCorrCellLearn0DaysToOwn{n}=[];
     allCorrCellLearnRand12DaysToOwn{n}=[];
        allCorrCellLearn12DaysToOwnMean{n}=[];
     allCorrCellLearn0DaysToOwnMean{n}=[]; 
     allCorrCellLearnRand12DaysToOwnMean{n}=[];
   end
     
  for n=1:length(foldersAll);
    cd(foldersAll{n});
    cd('corrBehaviorRunByRun')
load('learnIndivCorr.mat');
load('noLearnIndivCorr.mat');
load('learnIndivCorrRand.mat');

load('learnCorrToMean.mat');
load('noLearnCorrToMean.mat');
load('learnCorrToMeanRand.mat');

for m=1:11;
      A=mean(learnIndivCorr{m},1,'omitnan');
    B=mean(noLearnIndivCorr{m},1,'omitnan');
     C=mean(learnIndivCorrRand{m},1,'omitnan');
      D=mean(learnCorrToMean{m},1,'omitnan');
    E=mean(noLearnCorrToMean{m},1,'omitnan');
     F=mean(learnCorrToMeanRand{m},1,'omitnan');
   allCorrCellLearn12DaysToOwn{m}(end+1:end+length(A))=A;
  allCorrCellLearn0DaysToOwn{m}(end+1:end+length(B))=B;
  allCorrCellLearnRand12DaysToOwn{m}(end+1:end+length(C))=C;
    allCorrCellLearn12DaysToOwnMean{m}(end+1:end+length(D))=D;
  allCorrCellLearn0DaysToOwnMean{m}(end+1:end+length(E))=E;
  allCorrCellLearnRand12DaysToOwnMean{m}(end+1:end+length(F))=F;
end
  end
  cd(p);
  
  save('allCorrCellLearn12DaysToOwn.mat','allCorrCellLearn12DaysToOwn');
    save('allCorrCellLearn0DaysToOwn.mat','allCorrCellLearn0DaysToOwn');
    save('allCorrCellLearnRand12DaysToOwn.mat','allCorrCellLearnRand12DaysToOwn');
    
  save('allCorrCellLearn12DaysToOwnMean.mat','allCorrCellLearn12DaysToOwnMean');
    save('allCorrCellLearn0DaysToOwnMean.mat','allCorrCellLearn0DaysToOwnMean');
    save('allCorrCellLearnRand12DaysToOwnMean.mat','allCorrCellLearnRand12DaysToOwnMean');
    
  figure
  subplot(121)
    M=[];
  S=[];
  D=allCorrCellLearn12DaysToOwn;
  
  for n=1:length(D);
      M(n)=mean(D{n},'omitnan');
      S(n)=nanstd(D{n})/sqrt(length(D{n}));
  end
  hold on
  errorbar([1:1:length(D)],M,S);
  [~,p1]=corr([1:1:length(M)-1]',M(2:end)');

      M=[];
  S=[];
  D=allCorrCellLearn0DaysToOwn;
  
  for n=1:length(D);
      M(n)=mean(D{n},'omitnan');
      S(n)=nanstd(D{n})/sqrt(length(D{n}));
  end
  hold on
  errorbar([1:1:length(D)],M,S);
    [~,p2]=corr([1:1:length(M)-1]',M(2:end)');
    
     M=[];
  S=[];
  D=allCorrCellLearnRand12DaysToOwn;
  
  for n=1:length(D);
      M(n)=mean(D{n},'omitnan');
      S(n)=nanstd(D{n})/sqrt(length(D{n}));
  end
  hold on
  errorbar([1:1:length(D)],M,S);
    [~,p3]=corr([1:1:length(M)-1]',M(2:end)');
    

  sig=[];
  p=[];
  for n=1:11;
      [~,p(n)]=ttest2(allCorrCellLearn12DaysToOwn{n},allCorrCellLearn0DaysToOwn{n});
      
      if p(n)<0.05;
          sig(n)=1;
          
          hold on
          plot(n,mean(allCorrCellLearn12DaysToOwn{n},'omitnan')*1.05,'r*');
      else 
          sig(n)=0;
      end
  end
  
  ylim([0.24 0.4])
  legend('learn12','learn0','Location','Southeast');
    title(['learn vs nolearn indiv corr p1=',num2str(p1),' p2=',num2str(p2),' p3=',num2str(p3)])

  
    
    subplot(122)
    M=[];
  S=[];
  D=allCorrCellLearn12DaysToOwnMean;
  
  for n=1:length(D);
      M(n)=mean(D{n},'omitnan');
      S(n)=nanstd(D{n})/sqrt(length(D{n}));
  end
  hold on
  errorbar([1:1:length(D)],M,S);
  [~,p1]=corr([1:1:length(M)-1]',M(2:end)');

      M=[];
  S=[];
  D=allCorrCellLearn0DaysToOwnMean;
  
  for n=1:length(D);
      M(n)=mean(D{n},'omitnan');
      S(n)=nanstd(D{n})/sqrt(length(D{n}));
  end
  hold on
  errorbar([1:1:length(D)],M,S);
    [~,p2]=corr([1:1:length(M)-1]',M(2:end)');
    
    M=[];
  S=[];
  D=allCorrCellLearnRand12DaysToOwnMean;
  
  for n=1:length(D);
      M(n)=mean(D{n},'omitnan');
      S(n)=nanstd(D{n})/sqrt(length(D{n}));
  end
  hold on
  errorbar([1:1:length(D)],M,S);
    [~,p3]=corr([1:1:length(M)-1]',M(2:end)');


  sig=[];
  p=[];
  for n=1:11;
      [~,p(n)]=ttest2(allCorrCellLearn12DaysToOwnMean{n},allCorrCellLearn0DaysToOwnMean{n});
      
      if p(n)<0.05;
          sig(n)=1;
          
          hold on
          plot(n,mean(allCorrCellLearn12DaysToOwnMean{n},'omitnan')*1.05,'r*');
      else 
          sig(n)=0;
      end
  end
  

  legend('learn12','learn0','Location','Southeast');
    title(['learn vs nolearn corr mean p1=',num2str(p1),' p2=',num2str(p2),' p3=',num2str(p3)])
        saveas(gcf,'corr_allRunAllCellsLearnWithDays_toOwnRuns.fig');%in this figure, the left plot is for each run correlate to all other runs in their own category: learn runs correlate to other learn runs; unlearn runs clrrelate to other unlearn runs.
        %the second subplot is the runs in their own category correlate wo
        %their mean activity of that category. So for unlearned runs, since
        %there are less runs to make that mean activity, the correlation is
        %higher but this doesn't mean the unlearned runs correlate better.
        %I think it's better to use the correlation to other runs (subplot
        %121)

       
  
%% dfof_sig: see whether there is an incrased fluorescence during learning
load('foldersAll.mat');
dfofSigAllFOVs{1}=[];
dfofSigAllFOVs{2}=[];
dfofSigAllFOVs{3}=[];
dfofSigAllFOVs{4}=[];
dfofSigAllFOVs{5}=[];
dfofSigAllFOVs{6}=[];
dfofSigAllFOVs{7}=[];
dfofSigAllFOVs{8}=[];
dfofSigAllFOVs{9}=[];
dfofSigAllFOVs{10}=[];
dfofSigAllFOVs{11}=[];

p=pwd;

for n=1:length(foldersAll);
  cd(foldersAll{n});
  disp(n)
  
cd ..\
load('allDfof_sig.mat');
   for m=1:length(dfofSigAllFOVs);
      dfofSigAllFOVs{m}(end+1:end+size(allDfof_sig{m},1),:)=allDfof_sig{m};
   end
end

cd(p);

save('dfofSigAllFOVs.mat','dfofSigAllFOVs');

M=[]; %mean dfofSig of each cell each day

for n=1:length(dfofSigAllFOVs);
    for m=1:size(dfofSigAllFOVs{1},1);
        M(m,n)=mean(dfofSigAllFOVs{n}(m,:));
    end
end

%normalizd it based on the first day in new env
MNorm=[];
for n=1:size(M,1);
    MNorm(n,:)=M(n,:)/M(n,2);
end

meanMNorm=[];
semMNorm=[];
for n=1:size(MNorm,2);
    meanMNorm(1,n)=mean(MNorm(:,n));
    semMNorm(1,n)=std(MNorm(:,n))/sqrt(size(MNorm,1));
end
  

P=[]; %mean dfofSig of each cell each day

for n=1:length(dfofSigAllFOVs);
    for m=1:size(dfofSigAllFOVs{1},1);
        P(m,n)=max(dfofSigAllFOVs{n}(m,:));
    end
end

%normalizd it based on the first day in new env
PNorm=[];
for n=1:size(P,1);
    PNorm(n,:)=P(n,:)/P(n,2);
end

meanPNorm=[];
semPNorm=[];
for n=1:size(PNorm,2);
    meanPNorm(1,n)=mean(PNorm(:,n));
    semPNorm(1,n)=std(PNorm(:,n))/sqrt(size(PNorm,1));
end

figure,

subplot(221)
errorbar([1:1:length(meanMNorm)],meanMNorm,semMNorm)
[~,p]=corr([1:1:length(meanMNorm)-1]',meanMNorm(2:end)');
title(['mean dfofSig, p=',num2str(p)]);

subplot(222)
errorbar([1:1:length(meanPNorm)],meanPNorm,semPNorm)
[~,p]=corr([1:1:length(meanPNorm)-1]',meanPNorm(2:end)');
title(['peak dfofSig, p=',num2str(p)]);

%% dfof: see whether there is an incrased fluorescence during learning
load('foldersAll.mat');
dfofAllFOVs{1}=[];
dfofAllFOVs{2}=[];
dfofAllFOVs{3}=[];
dfofAllFOVs{4}=[];
dfofAllFOVs{5}=[];
dfofAllFOVs{6}=[];
dfofAllFOVs{7}=[];
dfofAllFOVs{8}=[];
dfofAllFOVs{9}=[];
dfofAllFOVs{10}=[];
dfofAllFOVs{11}=[];

p=pwd;

for n=1:length(foldersAll);
  cd(foldersAll{n});
  disp(n)
  
cd ..\
load('allDfof.mat');
   for m=1:length(dfofAllFOVs);
      dfofAllFOVs{m}(end+1:end+size(allDfof{m},1),:)=allDfof{m};
   end
end

cd(p);

save('dfofAllFOVs.mat','dfofAllFOVs');

M=[]; %mean dfofSig of each cell each day

for n=1:length(dfofAllFOVs);
    for m=1:size(dfofAllFOVs{1},1);
        M(m,n)=mean(dfofAllFOVs{n}(m,:));
    end
end

%normalizd it based on the first day in new env
MNorm=[];
for n=1:size(M,1);
    MNorm(n,:)=M(n,:)/M(n,2);
end

meanMNorm=[];
semMNorm=[];
for n=1:size(MNorm,2);
    meanMNorm(1,n)=mean(MNorm(:,n));
    semMNorm(1,n)=std(MNorm(:,n))/sqrt(size(MNorm,1));
end
  

P=[]; %mean dfofSig of each cell each day

for n=1:length(dfofAllFOVs);
    for m=1:size(dfofAllFOVs{1},1);
        P(m,n)=max(dfofAllFOVs{n}(m,:));
    end
end

%normalizd it based on the first day in new env
PNorm=[];
for n=1:size(P,1);
    PNorm(n,:)=P(n,:)/P(n,2);
end

meanPNorm=[];
semPNorm=[];
for n=1:size(PNorm,2);
    meanPNorm(1,n)=mean(PNorm(:,n));
    semPNorm(1,n)=std(PNorm(:,n))/sqrt(size(PNorm,1));
end


subplot(223)
errorbar([1:1:length(meanMNorm)],meanMNorm,semMNorm)
[~,p]=corr([1:1:length(meanMNorm)-1]',meanMNorm(2:end)');
title(['mean dfof, p=',num2str(p)]);

subplot(224)
errorbar([1:1:length(meanPNorm)],meanPNorm,semPNorm)
[~,p]=corr([1:1:length(meanPNorm)-1]',meanPNorm(2:end)');
title(['peak dfof, p=',num2str(p)]);

saveas(gcf,'dfofChangeOverDays.fig');

%% use raw dfof or raw dfof_sig
meanDfofAllFOVs=[];
meanDfof_sigAllFOVs=[];
maxDfofAllFOVs=[];
maxDfof_sigAllFOVs=[];
p=pwd;
for n=1:length(foldersAll);
    disp(n)
    cd(foldersAll{n});
    load('commonCells.mat');
    load('useFolders');
    meanDfof=[];%each column is one day, each row is one cell
meanDfof_sig=[];
maxDfof=[];
maxDfof_sig=[];
for m=1:length(useFolders);
    cd(useFolders{m});
    d=dir('dfof_*');
    for i=1:length(d);
        load(d(i).name);
    end
    cells=commonCells(:,m);
    D=dfof(:,cells)';
      S=dfof_sig(:,cells)';
      
      meanDfof(:,m)=mean(D,2);
      meanDfof_sig(:,m)=mean(S,2);
      maxDfof(:,m)=max(D,[],2);
      maxDfof_sig(:,m)=max(S,[],2);
end

cd(foldersAll{n});
save('meanDfof.mat','meanDfof');
save('meanDfof_sig.mat','meanDfof_sig');  
save('maxDfof.mat','maxDfof');
save('maxDfof_sig.mat','maxDfof_sig');   
    
meanDfofAllFOVs(end+1:end+size(meanDfof,1),:)=meanDfof;
meanDfof_sigAllFOVs(end+1:end+size(meanDfof_sig,1),:)=meanDfof_sig;
maxDfofAllFOVs(end+1:end+size(maxDfof,1),:)=maxDfof;
maxDfof_sigAllFOVs(end+1:end+size(maxDfof_sig,1),:)=maxDfof_sig;
end
    
cd(p);

meanDfofAllFOVsNorm=[];
for n=1:size(meanDfofAllFOVs,1);
    meanDfofAllFOVsNorm(n,:)=meanDfofAllFOVs(n,:)/meanDfofAllFOVs(n,2);
end

meanmeanDfofAllFOVsNorm=[];
semmeanDfofAllFOVsNorm=[];
for n=1:size(meanDfofAllFOVsNorm,2);
    meanmeanDfofAllFOVsNorm(1,n)=mean(meanDfofAllFOVsNorm(:,n));
    semmeanDfofAllFOVsNorm(1,n)=std(meanDfofAllFOVsNorm(:,n))/sqrt(size(meanDfofAllFOVsNorm,1));
end

figure
subplot(221)
errorbar([1:1:length(meanmeanDfofAllFOVsNorm)],meanmeanDfofAllFOVsNorm,semmeanDfofAllFOVsNorm)
[~,p]=corr([1:1:length(meanmeanDfofAllFOVsNorm)-1]',meanmeanDfofAllFOVsNorm(2:end)');
title(['mean raw dfof, p=',num2str(p)]);

meanDfof_sigAllFOVsNorm=[];
for n=1:size(meanDfof_sigAllFOVs,1);
    meanDfof_sigAllFOVsNorm(n,:)=meanDfof_sigAllFOVs(n,:)/meanDfof_sigAllFOVs(n,2);
end

meanmeanDfof_sigAllFOVsNorm=[];
semmeanDfof_sigAllFOVsNorm=[];
for n=1:size(meanDfof_sigAllFOVsNorm,2);
    meanmeanDfof_sigAllFOVsNorm(1,n)=mean(meanDfof_sigAllFOVsNorm(:,n));
    semmeanDfof_sigAllFOVsNorm(1,n)=std(meanDfof_sigAllFOVsNorm(:,n))/sqrt(size(meanDfof_sigAllFOVsNorm,1));
end

subplot(222)
errorbar([1:1:length(meanmeanDfof_sigAllFOVsNorm)],meanmeanDfof_sigAllFOVsNorm,semmeanDfof_sigAllFOVsNorm)
[~,p]=corr([1:1:length(meanmeanDfof_sigAllFOVsNorm)-1]',meanmeanDfof_sigAllFOVsNorm(2:end)');
title(['mean raw Dfof_sig, p=',num2str(p)]);

maxDfofAllFOVsNorm=[];
for n=1:size(maxDfofAllFOVs,1);
    maxDfofAllFOVsNorm(n,:)=maxDfofAllFOVs(n,:)/maxDfofAllFOVs(n,2);
end

meanmaxDfofAllFOVsNorm=[];
semmaxDfofAllFOVsNorm=[];
for n=1:size(maxDfofAllFOVsNorm,2);
    meanmaxDfofAllFOVsNorm(1,n)=mean(maxDfofAllFOVsNorm(:,n));
    semmaxDfofAllFOVsNorm(1,n)=std(maxDfofAllFOVsNorm(:,n))/sqrt(size(maxDfofAllFOVsNorm,1));
end

subplot(223)
errorbar([1:1:length(meanmaxDfofAllFOVsNorm)],meanmaxDfofAllFOVsNorm,semmaxDfofAllFOVsNorm)
[~,p]=corr([1:1:length(meanmaxDfofAllFOVsNorm)-1]',meanmaxDfofAllFOVsNorm(2:end)');
title(['max raw dfof, p=',num2str(p)]);

maxDfof_sigAllFOVsNorm=[];
for n=1:size(meanDfof_sigAllFOVs,1);
    maxDfof_sigAllFOVsNorm(n,:)=maxDfof_sigAllFOVs(n,:)/maxDfof_sigAllFOVs(n,2);
end

meanmaxDfof_sigAllFOVsNorm=[];
semmaxDfof_sigAllFOVsNorm=[];
for n=1:size(maxDfof_sigAllFOVsNorm,2);
    meanmaxDfof_sigAllFOVsNorm(1,n)=mean(maxDfof_sigAllFOVsNorm(:,n));
    semmaxDfof_sigAllFOVsNorm(1,n)=std(maxDfof_sigAllFOVsNorm(:,n))/sqrt(size(maxDfof_sigAllFOVsNorm,1));
end

subplot(224)
errorbar([1:1:length(meanmaxDfof_sigAllFOVsNorm)],meanmaxDfof_sigAllFOVsNorm,semmaxDfof_sigAllFOVsNorm)
[~,p]=corr([1:1:length(meanmaxDfof_sigAllFOVsNorm)-1]',meanmaxDfof_sigAllFOVsNorm(2:end)');
title(['max raw Dfof_sig, p=',num2str(p)]);


saveas(gcf,'rawDfofChangeOverDays.fig');


%% spatial information

%4 calculations
load('foldersAll.mat');
SI1=[];% spatial information, method 1 (mean lambda is mean fluorescence) using regular speed threshold. %each row is a cell and each column is a day
SI2=[];% spatial information,method 1  (mean lambda is sum of pi*lambdai) using 1cm/s speed thresh as speed threshold. %each row is a cell and each column is a day
SI3=[];% spatial information, new method 2 (mean lambda is mean fluorescence) using regular speed threshold. %each row is a cell and each column is a day
SI4=[];% spatial information,new method 2  (mean lambda is sum of pi*lambdai) using 1cm/s speed thresh as speed threshold. %each row is a cell and each column is a day

for n=1:length(foldersAll);
  load([foldersAll{n} '\' 'commonCells.mat']);
  load([foldersAll{n} '\' 'useFolders.mat']);
  disp(n)
  A1=[];
  A2=[];
  A3=[];
  A4=[];
  for m=1:size(commonCells,2);
       disp(m)
      a=commonCells(:,m);
      load([useFolders{m} '\spatialInfo_sig\SI.mat']);
      A1(:,m)=SI(a);
       load([useFolders{m} '\spatialInfo_sig_1cms\SI.mat']);
      A2(:,m)=SI(a);
      load([useFolders{m} '\spatialInfo_sig_New\SI.mat']);
       A3(:,m)=SI(a);
        load([useFolders{m} '\spatialInfo_sig_New_1cms\SI.mat']);
       A4(:,m)=SI(a);

  end
  SI1(end+1:end+size(A1,1),:)=A1;
  SI2(end+1:end+size(A2,1),:)=A2;
  SI3(end+1:end+size(A1,1),:)=A3;
  SI4(end+1:end+size(A2,1),:)=A4;
end


figure,
A=SI1;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));

subplot(221)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)','tail','right');
xlabel('days')
ylabel('spatial information');
title(['SI p=',num2str(p)]);

A=SI2;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));

subplot(222)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)','tail','right');
xlabel('days')
ylabel('spatial information 1cms');
title(['SI 1cms p=',num2str(p)]);

A=SI3;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));

subplot(223)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)','tail','right');
xlabel('days')
ylabel('spatial information new');
title(['SI new p=',num2str(p)]);

A=SI4;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));

subplot(224)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)','tail','right');
xlabel('days')
ylabel('spatial information new 1cms');
title(['SI new 1cms p=',num2str(p)]);

saveas(gcf,'spatialInformation.fig');
save('SI.mat','SI1','SI2','SI3','SI4');
%% fluorescence

%3 calculations
pp=pwd;
load('foldersAll.mat');
mF=[];% mean fluorescence. %each row is a cell and each column is a day
mDfof=[];% mean dfof, initial calculation,using 1cm/s speed thresh. %each row is a cell and each column is a day
mDfofSig=[];% mean dfof sig, new calculation,using 1cm/s speed thresh. %each row is a cell and each column is a day
mDfofMean=[]; %mean of mean dfof
mDfofSigMean=[]; %mean of mean dfof sig

for n=1:length(foldersAll);
  load([foldersAll{n} '\' 'commonCells.mat']);
  load([foldersAll{n} '\' 'useFolders.mat']);
  disp(n)
  A1=[];
  A2=[];
  A3=[];
  A4=[];
  A5=[];
  for m=1:size(commonCells,2);
       disp(m)
      a=commonCells(:,m);
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

%% calculations: Fluorescence after speed threshold
pp=pwd;
load('foldersAll.mat');
mFSP=[];% mean fluorescence. %each row is a cell and each column is a day
mDfofSP=[];% mean dfof, initial calculation,using 1cm/s speed thresh. %each row is a cell and each column is a day
mDfofSigSP=[];% mean dfof sig, new calculation,using 1cm/s speed thresh. %each row is a cell and each column is a day
mDfofMean=[]; %mean of mean dfof
mDfofSigMean=[]; %mean of mean dfof sig

for n=1:length(foldersAll);
  load([foldersAll{n} '\' 'commonCells.mat']);
  load([foldersAll{n} '\' 'useFolders.mat']);
  disp(n)
  A1=[];
  A2=[];
  A3=[];
  A4=[];
  A5=[];
  for m=1:size(commonCells,2);
       disp(m)
      a=commonCells(:,m);
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

%speed threshold
load('abfFake.mat');
speed=diff(abfFake.y);
[speedThreshold]= speedThreshold1D( 1,length(abfFake.t),abfFake );
close
fastEnough=[false ;speed>speedThreshold];
i=find(fastEnough);

     F=F(i,:)';
     dfof=dfof(i,:)';
     dfof_sig=dfof_sig(i,:)';
 
     dfofaveragesmooth=dfofaveragesmooth';
     dfofaveragesmooth_sig=dfofaveragesmooth_sig';
     
      A1(:,m)=mean(F(a,:),2);
       A2(:,m)=mean(dfof(a,:),2);
       A3(:,m)=mean(dfof_sig(a,:),2);
       A4(:,m)=mean(dfofaveragesmooth(a,:),2);
       A5(:,m)=mean( dfofaveragesmooth_sig(a,:),2);
  end
  mFSP(end+1:end+size(A1,1),:)=A1;
  mDfofSP(end+1:end+size(A2,1),:)=A2;
  mDfofSigSP(end+1:end+size(A3,1),:)=A3;
  mDfofMean(end+1:end+size(A2,1),:)=A4;
  mDfofSigMean(end+1:end+size(A3,1),:)=A5;
end


figure,
A=mFSP;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));

subplot(151)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean F');
title(['mean F p=',num2str(p)]);

A=mDfofSP;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));

subplot(152)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean dfof');
title(['mean dfof p=',num2str(p)]);

A=mDfofSigSP;
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
saveas(gcf,'allFSP.fig');
save('allFSP.mat','mFSP','mDfofSP','mDfofSigSP','mDfofMean','mDfofSigMean');


%% fluorescence of all cells

%3 calculations
pp=pwd;
load('foldersAll.mat');
mFAllCells={};% mean fluorescence. %each row is a cell and each column is a day
mDfofAllCells={};% mean dfof, initial calculation,using 1cm/s speed thresh. %each row is a cell and each column is a day
mDfofSigAllCells={};% mean dfof sig, new calculation,using 1cm/s speed thresh. %each row is a cell and each column is a day
mDfofMeanAllCells={}; %mean of mean dfof
mDfofSigMeanAllCells={}; %mean of mean dfof sig

for n=1:11;
    mFAllCells{n}=[];% mean fluorescence. %each row is a cell and each column is a day
mDfofAllCells{n}=[];% mean dfof, initial calculation,using 1cm/s speed thresh. %each row is a cell and each column is a day
mDfofSigAllCells{n}=[];% mean dfof sig, new calculation,using 1cm/s speed thresh. %each row is a cell and each column is a day
mDfofMeanAllCells{n}=[]; %mean of mean dfof
mDfofSigMeanAllCells{n}=[];%mean of mean dfof sig
end

for n=1:length(foldersAll);
%   load([foldersAll{n} '\' 'commonCells.mat']);
  load([foldersAll{n} '\' 'useFolders.mat']);
  disp(n)
  A1=[];
  A2=[];
  A3=[];
  A4=[];
  A5=[];
  for m=1:length(useFolders);
       disp(m)
%       a=commonCells(:,m);
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

load('roiIdxUse.mat');

     F=F(:,roiIdxUse)';
     dfof=dfof(:,roiIdxUse)';
     dfof_sig=dfof_sig(:,roiIdxUse)';
 
     dfofaveragesmooth=dfofaveragesmooth(:,roiIdxUse)';
     dfofaveragesmooth_sig=dfofaveragesmooth_sig(:,roiIdxUse)';
     
      A1=nanmean(F,2);
       A2=nanmean(dfof,2);
       A3=nanmean(dfof_sig,2);
       A4=nanmean(dfofaveragesmooth,2);
       A5=nanmean(dfofaveragesmooth_sig,2);
  
  mFAllCells{m}(end+1:end+length(A1),:)=A1;
  mDfofAllCells{m}(end+1:end+length(A2),:)=A2;
  mDfofSigAllCells{m}(end+1:end+length(A3),:)=A3;
  mDfofMeanAllCells{m}(end+1:end+length(A4),:)=A4;
  mDfofSigMeanAllCells{m}(end+1:end+length(A5),:)=A5;
end
end
cd(pp);
figure,
A=mFAllCells;
M=[];
S=[];
for n=1:length(A);
M(n)=mean(A{n});
S(n)=nansem(A{n},1);
end

subplot(151)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean F');
title(['mean F p=',num2str(p)]);

A=mDfofAllCells;
M=[];
S=[];
for n=1:length(A);
M(n)=mean(A{n});
S(n)=nansem(A{n},1);
end

subplot(152)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean dfof');
title(['mean dfof p=',num2str(p)]);

A=mDfofSigAllCells;
M=[];
S=[];
for n=1:length(A);
M(n)=mean(A{n});
S(n)=nansem(A{n},1);
end

subplot(153)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean dfof sig');
title(['mean dfof sig p=',num2str(p)]);

A=mDfofMeanAllCells;
M=[];
S=[];
for n=1:length(A);
M(n)=mean(A{n});
S(n)=nansem(A{n},1);
end

subplot(154)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean dfof mean');
title(['mean dfof mean p=',num2str(p)]);

A=mDfofSigMeanAllCells;
M=[];
S=[];
for n=1:length(A);
M(n)=mean(A{n});
S(n)=nansem(A{n},1);
end

subplot(155)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean dfof sig mean');
title(['mean dfof sig mean p=',num2str(p)]);
tightfig;

cd(pp);
%
saveas(gcf,'allFAllCells.fig');
save('allFAllCells.mat','mFAllCells','mDfofAllCells','mDfofSigAllCells','mDfofMeanAllCells','mDfofSigMeanAllCells');


%% fluorescence of all cells: AFTER SPEED THRESHOLD

%3 calculations
pp=pwd;
load('foldersAll.mat');
mFAllCellsSP={};% mean fluorescence. %each row is a cell and each column is a day
mDfofAllCellsSP={};% mean dfof, initial calculation,using 1cm/s speed thresh. %each row is a cell and each column is a day
mDfofSigAllCellsSP={};% mean dfof sig, new calculation,using 1cm/s speed thresh. %each row is a cell and each column is a day
mDfofMeanAllCells={}; %mean of mean dfof
mDfofSigMeanAllCells={}; %mean of mean dfof sig

for n=1:11;
    mFAllCellsSP{n}=[];% mean fluorescence. %each row is a cell and each column is a day
mDfofAllCellsSP{n}=[];% mean dfof, initial calculation,using 1cm/s speed thresh. %each row is a cell and each column is a day
mDfofSigAllCellsSP{n}=[];% mean dfof sig, new calculation,using 1cm/s speed thresh. %each row is a cell and each column is a day
mDfofMeanAllCells{n}=[]; %mean of mean dfof
mDfofSigMeanAllCells{n}=[];%mean of mean dfof sig
end

for n=1:length(foldersAll);
%   load([foldersAll{n} '\' 'commonCells.mat']);
  load([foldersAll{n} '\' 'useFolders.mat']);
  disp(n)
  A1=[];
  A2=[];
  A3=[];
  A4=[];
  A5=[];
  for m=1:length(useFolders);
       disp(m)
%       a=commonCells(:,m);
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

load('roiIdxUse.mat');
%speed threshold
load('abfFake.mat');
speed=diff(abfFake.y);
[speedThreshold]= speedThreshold1D( 1,length(abfFake.t),abfFake );
close
fastEnough=[false ;speed>speedThreshold];
i=find(fastEnough);

     F=F(i,roiIdxUse)';
     dfof=dfof(i,roiIdxUse)';
     dfof_sig=dfof_sig(i,roiIdxUse)';
 
     dfofaveragesmooth=dfofaveragesmooth(:,roiIdxUse)';
     dfofaveragesmooth_sig=dfofaveragesmooth_sig(:,roiIdxUse)';
     
      A1=nanmean(F,2);
       A2=nanmean(dfof,2);
       A3=nanmean(dfof_sig,2);
       A4=nanmean(dfofaveragesmooth,2);
       A5=nanmean(dfofaveragesmooth_sig,2);
  
  mFAllCellsSP{m}(end+1:end+length(A1),:)=A1;
  mDfofAllCellsSP{m}(end+1:end+length(A2),:)=A2;
  mDfofSigAllCellsSP{m}(end+1:end+length(A3),:)=A3;
  mDfofMeanAllCells{m}(end+1:end+length(A4),:)=A4;
  mDfofSigMeanAllCells{m}(end+1:end+length(A5),:)=A5;
end
end
cd(pp);
figure,
A=mFAllCellsSP;
M=[];
S=[];
for n=1:length(A);
M(n)=mean(A{n});
S(n)=nansem(A{n},1);
end

subplot(151)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean F');
title(['mean F p=',num2str(p)]);

A=mDfofAllCellsSP;
M=[];
S=[];
for n=1:length(A);
M(n)=mean(A{n});
S(n)=nansem(A{n},1);
end

subplot(152)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean dfof');
title(['mean dfof p=',num2str(p)]);

A=mDfofSigAllCellsSP;
M=[];
S=[];
for n=1:length(A);
M(n)=mean(A{n});
S(n)=nansem(A{n},1);
end

subplot(153)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean dfof sig');
title(['mean dfof sig p=',num2str(p)]);

A=mDfofMeanAllCells;
M=[];
S=[];
for n=1:length(A);
M(n)=mean(A{n});
S(n)=nansem(A{n},1);
end

subplot(154)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean dfof mean');
title(['mean dfof mean p=',num2str(p)]);

A=mDfofSigMeanAllCells;
M=[];
S=[];
for n=1:length(A);
M(n)=mean(A{n});
S(n)=nansem(A{n},1);
end

subplot(155)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean dfof sig mean');
title(['mean dfof sig mean p=',num2str(p)]);
tightfig;

cd(pp);
%
saveas(gcf,'allFAllCellsSP.fig');
save('allFAllCellsSP.mat','mFAllCellsSP','mDfofAllCellsSP','mDfofSigAllCellsSP','mDfofMeanAllCells','mDfofSigMeanAllCells');

figure,
A=mFAllCellsSP;
S=[];
for n=1:length(A);
S(n)=sum(A{n},1);
end

subplot(151)
plot(S)
% [r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('sum F');
title(['sum F']);

A=mDfofAllCellsSP;
S=[];
for n=1:length(A);
S(n)=sum(A{n},1);
end


subplot(152)
plot(S)
xlabel('days')
ylabel('sum dfof');
title(['sum dfof']);

A=mDfofSigAllCellsSP;
S=[];
for n=1:length(A);
S(n)=sum(A{n},1);
end
subplot(153)
plot(S)
xlabel('days')
ylabel('sum dfof sig');
title(['sum dfof sig']);

A=mDfofMeanAllCells;
S=[];
for n=1:length(A);
S(n)=sum(A{n},1);
end

subplot(154)
plot(S)

xlabel('days')
ylabel('sum dfof mean');
title(['sum dfof mean']);

A=mDfofSigMeanAllCells;
S=[];
for n=1:length(A);
S(n)=sum(A{n},1);
end

subplot(155)
plot(S)

xlabel('days')
ylabel('sum dfof sig mean');
title(['sum dfof sig mean']);
tightfig;

saveas(gcf,'allFAllCellsSP_sum.fig');



%% calculations
pp=pwd;
load('foldersAll.mat');
mFAllCellsSPFOV=[];% mean fluorescence. %each row is a cell and each column is a day
mDfofAllCellsSPFOV=[];% mean dfof, initial calculation,using 1cm/s speed thresh. %each row is a cell and each column is a day
mDfofSigAllCellsSPFOV=[];% mean dfof sig, new calculation,using 1cm/s speed thresh. %each row is a cell and each column is a day
mDfofMeanAllCellsFOV=[]; %mean of mean dfof
mDfofSigMeanAllCellsFOV=[]; %mean of mean dfof sig

mFAllCellsSPUCFOV=[];% mean fluorescence. %each row is a cell and each column is a day
mDfofAllCellsSPUCFOV=[];% mean dfof, initial calculation,using 1cm/s speed thresh. %each row is a cell and each column is a day
mDfofSigAllCellsSPUCFOV=[];% mean dfof sig, new calculation,using 1cm/s speed thresh. %each row is a cell and each column is a day
mDfofMeanAllCellsUCFOV=[]; %mean of mean dfof
mDfofSigMeanAllCellsCUFOV=[]; %mean of mean dfof sig


% for n=1:11;
%     mFAllCellsSPFOV{n}=[];% mean fluorescence. %each row is a cell and each column is a day
% mDfofAllCellsSPFOV{n}=[];% mean dfof, initial calculation,using 1cm/s speed thresh. %each row is a cell and each column is a day
% mDfofSigAllCellsSPFOV{n}=[];% mean dfof sig, new calculation,using 1cm/s speed thresh. %each row is a cell and each column is a day
% mDfofMeanAllCellsFOV{n}=[]; %mean of mean dfof
% mDfofSigMeanAllCellsFOV{n}=[];%mean of mean dfof sig
% end
% 
% 
% for n=1:11;
%     mFAllCellsSPUCFOV{n}=[];% mean fluorescence. %each row is a cell and each column is a day
% mDfofAllCellsSPUCFOV{n}=[];% mean dfof, initial calculation,using 1cm/s speed thresh. %each row is a cell and each column is a day
% mDfofSigAllCellsSPUCFOV{n}=[];% mean dfof sig, new calculation,using 1cm/s speed thresh. %each row is a cell and each column is a day
% mDfofMeanAllCellsUCFOV{n}=[]; %mean of mean dfof
% mDfofSigMeanAllCellsUCFOV{n}=[];%mean of mean dfof sig
% end



for n=1:length(foldersAll);
  load([foldersAll{n} '\' 'commonCells.mat']);
  load([foldersAll{n} '\' 'useFolders.mat']);
  disp(n)
  A1=[];
  A2=[];
  A3=[];
  A4=[];
  A5=[];
  
   A1UC=[];
  A2UC=[];
  A3UC=[];
  A4UC=[];
  A5UC=[];
  for m=1:length(useFolders);
       disp(m)
%       a=commonCells(:,m);
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

load('roiIdxUse.mat');
%speed threshold
load('abfFake.mat');
speed=diff(abfFake.y);
[speedThreshold]= speedThreshold1D( 1,length(abfFake.t),abfFake );
close
fastEnough=[false ;speed>speedThreshold];
i=find(fastEnough);

ucIdx=setdiff(find(roiIdxUse),commonCells(:,m));

     FAll=F(i,roiIdxUse)';
     dfofAll=dfof(i,roiIdxUse)';
     dfof_sigAll=dfof_sig(i,roiIdxUse)';
 
     dfofaveragesmoothAll=dfofaveragesmooth(:,roiIdxUse)';
     dfofaveragesmooth_sigAll=dfofaveragesmooth_sig(:,roiIdxUse)';
     
       FUC=F(i,ucIdx)';
     dfofUC=dfof(i,ucIdx)';
     dfof_sigUC=dfof_sig(i,ucIdx)';
 
     dfofaveragesmoothUC=dfofaveragesmooth(:,ucIdx)';
     dfofaveragesmooth_sigUC=dfofaveragesmooth_sig(:,ucIdx)';
     
      A1=nanmean(nanmean(FAll,2));
       A2=nanmean(nanmean(dfofAll,2));
       A3=nanmean(nanmean(dfof_sigAll,2));
       A4=nanmean(nanmean(dfofaveragesmoothAll,2));
       A5=nanmean(nanmean(dfofaveragesmooth_sigAll,2));
       
             A1UC=nanmean(nanmean(FUC,2));
       A2UC=nanmean(nanmean(dfofUC,2));
       A3UC=nanmean(nanmean(dfof_sigUC,2));
       A4UC=nanmean(nanmean(dfofaveragesmoothUC,2));
       A5UC=nanmean(nanmean(dfofaveragesmooth_sigUC,2));
  
  mFAllCellsSPFOV(n,m)=A1;
  mDfofAllCellsSPFOV(n,m)=A2;
  mDfofSigAllCellsSPFOV(n,m)=A3;
  mDfofMeanAllCellsFOV(n,m)=A4;
  mDfofSigMeanAllCellsFOV(n,m)=A5;
  
   mFAllCellsSPUCFOV(n,m)=A1UC;
  mDfofAllCellsSPUCFOV(n,m)=A2UC;
  mDfofSigAllCellsSPUCFOV(n,m)=A3UC;
  mDfofMeanAllCellsUCFOV(n,m)=A4UC;
  mDfofSigMeanAllCellsUCFOV(n,m)=A5UC;
  
end
end
cd(pp);

figure,
A=mFAllCellsSPFOV(2:end,:);

M=nanmean(A,1);
S=nansem(A,1);

subplot(151)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean F');
title(['mean F p=',num2str(p)]);

A=mDfofAllCellsSPFOV(2:end,:);

M=nanmean(A,1);
S=nansem(A,1);

subplot(152)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean dfof');
title(['mean dfof p=',num2str(p)]);

A=mDfofSigAllCellsSPFOV(2:end,:);
M=nanmean(A,1);
S=nansem(A,1);

subplot(153)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean dfof sig');
title(['mean dfof sig p=',num2str(p)]);

A=mDfofMeanAllCellsFOV(2:end,:);
M=nanmean(A,1);
S=nansem(A,1);

subplot(154)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean dfof mean');
title(['mean dfof mean p=',num2str(p)]);

A=mDfofSigMeanAllCellsFOV(2:end,:);
M=nanmean(A,1);
S=nansem(A,1);

subplot(155)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean dfof sig mean');
title(['mean dfof sig mean p=',num2str(p)]);
tightfig;

cd(pp);
%
saveas(gcf,'allFAllCellsSPMeanFOV.fig');
save('allFAllCellsSPFOV.mat','mFAllCellsSPFOV','mDfofAllCellsSPFOV','mDfofSigAllCellsSPFOV','mDfofMeanAllCellsFOV','mDfofSigMeanAllCellsFOV');

figure,
A=mFAllCellsSPUCFOV(2:end,:);

M=nanmean(A,1);
S=nansem(A,1);

subplot(151)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean F');
title(['mean F p=',num2str(p)]);

A=mDfofAllCellsSPUCFOV(2:end,:);

M=nanmean(A,1);
S=nansem(A,1);

subplot(152)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean dfof');
title(['mean dfof p=',num2str(p)]);

A=mDfofSigAllCellsSPUCFOV(2:end,:);
M=nanmean(A,1);
S=nansem(A,1);

subplot(153)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean dfof sig');
title(['mean dfof sig p=',num2str(p)]);

A=mDfofMeanAllCellsUCFOV(2:end,:);
M=nanmean(A,1);
S=nansem(A,1);

subplot(154)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean dfof mean');
title(['mean dfof mean p=',num2str(p)]);

A=mDfofSigMeanAllCellsUCFOV(2:end,:);
M=nanmean(A,1);
S=nansem(A,1);

subplot(155)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean dfof sig mean');
title(['mean dfof sig mean p=',num2str(p)]);
tightfig;

cd(pp);

saveas(gcf,'allFAllCellsSPMeanFOVUncommonCells.fig');
save('allFAllCellsSPUCFOV.mat','mFAllCellsSPUCFOV','mDfofAllCellsSPUCFOV','mDfofSigAllCellsSPUCFOV','mDfofMeanAllCellsUCFOV','mDfofSigMeanAllCellsUCFOV');



%% calculations
pp=pwd;
load('foldersAll.mat');
mFAllCellsSumSPFOV=[];% mean fluorescence. %each row is a cell and each column is a day
mDfofAllCellsSumSPFOV=[];% mean dfof, initial calculation,using 1cm/s speed thresh. %each row is a cell and each column is a day
mDfofSigAllCellsSumSPFOV=[];% mean dfof sig, new calculation,using 1cm/s speed thresh. %each row is a cell and each column is a day
mDfofMeanAllCellsSumFOV=[]; %mean of mean dfof
mDfofSigMeanAllCellsSumFOV=[]; %mean of mean dfof sig

mFAllCellsSumSPUCFOV=[];% mean fluorescence. %each row is a cell and each column is a day
mDfofAllCellsSumSPUCFOV=[];% mean dfof, initial calculation,using 1cm/s speed thresh. %each row is a cell and each column is a day
mDfofSigAllCellsSumSPUCFOV=[];% mean dfof sig, new calculation,using 1cm/s speed thresh. %each row is a cell and each column is a day
mDfofMeanAllCellsSumUCFOV=[]; %mean of mean dfof
mDfofSigMeanAllCellsSumUCFOV=[]; %mean of mean dfof sig


% for n=1:11;
%     mFAllCellsSPFOV{n}=[];% mean fluorescence. %each row is a cell and each column is a day
% mDfofAllCellsSPFOV{n}=[];% mean dfof, initial calculation,using 1cm/s speed thresh. %each row is a cell and each column is a day
% mDfofSigAllCellsSPFOV{n}=[];% mean dfof sig, new calculation,using 1cm/s speed thresh. %each row is a cell and each column is a day
% mDfofMeanAllCellsFOV{n}=[]; %mean of mean dfof
% mDfofSigMeanAllCellsFOV{n}=[];%mean of mean dfof sig
% end
% 
% 
% for n=1:11;
%     mFAllCellsSPUCFOV{n}=[];% mean fluorescence. %each row is a cell and each column is a day
% mDfofAllCellsSPUCFOV{n}=[];% mean dfof, initial calculation,using 1cm/s speed thresh. %each row is a cell and each column is a day
% mDfofSigAllCellsSPUCFOV{n}=[];% mean dfof sig, new calculation,using 1cm/s speed thresh. %each row is a cell and each column is a day
% mDfofMeanAllCellsUCFOV{n}=[]; %mean of mean dfof
% mDfofSigMeanAllCellsUCFOV{n}=[];%mean of mean dfof sig
% end



for n=1:length(foldersAll);
  load([foldersAll{n} '\' 'commonCells.mat']);
  load([foldersAll{n} '\' 'useFolders.mat']);
  disp(n)
  A1=[];
  A2=[];
  A3=[];
  A4=[];
  A5=[];
  
   A1UC=[];
  A2UC=[];
  A3UC=[];
  A4UC=[];
  A5UC=[];
  for m=1:length(useFolders);
       disp(m)
%       a=commonCells(:,m);
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

load('roiIdxUse.mat');
%speed threshold
load('abfFake.mat');
speed=diff(abfFake.y);
[speedThreshold]= speedThreshold1D( 1,length(abfFake.t),abfFake );
close
fastEnough=[false ;speed>speedThreshold];
i=find(fastEnough);

ucIdx=setdiff(find(roiIdxUse),commonCells(:,m));

     FAll=F(i,roiIdxUse)';
     dfofAll=dfof(i,roiIdxUse)';
     dfof_sigAll=dfof_sig(i,roiIdxUse)';
 
     dfofaveragesmoothAll=dfofaveragesmooth(:,roiIdxUse)';
     dfofaveragesmooth_sigAll=dfofaveragesmooth_sig(:,roiIdxUse)';
     
       FUC=F(i,ucIdx)';
     dfofUC=dfof(i,ucIdx)';
     dfof_sigUC=dfof_sig(i,ucIdx)';
 
     dfofaveragesmoothUC=dfofaveragesmooth(:,ucIdx)';
     dfofaveragesmooth_sigUC=dfofaveragesmooth_sig(:,ucIdx)';
     
      A1=nansum(nanmean(FAll,2));
       A2=nansum(nanmean(dfofAll,2));
       A3=nansum(nanmean(dfof_sigAll,2));
       A4=nansum(nanmean(dfofaveragesmoothAll,2));
       A5=nansum(nanmean(dfofaveragesmooth_sigAll,2));
       
             A1UC=nansum(nanmean(FUC,2));
       A2UC=nansum(nanmean(dfofUC,2));
       A3UC=nansum(nanmean(dfof_sigUC,2));
       A4UC=nansum(nanmean(dfofaveragesmoothUC,2));
       A5UC=nansum(nanmean(dfofaveragesmooth_sigUC,2));
  
  mFAllCellsSumSPFOV(n,m)=A1;
  mDfofAllCellsSumSPFOV(n,m)=A2;
  mDfofSigAllCellsSumSPFOV(n,m)=A3;
  mDfofMeanAllCellsSumFOV(n,m)=A4;
  mDfofSigMeanAllCellsSumFOV(n,m)=A5;
  
   mFAllCellsSumSPUCFOV(n,m)=A1UC;
  mDfofAllCellsSumSPUCFOV(n,m)=A2UC;
  mDfofSigAllCellsSumSPUCFOV(n,m)=A3UC;
  mDfofMeanAllCellsSumUCFOV(n,m)=A4UC;
  mDfofSigMeanAllCellsSumUCFOV(n,m)=A5UC;
  
end
end
cd(pp);

figure,
A=mFAllCellsSumSPFOV(2:end,:);

M=nanmean(A,1);
S=nansem(A,1);

subplot(151)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean F');
title(['mean F p=',num2str(p)]);

A=mDfofAllCellsSumSPFOV(2:end,:);

M=nanmean(A,1);
S=nansem(A,1);

subplot(152)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean dfof');
title(['mean dfof p=',num2str(p)]);

A=mDfofSigAllCellsSumSPFOV(2:end,:);
M=nanmean(A,1);
S=nansem(A,1);

subplot(153)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean dfof sig');
title(['mean dfof sig p=',num2str(p)]);

A=mDfofMeanAllCellsSumFOV(2:end,:);
M=nanmean(A,1);
S=nansem(A,1);

subplot(154)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean dfof mean');
title(['mean dfof mean p=',num2str(p)]);

A=mDfofSigMeanAllCellsSumFOV(2:end,:);
M=nanmean(A,1);
S=nansem(A,1);

subplot(155)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean dfof sig mean');
title(['mean dfof sig mean p=',num2str(p)]);
tightfig;

cd(pp);
%
saveas(gcf,'allFAllCellsSPMeanSumFOV.fig');
save('allFAllCellsSumSPFOV.mat','mFAllCellsSumSPFOV','mDfofAllCellsSumSPFOV','mDfofSigAllCellsSumSPFOV','mDfofMeanAllCellsSumFOV','mDfofSigMeanAllCellsSumFOV');

figure,
A=mFAllCellsSumSPUCFOV(2:end,:);

M=nanmean(A,1);
S=nansem(A,1);

subplot(151)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean F');
title(['mean F p=',num2str(p)]);

A=mDfofAllCellsSumSPUCFOV(2:end,:);

M=nanmean(A,1);
S=nansem(A,1);

subplot(152)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean dfof');
title(['mean dfof p=',num2str(p)]);

A=mDfofSigAllCellsSumSPUCFOV(2:end,:);
M=nanmean(A,1);
S=nansem(A,1);

subplot(153)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean dfof sig');
title(['mean dfof sig p=',num2str(p)]);

A=mDfofMeanAllCellsSumUCFOV(2:end,:);
M=nanmean(A,1);
S=nansem(A,1);

subplot(154)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean dfof mean');
title(['mean dfof mean p=',num2str(p)]);

A=mDfofSigMeanAllCellsSumUCFOV(2:end,:);
M=nanmean(A,1);
S=nansem(A,1);

subplot(155)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean dfof sig mean');
title(['mean dfof sig mean p=',num2str(p)]);
tightfig;

cd(pp);

saveas(gcf,'allFAllCellsSPMeanSumFOVUncommonCells.fig');
save('allFAllCellsSumSPUCFOV.mat','mFAllCellsSumSPUCFOV','mDfofAllCellsSumSPUCFOV','mDfofSigAllCellsSumSPUCFOV','mDfofMeanAllCellsSumUCFOV','mDfofSigMeanAllCellsSumUCFOV');


%% sig trans


pp=pwd;
load('foldersAll.mat');
mFreq=[];% mean frequency of sig trans per cm. %each row is a cell and each column is a day
mFreqRun=[];%mean frequency calcualted as mean of run by run
mDisCover=[];%distance covered by sig trans
mDisCoverRun=[];%distance covered by sig trans calcualted as mean of run by run
mAmp=[];% mean amplitudes of sig trans. %each row is a cell and each column is a day
mSigInteg=[];%significant transnet integrated by tune
mDur=[];% mean durations of sig trans. %each row is a cell and each column is a day
mSigIntegPerRun=[]; %integorated sig trans per run
for n=1:length(foldersAll);
  load([foldersAll{n} '\' 'commonCells.mat']);
  load([foldersAll{n} '\' 'useFolders.mat']);
  disp(n)
  
A1=[];
A2=[];
A3=[];
A4=[];
A5=[];
A6=[];
A7=[];
A8=[];

  for m=1:size(commonCells,2);
       disp(m)
      a=commonCells(:,m);
      cd(useFolders{m});
      
load('sigTransInfo.mat');
load('Freq.mat');
load('sigTransInteg.mat');
  load('sigTransIntegPerRun.mat');   
      A1(:,m)=Freq(a);
       A2(:,m)=FreqRun(a);
      A3(:,m)=distCovered(a);
       A4(:,m)=distCoveredRun(a);
        A6(:,m)=sigTransInteg(a);
       
      for i=1:length(a);
       A5(i,m)=mean(sigAmp{a(i)});
       
       aa=t_on_off{a(i)};
       A7(i,m)=mean(aa(:,3)-aa(:,2)+1)*0.1;%frame interval is 0.1 
      end
    A8(:,m)=sigTransIntegPerRun(a);    
  end

  
  mFreq(end+1:end+size(A1,1),:)=A1;
  mFreqRun(end+1:end+size(A2,1),:)=A2;
  mDisCover(end+1:end+size(A3,1),:)=A3;
   mDisCoverRun(end+1:end+size(A4,1),:)=A4;
    mAmp(end+1:end+size(A5,1),:)=A5;
    mSigInteg(end+1:end+size(A6,1),:)=A6;
    mDur(end+1:end+size(A7,1),:)=A7;
    mSigIntegPerRun(end+1:end+size(A8,1),:)=A8;
end

cd(pp);

figure,
A=mFreq;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));

subplot(241)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean Freq');
title(['mean Freq p=',num2str(p)]);

A=mFreqRun;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));

subplot(242)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean Freq run');
title(['mean Freq run p=',num2str(p)]);

A=mDisCover;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));

subplot(243)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean Dis cover');
title(['mean Dis cover p=',num2str(p)]);

A=mDisCoverRun;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));

subplot(244)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean Dis cover Run');
title(['mean Dis cover Run p=',num2str(p)]);

A=mAmp;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));

subplot(245)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean Amp');
title(['mean Amp p=',num2str(p)]);

A=mSigInteg;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));

subplot(246)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean SigInteg');
title(['mean SigInteg p=',num2str(p)]);

A=mDur;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));

subplot(247)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean Dur');
title(['mean Dur p=',num2str(p)]);

A=mSigIntegPerRun;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));

subplot(248)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean sig integ per run');
title(['mean sig integ per run p=',num2str(p)]);

saveas(gcf,'allSigTrans.fig');
save('allSigTrans.mat','mFreq','mFreqRun','mDisCover','mDisCoverRun','mAmp','mSigInteg','mDur','mSigIntegPerRun');

%% sig trans of all cells

pp=pwd;
load('foldersAll.mat');
mFreqAllCells=[];% mean frequency of sig trans per cm. %each row is a cell and each column is a day
mFreqRunAllCells=[];%mean frequency calcualted as mean of run by run
mDisCoverAllCells=[];%distance covered by sig trans
mDisCoverRunAllCells=[];%distance covered by sig trans calcualted as mean of run by run
mAmpAllCells=[];% mean amplitudes of sig trans. %each row is a cell and each column is a day
mSigIntegAllCells=[];%significant transnet integrated by tune
mDurAllCells=[];% mean durations of sig trans. %each row is a cell and each column is a day
mSigIntegPerRunAllCells=[]; %integorated sig trans per run

for n=1:11;
    mFreqAllCells{n}=[];% mean frequency of sig trans per cm. %each row is a cell and each column is a day
mFreqRunAllCells{n}=[];%mean frequency calcualted as mean of run by run
mDisCoverAllCells{n}=[];%distance covered by sig trans
mDisCoverRunAllCells{n}=[];%distance covered by sig trans calcualted as mean of run by run
mAmpAllCells{n}=[];% mean amplitudes of sig trans. %each row is a cell and each column is a day
mSigIntegAllCells{n}=[];%significant transnet integrated by tune
mDurAllCells{n}=[];% mean durations of sig trans. %each row is a cell and each column is a day
mSigIntegPerRunAllCells{n}=[]; %integorated sig trans per run
end

for n=1:length(foldersAll);
%   load([foldersAll{n} '\' 'commonCells.mat']);
  load([foldersAll{n} '\' 'useFolders.mat']);
  disp(n)
  

  for m=1:length(useFolders);
      
      A1=[];
A2=[];
A3=[];
A4=[];
A5=[];
A6=[];
A7=[];
A8=[];
       disp(m)
%       a=commonCells(:,m);
      cd(useFolders{m});
      
load('sigTransInfo.mat');
load('Freq.mat');
load('sigTransInteg.mat');
load('sigTransIntegPerRun.mat');
   load('roiIdxUse.mat');
   a=find(roiIdxUse);
   
      A1=Freq(a);
       A2=FreqRun(a);
      A3=distCovered(a);
       A4=distCoveredRun(a);
        A6=sigTransInteg(a);
       
      for i=1:length(a);
       A5(i,1)=mean(sigAmp{a(i)});
       
       if isempty(t_on_off{a(i)});
           A7(i,1)=nan;
       elseif t_on_off{a(i)}==0;
           A7(i,1)=nan;
       else
       aa=t_on_off{a(i)};
       A7(i,1)=mean(aa(:,3)-aa(:,2)+1)*0.1;%frame interval is 0.1        
       end
      end
  A8=sigTransIntegPerRun(a);   
  
  mFreqAllCells{m}(end+1:end+length(A1),:)=A1;
  mFreqRunAllCells{m}(end+1:end+length(A2),:)=A2;
  mDisCoverAllCells{m}(end+1:end+length(A3),:)=A3;
   mDisCoverRunAllCells{m}(end+1:end+length(A4),:)=A4;
    mAmpAllCells{m}(end+1:end+length(A5),:)=A5;
    mSigIntegAllCells{m}(end+1:end+length(A6),:)=A6;
    mDurAllCells{m}(end+1:end+length(A7),:)=A7;
    mSigIntegPerRunAllCells{m}(end+1:end+length(A8),:)=A8;
  end
end

cd(pp);

figure,
A=mFreqAllCells;
M=[];
S=[];
for n=1:length(A);
M(n)=nanmean(A{n});
S(n)=nansem(A{n},1);
end

subplot(241)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean Freq');
title(['mean Freq p=',num2str(p)]);

A=mFreqRunAllCells;
M=[];
S=[];
for n=1:length(A);
M(n)=nanmean(A{n});
S(n)=nansem(A{n},1);
end

subplot(242)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean Freq run');
title(['mean Freq run p=',num2str(p)]);

A=mDisCoverAllCells;
M=[];
S=[];
for n=1:length(A);
M(n)=nanmean(A{n});
S(n)=nansem(A{n},1);
end

subplot(243)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean Dis cover');
title(['mean Dis cover p=',num2str(p)]);

A=mDisCoverRunAllCells;
M=[];
S=[];
for n=1:length(A);
M(n)=nanmean(A{n});
S(n)=nansem(A{n},1);
end

subplot(244)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean Dis cover Run');
title(['mean Dis cover Run p=',num2str(p)]);

A=mAmpAllCells;
M=[];
S=[];
for n=1:length(A);
M(n)=nanmean(A{n});
S(n)=nansem(A{n},1);
end

subplot(245)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean Amp');
title(['mean Amp p=',num2str(p)]);

A=mSigIntegAllCells;
M=[];
S=[];
for n=1:length(A);
M(n)=nanmean(A{n});
S(n)=nansem(A{n},1);
end

subplot(246)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean SigInteg');
title(['mean SigInteg p=',num2str(p)]);

A=mDurAllCells;
M=[];
S=[];
for n=1:length(A);
M(n)=nanmean(A{n});
S(n)=nansem(A{n},1);
end

subplot(247)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean Dur');
title(['mean Dur p=',num2str(p)]);


A=mSigIntegPerRunAllCells;
M=[];
S=[];
for n=1:length(A);
M(n)=nanmean(A{n});
S(n)=nansem(A{n},1);
end

subplot(248)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean sig integ per run');
title(['mean sig integ per run p=',num2str(p)]);


saveas(gcf,'allSigTransAllCells.fig');
save('allSigTransAllCells.mat','mFreqAllCells','mFreqRunAllCells','mDisCoverAllCells','mDisCoverRunAllCells','mAmpAllCells','mSigIntegAllCells','mDurAllCells','mSigIntegPerRunAllCells');


%% sig trans_remove the one below speed threshold

pp=pwd;
load('foldersAll.mat');
mFreqSP=[];% mean frequency of sig trans per cm. %each row is a cell and each column is a day
mFreqRunSP=[];%mean frequency calcualted as mean of run by run
mDisCoverSP=[];%distance covered by sig trans
mDisCoverRunSP=[];%distance covered by sig trans calcualted as mean of run by run
mAmpSP=[];% mean amplitudes of sig trans. %each row is a cell and each column is a day
mSigIntegSP=[];%significant transnet integrated by tune
mDurSP=[];% mean durations of sig trans. %each row is a cell and each column is a day
mSigIntegPerRunSP=[]; %integorated sig trans per run
for n=1:length(foldersAll);
  load([foldersAll{n} '\' 'commonCells.mat']);
  load([foldersAll{n} '\' 'useFolders.mat']);
  disp(n)
  
A1=[];
A2=[];
A3=[];
A4=[];
A5=[];
A6=[];
A7=[];
A8=[];

  for m=1:size(commonCells,2);
       disp(m)
      a=commonCells(:,m);
      cd(useFolders{m});
      
load('sigTransInfo.mat');
load('Freq.mat');
load('sigTransIntegSP.mat');
  load('sigTransIntegPerRunSP.mat');   
      A1(:,m)=FreqSP(a);
       A2(:,m)=FreqRunSP(a);
      A3(:,m)=distCoveredSP(a);
       A4(:,m)=distCoveredRunSP(a);
        A6(:,m)=sigTransIntegSP(a);
       
      for i=1:length(a);
       A5(i,m)=mean(sigAmpSP{a(i)});
       
       aa=t_on_offSP{a(i)};
       if aa~=0;
       A7(i,m)=mean(aa(:,3)-aa(:,2)+1)*0.1;%frame interval is 0.1 
       end
      end
    A8(:,m)=sigTransIntegPerRunSP(a);    
  end

  
  mFreqSP(end+1:end+size(A1,1),:)=A1;
  mFreqRunSP(end+1:end+size(A2,1),:)=A2;
  mDisCoverSP(end+1:end+size(A3,1),:)=A3;
   mDisCoverRunSP(end+1:end+size(A4,1),:)=A4;
    mAmpSP(end+1:end+size(A5,1),:)=A5;
    mSigIntegSP(end+1:end+size(A6,1),:)=A6;
    mDurSP(end+1:end+size(A7,1),:)=A7;
    mSigIntegPerRunSP(end+1:end+size(A8,1),:)=A8;
end

cd(pp);

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

%% sig trans of all cells_remove the one below speed threshold

pp=pwd;
load('foldersAll.mat');
mFreqAllCellsSP=[];% mean frequency of sig trans per cm. %each row is a cell and each column is a day
mFreqRunAllCellsSP=[];%mean frequency calcualted as mean of run by run
mDisCoverAllCellsSP=[];%distance covered by sig trans
mDisCoverRunAllCellsSP=[];%distance covered by sig trans calcualted as mean of run by run
mAmpAllCellsSP=[];% mean amplitudes of sig trans. %each row is a cell and each column is a day
mSigIntegAllCellsSP=[];%significant transnet integrated by tune
mDurAllCellsSP=[];% mean durations of sig trans. %each row is a cell and each column is a day
mSigIntegPerRunAllCellsSP=[]; %integorated sig trans per run

for n=1:11;
    mFreqAllCellsSP{n}=[];% mean frequency of sig trans per cm. %each row is a cell and each column is a day
mFreqRunAllCellsSP{n}=[];%mean frequency calcualted as mean of run by run
mDisCoverAllCellsSP{n}=[];%distance covered by sig trans
mDisCoverRunAllCellsSP{n}=[];%distance covered by sig trans calcualted as mean of run by run
mAmpAllCellsSP{n}=[];% mean amplitudes of sig trans. %each row is a cell and each column is a day
mSigIntegAllCellsSP{n}=[];%significant transnet integrated by tune
mDurAllCellsSP{n}=[];% mean durations of sig trans. %each row is a cell and each column is a day
mSigIntegPerRunAllCellsSP{n}=[]; %integorated sig trans per run
end

for n=1:length(foldersAll);
%   load([foldersAll{n} '\' 'commonCells.mat']);
  load([foldersAll{n} '\' 'useFolders.mat']);
  disp(n)
  

  for m=1:length(useFolders);
      
      A1=[];
A2=[];
A3=[];
A4=[];
A5=[];
A6=[];
A7=[];
A8=[];
       disp(m)
%       a=commonCells(:,m);
      cd(useFolders{m});
      
load('sigTransInfo.mat');
load('Freq.mat');
load('sigTransIntegSP.mat');
load('sigTransIntegPerRunSP.mat');
   load('roiIdxUse.mat');
   a=find(roiIdxUse);
   
      A1=FreqSP(a);
       A2=FreqRunSP(a);
      A3=distCoveredSP(a);
       A4=distCoveredRunSP(a);
        A6=sigTransIntegSP(a);
       
      for i=1:length(a);
       A5(i,1)=mean(sigAmpSP{a(i)});
       
       if isempty(t_on_offSP{a(i)});
           A7(i,1)=nan;
       elseif t_on_offSP{a(i)}==0;
           A7(i,1)=nan;
       else
           
       aa=t_on_offSP{a(i)};
       if aa~=0;
       A7(i,1)=mean(aa(:,3)-aa(:,2)+1)*0.1;%frame interval is 0.1   
       end
       end
      end
  A8=sigTransIntegPerRunSP(a);   
  
  mFreqAllCellsSP{m}(end+1:end+length(A1),:)=A1;
  mFreqRunAllCellsSP{m}(end+1:end+length(A2),:)=A2;
  mDisCoverAllCellsSP{m}(end+1:end+length(A3),:)=A3;
   mDisCoverRunAllCellsSP{m}(end+1:end+length(A4),:)=A4;
    mAmpAllCellsSP{m}(end+1:end+length(A5),:)=A5;
    mSigIntegAllCellsSP{m}(end+1:end+length(A6),:)=A6;
    mDurAllCellsSP{m}(end+1:end+length(A7),:)=A7;
    mSigIntegPerRunAllCellsSP{m}(end+1:end+length(A8),:)=A8;
  end
end

cd(pp);

figure,
A=mFreqAllCellsSP;
M=[];
S=[];
for n=1:length(A);
M(n)=nanmean(A{n});
S(n)=nansem(A{n},1);
end

subplot(241)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean Freq');
title(['mean Freq p=',num2str(p)]);

A=mFreqRunAllCellsSP;
M=[];
S=[];
for n=1:length(A);
M(n)=nanmean(A{n});
S(n)=nansem(A{n},1);
end

subplot(242)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean Freq run');
title(['mean Freq run p=',num2str(p)]);

A=mDisCoverAllCellsSP;
M=[];
S=[];
for n=1:length(A);
M(n)=nanmean(A{n});
S(n)=nansem(A{n},1);
end

subplot(243)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean Dis cover');
title(['mean Dis cover p=',num2str(p)]);

A=mDisCoverRunAllCellsSP;
M=[];
S=[];
for n=1:length(A);
M(n)=nanmean(A{n});
S(n)=nansem(A{n},1);
end

subplot(244)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean Dis cover Run');
title(['mean Dis cover Run p=',num2str(p)]);

A=mAmpAllCellsSP;
M=[];
S=[];
for n=1:length(A);
M(n)=nanmean(A{n});
S(n)=nansem(A{n},1);
end

subplot(245)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean Amp');
title(['mean Amp p=',num2str(p)]);

A=mSigIntegAllCellsSP;
M=[];
S=[];
for n=1:length(A);
M(n)=nanmean(A{n});
S(n)=nansem(A{n},1);
end

subplot(246)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean SigInteg');
title(['mean SigInteg p=',num2str(p)]);

A=mDurAllCellsSP;
M=[];
S=[];
for n=1:length(A);
M(n)=nanmean(A{n});
S(n)=nansem(A{n},1);
end

subplot(247)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean Dur');
title(['mean Dur p=',num2str(p)]);


A=mSigIntegPerRunAllCellsSP;
M=[];
S=[];
for n=1:length(A);
M(n)=nanmean(A{n});
S(n)=nansem(A{n},1);
end

subplot(248)
errorbar([1:1:11],M,S)
[r,p]=corr([1:1:10]',M(2:end)');
xlabel('days')
ylabel('mean sig integ per run');
title(['mean sig integ per run p=',num2str(p)]);


saveas(gcf,'allSigTransAllCellsSP.fig');
save('allSigTransAllCellsSP.mat','mFreqAllCellsSP','mFreqRunAllCellsSP','mDisCoverAllCellsSP','mDisCoverRunAllCellsSP','mAmpAllCellsSP','mSigIntegAllCellsSP','mDurAllCellsSP','mSigIntegPerRunAllCellsSP');





%%
load('foldersAll.mat');

%speed scores
allSpeedScores=[];
p=pwd;

for n=1:length(foldersAll);
    disp(n)
    cd(foldersAll{n});
    load('commonCells.mat');
    load('useFolders');
    S=[];
    for m=1:length(useFolders);
        disp(m)
        a=commonCells(:,m);
        cd(useFolders{m});
        load('speed_dfof_sig\speed.mat');
        S(:,m)=speed.scores(a);
    end
    allSpeedScores(end+1:end+size(S,1),:)=S;
end
cd(p)
save('allSpeedScores.mat','allSpeedScores');
figure
errorbar([1:1:11],mean(allSpeedScores,1),nansem(allSpeedScores,1));
title('speed scores');
saveas(gcf,'speedScores.fig');

%% number of field over days
load('fieldDistri.mat');
NField={};
MField=[];
SField=[];
for n=1:length(fieldDistri);
    NField{n}=[];
    for m=1:length(fieldDistri{n});
        
NField{n}(m)=length(fieldDistri{n}{m});
    end
    MField(n,1)=nanmean(NField{n});
    SField(n,1)=nansem(NField{n},2);
end

figure,
errorbar([1:1:11],MField,SField)

load('fieldDistriUC.mat');
NFieldUC={};
MFieldUC=[];
SFieldUC=[];
for n=1:length(fieldDistriUC);
    NFieldUC{n}=[];
    for m=1:length(fieldDistriUC{n});
        
NFieldUC{n}(m)=length(fieldDistriUC{n}{m});
    end
    MFieldUC(n,1)=nanmean(NFieldUC{n});
    SFieldUC(n,1)=nansem(NFieldUC{n},2);
end

hold on
errorbar([1:1:11],MFieldUC,SFieldUC)
title('binsWithFieldsInCommonUnCommon');
legend('common','uncommon','Location','Northeast');
saveas(gcf,'binsWithFieldsInCommonUnCommon.fig');
save('NField.mat','NField');
save('MField.mat','MField');
save('SField.mat','SField');

save('NFieldUC.mat','NFieldUC');
save('MFieldUC.mat','MFieldUC');
save('SFieldUC.mat','SFieldUC');


%% FIELD DISTRIBUTION CORRELATION TO ENVIRONMENTS
load('fieldDistriEachBin.mat');
load('dfofSigNormAllFOVs.mat');

load('fieldDistriEachBinUC.mat');
load('corrAllFOVUncommon.mat');

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

corrTempRLC=[];
corrTempRC=[];
corrTempLC=[];
corrTempRLUC=[];
corrTempRUC=[];
corrTempLUC=[];

NCommon=size(dfofSigNormAllFOVs{1},1);
NUncommon=[];
for n=1:length(corrAllFOVUncommon.dfof{1});
    NUncommon(n)=length(corrAllFOVUncommon.dfof{1}{n});
end



for n=1:length(tempR);
    corrTempRLC(n)=corr(tempRL{n},fieldDistriEachBin(n,:)'/NCommon);
    corrTempRC(n)=corr(tempR{n},fieldDistriEachBin(n,:)'/NCommon);
    corrTempLC(n)=corr(tempL{n},fieldDistriEachBin(n,:)'/NCommon);
    corrTempRLUC(n)=corr(tempRL{n},fieldDistriEachBinUC(n,:)'/NUncommon(n));
    corrTempRUC(n)=corr(tempR{n},fieldDistriEachBinUC(n,:)'/NUncommon(n));
    corrTempLUC(n)=corr(tempL{n},fieldDistriEachBinUC(n,:)'/NUncommon(n));
end

figure,
subplot(131)
title('corr to temp RL');
plot(corrTempRLC,'r');
hold on
plot(corrTempRLUC,'k');
ylim([-0.5 0.4]);
title('L and R');
legend('common','uncommon','Location','Southeast');

subplot(132)
title('corr to temp R');
plot(corrTempRC,'r');
hold on
plot(corrTempRUC,'k');
ylim([-0.5 0.4]);
title('R');
legend('common','uncommon','Location','Southeast');

subplot(133)
title('corr to temp R');
plot(corrTempLC,'r');
hold on
plot(corrTempLUC,'k');
ylim([-0.5 0.4]);
title('L');
legend('common','uncommon','Location','Southeast');

saveas(gcf,'fieldCorrToEnv.fig');

save('environmentCorrelations.mat','corrTempRLC','corrTempLC','corrTempRC','corrTempRLUC','corrTempLUC','corrTempRUC');
save('NCommon.mat','NCommon');
save('NUncommon.mat','NUncommon');

%% correlation to template: individual cells
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

load('fieldDistri.mat');
fieldOnly={};%this is field only activity
for n=1:11;
    fieldOnly{n}=zeros(length(fieldDistri{1}),200);
end

for n=1:length(fieldDistri);
    for m=1:length(fieldDistri{n});
        a=fieldDistri{n}{m};
        if ~isnan(a);
        fieldOnly{n}(m,a)=1;
        end
    end
end

corrToEnvR=[];
corrToEnvL=[];
corrToEnvRL=[];

for n=1:length(fieldOnly);
    for m=1:size(fieldOnly{n},1);
        a=fieldOnly{n}(m,:)';
        corrToEnvR(n,m)=corr(tempR{n},a);
        corrToEnvL(n,m)=corr(tempL{n},a);
        corrToEnvRL(n,m)=corr(tempRL{n},a);
    end
end

save('fieldOnly.mat','fieldOnly');
save('corrToEnvR.mat','corrToEnvR');
save('corrToEnvRL.mat','corrToEnvRL');
save('corrToEnvL.mat','corrToEnvL');

load('fieldDistriUC.mat');
fieldOnlyUC={};%this is field only activity
for n=1:11;
    fieldOnlyUC{n}=zeros(length(fieldDistriUC{n}),200);
end

for n=1:length(fieldDistriUC);
    for m=1:length(fieldDistriUC{n});
        a=fieldDistriUC{n}{m};
        if ~isnan(a);
        fieldOnlyUC{n}(m,a)=1;
        end
    end
end

corrToEnvRUC={};
corrToEnvLUC={};
corrToEnvRLUC={};

for n=1:length(fieldOnlyUC);
    for m=1:size(fieldOnlyUC{n},1);
        a=fieldOnlyUC{n}(m,:)';
        corrToEnvRUC{n}(m)=corr(tempR{n},a);
        corrToEnvLUC{n}(m)=corr(tempL{n},a);
        corrToEnvRLUC{n}(m)=corr(tempRL{n},a);
    end
end

save('fieldOnlyUC.mat','fieldOnlyUC');
save('corrToEnvRUC.mat','corrToEnvRUC');
save('corrToEnvRLUC.mat','corrToEnvRLUC');
save('corrToEnvLUC.mat','corrToEnvLUC');

figure
subplot(131);
S=[];
M=[];
A=corrToEnvRLUC;
for n=1:length(A);
M(n)=nanmean(A{n},2);
S(n)=nansem(A{n},2);
end

A1=corrToEnvRL;
M1=nanmean(A1,2);
S1=nansem(A1,2);

errorbar([1:1:11],M1,S1);
hold on
errorbar([1:1:11],M,S);
title('RL');
ylim([-0.02 0.05]);
legend('common','uncommon','Location','Southeast');

subplot(132);
S=[];
M=[];
A=corrToEnvRUC;
for n=1:length(A);
M(n)=nanmean(A{n},2);
S(n)=nansem(A{n},2);
end

A1=corrToEnvR;
M1=nanmean(A1,2);
S1=nansem(A1,2);
errorbar([1:1:11],M1,S1);

hold on
errorbar([1:1:11],M,S);
title('R');
ylim([-0.02 0.05]);
legend('common','uncommon','Location','Southeast');

subplot(133);
S=[];
M=[];
A=corrToEnvLUC;
for n=1:length(A);
M(n)=nanmean(A{n},2);
S(n)=nansem(A{n},2);
end
A1=corrToEnvL;
M1=nanmean(A1,2);
S1=nansem(A1,2);
errorbar([1:1:11],M1,S1);

hold on
errorbar([1:1:11],M,S);
title('L');
ylim([-0.02 0.05]);
legend('common','uncommon','Location','Southeast');

saveas(gcf,'individualCellCorrToEnv.fig');


%% individual FOVs
corrTempLCAllFOVs=[];
corrTempRCAllFOVs=[];
corrTempRLCAllFOVs=[];

corrTempLUCAllFOVs=[];
corrTempRUCAllFOVs=[];
corrTempRLUCAllFOVs=[];
load('foldersAll.mat');
p=pwd;
for n=1:length(foldersAll);
    cd(foldersAll{n});
    load('fieldsCommonUncommon\environmentCorrelations.mat');
    corrTempLCAllFOVs(n,:)=corrTempLC;
    corrTempRCAllFOVs(n,:)=corrTempRC;
    corrTempRLCAllFOVs(n,:)=corrTempRLC;
    
    corrTempLUCAllFOVs(n,:)=corrTempLUC;
    corrTempRUCAllFOVs(n,:)=corrTempRUC;
    corrTempRLUCAllFOVs(n,:)=corrTempRLUC;
    
end
   
cd(p);

figure
subplot(131);
A=corrTempRLCAllFOVs;
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:11],M,S,'r');

A=corrTempRLUCAllFOVs;
M=nanmean(A,1);
S=nansem(A,1);
hold on
errorbar([1:1:11],M,S,'k');
title('RL')
legend('common','uncommon','Location','Southeast');

subplot(132);
A=corrTempRCAllFOVs;
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:11],M,S,'r');

A=corrTempRUCAllFOVs;
M=nanmean(A,1);
S=nansem(A,1);
hold on
errorbar([1:1:11],M,S,'k');
title('R')
legend('common','uncommon','Location','Southeast');

subplot(133);
A=corrTempLCAllFOVs;
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:11],M,S,'r');

A=corrTempLUCAllFOVs;
M=nanmean(A,1);
S=nansem(A,1);
hold on
errorbar([1:1:11],M,S,'k');
title('L')
legend('common','uncommon','Location','Southeast');
saveas(gcf,'individualFOVCorrToEnv.fig');
save('corrTempLCAllFOVs.mat','corrTempLCAllFOVs');
save('corrTempRCAllFOVs.mat','corrTempRCAllFOVs');
save('corrTempRLCAllFOVs.mat','corrTempRLCAllFOVs');

save('corrTempLUCAllFOVs.mat','corrTempLUCAllFOVs');
save('corrTempRUCAllFOVs.mat','corrTempRUCAllFOVs');
save('corrTempRLUCAllFOVs.mat','corrTempRLUCAllFOVs');

%% individual MICE
corrTempLCAllMice=[];
corrTempRCAllMice=[];
corrTempRLCAllMice=[];

corrTempLUCAllMice=[];
corrTempRUCAllMice=[];
corrTempRLUCAllMice=[];

mouseFolders={};
mouseFolders{1}='E:\learningAnalysis\IndividualMiceAnalysis\ID20210207';
mouseFolders{2}='E:\learningAnalysis\IndividualMiceAnalysis\ID20210519_1';
mouseFolders{3}='E:\learningAnalysis\IndividualMiceAnalysis\ID20210519_2';
mouseFolders{4}='E:\learningAnalysis\IndividualMiceAnalysis\ID20210811A';

save('mouseFolders.mat','mouseFolders');

p=pwd;
for n=1:length(mouseFolders);
    cd(mouseFolders{n});
    load('environmentCorrelations.mat');
    corrTempLCAllMice(n,:)=corrTempLC;
    corrTempRCAllMice(n,:)=corrTempRC;
    corrTempRLCAllMice(n,:)=corrTempRLC;
    
    corrTempLUCAllMice(n,:)=corrTempLUC;
    corrTempRUCAllMice(n,:)=corrTempRUC;
    corrTempRLUCAllMice(n,:)=corrTempRLUC;
    
end
   
cd(p);

figure
subplot(131);
A=corrTempRLCAllMice;
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:11],M,S,'r');

A=corrTempRLUCAllMice;
M=nanmean(A,1);
S=nansem(A,1);
hold on
errorbar([1:1:11],M,S,'k');
title('RL')
legend('common','uncommon','Location','Southeast');

subplot(132);
A=corrTempRCAllMice;
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:11],M,S,'r');

A=corrTempRUCAllMice;
M=nanmean(A,1);
S=nansem(A,1);
hold on
errorbar([1:1:11],M,S,'k');
title('R')
legend('common','uncommon','Location','Southeast');

subplot(133);
A=corrTempLCAllMice;
M=nanmean(A,1);
S=nansem(A,1);
errorbar([1:1:11],M,S,'r');

A=corrTempLUCAllMice;
M=nanmean(A,1);
S=nansem(A,1);
hold on
errorbar([1:1:11],M,S,'k');
title('L')
legend('common','uncommon','Location','Southeast');
saveas(gcf,'individualMiceCorrToEnv.fig');
save('corrTempLCAllMice.mat','corrTempLCAllMice');
save('corrTempRCAllMice.mat','corrTempRCAllMice');
save('corrTempRLCAllMice.mat','corrTempRLCAllMice');

save('corrTempLUCAllMice.mat','corrTempLUCAllMice');
save('corrTempRUCAllMice.mat','corrTempRUCAllMice');
save('corrTempRLUCAllMice.mat','corrTempRLUCAllMice');
%% gaining fields
load('fieldDistri.mat');
load('fieldDistriUC.mat');

fieldDistri_zeroOne={};
fieldDistriUC_zeroOne={};

for n=1:length(fieldDistri);
    fieldDistri_zeroOne{n}=zeros(length(fieldDistri{n}),200);
    fieldDistriUC_zeroOne{n}=zeros(length(fieldDistriUC{n}),200);
    for m=1:length(fieldDistri{n});
        i=fieldDistri{n}{m};
        if ~isnan(i);
        fieldDistri_zeroOne{n}(m,i)=1;
        end
    end
    for m=1:length(fieldDistriUC{n});
        i=fieldDistriUC{n}{m};
         if ~isnan(i);
        fieldDistriUC_zeroOne{n}(m,i)=1;
         end
    end
end

save('fieldDistri_zeroOne.mat','fieldDistri_zeroOne');
save('fieldDistriUC_zeroOne.mat','fieldDistriUC_zeroOne');

fieldDistriDiff=[];
b=[2 3];%day 1 and 2 in new env (2 and 3 in data)
a=[8 9 10 11];%AFter learning days
fieldDistriBeforeLearning=[];
fieldDistriAfterLearning=[];
for n=1:length(fieldDistri{1});
    beforeLearning=[];
    for i=1:length(b);
        beforeLearning(end+1,:)=fieldDistri_zeroOne{b(i)}(n,:);
    end
    afterLearning=[];
    for i=1:length(a);
        afterLearning(end+1,:)=fieldDistri_zeroOne{a(i)}(n,:);
    end
    fieldDistriDiff(n,:)=mean(afterLearning,1)-mean(beforeLearning);
    fieldDistriBeforeLearning(n,:)=mean(beforeLearning);
     fieldDistriAfterLearning(n,:)=mean(afterLearning);
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
%% activity correlatin of dfof to envirnoment
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

load('dfofSigNormAllFOVs.mat');
corrToEnvRDfof=[];
corrToEnvLDfof=[];
corrToEnvRLDfof=[];

for n=1:length(dfofSigNormAllFOVs);
    for m=1:size(dfofSigNormAllFOVs{n},1);
        a=dfofSigNormAllFOVs{n}(m,:)';
        corrToEnvRDfof(n,m)=corr(tempR{n},a);
        corrToEnvLDfof(n,m)=corr(tempL{n},a);
        corrToEnvRLDfof(n,m)=corr(tempRL{n},a);
    end
end

save('corrToEnvRDfof.mat','corrToEnvRDfof');
save('corrToEnvRLDfof.mat','corrToEnvRLDfof');
save('corrToEnvLDfof.mat','corrToEnvLDfof');

% 
% for n=1:length(fieldDistriUC);
%     for m=1:length(fieldDistriUC{n});
%         a=fieldDistriUC{n}{m};
%         if ~isnan(a);
%         fieldOnlyUC{n}(m,a)=1;
%         end
%     end
% end
% 
% corrToEnvRUC={};
% corrToEnvLUC={};
% corrToEnvRLUC={};
% 
% for n=1:length(fieldOnlyUC);
%     for m=1:size(fieldOnlyUC{n},1);
%         a=fieldOnlyUC{n}(m,:)';
%         corrToEnvRUC{n}(m)=corr(tempR{n},a);
%         corrToEnvLUC{n}(m)=corr(tempL{n},a);
%         corrToEnvRLUC{n}(m)=corr(tempRL{n},a);
%     end
% end

% save('fieldOnlyUC.mat','fieldOnlyUC');
% save('corrToEnvRUC.mat','corrToEnvRUC');
% save('corrToEnvRLUC.mat','corrToEnvRLUC');
% save('corrToEnvLUC.mat','corrToEnvLUC');


figure
subplot(131);
% S=[];
% M=[];
% A=corrToEnvRLUC;
% for n=1:length(A);
% M(n)=nanmean(A{n},2);
% S(n)=nansem(A{n},2);
% end

A1=corrToEnvRLDfof;
M1=nanmean(A1,2);
S1=nansem(A1,2);

errorbar([1:1:11],M1,S1);
% hold on
% errorbar([1:1:11],M,S);
title('RL');
ylim([-0.02 0.05]);
% legend('common','uncommon','Location','Southeast');

subplot(132);
% S=[];
% M=[];
% A=corrToEnvRUC;
% for n=1:length(A);
% M(n)=nanmean(A{n},2);
% S(n)=nansem(A{n},2);
% end

A1=corrToEnvRDfof;
M1=nanmean(A1,2);
S1=nansem(A1,2);
errorbar([1:1:11],M1,S1);

% hold on
% errorbar([1:1:11],M,S);
title('R');
ylim([-0.02 0.05]);
% legend('common','uncommon','Location','Southeast');

subplot(133);
% S=[];
% M=[];
% A=corrToEnvLUC;
% for n=1:length(A);
% M(n)=nanmean(A{n},2);
% S(n)=nansem(A{n},2);
% end
A1=corrToEnvLDfof;
M1=nanmean(A1,2);
S1=nansem(A1,2);
errorbar([1:1:11],M1,S1);

% hold on
% errorbar([1:1:11],M,S);
title('L');
ylim([-0.02 0.05]);
% legend('common','uncommon','Location','Southeast');

saveas(gcf,'individualCellCorrToEnv_dfof.fig');

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
