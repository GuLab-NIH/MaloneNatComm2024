load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\indicesAllNewCueThresh.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\corrAllFOVs.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\dfofSigNormAllFOVs.mat');

i1=indicesAllNewCueThresh.CueRIdx;
i2=indicesAllNewCueThresh.CueLIdx;
i=[i1;i2];
corrAllCueRL.dfof{1}=corrAllFOVs.dfof{1}(i,:);
corrAllCueRL.dfof{2}=corrAllFOVs.dfof{2}(i,:);
corrAllCueRL.dfof{3}=corrAllFOVs.dfof{3}(i,:);

corrAllCueRL.dfofSig{1}=corrAllFOVs.dfofSig{1}(i,:);
corrAllCueRL.dfofSig{2}=corrAllFOVs.dfofSig{2}(i,:);


corrAllCueRL.dfofSig{3}=corrAllFOVs.dfofSig{3}(i,:);

corrAllCueRL.dfofMeanToMean=mean(corrAllCueRL.dfof{1},1);
corrAllCueRL.dfofMeanToNext=mean(corrAllCueRL.dfof{2},1);
corrAllCueRL.dfofMeanToOthers=mean(corrAllCueRL.dfof{3},1);

corrAllCueRL.dfofSigMeanToMean=mean(corrAllCueRL.dfofSig{1},1);
corrAllCueRL.dfofSigMeanToNext=mean(corrAllCueRL.dfofSig{2},1);
corrAllCueRL.dfofSigMeanToOthers=mean(corrAllCueRL.dfofSig{3},1);

save('corrAllCueRL.mat','corrAllCueRL');
        
figure
subplot(141)
title('dfof');
plot(corrAllCueRL.dfofMeanToMean);
hold on
plot(corrAllCueRL.dfofMeanToNext);
hold on
plot(corrAllCueRL.dfofMeanToOthers);
legend('mean','next','others','Location','SouthEast');

subplot(142)
title('dfof sig');
plot(corrAllCueRL.dfofSigMeanToMean);
hold on
plot(corrAllCueRL.dfofSigMeanToNext);
hold on
plot(corrAllCueRL.dfofSigMeanToOthers);
legend('mean','next','others','Location','SouthEast');

%only use dfofSig Mean To Mean and do shuffle
A=corrAllCueRL.dfofSig{1};
MTempRL=corrAllCueRL.dfofSigMeanToMean;
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

title(['sig meanToMean,p=',num2str(p)]);

sigToNDay1=[];%significance compared to day1 in new env
day1=1;%day1 is the first column of "A".
A=corrAllCueRL.dfofSig{1};

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
saveas(gcf,'corrMeanCueRL.fig');

%% PLOTTING MATRIX CORRELATION: FIRST RIGHT, SECOND LEFT
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\cueCellRNewCueThresh\dfofSigNormCueRSort_Max1.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\cueCellLNewCueThresh\dfofSigNormCueLSort_Max1.mat');

dfofSigNormCueRLSort_Max1={};
for n=1:length(dfofSigNormCueRSort_Max1);
    dfofSigNormCueRLSort_Max1{n}=[dfofSigNormCueRSort_Max1{n};dfofSigNormCueLSort_Max1{n}];
end
save('dfofSigNormCueRLSort_Max1.mat','dfofSigNormCueRLSort_Max1');
matrixCorr=[]; %matrix n to matrix n+1
a=dfofSigNormCueRLSort_Max1;
for n=1:length(a)-1;
    matrixCorr(n,1)=corr2(a{n},a{n+1});
end
save('matrixCorr.mat','matrixCorr');

figure,plot(matrixCorr,'r');
hold on
plot(matrixCorr,'r.','MarkerSize',10);
[r,p]=corr([1:1:9]',matrixCorr(2:end));
title(['matrix corr p=',num2str(p)]);

saveas(gcf,'matrixCorr_max1.fig');

%% individual corr
load('dfofSigNormCueRLSort_Max1.mat');
A=dfofSigNormCueRLSort_Max1;
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
