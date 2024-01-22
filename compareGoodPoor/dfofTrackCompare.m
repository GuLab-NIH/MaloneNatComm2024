
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\dfofTrack\meanAllDfof.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\dfofTrack\semAllDfof.mat');
meanG=meanAllDfof;
semG=semAllDfof;

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\dfofTrack\meanAllDfof.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\dfofTrack\semAllDfof.mat');
meanB=meanAllDfof;
semB=semAllDfof;

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\dfofTrack\meanAllDfof_sig.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\dfofTrack\semAllDfof_sig.mat');
meanGS=meanAllDfof_sig;
semGS=semAllDfof_sig;

load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\dfofTrack\meanAllDfof_sig.mat');
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\dfofTrack\semAllDfof_sig.mat');
meanBS=meanAllDfof_sig;
semBS=semAllDfof_sig;

load('E:\ID20211030\20211123Loc1Old\loc1Old\pcaica\cueAnalysisNew_sig\tempL.mat');
load('E:\ID20211030\20211123Loc1Old\loc1Old\pcaica\cueAnalysisNew_sig\tempR.mat');
tempLO=tempL;
tempRO=tempR;

load('E:\ID20211030\20211123New\loc1New\pcaica\cueAnalysisNew_sig\tempL.mat');
load('E:\ID20211030\20211123New\loc1New\pcaica\cueAnalysisNew_sig\tempR.mat');

%

figure
for n=1:11;
    subplot(4,3,n);
    a=meanG(n,:);
    b=semG(n,:);
    c=meanB(n,:);
    d=semB(n,:);
 
    if n==1;
        hold on
        plot([1:1:200],tempLO*max(a),'k');
        hold on
        plot([1:1:200],tempRO*max(a),'k');
        hold on
        line([7 7],[0 max(a)],'Color','b');
        hold on
        line([107 107],[0 max(a)],'Color','b');
    else
         hold on
        plot([1:1:200],tempL*max(a),'k');
        hold on
        plot([1:1:200],tempR*max(a),'k');
        hold on
        line([178 178],[0 max(a)],'Color','b');
    end
    hold on
    errorbar([1:1:length(a)],a,b,'g');
    hold on
    errorbar([1:1:length(a)],c,d,'m');
    xlim([1 200]);
    title(['day ',num2str(n-1)]);
end

saveas(gcf,'allCellsCompareDfof.fig');


figure
for n=1:11;
    subplot(4,3,n);
    a=meanGS(n,:);
    b=semGS(n,:);
    c=meanBS(n,:);
    d=semBS(n,:);
 
    if n==1;
        hold on
        plot([1:1:200],tempLO*max(a),'k');
        hold on
        plot([1:1:200],tempRO*max(a),'k');
        hold on
        line([7 7],[0 max(a)],'Color','b');
        hold on
        line([107 107],[0 max(a)],'Color','b');
    else
         hold on
        plot([1:1:200],tempL*max(a),'k');
        hold on
        plot([1:1:200],tempR*max(a),'k');
        hold on
        line([178 178],[0 max(a)],'Color','b');
    end
    hold on
    errorbar([1:1:length(a)],a,b,'g');
    hold on
    errorbar([1:1:length(a)],c,d,'m');
    xlim([1 200]);
    title(['day ',num2str(n-1)]);
end

saveas(gcf,'allCellsCompareDfof_sig.fig');

%
diffGB=[];
for n=1:size(meanG,1);
    diffGB(n,:)=(meanG(n,:)-meanB(n,:));
end

save('diffGB.mat','diffGB');

figure,
for n=1:11;
    subplot(4,3,n);
    a=diffGB(n,:);
    
    if n==1;
        hold on
        plot([1:1:200],tempLO*max(a),'k');
        hold on
        plot([1:1:200],tempRO*max(a),'k');
        hold on
        line([7 7],[0 max(a)],'Color','b');
        hold on
        line([107 107],[0 max(a)],'Color','b');
    else
         hold on
        plot([1:1:200],tempL*max(a),'k');
        hold on
        plot([1:1:200],tempR*max(a),'k');
        hold on
        line([178 178],[0 max(a)],'Color','b');
    end
    hold on
     plot([1:1:length(a)],a,'r');
    xlim([1 200]);
    title(['day ',num2str(n-1)]);
end
saveas(gcf,'allCellsCompareDfof_good_poor.fig');


%
diffGBS=[];
for n=1:size(meanGS,1);
    diffGBS(n,:)=(meanGS(n,:)-meanBS(n,:));
end

save('diffGBS.mat','diffGBS');

figure,
for n=1:11;
    subplot(4,3,n);
    a=diffGBS(n,:);
    
    if n==1;
        hold on
        plot([1:1:200],tempLO*max(a),'k');
        hold on
        plot([1:1:200],tempRO*max(a),'k');
        hold on
        line([7 7],[0 max(a)],'Color','b');
        hold on
        line([107 107],[0 max(a)],'Color','b');
    else
         hold on
        plot([1:1:200],tempL*max(a),'k');
        hold on
        plot([1:1:200],tempR*max(a),'k');
        hold on
        line([178 178],[0 max(a)],'Color','b');
    end
    hold on
     plot([1:1:length(a)],a,'r');
    xlim([1 200]);
    title(['day ',num2str(n-1)]);
end
saveas(gcf,'allCellsCompareDfof_sig_good_poor.fig');



%%

load('E:\ID20211030\20211123Loc1Old\loc1Old\pcaica\cueAnalysisNew_sig\tempL.mat');
load('E:\ID20211030\20211123Loc1Old\loc1Old\pcaica\cueAnalysisNew_sig\tempR.mat');
tempLO=tempL;
tempRO=tempR;
%EXPAND CUE + 2 BINS BOTH SIDES
bin=2;
r=contiguous(tempLO,1);
r=r{1,2};
for n=1:size(r,1)
    a=r(n,:);
    a(1)=a(1)-bin;
    a(2)=a(2)+bin;
    tempLO(a)=1;
end
r=contiguous(tempRO,1);
r=r{1,2};
for n=1:size(r,1)
    a=r(n,:);
    a(1)=a(1)-bin;
    a(2)=a(2)+bin;
    tempRO(a)=1;
end

%new env

load('E:\ID20211030\20211123New\loc1New\pcaica\cueAnalysisNew_sig\tempL.mat');
load('E:\ID20211030\20211123New\loc1New\pcaica\cueAnalysisNew_sig\tempR.mat');


%EXPAND CUE + 2 BINS BOTH SIDES
r=contiguous(tempL,1);
r=r{1,2};
for n=1:size(r,1)
    a=r(n,:);
    a(1)=a(1)-bin;
    a(2)=a(2)+bin;
    tempL(a)=1;
end
r=contiguous(tempR,1);
r=r{1,2};
for n=1:size(r,1)
    a=r(n,:);
    a(1)=a(1)-bin;
    a(2)=a(2)+bin;
    tempR(a)=1;
end


tempO=tempLO+tempRO;
inCueO=find(tempO>0);
outCueO=find(tempO==0);

inCueO=inCueO(inCueO<=196);
outCueO=outCueO(outCueO<=196); %for old env, since there are two reward, it is hard to remove all bin data after reward

% beforeRewardO=[4;5;6;7;104;105;106;107];
% inCueO=[inCueO;beforeRewardO];
% outCueO=setdiff(outCueO,beforeRewardO);

%for the new env,remove bins after reward

temp=tempL+tempR;
temp=temp(1:178);%remove everything after reward

inCue=find(temp>0);
outCue=find(temp==0);

inCue=inCue(inCue<=196);
outCue=outCue(outCue<=196);

beforeReward=[175;176;177;178];
% inCue=[inCue;beforeReward];
% outCue=setdiff(outCue,beforeReward);
% outCue=outCue(outCue<=178);
% %%%%%%%%%
% beforeReward=[175;176;177];
inCue=[inCue;beforeReward];
outCue=setdiff(outCue,[beforeReward;178]);
% %%%%%%%%%
%%
load('diffGB.mat');

% %use diffGB
inCueBinDiff={};%each bin, the above was each cue
outCueBinDiff={};

for n=1:11;
    inCueBinDiff{n}=[];
    outCueBinDiff{n}=[];
    if n==1;
        inCueBinDiff{n}=diffGB(n,inCueO);
        outCueBinDiff{n}=diffGB(n,outCueO);
    else
         inCueBinDiff{n}=diffGB(n,inCue);
        outCueBinDiff{n}=diffGB(n,outCue);
    end
end
inCueNoR=setdiff(inCue,beforeReward);

save
%use diffGB
inCueNoRBinDiff={};%each bin, the above was each cue
beforeRBinDiff={};

for n=1:11;
    inCueNoRBinDiff{n}=[];
    beforeRBinDiff{n}=[];

         inCueNoRBinDiff{n}=diffGB(n,inCueNoR);
        beforeRBinDiff{n}=diffGB(n,beforeReward);
  
end

M=[];
E=[];
for n=1:11;
    M(n,3)=mean(beforeRBinDiff{n});
    M(n,2)=mean(inCueNoRBinDiff{n});
    M(n,1)=mean(outCueBinDiff{n});
    E(n,3)=nansem(beforeRBinDiff{n},2);
    E(n,2)=nansem(inCueNoRBinDiff{n},2);
    E(n,1)=nansem(outCueBinDiff{n},2);
end

sig=[];

for n=1:11;
    [~,sig(n,1)]=ttest2(outCueBinDiff{n},inCueNoRBinDiff{n});
    [~,sig(n,2)]=ttest2(outCueBinDiff{n},beforeRBinDiff{n});
end
 errorbarGroups(M(2:end,:),E(2:end,:),0)
 title('remove bins after reward dfof')
 saveas(gcf,'rewardInCueOutCueNewEnvTogether_removeBinAfterReward');
 
  sigAll=[];
 for n=1:11;
     [~,sigAll(n,1)]=ttest2(beforeRBinDiff{n},inCueNoRBinDiff{n});
     [~,sigAll(n,2)]=ttest2(outCueBinDiff{n},inCueNoRBinDiff{n});
 end
 
 
 
 %%
load('diffGBS.mat');

% %use diffGB
inCueBinDiff={};%each bin, the above was each cue
outCueBinDiff={};

for n=1:11;
    inCueBinDiff{n}=[];
    outCueBinDiff{n}=[];
    if n==1;
        inCueBinDiff{n}=diffGBS(n,inCueO);
        outCueBinDiff{n}=diffGBS(n,outCueO);
    else
         inCueBinDiff{n}=diffGBS(n,inCue);
        outCueBinDiff{n}=diffGBS(n,outCue);
    end
end
inCueNoR=setdiff(inCue,beforeReward);

save('zones.mat','inCueNoR','beforeReward','outCue')
%use diffGB
inCueNoRBinDiff={};%each bin, the above was each cue
beforeRBinDiff={};

for n=1:11;
    inCueNoRBinDiff{n}=[];
    beforeRBinDiff{n}=[];

         inCueNoRBinDiff{n}=diffGBS(n,inCueNoR);
        beforeRBinDiff{n}=diffGBS(n,beforeReward);
  
end

M=[];
E=[];
for n=1:11;
    M(n,3)=mean(beforeRBinDiff{n});
    M(n,2)=mean(inCueNoRBinDiff{n});
    M(n,1)=mean(outCueBinDiff{n});
    E(n,3)=nansem(beforeRBinDiff{n},2);
    E(n,2)=nansem(inCueNoRBinDiff{n},2);
    E(n,1)=nansem(outCueBinDiff{n},2);
end

sig=[];

for n=1:11;
    [~,sig(n,1)]=ttest2(outCueBinDiff{n},inCueNoRBinDiff{n});
    [~,sig(n,2)]=ttest2(outCueBinDiff{n},beforeRBinDiff{n});
end
 errorbarGroups(M(2:end,:),E(2:end,:),0)
 title('remove bins after reward dfof sig')
%  ylim([-0.028 0.015])
%  breakyaxis([-0.0228 -0.0052]);
%for paper data
for n=1:length(outCueBinDiff);
    outCueBinDiff{n}=outCueBinDiff{n}';
    inCueNoRBinDiff{n}=inCueNoRBinDiff{n}';
    beforeRBinDiff{n}=beforeRBinDiff{n}';
end

save('outCueBinDiff.mat','outCueBinDiff');
save('inCueNoRBinDiff.mat','inCueNoRBinDiff');
save('beforeRBinDiff.mat','beforeRBinDiff');
%change it back
for n=1:length(outCueBinDiff);
    outCueBinDiff{n}=outCueBinDiff{n}';
    inCueNoRBinDiff{n}=inCueNoRBinDiff{n}';
    beforeRBinDiff{n}=beforeRBinDiff{n}';
end
 saveas(gcf,'rewardInCueOutCueNewEnvTogether_SIG_removeBinAfterReward');
 save('sig.mat','sig')
  sigAll=[];
 for n=1:11;
     [~,sigAll(n,1)]=ttest2(beforeRBinDiff{n},inCueNoRBinDiff{n});
     [~,sigAll(n,2)]=ttest2(outCueBinDiff{n},inCueNoRBinDiff{n});
 end

 %put all days together
 
 beforeRBinDiffAll=cell2mat(beforeRBinDiff(2:end));
 outCueBinDiffAll=cell2mat(outCueBinDiff(2:end));
 inCueNoRBinDiffAll=cell2mat(inCueNoRBinDiff(2:end));
 M=[];
 M(1)=nanmean(outCueBinDiffAll);
  M(2)=nanmean(inCueNoRBinDiffAll);
   M(3)=nanmean(beforeRBinDiffAll);
    E=[];
 E(1)=nansem(outCueBinDiffAll,2);
  E(2)=nansem(inCueNoRBinDiffAll,2);
   E(3)=nansem(beforeRBinDiffAll,2);
figure
 bar([1:1:3],M)
 hold on
 errorbar([1:1:3],M,E,'.')
 
 [~,sigAllDay1]=ttest2( outCueBinDiffAll, inCueNoRBinDiffAll);
 [~,sigAllDay2]=ttest2( outCueBinDiffAll,  beforeRBinDiffAll);
 outCueBinDiffAll=outCueBinDiffAll';
 inCueNoRBinDiffAll=inCueNoRBinDiffAll';
 beforeRBinDiffAll=beforeRBinDiffAll';
 
 save('sigAllDay1.mat','sigAllDay1')
  save('sigAllDay2.mat','sigAllDay2')
  save('outCueBinDiffAll.mat','outCueBinDiffAll');
  save('inCueNoRBinDiffAll.mat','inCueNoRBinDiffAll');
  save('beforeRBinDiffAll.mat','beforeRBinDiffAll');
 saveas(gcf,'AllDays_rewardInCueOutCueNewEnvTogether_SIG_removeBinAfterReward.fig')
 
 %% see whether reward activity was larger than zero
 [r,p]=ttest(beforeRBinDiffAll,0)
%% 
  %% separately plot good and poor in bar graph sig trans
 
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\dfofTrack\allCellsDfof_sig.mat');
G=allCellsDfof_sig;
 
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\dfofTrack\allCellsDfof_sig.mat');
B=allCellsDfof_sig;

GG={};
BB={};
for n=1:11;
    GG{n}=[];
    BB{n}=[];
GG{n}(:,1)=nanmean(G{n}(:,outCue),2);
GG{n}(:,2)=nanmean(G{n}(:,inCueNoR),2);
GG{n}(:,3)=nanmean(G{n}(:,beforeReward),2);

BB{n}(:,1)=nanmean(B{n}(:,outCue),2);
BB{n}(:,2)=nanmean(B{n}(:,inCueNoR),2);
BB{n}(:,3)=nanmean(B{n}(:,beforeReward),2);
end

M=[];
E=[];
for n=1:11;
    M(n,:)=nanmean(GG{n},1);
     E(n,:)=nansem(GG{n},1);
end

figure,
subplot(211)
 errorbarGroups(M(2:end,:),E(2:end,:),1)
 title('DFOFSIG  good')
  ylim([0 0.05])

 M=[];
E=[];
for n=1:11;
    M(n,:)=nanmean(BB{n},1);
     E(n,:)=nansem(BB{n},1);
end

subplot(212)
 errorbarGroups(M(2:end,:),E(2:end,:),1)
 title('DFOFSIG poor')
 ylim([0 0.05])
 saveas(gcf,'dfofSig_goodAndPoor');


