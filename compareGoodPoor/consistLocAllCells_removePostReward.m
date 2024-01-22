% load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\runingLocationConsistency\meanAllConsLoc.mat');
% load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\runingLocationConsistency\semAllConsLoc.mat');
% meanG=meanAllConsLoc;
% semG=semAllConsLoc;
% 
% load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\runingLocationConsistency\meanAllConsLoc.mat');
% load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\runingLocationConsistency\semAllConsLoc.mat');
% meanB=meanAllConsLoc;
% semB=semAllConsLoc;
% 
% 
% load('E:\ID20211030\20211123Loc1Old\loc1Old\pcaica\cueAnalysisNew_sig\tempL.mat');
% load('E:\ID20211030\20211123Loc1Old\loc1Old\pcaica\cueAnalysisNew_sig\tempR.mat');
% tempLO=tempL;
% tempRO=tempR;
% 
% load('E:\ID20211030\20211123New\loc1New\pcaica\cueAnalysisNew_sig\tempL.mat');
% load('E:\ID20211030\20211123New\loc1New\pcaica\cueAnalysisNew_sig\tempR.mat');
% 
% %% 
% 
% figure
% for n=1:11;
%     subplot(4,3,n);
%     a=meanG(n,:);
%     b=semG(n,:);
%     c=meanB(n,:);
%     d=semB(n,:);
%  
%     if n==1;
%         hold on
%         plot([1:1:200],tempLO*max(a),'k');
%         hold on
%         plot([1:1:200],tempRO*max(a),'k');
%         hold on
%         line([7 7],[0 max(a)],'Color','b');
%         hold on
%         line([107 107],[0 max(a)],'Color','b');
%     else
%          hold on
%         plot([1:1:200],tempL*max(a),'k');
%         hold on
%         plot([1:1:200],tempR*max(a),'k');
%         hold on
%         line([178 178],[0 max(a)],'Color','b');
%     end
%     hold on
%     errorbar([1:1:length(a)],a,b,'g');
%     hold on
%     errorbar([1:1:length(a)],c,d,'m');
%     xlim([1 200]);
%     title(['day ',num2str(n-1)]);
% end
% 
% saveas(gcf,'allCellsCompare.fig');
% 
% %% compare cue zone
% diffGB=[];
% for n=1:size(meanG,1);
%     diffGB(n,:)=(meanG(n,:)-meanB(n,:));
% end
% 
% save('diffGB.mat','diffGB');
% 
% figure,
% for n=1:11;
%     subplot(4,3,n);
%     a=diffGB(n,:);
%     
%     if n==1;
%         hold on
%         plot([1:1:200],tempLO*max(a),'k');
%         hold on
%         plot([1:1:200],tempRO*max(a),'k');
%         hold on
%         line([7 7],[0 max(a)],'Color','b');
%         hold on
%         line([107 107],[0 max(a)],'Color','b');
%     else
%          hold on
%         plot([1:1:200],tempL*max(a),'k');
%         hold on
%         plot([1:1:200],tempR*max(a),'k');
%         hold on
%         line([178 178],[0 max(a)],'Color','b');
%     end
%     hold on
%      plot([1:1:length(a)],a,'r');
%     xlim([1 200]);
%     title(['day ',num2str(n-1)]);
% end
% saveas(gcf,'allCellsCompare_good_poor.fig');




%% within and outside consistency differences: all points within and outside cue and before reward
load('E:\ID20211030\20211123Loc1Old\loc1Old\pcaica\cueAnalysisNew_sig\tempL.mat');
load('E:\ID20211030\20211123Loc1Old\loc1Old\pcaica\cueAnalysisNew_sig\tempR.mat');
tempLO=tempL;
tempRO=tempR;

%new env

load('E:\ID20211030\20211123New\loc1New\pcaica\cueAnalysisNew_sig\tempL.mat');
load('E:\ID20211030\20211123New\loc1New\pcaica\cueAnalysisNew_sig\tempR.mat');

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
inCue=[inCue;beforeReward];
outCue=setdiff(outCue,beforeReward);
save('beforeReward.mat','beforeReward')
save('outCue.mat','outCue')
save('outCueO.mat','outCueO')
% load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\runingLocationConsistency\allConsLoc.mat');
% G=allConsLoc;
% 
% load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\runingLocationConsistency\allConsLoc.mat');
% B=allConsLoc;
% 
% inCueCorrG={};
% inCueCorrMG=[];
% inCueCorrSEMG=[];
% 
% inCueCorrB={};
% inCueCorrMB=[];
% inCueCorrSEMB=[];
% 
% outCueCorrG={};
% outCueCorrMG=[];
% outCueCorrSEMG=[];
% 
% outCueCorrB={};
% outCueCorrMB=[];
% outCueCorrSEMB=[];
% 
% for n=1:length(G);
%     if n==1;
%         %familar env
%     inCueCorrG{n}=G{n}(:,inCueO);
%     g=reshape(inCueCorrG{n},[1 size(inCueCorrG{n},1)*size(inCueCorrG{n},2)]);
%     inCueCorrMG(n)=nanmean(g,2);
%     inCueCorrSEMG(n)=nansem(g,2);
%     
%      inCueCorrB{n}=B{n}(:,inCueO);
%     b=reshape(inCueCorrB{n},[1 size(inCueCorrB{n},1)*size(inCueCorrB{n},2)]);
%     inCueCorrMB(n)=nanmean(b,2);
%     inCueCorrSEMB(n)=nansem(b,2);
%     
%      outCueCorrG{n}=G{n}(:,outCueO);
%     g=reshape(outCueCorrG{n},[1 size(outCueCorrG{n},1)*size(outCueCorrG{n},2)]);
%     outCueCorrMG(n)=nanmean(g,2);
%     outCueCorrSEMG(n)=nansem(g,2);
%     
%      outCueCorrB{n}=B{n}(:,outCueO);
%     b=reshape(outCueCorrB{n},[1 size(outCueCorrB{n},1)*size(outCueCorrB{n},2)]);
%     outCueCorrMB(n)=nanmean(b,2);
%     outCueCorrSEMB(n)=nansem(b,2);
%     else
%                 %novel env
%   inCueCorrG{n}=G{n}(:,inCue);
%     g=reshape(inCueCorrG{n},[1 size(inCueCorrG{n},1)*size(inCueCorrG{n},2)]);
%     inCueCorrMG(n)=nanmean(g,2);
%     inCueCorrSEMG(n)=nansem(g,2);
%     
%      inCueCorrB{n}=B{n}(:,inCue);
%     b=reshape(inCueCorrB{n},[1 size(inCueCorrB{n},1)*size(inCueCorrB{n},2)]);
%     inCueCorrMB(n)=nanmean(b,2);
%     inCueCorrSEMB(n)=nansem(b,2);
%     
%      outCueCorrG{n}=G{n}(:,outCue);
%     g=reshape(outCueCorrG{n},[1 size(outCueCorrG{n},1)*size(outCueCorrG{n},2)]);
%     outCueCorrMG(n)=nanmean(g,2);
%     outCueCorrSEMG(n)=nansem(g,2);
%     
%      outCueCorrB{n}=B{n}(:,outCue);
%     b=reshape(outCueCorrB{n},[1 size(outCueCorrB{n},1)*size(outCueCorrB{n},2)]);
%     outCueCorrMB(n)=nanmean(b,2);
%     outCueCorrSEMB(n)=nansem(b,2);
%     end
% end
% 
% %%
% 
% figure, 
% 
% for n=1:11;
%     A=[inCueCorrMG(n) inCueCorrMB(n) outCueCorrMG(n) outCueCorrMB(n)];
%      B=[inCueCorrSEMG(n) inCueCorrSEMB(n) outCueCorrSEMG(n) outCueCorrSEMB(n)]; 
%      subplot(3,4,n);
%      bar([1:1:4],A);
%      hold on
%      errorbar([1:1:4],A,B,'.');
%      title(['day',num2str(n-1)]);
% end
% 
% saveas(gcf,'inOutCueDiff.fig');%error bar is all values within the cue
% 
% 
% 
% %% calculate for each cue
% 
% load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\runingLocationConsistency\allConsLoc.mat');
% G=allConsLoc;
% 
% load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\runingLocationConsistency\allConsLoc.mat');
% B=allConsLoc;
% 
% 
% load('E:\ID20211030\20211123Loc1Old\loc1Old\pcaica\cueAnalysisNew_sig\tempL.mat');
% load('E:\ID20211030\20211123Loc1Old\loc1Old\pcaica\cueAnalysisNew_sig\tempR.mat');
% tempLO=tempL;
% tempRO=tempR;
% tempO=tempLO+tempRO;
% tempO(beforeRewardO)=1;
% tempO=tempO(1:196);
% 
% load('E:\ID20211030\20211123New\loc1New\pcaica\cueAnalysisNew_sig\tempL.mat');
% load('E:\ID20211030\20211123New\loc1New\pcaica\cueAnalysisNew_sig\tempR.mat');
% temp=tempL+tempR;
% temp(beforeReward)=1;
% temp=temp(1:196);
% 
% r=contiguous(tempO,[1]);
% r=r{1,2};
% tempOInCueP={};%PIECES for in cue areas old env
% for n=1:size(r,1);
% tempOInCueP{n}=[r(n,1):1:r(n,2)];
% end
% 
% r=contiguous(temp,[1]);
% r=r{1,2};
% tempInCueP={};%PIECES for in cue areas new env
% for n=1:size(r,1);
% tempInCueP{n}=[r(n,1):1:r(n,2)];
% end
% 
% 
% r=contiguous(tempO,[0]);
% r=r{1,2};
% tempOOutCueP={};%PIECES for in cue areas old env
% for n=1:size(r,1);
% tempOOutCueP{n}=[r(n,1):1:r(n,2)];
% end
% 
% r=contiguous(temp,[0]);
% r=r{1,2};
% tempOutCueP={};%PIECES for in cue areas new env
% for n=1:size(r,1);
% tempOutCueP{n}=[r(n,1):1:r(n,2)];
% end
% 
% %for each day, we collect n number of period to compare activity in good
% %and poor performers
% 
% inCuePDiff=[];
% outCuePDiff=[];
% 
% for n=1:11;
%     inCuePDiff{n}=[];
%     outCuePDiff{n}=[];
%     if n==1;
%         for m=1:length(tempOInCueP);
%             g=G{n}(:,tempOInCueP{m});
%             inCuePDiff{n}(m,1)=nanmean(nanmean(g));
%             b=B{n}(:,tempOInCueP{m});
%             inCuePDiff{n}(m,2)=nanmean(nanmean(b));
%         end
%         
%         for m=1:length(tempOOutCueP);
%             g=G{n}(:,tempOOutCueP{m});
%             outCuePDiff{n}(m,1)=nanmean(nanmean(g));
%             b=B{n}(:,tempOOutCueP{m});
%             outCuePDiff{n}(m,2)=nanmean(nanmean(b));
%         end
%     else
%          for m=1:length(tempInCueP);
%             g=G{n}(:,tempInCueP{m});
%             inCuePDiff{n}(m,1)=nanmean(nanmean(g));
%             b=B{n}(:,tempInCueP{m});
%             inCuePDiff{n}(m,2)=nanmean(nanmean(b));
%         end
%         
%         for m=1:length(tempOutCueP);
%             g=G{n}(:,tempOutCueP{m});
%             outCuePDiff{n}(m,1)=nanmean(nanmean(g));
%             b=B{n}(:,tempOutCueP{m});
%             outCuePDiff{n}(m,2)=nanmean(nanmean(b));
%         end
%     end
%     
% end
    %%        
        
%the error bar is different cues: compare the mean at cue by cue: less
%significant
% figure
% for n=1:length(inCuePDiff);
%     
%     subplot(3,4,n)
%    in=inCuePDiff{n}(:,1)-inCuePDiff{n}(:,2);
%    out=outCuePDiff{n}(:,1)-outCuePDiff{n}(:,2);
%    [r,p]=ttest2(in,out);
%    bar([1 2],[mean(in) mean(out)]);
%    hold on
%    errorbar([1 2],[mean(in) mean(out)],[nansem(in,1) nansem(in,1)],'.');
%    title(['day',num2str(n-1),' p=',num2str(p)]);
% end
% 
% saveas(gcf,'diffInCueOutCue_errorBarIsEachPeriod.fig');
% 
%% compare difference: difference within cue bins and outside of cue bins: errorbar is each bin. Use the above data
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

% save('inCueO.mat','inCueO');
% save('outCueO.mat','outCueO');
% save('inCue.mat','inCue');
% save('outCue.mat','outCue');
% 
%     
% figure
% for n=1:length(inCueBinDiff);
%     
%     subplot(3,4,n)
%    in=inCueBinDiff{n};
%    out=outCueBinDiff{n};
%    [r,p]=ttest2(in,out);
%    bar([1 2],[mean(in) mean(out)]);
%    hold on
%    errorbar([1 2],[mean(in) mean(out)],[nansem(in,2) nansem(out,2)],'.');
%    title(['day',num2str(n-1),' p=',num2str(p)]);
% end
% 
% saveas(gcf,'diffInCueOutCue_errorBarIsAllBins.fig');
% 
%% plotting above together (no old env)

% M=[];
% E=[];
% for n=1:11;
%     M(n,1)=mean(inCueBinDiff{n});
%     M(n,2)=mean(outCueBinDiff{n});
%     E(n,1)=nansem(inCueBinDiff{n},2);
%     E(n,2)=nansem(outCueBinDiff{n},2);
% end
% 
%  errorbarGroups(M(2:end,:),E(2:end,:))
%  
%  saveas(gcf,'inOutCueNewEnvTogether');
%% separating reward zone and cues
inCueNoR=setdiff(inCue,beforeReward);
save('inCueNoR.mat','inCueNoR')
%use diffGB
inCueNoRBinDiff={};%each bin, the above was each cue
beforeRBinDiff={};

for n=1:11;
    inCueNoRBinDiff{n}=[];
    beforeRBinDiff{n}=[];

         inCueNoRBinDiff{n}=diffGB(n,inCueNoR);
        beforeRBinDiff{n}=diffGB(n,beforeReward);
  
end

save('beforeRBinDiff.mat','beforeRBinDiff')
save('inCueNoRBinDiff.mat','inCueNoRBinDiff');
save('outCueBinDiff.mat','outCueBinDiff')
% figure
% for n=1:length(inCueNoRBinDiff);
%     
%     subplot(3,4,n)
%    in=inCueNoRBinDiff{n};
%    out=beforeRBinDiff{n};
%    [r,p]=ttest2(in,out);
%    bar([1 2],[mean(in) mean(out)]);
%    hold on
%    errorbar([1 2],[mean(in) mean(out)],[nansem(in,2) nansem(out,2)],'.');
%    title(['day',num2str(n),' p=',num2str(p)]);
% end
% 
% saveas(gcf,'diffCueAndReward_errorBarIsAllBins.fig');
% 
% M=[];
% E=[];
% for n=1:11;
%     M(n,1)=mean(inCueNoRBinDiff{n});
%     M(n,2)=mean(beforeRBinDiff{n});
%     E(n,1)=nansem(inCueNoRBinDiff{n},2);
%     E(n,2)=nansem(beforeRBinDiff{n},2);
% end
% 
%  errorbarGroups(M(2:end,:),E(2:end,:),0)
%  saveas(gcf,'diffCueAndRewardNewEnvTogether');
 
 %% putting inCueNoReward, outCue, reward together
%  figure
% for n=1:length(inCueNoRBinDiff);
%     
%    subplot(3,4,n)
%    r=beforeRBinDiff{n};
%    in=inCueNoRBinDiff{n};
%    out=outCueBinDiff{n};
%    [~,p]=ttest2(in,out);
%    bar([1 2 3],[mean(r) mean(in) mean(out)]);
%    hold on
%    errorbar([1 2 3],[mean(r) mean(in) mean(out)],[nansem(r,2) nansem(in,2) nansem(out,2)],'.');
%    title(['day',num2str(n-1),' p=',num2str(p)]);
% end
% 
% saveas(gcf,'rewardInCueOutCue_errorBarIsAllBins.fig');

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

save('sig.mat','sig')
 errorbarGroups(M(2:end,:),E(2:end,:),0)
 title('remove bins after reward')
 saveas(gcf,'rewardInCueOutCueNewEnvTogether_removeBinAfterReward');
 
 %% siginficance in this plot
 sigAll=[];
 for n=1:11;
     [~,sigAll(n,1)]=ttest2(beforeRBinDiff{n},inCueNoRBinDiff{n});
     [~,sigAll(n,2)]=ttest2(outCueBinDiff{n},inCueNoRBinDiff{n});
 end
 
%% putting all days together

%put all days together
 
 beforeRBinDiffAll=cell2mat(beforeRBinDiff(2:11));
 outCueBinDiffAll=cell2mat(outCueBinDiff(2:11));
 inCueNoRBinDiffAll=cell2mat(inCueNoRBinDiff(2:11));
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
 save('sigAllDay1.mat','sigAllDay1')
  save('sigAllDay2.mat','sigAllDay2')
 saveas(gcf,'AllDays_rewardInCueOutCueNewEnvTogether_removeBinAfterReward.fig')

 save('outCueBinDiffAll.mat','outCueBinDiffAll')
  save('inCueNoRBinDiffAll.mat','inCueNoRBinDiffAll')
  save('beforeRBinDiffAll.mat','beforeRBinDiffAll')
%%
  %% separately plot good and poor in bar graph sig trans
 
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\runingLocationConsistency\allConsLoc.mat');
G=allConsLoc;
 
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\runingLocationConsistency\allConsLoc.mat');
B=allConsLoc;

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
 title('RBR  good')
 ylim([0 0.25])

 M=[];
E=[];
for n=1:11;
    M(n,:)=nanmean(BB{n},1);
     E(n,:)=nansem(BB{n},1);
end

subplot(212)
 errorbarGroups(M(2:end,:),E(2:end,:),1)
 title('RBR poor')
 ylim([0 0.25])
 saveas(gcf,'RBR_goodAndPoor_removePostReward');

%  %% calculate the fractional changes
%  
%  %% compare cue zone
% diffGBF=[];
% for n=1:size(meanG,1);
%     diffGBF(n,:)=(meanG(n,:)-meanB(n,:))./meanB(n,:);
% end
% 
% save('diffGBF.mat','diffGBF');
% 
% figure,
% for n=1:11;
%     subplot(4,3,n);
%     a=diffGBF(n,:);
%     
%     if n==1;
%         hold on
%         plot([1:1:200],tempLO*max(a),'k');
%         hold on
%         plot([1:1:200],tempRO*max(a),'k');
%         hold on
%         line([7 7],[0 max(a)],'Color','b');
%         hold on
%         line([107 107],[0 max(a)],'Color','b');
%     else
%          hold on
%         plot([1:1:200],tempL*max(a),'k');
%         hold on
%         plot([1:1:200],tempR*max(a),'k');
%         hold on
%         line([178 178],[0 max(a)],'Color','b');
%     end
%     hold on
%      plot([1:1:length(a)],a,'r');
%     xlim([1 200]);
%     title(['day ',num2str(n-1)]);
% end
% saveas(gcf,'allCellsCompare_good_poorFraction.fig');
% 
%  
