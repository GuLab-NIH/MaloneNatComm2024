cd ..\
load('foldersAllIndvFOV.mat');
cd('runingLocationConsistency');

load('E:\ID20211030\20211123Loc1Old\loc1Old\pcaica\cueAnalysisNew_sig\tempL.mat');
load('E:\ID20211030\20211123Loc1Old\loc1Old\pcaica\cueAnalysisNew_sig\tempR.mat');
tempLO=tempL;
tempRO=tempR;

load('E:\ID20211030\20211123New\loc1New\pcaica\cueAnalysisNew_sig\tempL.mat');
load('E:\ID20211030\20211123New\loc1New\pcaica\cueAnalysisNew_sig\tempR.mat');

%% all cells

allConsLoc={};%each cell is one day
for n=1:11;
    allConsLoc{n}=[];
end
p=pwd;
for n=1:length(foldersAllIndvFOV);
    disp(n)
    for m=1:length(foldersAllIndvFOV{n});
        cd(foldersAllIndvFOV{n}{m})
        load('RunByRun_dfof\corrInfoLocation.mat');
        c=corrInfoLocation.toOthersMean;
        allConsLoc{m}(end+1:end+size(c,1),:)=c;
    end
end
cd(p);

save('allConsLoc.mat','allConsLoc');

%% 
meanAllConsLoc=[];
semAllConsLoc=[];

for n=1:length(allConsLoc);
meanAllConsLoc(n,:)=nanmean(allConsLoc{n},1);
semAllConsLoc(n,:)=nansem(allConsLoc{n},1); 
end

save('meanAllConsLoc.mat','meanAllConsLoc');
save('semAllConsLoc.mat','semAllConsLoc');

%% plot

figure
for n=1:11;
    subplot(4,3,n);
    a=meanAllConsLoc(n,:);
    b=semAllConsLoc(n,:);
    errorbar([1:1:length(a)],a,b,'k');
    if n==1;
        hold on
        plot([1:1:200],tempLO*max(a),'m');
        hold on
        plot([1:1:200],tempRO*max(a),'g');
        hold on
        line([7 7],[0 max(a)],'Color','r');
        hold on
        line([107 107],[0 max(a)],'Color','r');
    else
         hold on
        plot([1:1:200],tempL*max(a),'m');
        hold on
        plot([1:1:200],tempR*max(a),'g');
        hold on
        line([178 178],[0 max(a)],'Color','r');
    end
    xlim([1 200]);
    title(['day ',num2str(n)]);
end

saveas(gcf,'allCells.fig');
       
        
    
    