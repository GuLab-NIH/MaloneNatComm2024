cd ..\
load('foldersAllIndvFOV.mat');
cd('dfofTrack');

%% 

p=pwd;

folders=foldersAllIndvFOV;
folders=folders([1 2 3 4 5 6 7 8 9]);

allCellsDfof={};
allCellsDfof_sig={};

for n=1:11;
    allCellsDfof{n}=[];
    allCellsDfof_sig{n}=[];
end

% M=[];

for n=1:length(folders);
    disp(n)
    for m=1:length(folders{n});
        disp(m)
        cd(folders{n}{m});
        clear dfofaveragesmooth
        clear dfofaveragesmooth_sig
        d=dir('dfofaveragesmooth_*.mat');
        for i=1:length(d);
            load(d(i).name);
        end
        dfofaveragesmooth=dfofaveragesmooth';
%         M(n,m)=mean(dfofaveragesmooth,'All');
        dfofaveragesmooth_sig=dfofaveragesmooth_sig';
        allCellsDfof{m}(end+1:end+size(dfofaveragesmooth,1),:)=dfofaveragesmooth;
        allCellsDfof_sig{m}(end+1:end+size(dfofaveragesmooth_sig,1),:)=dfofaveragesmooth_sig;
    end
end
cd(p)


save('allCellsDfof.mat','allCellsDfof')
save('allCellsDfof_sig.mat','allCellsDfof_sig')

%%

meanAllDfof=[];
semAllDfof=[];
meanAllDfof_sig=[];
semAllDfof_sig=[];
for n=1:length(allCellsDfof);
    meanAllDfof(n,:)=nanmean(allCellsDfof{n},1);
     semAllDfof(n,:)=nansem(allCellsDfof{n},1);
       meanAllDfof_sig(n,:)=nanmean(allCellsDfof_sig{n},1);
     semAllDfof_sig(n,:)=nansem(allCellsDfof_sig{n},1);
end

load('E:\ID20211030\20211123Loc1Old\loc1Old\pcaica\cueAnalysisNew_sig\tempL.mat');
load('E:\ID20211030\20211123Loc1Old\loc1Old\pcaica\cueAnalysisNew_sig\tempR.mat');
tempLO=tempL;
tempRO=tempR;

load('E:\ID20211030\20211123New\loc1New\pcaica\cueAnalysisNew_sig\tempL.mat');
load('E:\ID20211030\20211123New\loc1New\pcaica\cueAnalysisNew_sig\tempR.mat');

figure
for n=1:11;
    subplot(4,3,n);
    a=meanAllDfof(n,:);
    b=semAllDfof(n,:);
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

saveas(gcf,'dfofAllCells.fig')

figure
for n=1:11;
    subplot(4,3,n);
    a=meanAllDfof_sig(n,:);
    b=semAllDfof_sig(n,:);
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

saveas(gcf,'dfof_sigAllCells.fig')

save('meanAllDfof.mat','meanAllDfof');
save('semAllDfof.mat','semAllDfof');
save('meanAllDfof_sig.mat','meanAllDfof_sig');
save('semAllDfof_sig.mat','semAllDfof_sig');