p=pwd;
cd ..\
load('foldersAll.mat');
cd(p);

errorsCommon={};
for n=1:11;
    errorsCommon{n}=[];
end
meanErrorsCommon=[];
perCommon=[];
for n=1:length(foldersAll);
    cd(foldersAll{n});
    cd('decodingCommonCells');
    load('error.mat');
    load('meanErrors.mat');
    load('perc.mat');

for m=1:length(error);
    errorsCommon{m}(end+1:end+length(error{m}))=error{m};
end
meanErrorsCommon(n,:)=meanErrors;
perCommon(n,:)=perc;
end

cd(p);

figure,
subplot(121)
errorbar([1:1:11],mean(meanErrorsCommon,1),nansem(meanErrorsCommon,1))
title('error');
subplot(122)
errorbar([1:1:11],mean(perCommon,1),nansem(perCommon,1))
title('percentage');
%%
p=pwd;
cd ..\
load('foldersAll.mat');
cd(p);

% errorsAll={};
% for n=1:11;
%     errorsAll{n}=[];
% end
meanErrorsAll=[];
perAll=[];
for n=1:length(foldersAll);
    cd(foldersAll{n});
    cd('decodingAllCells');
    load('error.mat');
    load('meanErrors.mat');
    load('perc.mat');

% for m=1:length(error);
%     errorsAll{m}(end+1:end+length(error{m}))=error{m};
% end
meanErrorsAll(n,:)=meanErrors;
perAll(n,:)=perc;
end

cd(p);

figure,
subplot(121)
errorbar([1:1:11],mean(meanErrorsAll,1),nansem(meanErrorsAll,1))
title('error');
subplot(122)
errorbar([1:1:11],mean(perAll,1),nansem(perAll,1))
title('percentage');


%% find number of cells
p=pwd;
cd ..\
load('foldersAllIndvFOV.mat');
cd(p);
NCells=[];
for n=1:length(foldersAllIndvFOV);
    for m=1:length(foldersAllIndvFOV{n})
   
    load([foldersAllIndvFOV{n}{m} '\roiIdxUse.mat']);
    NCells(end+1)=length(find(roiIdxUse));
end
end

min(NCells)

%%
p=pwd;
cd ..\
load('foldersAll.mat');
cd(p);

medianErrorsAll=[];
% for n=1:11;
%     errorsAll{n}=[];
% end
meanErrorsAll=[];
perAll=[];
% use=[1:1:22];
% use=use(1:end~=6);
for n=1:length(foldersAll)
    cd(foldersAll{n});
    cd('decoding50Cells');
    load('error.mat');
    load('meanErrors.mat');
    load('perc.mat');
M=[];
MPer=[];
for m=1:length(error);
    for mm=1:2;
    M(mm,m)=median(error{m}{mm});
    MPer(mm,m)=length(find(error{m}{mm}<=1))/length(error{m}{mm});
    end
end
medianErrorsAll(n,:)=mean(M,1);
meanErrorsAll(n,:)=meanErrors;
perAll(n,:)=mean(perc,1);
perAllAdjust(n,:)=mean(MPer,1);
end

cd(p);

figure,

subplot(141)
errorbar([1:1:11],mean(medianErrorsAll,1),nansem(medianErrorsAll,1))
title('medianError');

subplot(142)
errorbar([1:1:11],mean(meanErrorsAll,1),nansem(meanErrorsAll,1))
title('error');
subplot(143)
errorbar([1:1:11],mean(perAll,1),nansem(perAll,1))
title('percentage');
subplot(144)
errorbar([1:1:11],mean(perAllAdjust,1),nansem(perAllAdjust,1))
title('percentage diff <=1');

%% only consider error between a limited range of bins because the middles are more reliable. Based on this paper: https://www.sciencedirect.com/science/article/pii/S0896627321006929: the both ends of the track have uneven activity due to speed and also reward (our reward is around 890cm
p=pwd;
cd ..\
load('foldersAll.mat');
cd(p);

startBin=[5 10 15 20 25 30 35 40];
endBin=200-startBin;
% use=[1:1:22];
% use=use(1:end~=6);
perAllAdjust=[];
errorLessBinsAll={};
perDecoded={};
for fov=1:length(foldersAll)
    cd(foldersAll{fov});
    cd('decoding50Cells');
    d=dir('decd_*.mat');
    for m=1:length(startBin);
    errorS=[];
    perS=[];
    for n=1:length(d);
    load(d(n).name);
    iUse=find(decd.testYBins>=startBin(m)&decd.testYBins<=endBin(m));
    useE=decd.error(iUse);
    errorS(1,n)=mean(useE);
    perS(1,n)=length(find(useE<=2))/length(useE);
    end
    errorS=(errorS(1:2:21)+errorS(2:2:22))/2;
    perS=(perS(1:2:21)+perS(2:2:22))/2;
errorLessBinsAll{m}(fov,:)=errorS;
perDecoded{m}(fov,:)=perS;
    end
end
cd(p)
figure
for n=1:length(startBin);
    subplot(2,length(startBin)/2,n);
errorbar([1:1:11],mean(errorLessBinsAll{n}*5,1),nansem(errorLessBinsAll{n}*5,1))
title(['error',num2str(startBin(n))]);
end

saveas(gcf,'errors50cells.fig');
figure
for n=1:length(startBin);
    subplot(2,length(startBin)/2,n);
errorbar([1:1:11],mean(perDecoded{n},1),nansem(perDecoded{n},1))
title(['per',num2str(startBin(n))]);
end
saveas(gcf,'percentEncoded50cells.fig');
save('rand50Cells.mat','errorLessBinsAll','perDecoded');

%% try the same method on all cell data

p=pwd;
cd ..\
load('foldersAll.mat');
cd(p);

startBin=[5 10 15 20 25 30 35 40];
endBin=200-startBin;
% use=[1:1:22];
% use=use(1:end~=6);
perAllAdjust=[];
errorLessBinsAll={};
perDecoded={};
for fov=1:length(foldersAll)
    cd(foldersAll{fov});
    cd('decodingAllCells');
    d=dir('decd_*.mat');
    for m=1:length(startBin);
    errorS=[];
    perS=[];
    for n=1:length(d);
    load(d(n).name);
    iUse=find(decd.testYBins>=startBin(m)&decd.testYBins<=endBin(m));
    useE=decd.error(iUse);
    errorS(1,n)=mean(useE);
    perS(1,n)=length(find(useE<=2))/length(useE);
    end
%     errorS=(errorS(1:2:21)+errorS(2:2:22))/2;
%     perS=(perS(1:2:21)+perS(2:2:22))/2;
errorLessBinsAll{m}(fov,:)=errorS;
perDecoded{m}(fov,:)=perS;
    end
end
cd(p)
figure
for n=1:length(startBin);
    subplot(2,length(startBin)/2,n);
errorbar([1:1:11],mean(errorLessBinsAll{n}(2:end,:)*5,1),nansem(errorLessBinsAll{n}(2:end,:)*5,1))
title(['error',num2str(startBin(n))]);
end
saveas(gcf,'errorsAllcells.fig');
figure
for n=1:length(startBin);
    subplot(2,length(startBin)/2,n);
errorbar([1:1:11],mean(perDecoded{n}(2:end,:),1),nansem(perDecoded{n}(2:end,:),1))
title(['per',num2str(startBin(n))]);
end

saveas(gcf,'percentEncodedAllcells.fig');

save('allCells.mat','errorLessBinsAll','perDecoded');

%% only consider error between 150 and 850cm
p=pwd;
cd ..\
load('foldersAll.mat');
cd(p);

startBin=[5 10 15 20 25 30 35 40];
endBin=200-startBin;
% use=[1:1:22];
% use=use(1:end~=6);
perAllAdjust=[];
errorLessBinsAll={};
perDecoded={};
for fov=1:length(foldersAll)
    cd(foldersAll{fov});
    cd('decoding50CellsDown3');
    d=dir('decd_*.mat');
    for m=1:length(startBin);
    errorS=[];
    perS=[];
    for n=1:length(d);
    load(d(n).name);
    iUse=find(decd.testYBins>=startBin(m)&decd.testYBins<=endBin(m));
    useE=decd.error(iUse);
    errorS(1,n)=mean(useE);
    perS(1,n)=length(find(useE<=2))/length(useE);
    end
    errorS=(errorS(1:2:21)+errorS(2:2:22))/2;
    perS=(perS(1:2:21)+perS(2:2:22))/2;
errorLessBinsAll{m}(fov,:)=errorS;
perDecoded{m}(fov,:)=perS;
    end
end
cd(p)
figure
for n=1:length(startBin);
    subplot(2,length(startBin)/2,n);
errorbar([1:1:11],mean(errorLessBinsAll{n}*5,1),nansem(errorLessBinsAll{n}*5,1))
title(['error',num2str(startBin(n))]);
end
saveas(gcf,'errors50cellsDownSample3.fig');
figure
for n=1:length(startBin);
    subplot(2,length(startBin)/2,n);
errorbar([1:1:11],mean(perDecoded{n},1),nansem(perDecoded{n},1))
title(['per',num2str(startBin(n))]);
end
saveas(gcf,'percentEncoded50cellsDownSample3.fig');

save('rand50CellsDownSample3.mat','errorLessBinsAll','perDecoded');

%% try the same method on all cell data

p=pwd;
cd ..\
load('foldersAll.mat');
cd(p);

startBin=[5 10 15 20 25 30 35 40];
endBin=200-startBin;
% use=[1:1:22];
% use=use(1:end~=6);
perAllAdjust=[];
errorLessBinsAll={};
perDecoded={};
for fov=1:length(foldersAll)
    cd(foldersAll{fov});
    cd('decodingAllCellsDown3');
    d=dir('decd_*.mat');
    for m=1:length(startBin);
    errorS=[];
    perS=[];
    for n=1:length(d);
    load(d(n).name);
    iUse=find(decd.testYBins>=startBin(m)&decd.testYBins<=endBin(m));
    useE=decd.error(iUse);
    errorS(1,n)=mean(useE);
    perS(1,n)=length(find(useE<=2))/length(useE);
    end
%     errorS=(errorS(1:2:21)+errorS(2:2:22))/2;
%     perS=(perS(1:2:21)+perS(2:2:22))/2;
errorLessBinsAll{m}(fov,:)=errorS;
perDecoded{m}(fov,:)=perS;
    end
end
cd(p)
figure
for n=1:length(startBin);
    subplot(2,length(startBin)/2,n);
errorbar([1:1:11],mean(errorLessBinsAll{n}(2:end,:)*5,1),nansem(errorLessBinsAll{n}(2:end,:)*5,1))
title(['error',num2str(startBin(n))]);
end
saveas(gcf,'errorsAllcellsDownSample3.fig');

figure
for n=1:length(startBin);
    subplot(2,length(startBin)/2,n);
errorbar([1:1:11],mean(perDecoded{n}(2:end,:),1),nansem(perDecoded{n}(2:end,:),1))
title(['per',num2str(startBin(n))]);
end
saveas(gcf,'percentEncodedAllcellsDownSample3.fig');
save('allCellsDownSample3.mat','errorLessBinsAll','perDecoded');
%% only consider error between a limited range of bins because the middles are more reliable. Based on this paper: https://www.sciencedirect.com/science/article/pii/S0896627321006929: the both ends of the track have uneven activity due to speed and also reward (our reward is around 890cm
p=pwd;
cd ..\
load('foldersAll.mat');
cd(p);

startBin=[5 10 15 20 25 30 35 40];
endBin=200-startBin;
% use=[1:1:22];
% use=use(1:end~=6);
perAllAdjust=[];
errorLessBinsAll={};
errorLessBinsAllNorm={};
perDecoded={};
for fov=1:length(foldersAll)
    cd(foldersAll{fov});
    cd('decodingCommonCells');
    d=dir('decd_*.mat');
    for m=1:length(startBin);
    errorS=[];
    perS=[];
    for n=1:length(d);
    load(d(n).name);
    iUse=find(decd.testYBins>=startBin(m)&decd.testYBins<=endBin(m));
    useE=decd.error(iUse);
    errorS(1,n)=mean(useE);
    perS(1,n)=length(find(useE<=2))/length(useE);
    end
%     errorS=(errorS(1:2:21)+errorS(2:2:22))/2;
%     perS=(perS(1:2:21)+perS(2:2:22))/2;
errorLessBinsAll{m}(fov,:)=errorS;
errorLessBinsAllNorm{m}(fov,:)=errorS/errorS(2);
perDecoded{m}(fov,:)=perS;
    end
end
cd(p)
figure
for n=1:length(startBin);
    subplot(2,length(startBin)/2,n);
errorbar([1:1:11],mean(errorLessBinsAll{n}*5,1),nansem(errorLessBinsAll{n}*5,1))
title(['error',num2str(startBin(n))]);
end

saveas(gcf,'errorsCommonCells.fig');
figure
for n=1:length(startBin);
    subplot(2,length(startBin)/2,n);
errorbar([1:1:11],mean(perDecoded{n},1),nansem(perDecoded{n},1))
title(['per',num2str(startBin(n))]);
end
saveas(gcf,'percentEncodedCommonCells.fig');
save('commonCells.mat','errorLessBinsAll','perDecoded');

