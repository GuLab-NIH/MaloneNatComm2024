
figure;
iB=[3 7 8 10];
iG=setdiff([1:1:15],iB);
% save('iB.mat','iB');
% save('iG.mat','iG');

%% Grid: Matrix correlation %2 figures

load('E:\learningAnalysis\behaviorVSActivity\MatrixCorrelation\deceNTDiffCM\indivMiceFolders.mat');
nData=9;%only look at new FOVs, 
corrM=[];

for n=1:length(indivMiceFolders);  
    load([indivMiceFolders{n} '\gridCellsNewCueThresh\matrixCorrMax1.mat']);
    corrM(end+1:end+nData)=matrixCorr(2:2+nData-1);
end

% save('corrM.mat','corrM');

NPerMouse=9;

idxB=[];
for n=1:length(iB);
    a=[NPerMouse*(iB(n)-1)+1:1:NPerMouse*iB(n)];
    idxB(end+1:end+length(a))=a;
end

idxG=[];
for n=1:length(iG);
    a=[NPerMouse*(iG(n)-1)+1:1:NPerMouse*iG(n)];
    idxG(end+1:end+length(a))=a;
end

load('E:\learningAnalysis\behaviorVSActivity\MatrixCorrelation\deceNTDiffCM\dece12.mat');

subplot(3,6,1),
a=dece12{4};
b=corrM;
plot(a(idxG),b(idxG),'g.')
hold on
plot(a(idxB),b(idxB),'m.')
[r,p]=corr(a',b');
title(['Grid matrix slow ',num2str([r p])]);

[p]=polyfit(a,b,1);
x=[0.9*min(a) 1.1*max(a)];
y=[p(1)*0.9*min(a)+p(2) p(1)*1.1*max(a)+p(2)];
hold on
line(x,y)

load('E:\learningAnalysis\behaviorVSActivity\MatrixCorrelation\lickingDiffCM\lick12.mat');

subplot(3,6,7),
a=lick12{1};
b=corrM;
plot(a(idxG),b(idxG),'g.')
hold on
plot(a(idxB),b(idxB),'m.')
[r,p]=corr(a',b');
title(['Grid matrix lick ',num2str([r p])]);
[p]=polyfit(a,b,1);
x=[0.9*min(a) 1.1*max(a)];
y=[p(1)*0.9*min(a)+p(2) p(1)*1.1*max(a)+p(2)];
hold on
line(x,y)
%% Grid:individual cell mean activity correlation %2 figures

nData=9;%only look at new FOVs, 
corrCA=[];

for n=1:length(indivMiceFolders);
    load([indivMiceFolders{n} '\gridCellsNewCueThresh\corrIndivi.mat']);
    corrCA(end+1:end+nData)=mean(corrIndivi(:,2:2+nData-1),1);
end


subplot(3,6,2),
a=dece12{4};
b=corrCA;
plot(a(idxG),b(idxG),'g.')
hold on
plot(a(idxB),b(idxB),'m.')
[r,p]=corr(a',b');
title(['Grid indivCell slow ',num2str([r p])]);
[p]=polyfit(a,b,1);
x=[0.9*min(a) 1.1*max(a)];
y=[p(1)*0.9*min(a)+p(2) p(1)*1.1*max(a)+p(2)];
hold on
line(x,y)

subplot(3,6,8),
a=lick12{1};
b=corrCA;
plot(a(idxG),b(idxG),'g.')
hold on
plot(a(idxB),b(idxB),'m.')
[r,p]=corr(a',b');
title(['Grid indivCell lick ',num2str([r p])]);
[p]=polyfit(a,b,1);
x=[0.9*min(a) 1.1*max(a)];
y=[p(1)*0.9*min(a)+p(2) p(1)*1.1*max(a)+p(2)];
hold on
line(x,y)

%% Grid: RBR activity correlation %2 figures
load('E:\learningAnalysis\behaviorVSActivity\RunByRun\Final_LickingAllOtherActvParams_Grid\corrGrid.mat');
b=corrGrid{3};

NPerMouse=10;

idxB=[];
for n=1:length(iB);
    a=[NPerMouse*(iB(n)-1)+1:1:NPerMouse*iB(n)];
    idxB(end+1:end+length(a))=a;
end

idxG=[];
for n=1:length(iG);
    a=[NPerMouse*(iG(n)-1)+1:1:NPerMouse*iG(n)];
    idxG(end+1:end+length(a))=a;
end


load('E:\learningAnalysis\behaviorVSActivity\RunByRun\Final_NT_decelerationInTimDiffcm\Dece.mat');
subplot(3,6,3);
a=Dece{4};
plot(a(idxG),b(idxG),'g.')
hold on
plot(a(idxB),b(idxB),'m.')
[r,p]=corr(a',b');
title(['Grid RBRCorr slow ',num2str([r p])]);
[p]=polyfit(a,b,1);
x=[0.9*min(a) 1.1*max(a)];
y=[p(1)*0.9*min(a)+p(2) p(1)*1.1*max(a)+p(2)];
hold on
line(x,y)

load('E:\learningAnalysis\behaviorVSActivity\RunByRun\Final_LickingAllDiffcm\Lick.mat');
subplot(3,6,9);
a=Lick{1};
plot(a(idxG),b(idxG),'g.')
hold on
plot(a(idxB),b(idxB),'m.')
[r,p]=corr(a',b');
title(['Grid RBRCorr lick ',num2str([r p])]);
[p]=polyfit(a,b,1);
x=[0.9*min(a) 1.1*max(a)];
y=[p(1)*0.9*min(a)+p(2) p(1)*1.1*max(a)+p(2)];
hold on
line(x,y)
%% Cue: Matrix correlation %2 figures

load('E:\learningAnalysis\behaviorVSActivity\MatrixCorrelation\deceNTDiffCM\indivMiceFolders.mat');
load('E:\learningAnalysis\behaviorVSActivity\RunByRun\corrAllCueRL\idxUse.mat');
indivMiceFoldersUse=indivMiceFolders(find(idxUse));
nData=9;%only look at new FOVs, 
corrM=[];

for n=1:length(indivMiceFoldersUse);  
    load([indivMiceFoldersUse{n} '\cueCellRLNewCueThresh\matrixCorrMax1.mat']);
    corrM(end+1:end+nData)=matrixCorr(2:2+nData-1);
end

% save('corrM.mat','corrM');

NPerMouse=9;

   %need new iB and iG because some mice have no cue cells 
iBUse=[3 6 7];
iGUse=setdiff([1:1:11],iBUse);

idxB=[];
for n=1:length(iBUse);
    a=[NPerMouse*(iBUse(n)-1)+1:1:NPerMouse*iBUse(n)];
    idxB(end+1:end+length(a))=a;
end

idxG=[];
for n=1:length(iGUse);
    a=[NPerMouse*(iGUse(n)-1)+1:1:NPerMouse*iGUse(n)];
    idxG(end+1:end+length(a))=a;
end

useData=[];
for n=1:length(idxUse);
    if idxUse(n)==1;
        useData(end+1:end+NPerMouse)=ones(1,NPerMouse);
    else
        useData(end+1:end+NPerMouse)=zeros(1,NPerMouse);
    end
end
      
load('E:\learningAnalysis\behaviorVSActivity\MatrixCorrelation\deceNTDiffCM\dece12.mat');

subplot(3,6,4),
a=dece12{4}(find(useData));
b=corrM;
plot(a(idxG),b(idxG),'g.')
hold on
plot(a(idxB),b(idxB),'m.')
[r,p]=corr(a',b');
title(['Cue matrix slow ',num2str([r p])]);
[p]=polyfit(a,b,1);
x=[0.9*min(a) 1.1*max(a)];
y=[p(1)*0.9*min(a)+p(2) p(1)*1.1*max(a)+p(2)];
hold on
line(x,y)

load('E:\learningAnalysis\behaviorVSActivity\MatrixCorrelation\lickingDiffCM\lick12.mat');

subplot(3,6,10),
a=lick12{1}(find(useData));
b=corrM;
plot(a(idxG),b(idxG),'g.')
hold on
plot(a(idxB),b(idxB),'m.')
[r,p]=corr(a',b');
title(['Cue matrix lick ',num2str([r p])]);
[p]=polyfit(a,b,1);
x=[0.9*min(a) 1.1*max(a)];
y=[p(1)*0.9*min(a)+p(2) p(1)*1.1*max(a)+p(2)];
hold on
line(x,y)
%% Cue: individual cell correlation %2 figures

load('E:\learningAnalysis\behaviorVSActivity\MatrixCorrelation\deceNTDiffCM\indivMiceFolders.mat');
load('E:\learningAnalysis\behaviorVSActivity\RunByRun\corrAllCueRL\idxUse.mat');
indivMiceFoldersUse=indivMiceFolders(find(idxUse));
nData=9;%only look at new FOVs, 
corrCA=[];

for n=1:length(indivMiceFoldersUse);
    load([indivMiceFoldersUse{n} '\corrIndivi.mat']);
    corrCA(end+1:end+nData)=mean(corrIndivi(:,2:2+nData-1),1);
end

% save('corrM.mat','corrM');

NPerMouse=9;

   %need new iB and iG because some mice have no cue cells 
iBUse=[3 6 7];
iGUse=setdiff([1:1:11],iBUse);

idxB=[];
for n=1:length(iBUse);
    a=[NPerMouse*(iBUse(n)-1)+1:1:NPerMouse*iBUse(n)];
    idxB(end+1:end+length(a))=a;
end

idxG=[];
for n=1:length(iGUse);
    a=[NPerMouse*(iGUse(n)-1)+1:1:NPerMouse*iGUse(n)];
    idxG(end+1:end+length(a))=a;
end

useData=[];
for n=1:length(idxUse);
    if idxUse(n)==1;
        useData(end+1:end+NPerMouse)=ones(1,NPerMouse);
    else
        useData(end+1:end+NPerMouse)=zeros(1,NPerMouse);
    end
end
      
load('E:\learningAnalysis\behaviorVSActivity\MatrixCorrelation\deceNTDiffCM\dece12.mat');

subplot(3,6,5),
a=dece12{4}(find(useData));
b=corrCA;
plot(a(idxG),b(idxG),'g.')
hold on
plot(a(idxB),b(idxB),'m.')
[r,p]=corr(a',b');
title(['Cue indivCell slow ',num2str([r p])]);
[p]=polyfit(a,b,1);
x=[0.9*min(a) 1.1*max(a)];
y=[p(1)*0.9*min(a)+p(2) p(1)*1.1*max(a)+p(2)];
hold on
line(x,y)
load('E:\learningAnalysis\behaviorVSActivity\MatrixCorrelation\lickingDiffCM\lick12.mat');

subplot(3,6,11),
a=lick12{1}(find(useData));
b=corrCA;
plot(a(idxG),b(idxG),'g.')
hold on
plot(a(idxB),b(idxB),'m.')
[r,p]=corr(a',b');
title(['Cue indivCell lick ',num2str([r p])]);
[p]=polyfit(a,b,1);
x=[0.9*min(a) 1.1*max(a)];
y=[p(1)*0.9*min(a)+p(2) p(1)*1.1*max(a)+p(2)];
hold on
line(x,y)
%% Cue: RBR activity correlation %2 figures
load('E:\learningAnalysis\behaviorVSActivity\RunByRun\corrAllCueRL\corrCueRL.mat');
b=corrCueRL{3};
load('E:\learningAnalysis\behaviorVSActivity\RunByRun\corrAllCueRL\idxUse.mat');

NPerMouse=10;
   %need new iB and iG because some mice have no cue cells 
iBUse=[3 6 7];
iGUse=setdiff([1:1:11],iBUse);

idxB=[];
for n=1:length(iBUse);
    a=[NPerMouse*(iBUse(n)-1)+1:1:NPerMouse*iBUse(n)];
    idxB(end+1:end+length(a))=a;
end

idxG=[];
for n=1:length(iGUse);
    a=[NPerMouse*(iGUse(n)-1)+1:1:NPerMouse*iGUse(n)];
    idxG(end+1:end+length(a))=a;
end

useData=[];
for n=1:length(idxUse);
    if idxUse(n)==1;
        useData(end+1:end+NPerMouse)=ones(1,NPerMouse);
    else
        useData(end+1:end+NPerMouse)=zeros(1,NPerMouse);
    end
end
        

load('E:\learningAnalysis\behaviorVSActivity\RunByRun\Final_NT_decelerationInTimDiffcm\Dece.mat');
subplot(3,6,6);
a=Dece{4}(find(useData));

plot(a(idxG),b(idxG),'g.')
hold on
plot(a(idxB),b(idxB),'m.')
[r,p]=corr(a',b');
title(['RBRCorr slow ',num2str([r p])]);
[p]=polyfit(a,b,1);
x=[0.9*min(a) 1.1*max(a)];
y=[p(1)*0.9*min(a)+p(2) p(1)*1.1*max(a)+p(2)];
hold on
line(x,y)
load('E:\learningAnalysis\behaviorVSActivity\RunByRun\Final_LickingAllDiffcm\Lick.mat');
subplot(3,6,12);
a=Lick{1}(find(useData));
plot(a(idxG),b(idxG),'g.')
hold on
plot(a(idxB),b(idxB),'m.')
[r,p]=corr(a',b');
title(['Cue RBRCorr lick ',num2str([r p])]);
[p]=polyfit(a,b,1);
x=[0.9*min(a) 1.1*max(a)];
y=[p(1)*0.9*min(a)+p(2) p(1)*1.1*max(a)+p(2)];
hold on
line(x,y)
%% cue with ranking and licking slowing: no individual days

iBUse=[3 6 7];
iGUse=setdiff([1:1:11],iBUse);

load('E:\learningAnalysis\behaviorVSActivity\RunByRun\corrAllCueRL\corrCueRLIndv.mat');
load('E:\learningAnalysis\behaviorVSActivity\RunByRun\Final_LickingAllOtherActvParams_Cue\idxUse.mat');
load('E:\learningAnalysis\behaviorVSActivity\RunByRun\Final_NT_decelerationInTimDiffcm\deceIndv15.mat');
L=corrCueRLIndv{3};
L=mean(L,2);
a=deceIndv15{4}(find(idxUse),:);
a=mean(a,2);
subplot(3,6,13)
plot(a(iGUse),L(iGUse),'g.','MarkerSize',10);
hold on
plot(a(iBUse),L(iBUse),'m.','MarkerSize',10);
[r,p]=corr(L,a,'Tail','right');
% [r,p]=corr(L,rrr,'Tail','right','Type','Spearman')
title(['corr cue slow ',num2str([r p])]);
[p]=polyfit(a,L,1);
x=[0.9*min(a) 1.1*max(a)];
y=[p(1)*0.9*min(a)+p(2) p(1)*1.1*max(a)+p(2)];
hold on
line(x,y)

load('E:\learningAnalysis\behaviorVSActivity\RunByRun\corrAllCueRL\corrCueRLIndv.mat');
load('E:\learningAnalysis\behaviorVSActivity\RunByRun\Final_LickingAllOtherActvParams_Cue\idxUse.mat');
load('E:\learningAnalysis\behaviorVSActivity\RunByRun\Final_LickingAllOtherActvParams\lickIndv15.mat');
L=corrCueRLIndv{3};
L=mean(L,2);
a=lickIndv15{1}(find(idxUse),:);
a=mean(a,2);
subplot(3,6,14)
plot(a(iGUse),L(iGUse),'g.','MarkerSize',10);
hold on
plot(a(iBUse),L(iBUse),'m.','MarkerSize',10);
[r,p]=corr(L,a,'Tail','right');
% [r,p]=corr(L,rrr,'Tail','right','Type','Spearman')
title(['corr cue lick ',num2str([r p])]);
[p]=polyfit(a,L,1);
x=[0.9*min(a) 1.1*max(a)];
y=[p(1)*0.9*min(a)+p(2) p(1)*1.1*max(a)+p(2)];
hold on
line(x,y)

load('E:\learningAnalysis\behaviorVSActivity\RunByRun\corrAllCueRL\corrCueRLIndv.mat');
load('E:\learningAnalysis\behaviorVSActivity\RunByRun\Final_LickingAllOtherActvParams_Cue\idxUse.mat');
load('E:\learningAnalysis\behaviorVSActivity\RunByRun\Param20_newRanking\rr.mat');
L=corrCueRLIndv{3};
L=mean(L,2);
rrr=rr(find(idxUse));
subplot(3,6,15)
plot(rrr(iGUse),L(iGUse),'g.','MarkerSize',10);
hold on
plot(rrr(iBUse),L(iBUse),'m.','MarkerSize',10);

[r,p]=corr(L,rrr,'Tail','right');
% [r,p]=corr(L,rrr,'Tail','right','Type','Spearman')
title(['corr cue ranking ',num2str([r p])]);
[p]=polyfit(rrr,L,1);
x=[0.9*min(rrr) 1.1*max(rrr)];
y=[p(1)*0.9*min(rrr)+p(2) p(1)*1.1*max(rrr)+p(2)];
hold on
line(x,y)

%%
%
tightfig;
saveas(gcf,'allActivityGridCuePSPL.fig')
