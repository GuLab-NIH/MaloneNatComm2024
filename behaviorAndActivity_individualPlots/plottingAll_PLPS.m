figure;
iB=[3 7 8 10];
iG=setdiff([1:1:15],iB);

%% Matrix correlation %2 figures
load('E:\learningAnalysis\behaviorVSActivity\MatrixCorrelation\deceNTDiffCM\corrM.mat');
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
load('E:\learningAnalysis\behaviorVSActivity\MatrixCorrelation\lickingDiffCM\lick12.mat');

subplot(2,3,1),
a1=(dece12{4}-50)/50;
a2=lick12{1};
a=(a1+a2)/2;
b=corrM;
plot(a(idxG),b(idxG),'g.')
hold on
plot(a(idxB),b(idxB),'m.')
[r,p]=corr(a',b');
title(['matrix slow ',num2str([r p])]);
[p]=polyfit(a,b,1);

x=[0.9*min(a) 1.1*max(a)];
y=[p(1)*0.9*min(a)+p(2) p(1)*1.1*max(a)+p(2)];
hold on
line(x,y)

%% individual cell mean activity correlation %2 figures
load('E:\learningAnalysis\behaviorVSActivity\indivActivityCorrelation\deceDiffCM\corrCA.mat');

subplot(2,3,2),
a1=(dece12{4}-50)/50;
a2=lick12{1};
a=(a1+a2)/2;
b=corrCA;
plot(a(idxG),b(idxG),'g.')
hold on
plot(a(idxB),b(idxB),'m.')
[r,p]=corr(a',b');
title(['indivCell slow ',num2str([r p])]);
[p]=polyfit(a,b,1);

x=[0.9*min(a) 1.1*max(a)];
y=[p(1)*0.9*min(a)+p(2) p(1)*1.1*max(a)+p(2)];
hold on
line(x,y)

%% field numbers %2 figures
load('E:\learningAnalysis\behaviorVSActivity\numberOfField\NFieldsCommonCells.mat');
b1=NFieldsCommonCells;
b=[];
for n=1:size(b1,1);
    b(end+1:end+size(b1,2))=b1(n,:);
end
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

load('E:\learningAnalysis\behaviorVSActivity\RunByRun\Final_NT_decelerationInTimDiffcm\deceIndv15.mat');
load('E:\learningAnalysis\behaviorVSActivity\RunByRun\Final_LickingAllOtherActvParams\lickIndv15.mat');

subplot(2,3,3),
a1=deceIndv15{4};
a=[];
for n=1:size(a1,1);
    a(end+1:end+size(a1,2))=a1(n,:);
end
a1=(a-50)/50;

a2=lickIndv15{1};
a=[];
for n=1:size(a2,1);
    a(end+1:end+size(a2,2))=a2(n,:);
end
a2=a;

a=(a1+a2)/2;

plot(a(idxG),b(idxG),'g.')
hold on
plot(a(idxB),b(idxB),'m.')
[r,p]=corr(a',b');
title(['nfield slow ',num2str([r p])]);
[p]=polyfit(a,b,1);

x=[0.9*min(a) 1.1*max(a)];
y=[p(1)*0.9*min(a)+p(2) p(1)*1.1*max(a)+p(2)];
hold on
line(x,y)

%% RBR activity correlation %2 figures
load('E:\learningAnalysis\behaviorVSActivity\RunByRun\Final_LickingAllDiffcm\corrAll.mat');
b=corrAll{3};

load('E:\learningAnalysis\behaviorVSActivity\RunByRun\Final_NT_decelerationInTimDiffcm\Dece.mat');
load('E:\learningAnalysis\behaviorVSActivity\RunByRun\Final_LickingAllDiffcm\Lick.mat');

subplot(2,3,4);
a1=Dece{4};
a2=Lick{1};
a=((a1-50)/50+a2)/2;
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

%% Decoding_error %2 figures
load('E:\learningAnalysis\behaviorVSActivity\decoding\error50.mat');
b1=error50;
b=[];
for n=1:size(b1,1);
    b(end+1:end+size(b1,2))=b1(n,:);
end
load('E:\learningAnalysis\behaviorVSActivity\RunByRun\Final_NT_decelerationInTimDiffcm\deceIndv15.mat');
load('E:\learningAnalysis\behaviorVSActivity\RunByRun\Final_LickingAllOtherActvParams\lickIndv15.mat');

subplot(2,3,5),
a1=deceIndv15{4};
a=[];
for n=1:size(a1,1);
    a(end+1:end+size(a1,2))=a1(n,:);
end
a1=(a-50)/50;

a2=lickIndv15{1};
a=[];
for n=1:size(a2,1);
    a(end+1:end+size(a2,2))=a2(n,:);
end
a2=a;

a=(a1+a2)/2;

plot(a(idxG),b(idxG),'g.')
hold on
plot(a(idxB),b(idxB),'m.')
[r,p]=corr(a',b');
title(['deco error slow ',num2str([r p])]);
[p]=polyfit(a,b,1);
x=[0.9*min(a) 1.1*max(a)];
y=[p(1)*0.9*min(a)+p(2) p(1)*1.1*max(a)+p(2)];
hold on
line(x,y)

%% Decoding_percent track decoded %2 figures
load('E:\learningAnalysis\behaviorVSActivity\decoding\perc50.mat');
b1=perc50;
b=[];
for n=1:size(b1,1);
    b(end+1:end+size(b1,2))=b1(n,:);
end

load('E:\learningAnalysis\behaviorVSActivity\RunByRun\Final_NT_decelerationInTimDiffcm\deceIndv15.mat');
load('E:\learningAnalysis\behaviorVSActivity\RunByRun\Final_LickingAllOtherActvParams\lickIndv15.mat');

subplot(2,3,6),
a1=deceIndv15{4};
a=[];
for n=1:size(a1,1);
    a(end+1:end+size(a1,2))=a1(n,:);
end
a1=(a-50)/50;

a2=lickIndv15{1};
a=[];
for n=1:size(a2,1);
    a(end+1:end+size(a2,2))=a2(n,:);
end
a2=a;

a=(a1+a2)/2;
plot(a(idxG),b(idxG),'g.')
hold on
plot(a(idxB),b(idxB),'m.')
[r,p]=corr(a',b');
title(['deco perc slow ',num2str([r p])]);
[p]=polyfit(a,b,1);
x=[0.9*min(a) 1.1*max(a)];
y=[p(1)*0.9*min(a)+p(2) p(1)*1.1*max(a)+p(2)];
hold on
line(x,y)
%%
tightfig;
saveas(gcf,'allActivityPSPLCombined.fig')