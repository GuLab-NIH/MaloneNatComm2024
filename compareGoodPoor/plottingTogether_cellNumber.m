%% cell numbers without normalization
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\cellNumbersEachFOV\cellNumbers.mat');
G=cellNumbers(2:end,:);%remove the first mouse because it was 40x
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\cellNumbersEachFOV\cellNumbers.mat');
B=cellNumbers;

A=G;
M=mean(A,1);
[r,p]=corr([1:1:size(A,2)-1]',M(2:end)')
S=std(A,1)/sqrt(size(A,1));

figure,
subplot(421)
hold on,errorbar(1:1:length(M),M,S,'g-')
hold on, plot(1:1:length(M),M,'g.','MarkerSize',10)

A=B;
M=mean(A,1);
[r,p]=corr([1:1:size(A,2)-1]',M(2:end)')
S=std(A,1)/sqrt(size(A,1));
hold on,errorbar(1:1:length(M),M,S,'m-')
hold on, plot(1:1:length(M),M,'m.','MarkerSize',10)
xlim([0.7 11.3]);
title('cell numbers no normalization')

sig=[];
for n=1:length(M);
    [sig(n),~]=ttest2(G(:,n),B(:,n));
end
[pAnova,pMC] = anovaRM2W(G,B)
% P=0.0049
%% plot just the first two days

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\cellNumbersEachFOV\cellNumbers.mat');
G=cellNumbers;% no need to remove remove the first mouse because it was 40x
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\cellNumbersEachFOV\cellNumbers.mat');
B=cellNumbers;
%add the last mouse 1105_3
G(23,1)=279;
G(23,2)=328;

subplot(422)
hold on
bar([1 2 4 5],[mean(G(:,1)) mean(G(:,2)) mean(B(:,1)) mean(B(:,2))])
hold on
errorbar([1 2 4 5],[mean(G(:,1)) mean(G(:,2)) mean(B(:,1)) mean(B(:,2))],[nansem(G(:,1),1) nansem(G(:,2),1) nansem(B(:,1),1) nansem(B(:,2),1)],'k.')
for n=1:size(G,1);
    hold on
    plot([1 2],G(n,[1 2]),'Color',[0.7 0.7 0.7])
    hold on
    plot([1 2],G(n,[1 2]),'.','MarkerSize',10,'Color',[0.7 0.7 0.7])
end
for n=1:size(B,1);
    hold on
    plot([4 5],B(n,[1 2]),'Color',[0.7 0.7 0.7])
    hold on
    plot([4 5],B(n,[1 2]),'.','MarkerSize',10,'Color',[0.7 0.7 0.7])
end

[r1,p1]=ttest(G(:,1),G(:,2))%p1=0.0021
[r2,p2]=ttest(B(:,1),B(:,2))%p2=0.0128

%% plot N cell as a function of behavior
load('E:\learningAnalysis\behaviorVSActivity\RunByRun\Param18_Ranking\rr.mat');

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\cellNumbersEachFOV\cellNumbers.mat');
G=cellNumbers;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\cellNumbersEachFOV\cellNumbers.mat');
B=cellNumbers;

idxG={};%indices of what FOV for what cell
idxG{1}=[1];%remove the first mouse
idxG{2}=[2 3 4];
idxG{3}=[5 6];
idxG{4}=[7 8 9];
idxG{5}=[10 11];
idxG{6}=[12 13];
idxG{7}=[14 15];
idxG{8}=[16];
idxG{9}=[17 18 19];
idxG{10}=[20 21 22];
% idxG{11}=[23 24 25];

idxB={};
idxB{1}=[1 2];
idxB{2}=[3 4];
idxB{3}=[5 6 7];
idxB{4}=[8 9];

GN=[];
for n=1:length(idxG);
    GN(n)=mean2(G(idxG{n},:));
end

BN=[];
for n=1:length(idxB);
    BN(n)=mean2(B(idxB{n},:));
end

N=[];
N([1 2])=GN(1:2);
N(3)=BN(3);
N(4:6)=GN(3:5);
N(7:8)=BN(2:3);
N(9)=GN(6);
N(10)=BN(4);
N(11:14)=GN(7:10);
N15=[279 328 307 319 264 309 285];%manually input the number of active cells for mouse 15 before the imaging got bad
N(15)=mean(N15(1:6));
% [r,p]=corr(N(2:end)',rr(2:end-1));

%use bahavioral ranking: see whether better mice had less number of active
%cells
load('E:\learningAnalysis\behaviorVSActivity\RunByRun\Param20_newRanking\rr.mat');
subplot(423)
iB=[3 7 8 10];
iG=setdiff([2:1:14],iB);
plot(rr(iB),N(iB),'m.');%NOT INCLUDING 15
hold on
plot(rr(iG),N(iG),'g.');
[r,p]=corr(rr(2:end-1),N(2:end-1)');
title(['mean N cell vs behavior ranking, p=',num2str([r p])]);
[p]=polyfit(rr,N,1);

x=[0.9*min(rr) 1.1*max(rr)];
y=[p(1)*0.9*min(rr)+p(2) p(1)*1.1*max(rr)+p(2)];
hold on
line(x,y)

load('E:\learningAnalysis\behaviorVSActivity\RunByRun\Final_NT_decelerationInTimDiffcm\deceIndv15.mat');
a1=deceIndv15{4};
a1=mean(a1,2);
subplot(425)
plot(a1(iB),N(iB),'m.');
hold on
plot(a1(iG),N(iG),'g.');
[r,p]=corr(mean(a1(2:14,:),2),N(2:end-1)');
title(['mean N cell vs slowing, p=',num2str([r p])]);

[p]=polyfit(a1,N,1);

x=[0.9*min(a1) 1.1*max(a1)];
y=[p(1)*0.9*min(a1)+p(2) p(1)*1.1*max(a1)+p(2)];
hold on
line(x,y)


load('E:\learningAnalysis\behaviorVSActivity\RunByRun\Final_LickingAllOtherActvParams\lickIndv15.mat');
a1=lickIndv15{1};
a1=mean(a1,2);
subplot(427)
plot(a1(iB),N(iB),'m.');
hold on
plot(a1(iG),N(iG),'g.');
[r,p]=corr(mean(a1(2:14,:),2),N(2:end-1)');
title(['mean N cell vs licking, p=',num2str([r p])]);
[p]=polyfit(a1,N,1);

x=[0.9*min(a1) 1.1*max(a1)];
y=[p(1)*0.9*min(a1)+p(2) p(1)*1.1*max(a1)+p(2)];
hold on
line(x,y)

%the fractional increase of number of cells with behavior
idxG={};%indices of what FOV for what cell
idxG{1}=[1];%remove the first mouse
idxG{2}=[2 3 4];
idxG{3}=[5 6];
idxG{4}=[7 8 9];
idxG{5}=[10 11];
idxG{6}=[12 13];
idxG{7}=[14 15];
idxG{8}=[16];
idxG{9}=[17 18 19];
idxG{10}=[20 21 22];
% idxG{11}=[23];%since this is only for the first two days, so can include the last mouse and first

GN=[];
for n=1:length(idxG);
    for m=1:length(idxG{n});
    k=(G(idxG{n}(m),2)-G(idxG{n}(m),1))/G(idxG{n}(m),1);
    GN(n)=mean(k);
    end
end

BN=[];
for n=1:length(idxB);
    for m=1:length(idxB{n});
    k=(B(idxB{n}(m),2)-B(idxB{n}(m),1))/B(idxB{n}(m),1);
    BN(n)=mean(k);
end
end

N=[];
N([1 2])=GN(1:2);
N(3)=BN(3);
N(4:6)=GN(3:5);
N(7:8)=BN(2:3);
N(9)=GN(6);
N(10)=BN(4);
N(11:14)=GN(7:10);
N(15)=0.1756;%manually put in, because G doesn't include mouse 15
subplot(424)
plot(rr,N,'r.');
[r,p]=corr(rr,N');
title(['percent cell increase vs behavior, p=',num2str(p)]);

load('E:\learningAnalysis\behaviorVSActivity\RunByRun\Final_NT_decelerationInTimDiffcm\deceIndv15.mat');
a1=deceIndv15{4};
subplot(426)
plot(mean(a1,2),N,'r.');
[r,p]=corr(mean(a1,2),N');
title(['percent cell increase vs slowing, p=',num2str(p)]);

load('E:\learningAnalysis\behaviorVSActivity\RunByRun\Final_LickingAllOtherActvParams\lickIndv15.mat');
a1=lickIndv15{1};
subplot(428)
plot(mean(a1,2),N,'r.');
[r,p]=corr(mean(a1,2),N');
title(['percent cell increase vs licking, p=',num2str(p)]);



saveas(gcf,'cellNumbersNoNormalizeCompareOneLastOldDay.fig');
% saveas(gcf,'cellNumbersNoNormalizeCompareOneLastOldDay.eps');



%%
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\cellNumbersEachFOV\cellNumbersNorm.mat');
G=cellNumbersNorm;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\cellNumbersEachFOV\cellNumbersNorm.mat');
B=cellNumbersNorm;

%in this data, since the G data the first FOV (old) is somehow very low, we
%used fov2 (new in env) to normalize all numbers
% 
% G=G(1:end~=9,:);
A=G;
M=mean(A,1);
[r,p]=corr([1:1:size(A,2)-1]',M(2:end)')
S=std(A,1)/sqrt(size(A,1));

figure,
subplot(221)
hold on,errorbar(1:1:length(M),M,S,'g-')
hold on, plot(1:1:length(M),M,'g.','MarkerSize',10)

A=B;
M=mean(A,1);
[r,p]=corr([1:1:size(A,2)-1]',M(2:end)')
S=std(A,1)/sqrt(size(A,1));
hold on,errorbar(1:1:length(M),M,S,'m-')
hold on, plot(1:1:length(M),M,'m.','MarkerSize',10)
xlim([0.7 11.3]);
title('cell numbers normalized by old day')

sig=[];
for n=2:length(M);
    [sig(n),~]=ttest2(G(:,n),B(:,n));
end



%% cell numbers norm old to mean std to 1
load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\cellNumbersEachFOV\cellNumbersNorm2.mat');
G=cellNumbersNorm2;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\cellNumbersEachFOV\cellNumbersNorm2.mat');
B=cellNumbersNorm2;

A=G;
M=mean(A,1);
[r,p]=corr([1:1:size(A,2)-1]',M(2:end)')
S=std(A,1)/sqrt(size(A,1));


subplot(222)
hold on,errorbar(1:1:length(M),M,S,'g-')
hold on, plot(1:1:length(M),M,'g.','MarkerSize',10)

A=B;
M=mean(A,1);
[r,p]=corr([1:1:size(A,2)-1]',M(2:end)')
S=std(A,1)/sqrt(size(A,1));
hold on,errorbar(1:1:length(M),M,S,'m-')
hold on, plot(1:1:length(M),M,'m.','MarkerSize',10)
xlim([0.7 11.3]);
title('cell numbers normalized by mean std old day')


%% plotting the new cell numbers (the old fov is the averaged one from multiple sessions)

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\cellNumbersEachFOV\newCellNumbersNorm.mat');
G=newCellNumbersNorm;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\cellNumbersEachFOV\newCellNumbersNorm.mat');
B=newCellNumbersNorm;

A=G;
M=mean(A,1);
[r,p]=corr([1:1:size(A,2)-1]',M(2:end)')
S=std(A,1)/sqrt(size(A,1));

subplot(223)
hold on,errorbar(1:1:length(M),M,S,'g-')
hold on, plot(1:1:length(M),M,'g.','MarkerSize',10)

A=B;
M=mean(A,1);
[r,p]=corr([1:1:size(A,2)-1]',M(2:end)')
S=std(A,1)/sqrt(size(A,1));
hold on,errorbar(1:1:length(M),M,S,'m-')
hold on, plot(1:1:length(M),M,'m.','MarkerSize',10)
title('cell numbers with old mean normalized by old mean')
xlim([0.7 11.3]);

sig=[];
for n=2:length(M);
    [sig(n),~]=ttest2(G(:,n),B(:,n));
end

%% plotting the new cell numbers (the old fov is the averaged one from multiple sessions)

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\cellNumbersEachFOV\newCellNumbersNorm2.mat');
G=newCellNumbersNorm2;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\cellNumbersEachFOV\newCellNumbersNorm2.mat');
B=newCellNumbersNorm2;

A=G;
M=mean(A,1);
[r,p]=corr([1:1:size(A,2)-1]',M(2:end)')
S=std(A,1)/sqrt(size(A,1));

subplot(224)
hold on,errorbar(1:1:length(M),M,S,'g-')
hold on, plot(1:1:length(M),M,'g.','MarkerSize',10)

A=B;
M=mean(A,1);
[r,p]=corr([1:1:size(A,2)-1]',M(2:end)')
S=std(A,1)/sqrt(size(A,1));
hold on,errorbar(1:1:length(M),M,S,'m-')
hold on, plot(1:1:length(M),M,'m.','MarkerSize',10)
title('cell numbers with old mean normalized by old mean std')
xlim([0.7 11.3]);

sig=[];
for n=2:length(M);
    [sig(n),~]=ttest2(G(:,n),B(:,n));
end

saveas(gcf,'cellNumbersCompare.fig');
saveas(gcf,'cellNumbersCompare.eps');

%% plotting the new cell numbers (the old fov is the averaged one from multiple sessions)
%only the first FOV


load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\cellNumbersEachFOV\cellNumbersNormFOV1.mat');
G=cellNumbersNormFOV1;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\cellNumbersEachFOV\cellNumbersNormFOV1.mat');
B=cellNumbersNormFOV1;

A=G;
M=mean(A,1);
[r,p]=corr([1:1:size(A,2)-1]',M(2:end)')
S=std(A,1)/sqrt(size(A,1));
figure
subplot(221)
hold on,errorbar(1:1:length(M),M,S,'g-')
hold on, plot(1:1:length(M),M,'g.','MarkerSize',10)

A=B;
M=mean(A,1);
[r,p]=corr([1:1:size(A,2)-1]',M(2:end)')
S=std(A,1)/sqrt(size(A,1));
hold on,errorbar(1:1:length(M),M,S,'m-')
hold on, plot(1:1:length(M),M,'m.','MarkerSize',10)
title('cell numbers with old mean normalized by old mean fov1')
xlim([0.7 11.3]);

sig=[];
for n=2:length(M);
    [sig(n),~]=ttest2(G(:,n),B(:,n));
end


% saveas(gcf,'cellNumbersCompareFOV1.fig');
% saveas(gcf,'cellNumbersCompareFOV1.eps');


%% plotting the new cell numbers (the old fov is the averaged one from multiple sessions)
%only the first FOV


load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\cellNumbersEachFOV\cellNumbersNorm2FOV1.mat');
G=cellNumbersNorm2FOV1;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\cellNumbersEachFOV\cellNumbersNorm2FOV1.mat');
B=cellNumbersNorm2FOV1;

A=G;
M=mean(A,1);
[r,p]=corr([1:1:size(A,2)-1]',M(2:end)')
S=std(A,1)/sqrt(size(A,1));

subplot(222)
hold on,errorbar(1:1:length(M),M,S,'g-')
hold on, plot(1:1:length(M),M,'g.','MarkerSize',10)

A=B;
M=mean(A,1);
[r,p]=corr([1:1:size(A,2)-1]',M(2:end)')
S=std(A,1)/sqrt(size(A,1));
hold on,errorbar(1:1:length(M),M,S,'m-')
hold on, plot(1:1:length(M),M,'m.','MarkerSize',10)
title('cell numbers with old mean normalized by old mean std fov1')
xlim([0.7 11.3]);

sig=[];
for n=2:length(M);
    [sig(n),~]=ttest2(G(:,n),B(:,n));
end


%% plotting the new cell numbers (the old fov is the averaged one from multiple sessions)
%only the first FOV


load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\cellNumbersEachFOV\newCellNumbersNormFOV1.mat');
G=newCellNumbersNormFOV1;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\cellNumbersEachFOV\newCellNumbersNormFOV1.mat');
B=newCellNumbersNormFOV1;

A=G;
M=mean(A,1);
[r,p]=corr([1:1:size(A,2)-1]',M(2:end)')
S=std(A,1)/sqrt(size(A,1));

subplot(223)
hold on,errorbar(1:1:length(M),M,S,'g-')
hold on, plot(1:1:length(M),M,'g.','MarkerSize',10)

A=B;
M=mean(A,1);
[r,p]=corr([1:1:size(A,2)-1]',M(2:end)')
S=std(A,1)/sqrt(size(A,1));
hold on,errorbar(1:1:length(M),M,S,'m-')
hold on, plot(1:1:length(M),M,'m.','MarkerSize',10)
title('cell numbers with old mean normalized by old mean fov1')
xlim([0.7 11.3]);

sig=[];
for n=2:length(M);
    [sig(n),~]=ttest2(G(:,n),B(:,n));
end


% saveas(gcf,'cellNumbersCompareFOV1.fig');
% saveas(gcf,'cellNumbersCompareFOV1.eps');


%% plotting the new cell numbers (the old fov is the averaged one from multiple sessions)
%only the first FOV


load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\cellNumbersEachFOV\newCellNumbersNorm2FOV1.mat');
G=newCellNumbersNorm2FOV1;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\cellNumbersEachFOV\newCellNumbersNorm2FOV1.mat');
B=newCellNumbersNorm2FOV1;

A=G;
M=mean(A,1);
[r,p]=corr([1:1:size(A,2)-1]',M(2:end)')
S=std(A,1)/sqrt(size(A,1));

subplot(224)
hold on,errorbar(1:1:length(M),M,S,'g-')
hold on, plot(1:1:length(M),M,'g.','MarkerSize',10)

A=B;
M=mean(A,1);
[r,p]=corr([1:1:size(A,2)-1]',M(2:end)')
S=std(A,1)/sqrt(size(A,1));
hold on,errorbar(1:1:length(M),M,S,'m-')
hold on, plot(1:1:length(M),M,'m.','MarkerSize',10)
title('cell numbers with old mean normalized by old mean std fov1')
xlim([0.7 11.3]);

sig=[];
for n=2:length(M);
    [sig(n),~]=ttest2(G(:,n),B(:,n));
end


saveas(gcf,'cellNumbersCompareFOV1.fig');
saveas(gcf,'cellNumbersCompareFOV1.eps');


%% plot a bar graph

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\cellNumbersEachFOV\newCellNumbersNorm.mat');
G=newCellNumbersNorm;
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\cellNumbersEachFOV\newCellNumbersNorm.mat');
B=newCellNumbersNorm;

figure
A=G;
M=mean(A,1);
S=std(A,1)/sqrt(size(A,1));
bar(M(1:2))

[r,p]=ttest(A(:,1),A(:,2))
%% Whether the fractional changes in each FOV in good and bad performers are different

load('E:\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\cellNumbersEachFOV\cellNumbers.mat');
G=cellNumbers;% no need to remove remove the first mouse because it was 40x
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\cellNumbersEachFOV\cellNumbers.mat');
B=cellNumbers;
%add the last mouse 1105_3
G(23,1)=279;
G(23,2)=328;

GI=[];
for n=1:size(G,1);
    GI(n)=(G(n,2)-G(n,1))/G(n,1);
end

BI=[];
for n=1:size(B,1);
    BI(n)=(B(n,2)-B(n,1))/B(n,1);
end

[r,p]=ttest2(GI,BI);

figure,
bar([1 2],[mean(GI) mean(BI)]);
hold on
errorbar([1 2],[mean(GI) mean(BI)],[nansem(GI,2) nansem(BI,2)],'.')

title(['fractional change p=',num2str(p)]);

saveas(gcf,'fractionalChanges.fig');