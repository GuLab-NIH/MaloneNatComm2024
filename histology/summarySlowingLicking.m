%mice are in this order:

%(different from here: Z:\SNM\labMembers\LC\Lian
%Cui\behavioralParamsByYG_NT\summary.mat')

%same as here: Z:\SNM\labMembers\LC\Lian Cui\IHC\IHC Matlab\data statistic graph
%rearranged the order in "summary" so that the order can be identical to
%the histology data

% %mouse order
% (1) 302
% (2) 303
% (3) 304
% (4) 305
% (5) 306
% (6) 307
% (7) 401
% (8) 402
% (9) 403
% (10) 405
% (11) LC01
% (12) LC02
% (13) LC03
% (14) LC04
% (15) LC05
% (16) LC06
% (17) LC07
% (18) LC08
% (19) LC09
% (20) LC10

last6=[];%this includes all track (batch 3 and 4 (both are 10 meter) and also 2 (4 meter))
%first column is slowing
%second column is licking
%copy and rearrange values from here:Z:\SNM\labMembers\LC\Lian
%Cui\behavioralParamsByYG_NT\summary.mat'

save('last6.mat','last6');

last7=[];%this doesn't include 4m track (batch 3 and 4 (both are 10 meter))
%first column is slowing
%second column is licking
%copy and rearrange values from here:Z:\SNM\labMembers\LC\Lian
%Cui\behavioralParamsByYG_NT\summary.mat'
%4m data are nans

last7(7:10,1:2)=nan(4,2);
save('last7.mat','last7');

%% compare behavior of IHC mice with imaging mice
%imaging mice
load('lickIndv15.mat');
L=lickIndv15{1}*100;
load('deceIndv15.mat');
D=deceIndv15{4};
badIdx=[3 7 8 10];
goodIdx=setdiff([1:1:15],badIdx);
LB=L(badIdx,:); %LICK OF BAD: IMAGING MICE
LG=L(goodIdx,:); %LICK OF GOOD: IMAGING MICE
DB=D(badIdx,:); %DECELERATION OF BAD: IMAGING MICE
DG=D(goodIdx,:); %DECELERATION OF GOOD: IMAGING MICE



%(1) %no 4m track, last 7 days, threshold: all 10 days mean of imaging mice
load('last7.mat');
%GOOD BAD LEARNERS IDENTIFIED BASED ON THE ABOVE THRESH (Z:\SNM\labMembers\LC\Lian
%Cui\behavioralParamsByYG_NT\summary.mat')
idxGoodF=[2 3 13 17 18 19];
idxGoodN=[1 4 5 14 15 20];
idxGoodH=[idxGoodF idxGoodN];

idxBadF=[6 12];
idxBadN=[11 16];
idxBadH=[idxBadF idxBadN];

LBH=last7(idxBadH,2)*100;%LICK OF BAD: HISTOLOGY MICE
LGH=last7(idxGoodH,2)*100;%LICK OF GOOD: HISTOLOGY MICE

DBH=last7(idxBadH,1);%DECELERATION OF BAD: HISTOLOGY MICE
DGH=last7(idxGoodH,1);%DECELERATION OF GOOD: HISTOLOGY MICE

%plotting imaging and histology mice together
day=1;
%IMAGING MICE: MEAN LICKING ACROSS 10 DAYS)
LBI=mean(LB(:,day:end),2);
LGI=mean(LG(:,day:end),2);
DBI=mean(DB(:,day:end),2);
DGI=mean(DG(:,day:end),2);

M=[];
S=[];

A=LBI;
M(1)=mean(A);
S(1)=std(A)/sqrt(length(A));

A=LGI;
M(2)=mean(A);
S(2)=std(A)/sqrt(length(A));

A=LBH;
M(3)=mean(A);
S(3)=std(A)/sqrt(length(A));

A=LGH;
M(4)=mean(A);
S(4)=std(A)/sqrt(length(A));


A=DBI;
M(5)=mean(A);
S(5)=std(A)/sqrt(length(A));

A=DGI;
M(6)=mean(A);
S(6)=std(A)/sqrt(length(A));

A=DBH;
M(7)=mean(A);
S(7)=std(A)/sqrt(length(A));

A=DGH;
M(8)=mean(A);
S(8)=std(A)/sqrt(length(A));

N=1;
[SIG(1,N),~]=ttest2(LBI,LBH);
[SIG(2,N),~]=ttest2(LGI,LGH);
[SIG(3,N),~]=ttest2(DBI,DBH);
[SIG(4,N),~]=ttest2(DGI,DGH);
[SIG(5,N),~]=ttest2(LBI,LGI);
[SIG(6,N),~]=ttest2(LBH,LGH);
[SIG(7,N),~]=ttest2(DBI,DGI);
[SIG(8,N),~]=ttest2(DBH,DGH);
figure
subplot(411)
bar(M);
hold on
errorbar([1:1:8],M,S,'r.')
names={'LBI';'LGI';'LBH';'LGH';'DBI';'DGI';'DBH';'DGH'};
set(gca,'xtick',[1:8],'xticklabel',names)
title('no 4m, 10 day thresh');

%%
%(2) %no 4m track, last 7 days, threshold: the last 7 mean of imaging mice
load('last7.mat');
%GOOD BAD LEARNERS IDENTIFIED BASED ON THE ABOVE THRESH (Z:\SNM\labMembers\LC\Lian
%Cui\behavioralParamsByYG_NT\summary.mat')
idxGoodF=[2 3 13 19];
idxGoodN=[1 4 5 14 15];
idxGoodH=[idxGoodF idxGoodN];

idxBadF=[6 12 17 18];
idxBadN=[11 16 20];
idxBadH=[idxBadF idxBadN];

LBH=last7(idxBadH,2)*100;%LICK OF BAD: HISTOLOGY MICE
LGH=last7(idxGoodH,2)*100;%LICK OF GOOD: HISTOLOGY MICE

DBH=last7(idxBadH,1);%DECELERATION OF BAD: HISTOLOGY MICE
DGH=last7(idxGoodH,1);%DECELERATION OF GOOD: HISTOLOGY MICE

%plotting imaging and histology mice together
day=4;
%IMAGING MICE: MEAN LICKING ACROSS 10 DAYS)
LBI=mean(LB(:,day:end),2);
LGI=mean(LG(:,day:end),2);
DBI=mean(DB(:,day:end),2);
DGI=mean(DG(:,day:end),2);

M=[];
S=[];

A=LBI;
M(1)=mean(A);
S(1)=std(A)/sqrt(length(A));

A=LGI;
M(2)=mean(A);
S(2)=std(A)/sqrt(length(A));

A=LBH;
M(3)=mean(A);
S(3)=std(A)/sqrt(length(A));

A=LGH;
M(4)=mean(A);
S(4)=std(A)/sqrt(length(A));


A=DBI;
M(5)=mean(A);
S(5)=std(A)/sqrt(length(A));

A=DGI;
M(6)=mean(A);
S(6)=std(A)/sqrt(length(A));

A=DBH;
M(7)=mean(A);
S(7)=std(A)/sqrt(length(A));

A=DGH;
M(8)=mean(A);
S(8)=std(A)/sqrt(length(A));

N=2;
[SIG(1,N),~]=ttest2(LBI,LBH);
[SIG(2,N),~]=ttest2(LGI,LGH);
[SIG(3,N),~]=ttest2(DBI,DBH);
[SIG(4,N),~]=ttest2(DGI,DGH);
[SIG(5,N),~]=ttest2(LBI,LGI);
[SIG(6,N),~]=ttest2(LBH,LGH);
[SIG(7,N),~]=ttest2(DBI,DGI);
[SIG(8,N),~]=ttest2(DBH,DGH);

subplot(412)
bar(M);
hold on
errorbar([1:1:8],M,S,'r.')
names={'LBI';'LGI';'LBH';'LGH';'DBI';'DGI';'DBH';'DGH'};
set(gca,'xtick',[1:8],'xticklabel',names)
title('no 4m, 7 day thresh');

%%
%(3) %all tracks (including 4m), last 6 days, threshold: all 10 days mean of imaging mice
load('last6.mat');
%GOOD BAD LEARNERS IDENTIFIED BASED ON THE ABOVE THRESH (Z:\SNM\labMembers\LC\Lian
%Cui\behavioralParamsByYG_NT\summary.mat')
idxGoodF=[9 2 3 6 13 17 19];
idxGoodN=[8 4 5 14 15 20];
idxGoodH=[idxGoodF idxGoodN];

idxBadF=[7 12 18];
idxBadN=[10 1 11 16];
idxBadH=[idxBadF idxBadN];

LBH=last6(idxBadH,2)*100;%LICK OF BAD: HISTOLOGY MICE
LGH=last6(idxGoodH,2)*100;%LICK OF GOOD: HISTOLOGY MICE

DBH=last6(idxBadH,1);%DECELERATION OF BAD: HISTOLOGY MICE
DGH=last6(idxGoodH,1);%DECELERATION OF GOOD: HISTOLOGY MICE

%plotting imaging and histology mice together
day=1;
%IMAGING MICE: MEAN LICKING ACROSS 10 DAYS)
LBI=mean(LB(:,day:end),2);
LGI=mean(LG(:,day:end),2);
DBI=mean(DB(:,day:end),2);
DGI=mean(DG(:,day:end),2);

M=[];
S=[];

A=LBI;
M(1)=mean(A);
S(1)=std(A)/sqrt(length(A));

A=LGI;
M(2)=mean(A);
S(2)=std(A)/sqrt(length(A));

A=LBH;
M(3)=mean(A);
S(3)=std(A)/sqrt(length(A));

A=LGH;
M(4)=mean(A);
S(4)=std(A)/sqrt(length(A));


A=DBI;
M(5)=mean(A);
S(5)=std(A)/sqrt(length(A));

A=DGI;
M(6)=mean(A);
S(6)=std(A)/sqrt(length(A));

A=DBH;
M(7)=mean(A);
S(7)=std(A)/sqrt(length(A));

A=DGH;
M(8)=mean(A);
S(8)=std(A)/sqrt(length(A));

N=3;
[SIG(1,N),~]=ttest2(LBI,LBH);
[SIG(2,N),~]=ttest2(LGI,LGH);
[SIG(3,N),~]=ttest2(DBI,DBH);
[SIG(4,N),~]=ttest2(DGI,DGH);
[SIG(5,N),~]=ttest2(LBI,LGI);
[SIG(6,N),~]=ttest2(LBH,LGH);
[SIG(7,N),~]=ttest2(DBI,DGI);
[SIG(8,N),~]=ttest2(DBH,DGH);

subplot(413)
bar(M);
hold on
errorbar([1:1:8],M,S,'r.')
names={'LBI';'LGI';'LBH';'LGH';'DBI';'DGI';'DBH';'DGH'};
set(gca,'xtick',[1:8],'xticklabel',names)
title('all, 10 day thresh');

%%
%(4) %all tracks (including 4m), last 6 days, threshold: the last 6 days mean of imaging mice
load('last6.mat');
%GOOD BAD LEARNERS IDENTIFIED BASED ON THE ABOVE THRESH (Z:\SNM\labMembers\LC\Lian
%Cui\behavioralParamsByYG_NT\summary.mat')
idxGoodF=[2 3 13 17 19];
idxGoodN=[8 4 5 14 15 20];
idxGoodH=[idxGoodF idxGoodN];

idxBadF=[7 9 6 12 18];
idxBadN=[10 1 11 16];
idxBadH=[idxBadF idxBadN];

LBH=last6(idxBadH,2)*100;%LICK OF BAD: HISTOLOGY MICE
LGH=last6(idxGoodH,2)*100;%LICK OF GOOD: HISTOLOGY MICE

DBH=last6(idxBadH,1);%DECELERATION OF BAD: HISTOLOGY MICE
DGH=last6(idxGoodH,1);%DECELERATION OF GOOD: HISTOLOGY MICE

%plotting imaging and histology mice together
day=5;
%IMAGING MICE: MEAN LICKING ACROSS 10 DAYS)
LBI=mean(LB(:,day:end),2);
LGI=mean(LG(:,day:end),2);
DBI=mean(DB(:,day:end),2);
DGI=mean(DG(:,day:end),2);

M=[];
S=[];

A=LBI;
M(1)=mean(A);
S(1)=std(A)/sqrt(length(A));

A=LGI;
M(2)=mean(A);
S(2)=std(A)/sqrt(length(A));

A=LBH;
M(3)=mean(A);
S(3)=std(A)/sqrt(length(A));

A=LGH;
M(4)=mean(A);
S(4)=std(A)/sqrt(length(A));


A=DBI;
M(5)=mean(A);
S(5)=std(A)/sqrt(length(A));

A=DGI;
M(6)=mean(A);
S(6)=std(A)/sqrt(length(A));

A=DBH;
M(7)=mean(A);
S(7)=std(A)/sqrt(length(A));

A=DGH;
M(8)=mean(A);
S(8)=std(A)/sqrt(length(A));

N=4;
[SIG(1,N),~]=ttest2(LBI,LBH);
[SIG(2,N),~]=ttest2(LGI,LGH);
[SIG(3,N),~]=ttest2(DBI,DBH);
[SIG(4,N),~]=ttest2(DGI,DGH);
[SIG(5,N),~]=ttest2(LBI,LGI);
[SIG(6,N),~]=ttest2(LBH,LGH);
[SIG(7,N),~]=ttest2(DBI,DGI);
[SIG(8,N),~]=ttest2(DBH,DGH);

subplot(414)
bar(M);
hold on
errorbar([1:1:8],M,S,'r.')
names={'LBI';'LGI';'LBH';'LGH';'DBI';'DGI';'DBH';'DGH'};
set(gca,'xtick',[1:8],'xticklabel',names)
title('all, 6 day thresh');


saveas(gcf,'compareBehavior.fig');

save('SIG.mat','SIG');