% Age and gender of all 15 mice

%7 females and 8 males. For bad learners, there are 2 males and 2 females
% 
%below: the experiment start means the date of the last day in old env

% 1118 f (DOB 8/8/20, experiment start 12/24/20: 4.5 months)  138 days
% 
% 0206 m (DOB 8/9/20, experiment 03/06/21: 7 months) 209 days
% 
% 0207 f (DOB 8/8/20, experiment 03/23/21: 7.5 months) 227 days
% 
% 0208 m (DOB 8/8/20, experiment 03/05/21: 7 months) 209 days
% 
% 0209 m (DOB 8/8/20, experiment 03/12/21: 7 months)  216 days
% 
% 0413 m (DOB 12/10/20, experiment 06/03/21: 6 months) 175 days
% 
% 0519_1 m (DOB 1/21/21, experiment 06/17/21: 5 months) 147 days
% 
% 0519_2 m (DOB 1/21/21, experiment 06/12/21: 5 months)  142 days
% 
% 0802_S5E2 m (DOB 2/8/21, experiment 10/3/21: 8 months)  237 days
% 
% 0811A f (DOB 5/25/21, experiment 11/5/21: 5.5 months)  164 days
% 
% 0811B f (DOB 5/27/21, experiment 11/29/21: 6 months)  186 days
% 
% 1029 m (DOB 7/18/21, experiment 11/30/21: 4.5 months)  135 days
% 
% 1030 f (DOB 7/18/21, experiment 11/23/21: 4 months) 128 days
% 
% 1105_2 f (DOB 7/18/21, experiment 11/21/21: 4 months)  126 days
% 
% 1105_3 f (DOB 7/14/21, experiment 12/13/21: 4 months) 152 days
% 
%
% Did not include this because unimageable after 4-5 days
% 
% 0810 m (DOB 5/27/21, experiment 11/30/21: 6 months) 

%%
ageGood=[138 209 209 216 175 237 186 135 128 126 152];
ageBad=[227 147 142 164];

figure
bar([1 2],[mean(ageGood) mean(ageBad)]);
hold on
errorbar([1 2],[mean(ageGood) mean(ageBad)],[nansem(ageGood) nansem(ageBad)],'k.');
xlim([0 3]);
ylabel('days');
[r,p]=ttest2(ageGood,ageBad)
title(['age p=',num2str(p)]);
saveas(gcf,'age.fig');
save('ageGood.mat','ageGood');
save('ageBad.mat','ageBad');

%% plot performace and age
age=[138 209 227 209 216 175 147 142 237 164 186 135 128 126 152];

   %see whether this correlates with the new env behavior
   load('E:\learningAnalysis\behaviorVSActivity\RunByRun\Final_LickingAllOtherActvParams\lickIndv15.mat');
   LN=lickIndv15{1};
LN=mean(LN,2);

   load('E:\learningAnalysis\behaviorVSActivity\RunByRun\Final_NT_decelerationInTimDiffcm\deceIndv15.mat');
  SN=deceIndv15{4};%150CM before reward
SN=mean(SN,2);

figure
iB=[3 7 8 10];
iG=setdiff([1:1:15],iB);
subplot(121);
plot(age(iB),LN(iB),'m.','MarkerSize',15);
hold on
plot(age(iG),LN(iG),'g.','MarkerSize',15);

[r,p]=corr(age',LN,'tail','right');
title(['age lick p=',num2str(p)]);

subplot(122);
plot(age(iB),SN(iB),'m.','MarkerSize',15);
hold on
plot(age(iG),SN(iG),'g.','MarkerSize',15);

[r,p]=corr(age',SN,'tail','right');
title(['age slow p=',num2str(p)]);

saveas(gcf,'ageLickSlow.fig');

%% gender
iF=[1 3 10 11 13 14 15];
iM=setdiff([1:1:15],iF);
   load('E:\learningAnalysis\behaviorVSActivity\RunByRun\Final_LickingAllOtherActvParams\lickIndv15.mat');
   LN=lickIndv15{1};
LN=mean(LN,2);

   load('E:\learningAnalysis\behaviorVSActivity\RunByRun\Final_NT_decelerationInTimDiffcm\deceIndv15.mat');
  SN=deceIndv15{4};%150CM before reward
SN=mean(SN,2);

LNF=LN(iF);
LNM=LN(iM);

SNF=SN(iF);
SNM=SN(iM);

figure
A=LNF;
B=LNM;
M=[];
S=[];
M(1)=mean(A);
M(2)=mean(B);
S(1)=nansem(A,1);
S(2)=nansem(B,1);

subplot(121)
bar([1 2],M);
hold on
errorbar([1 2],M,S,'.');
names={'female','male'};
set(gca,'xtick',[1 2],'xticklabel',names);
[r,p]=ttest2(A,B);
title(['lick p=',num2str(p)]);


A=SNF;
B=SNM;
M=[];
S=[];
M(1)=mean(A);
M(2)=mean(B);
S(1)=nansem(A,1);
S(2)=nansem(B,1);

subplot(122)
bar([1 2],M);
hold on
errorbar([1 2],M,S,'.');
names={'female','male'};
set(gca,'xtick',[1 2],'xticklabel',names);
[r,p]=ttest2(A,B);
title(['slow p=',num2str(p)]);
saveas(gcf,'genderLickSlow.fig');

