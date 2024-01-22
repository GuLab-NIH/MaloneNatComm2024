load("ChR2_6P7G-5CM2MS-PS.mat");
poor=[1 4 6 7 10 11 14 16];%28, 34, 37, 38, 35, 33, 45, 44
poorChR2=CHR2(poor);

C=CHR2(poor);
C1=cell2mat(C);% all mice together
M=nanmean(C1,1);
S=nansem(C1,1);
figure,
subplot(2,6,1) %per sample CS
bar([1:1:3],M(1:3),'FaceColor','w');
hold on
errorbar([1:1:3],M(1:3),S(1:3),'k.')

ylim([0 100])

for n=1:size(C1,1);
    hold on
    plot(C1(n,[1 2 3]),'Color','m','LineWidth',0.5)
end

[h1,p1] = ttest(C1(:,1),C1(:,2));
[h2,p2] = ttest(C1(:,1),C1(:,3));
[h3,p3] = ttest(C1(:,2),C1(:,3));

[isSignificant,adjusted_pvals,~]= bonferroni_holm([p1 p2 p3],0.05); % adjusted_pvals are the corrected p values
%adjusted_pvals
title(['PoorChR2CS',num2str(adjusted_pvals([1 2 3]))]);

% title(['PoorChR2 CS p',num2str([p1 p3])]);

subplot(2,6,3) %per sample RS
bar([1:1:3],M(4:6),'FaceColor','w');
hold on
errorbar([1:1:3],M(4:6),S(4:6),'k.')
ylim([0 100])


for n=1:size(C1,1);
    hold on
    plot(C1(n,[4 5 6]),'Color','m','LineWidth',0.5)
end

[h1,p1] = ttest(C1(:,4),C1(:,5));
[h2,p2] = ttest(C1(:,4),C1(:,6));
[h3,p3] = ttest(C1(:,5),C1(:,6));

% Use bonferroni to correct for the multiple comparison
[isSignificant,adjusted_pvals,~]= bonferroni_holm([p1 p2 p3],0.05); % adjusted_pvals are the corrected p values
%adjusted_pvals
title(['PoorChR2RS',num2str(adjusted_pvals([1 2 3]))]);
% title(['PoorChR2RS',num2str([p1 p3])]);

C2=[];% all mice together
C2(:,[1 2])=C1(:,[2 3])./C1(:,1);
C2(:,[3 4])=C1(:,[5 6])./C1(:,4);
C2(C2==inf)=NaN;
M=nanmean(C2,1);
S=nansem(C2,1);

subplot(2,6,5) %change amplitude per sample
bar([1 2],M([1 3]),'FaceColor','w');
hold on
errorbar([1 2],M([1 3]),S([1 3]),'k.')
for n=1:size(C2,1);
hold on
plot(C2(n,[1 3]),'Color','m','LineWidth',0.5)
end
[r,p]=ttest(C2(:,1),C2(:,3));
title(['change amp Sti',num2str(p)])

ylim([0.3 2.1])
save('poorIndivSession.mat','C1','C2')
% individual mice
% 
CHR2Mouse=[];
for n=1:length(CHR2);
    CHR2Mouse(n,:)=nanmean(CHR2{n},1);
end
% for n=1:length(CHR2);
%     A=CHR2{n};
%     A(:,[1 2 3])= A(:,[1 2 3])./A(:,1);
%      A(:,[4 5 6])= A(:,[4 5 6])./A(:,4);
%     CHR2Mouse(n,:)=nanmean(A,1);
% end
C1=CHR2Mouse(poor,:);% Poor mice
M=nanmean(C1,1);
S=nansem(C1,1);

subplot(2,6,2) %per mouse CS
bar([1:1:3],M(1:3),'FaceColor','w');
hold on
errorbar([1:1:3],M(1:3),S(1:3),'k.');
ylim([0 100])

for n=1:size(C1,1);
    hold on
    plot(C1(n,[1 2 3]),'Color','m','LineWidth',0.5)
end
% ylim([20 100])
[h1,p1] = ttest(C1(:,1),C1(:,2));
[h2,p2] = ttest(C1(:,1),C1(:,3));
[h3,p3] = ttest(C1(:,2),C1(:,3));

% Use bonferroni to correct for the multiple comparison
[isSignificant,adjusted_pvals,~]= bonferroni_holm([p1 p2 p3],0.05); % adjusted_pvals are the corrected p values
%adjusted_pvals
title(['PoorChR2CS',num2str(adjusted_pvals([1 2 3]))]);

subplot(2,6,4) %per mouse CS
bar([1:1:3],M(4:6),'FaceColor','w');
hold on
errorbar([1:1:3],M(4:6),S(4:6),'k.');

for n=1:size(C1,1);
    hold on
    plot(C1(n,[4 5 6]),'Color','m','LineWidth',0.5)
end
% ylim([20 100])
[h1,p1] = ttest(C1(:,4),C1(:,5));
[h2,p2] = ttest(C1(:,4),C1(:,6));
[h3,p3] = ttest(C1(:,5),C1(:,6));
ylim([0 100])

% Use bonferroni to correct for the multiple comparison
[isSignificant,adjusted_pvals,~]= bonferroni_holm([p1 p2 p3],0.05); % adjusted_pvals are the corrected p values
%adjusted_pvals
title(['PoorChR2RS',num2str(adjusted_pvals([1 2 3]))]);

subplot(2,6,6);%per mouse change amplitude

C2=[];% 
C2(:,[1 2])=C1(:,[2 3])./C1(:,1);
C2(:,[3 4])=C1(:,[5 6])./C1(:,4);
C2(C2==inf)=NaN;
M=nanmean(C2,1);
S=nansem(C2,1);

bar([1 2],M([1 3]),'FaceColor','w');
hold on
errorbar([1 2],M([1 3]),S([1 3]),'k.')
for n=1:size(C2,1);
hold on
plot(C2(n,[1 3]),'Color','m','LineWidth',0.5)
end
[r,p]=ttest(C2(:,1),C2(:,3));
title(['change amp Sti',num2str(p)])
ylim([0.3 2.1])
save('poorIndivMice.mat','C1','C2')
%%

load("Control_5P5G-5cm2ms-PS.mat");
poor=[6 7 8 9 10];% 40, 43, 31, 30, 42

C=CONTROL(poor);
C1=cell2mat(C);% all mice together
M=nanmean(C1,1);
S=nansem(C1,1);

subplot(2,6,7) %per sample CS
bar([1:1:3],M(1:3),'FaceColor','w');
hold on
errorbar([1:1:3],M(1:3),S(1:3),'k.')

ylim([0 100])

for n=1:size(C1,1);
    hold on
    plot(C1(n,[1 2 3]),'Color','m','LineWidth',0.5)
end

[h1,p1] = ttest(C1(:,1),C1(:,2));
[h2,p2] = ttest(C1(:,1),C1(:,3));
[h3,p3] = ttest(C1(:,2),C1(:,3));

[isSignificant,adjusted_pvals,~]= bonferroni_holm([p1 p2 p3],0.05); % adjusted_pvals are the corrected p values
%adjusted_pvals
title(['PoorCtlCS',num2str(adjusted_pvals([1 2 3]))]);

subplot(2,6,9) %per sample RS
bar([1:1:3],M(4:6),'FaceColor','w');
hold on
errorbar([1:1:3],M(4:6),S(4:6),'k.')
ylim([0 100])


for n=1:size(C1,1);
    hold on
    plot(C1(n,[4 5 6]),'Color','m','LineWidth',0.5)
end

[h1,p1] = ttest(C1(:,4),C1(:,5));
[h2,p2] = ttest(C1(:,4),C1(:,6));
[h3,p3] = ttest(C1(:,5),C1(:,6));

% Use bonferroni to correct for the multiple comparison
[isSignificant,adjusted_pvals,~]= bonferroni_holm([p1 p2 p3],0.05); % adjusted_pvals are the corrected p values
%adjusted_pvals
title(['PoorCtlRS',num2str(adjusted_pvals([1 2 3]))]);

C2=[];% all mice together
C2(:,[1 2])=C1(:,[2 3])./C1(:,1);
C2(:,[3 4])=C1(:,[5 6])./C1(:,4);
C2(C2==inf)=NaN;
M=nanmean(C2,1);
S=nansem(C2,1);

subplot(2,6,11) %change amplitude per sample
bar([1 2],M([1 3]),'FaceColor','w');
hold on
errorbar([1 2],M([1 3]),S([1 3]),'k.')
for n=1:size(C2,1);
hold on
plot(C2(n,[1 3]),'Color','m','LineWidth',0.5)
end
[r,p]=ttest(C2(:,1),C2(:,3));
title(['change amp Sti',num2str(p)])
ylim([0.3 2.1])
save('poorIndivSessionGFP.mat','C1','C2')

% individual mice

CONTROLMouse=[];
for n=1:length(CONTROL);
    CONTROLMouse(n,:)=nanmean(CONTROL{n},1);
end

% for n=1:length(CONTROL);
%     A=CONTROL{n};
%     A(:,[1 2 3])= A(:,[1 2 3])./A(:,1);
%      A(:,[4 5 6])= A(:,[4 5 6])./A(:,4);
%     CONTROLMouse(n,:)=nanmean(A,1);
% end

C1=CONTROLMouse(poor,:);% Poor mice

M=nanmean(C1,1);
S=nansem(C1,1);

subplot(2,6,8) %per mouse CS
bar([1:1:3],M(1:3),'FaceColor','w');
hold on
errorbar([1:1:3],M(1:3),S(1:3),'k.');
ylim([0 100])

for n=1:size(C1,1);
    hold on
    plot(C1(n,[1 2 3]),'Color','m','LineWidth',0.5)
end
% ylim([20 100])
[h1,p1] = ttest(C1(:,1),C1(:,2));
[h2,p2] = ttest(C1(:,1),C1(:,3));
[h3,p3] = ttest(C1(:,2),C1(:,3));

% Use bonferroni to correct for the multiple comparison
[isSignificant,adjusted_pvals,~]= bonferroni_holm([p1 p2 p3],0.05); % adjusted_pvals are the corrected p values
%adjusted_pvals
title(['PoorCtlCS',num2str(adjusted_pvals([1 2 3]))]);

subplot(2,6,10) %per mouse CS
bar([1:1:3],M(4:6),'FaceColor','w');
hold on
errorbar([1:1:3],M(4:6),S(4:6),'k.');

for n=1:size(C1,1);
    hold on
    plot(C1(n,[4 5 6]),'Color','m','LineWidth',0.5)
end
% ylim([20 100])
[h1,p1] = ttest(C1(:,4),C1(:,5));
[h2,p2] = ttest(C1(:,4),C1(:,6));
[h3,p3] = ttest(C1(:,5),C1(:,6));
ylim([0 100])

% Use bonferroni to correct for the multiple comparison
[isSignificant,adjusted_pvals,~]= bonferroni_holm([p1 p2 p3],0.05); % adjusted_pvals are the corrected p values
%adjusted_pvals
title(['PoorCtlRS',num2str(adjusted_pvals([1 2 3]))]);

subplot(2,6,12);%per mouse change amplitude

C2=[];% 
C2(:,[1 2])=C1(:,[2 3])./C1(:,1);
C2(:,[3 4])=C1(:,[5 6])./C1(:,4);
C2(C2==inf)=NaN;
M=nanmean(C2,1);
S=nansem(C2,1);

bar([1 2],M([1 3]),'FaceColor','w');
hold on
errorbar([1 2],M([1 3]),S([1 3]),'k.')
for n=1:size(C2,1);
hold on
plot(C2(n,[1 3]),'Color','m','LineWidth',0.5)
end
[r,p]=ttest(C2(:,1),C2(:,3));
title(['change amp Sti',num2str(p)])
ylim([0.3 2.1])
save('poorIndivMiceGFP.mat','C1','C2')
saveas(gcf,'forPaperPoor.fig')