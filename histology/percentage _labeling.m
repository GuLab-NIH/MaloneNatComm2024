%% This script runs percentage of labeling
GFP = [80.23255814	80.859375	24.03100775	29.80769231;48.23529412	68.7150838	15.29411765	30];
ChR2 = [60.6741573	35.7615894	3.370786517	1.492537313;46.42857143	36.79245283	5.952380952	3.521126761];

Reelin_GFP = nanmean(GFP(:,1));
Reelin_ChR2 = nanmean(ChR2(:,1));
Cal_GFP = nanmean(GFP(:,3));
Cal_ChR2 = nanmean(ChR2(:,3));

GFP_Reelin = nanmean(GFP(:,2));
ChR2_Reelin = nanmean(ChR2(:,2));
GFP_Cal = nanmean(GFP(:,4));
ChR2_Cal = nanmean(ChR2(:,4));

Reelin_GFP_sem = std(GFP(:,1))/sqrt(length(GFP(:,1)));
Reelin_ChR2_sem = std(ChR2(:,1))/sqrt(length(ChR2(:,1)));
Cal_GFP_sem = std(GFP(:,3))/sqrt(length(GFP(:,3)));
Cal_ChR2_sem = std(ChR2(:,3))/sqrt(length(ChR2(:,3)));

GFP_Reelin_sem = std(GFP(:,2))/sqrt(length(GFP(:,2)));
ChR2_Reelin_sem = std(ChR2(:,2))/sqrt(length(ChR2(:,2)));
GFP_Cal_sem = std(GFP(:,4))/sqrt(length(GFP(:,4)));
ChR2_Cal_sem = std(ChR2(:,4))/sqrt(length(ChR2(:,4)));


figure;

subplot(2,1,1)
x1 = 1;
y1 = Reelin_GFP;
e1 = Reelin_GFP_sem;

x2 = 2;
y2 = Cal_GFP;
e2 = Cal_GFP_sem;

x3 = 3.5;
y3 = GFP_Reelin;
e3 = GFP_Reelin_sem;

x4 = 4.5;
y4 = GFP_Cal;
e4 = GFP_Cal_sem;

b = bar(x1,y1,'FaceColor','flat','EdgeColor','flat');
hold on
er = errorbar(x1,y1,e1);
er.Color = [0 0 0];                            
er.LineStyle = 'none';
b.CData = [1, 1, 1];
b.CData = [0.5,0.5,0.5];
hold on

b = bar(x2,y2,'FaceColor','flat','EdgeColor','flat');
b.CData = [1, 1, 1];
b.CData = [0.5,0.5,0.5];
hold on
er = errorbar(x2,y2,e2);
er.Color = [0 0 0];                            
er.LineStyle = 'none';
hold on

b = bar(x3,y3,'FaceColor','flat','EdgeColor','flat');
b.CData = [1, 1, 1];
b.CData = [0.5,0.5,0.5];
hold on
er = errorbar(x3,y3,e3);
er.Color = [0 0 0];                            
er.LineStyle = 'none';
hold on

b = bar(x4,y4,'FaceColor','flat','EdgeColor','flat');
b.CData = [1, 1, 1];
b.CData = [0.5,0.5,0.5];
hold on
er = errorbar(x4,y4,e4);
er.Color = [0 0 0];                            
er.LineStyle = 'none';
box off
xticks([1 2 3.5 4.5])
xticklabels({'R+/GFP+','C+/GFP+','GFP+/R+','GFP+/C+'});

% plot chr2
subplot(2,1,2)
x1 = 1;
y1 = Reelin_ChR2;
e1 = Reelin_ChR2_sem;

x2 = 2;
y2 = Cal_ChR2;
e2 = Cal_ChR2_sem;

x3 = 3.5;
y3 = ChR2_Reelin;
e3 = ChR2_Reelin_sem;

x4 = 4.5;
y4 = ChR2_Cal;
e4 = ChR2_Cal_sem;

b = bar(x1,y1,'FaceColor','flat','EdgeColor','flat');
hold on
er = errorbar(x1,y1,e1);
er.Color = [0 0 0];                            
er.LineStyle = 'none';
b.CData = [1, 1, 1];
b.CData = [0.5,0.5,0.5];
hold on

b = bar(x2,y2,'FaceColor','flat','EdgeColor','flat');
b.CData = [1, 1, 1];
b.CData = [0.5,0.5,0.5];
hold on
er = errorbar(x2,y2,e2);
er.Color = [0 0 0];                            
er.LineStyle = 'none';
hold on

b = bar(x3,y3,'FaceColor','flat','EdgeColor','flat');
b.CData = [1, 1, 1];
b.CData = [0.5,0.5,0.5];
hold on
er = errorbar(x3,y3,e3);
er.Color = [0 0 0];                            
er.LineStyle = 'none';
hold on

b = bar(x4,y4,'FaceColor','flat','EdgeColor','flat');
b.CData = [1, 1, 1];
b.CData = [0.5,0.5,0.5];
hold on
er = errorbar(x4,y4,e4);
er.Color = [0 0 0];                            
er.LineStyle = 'none';
box off
xticks([1 2 3.5 4.5])
xticklabels({'R+/GFP+','C+/GFP+','GFP+/R+','GFP+/C+'});

ylim([0 80])
yticks([0 40 80])
yticklabels({'0','40','80'});

%% add individual data points
% to load data('Z:\SNM\labMembers\NT\paper\histology_cell type optogenetics\histology.mat')
% open the figure
hold on
subplot(2,1,1)

Reelin_GFP_data = GFP(:,1);
Cal_GFP_data = GFP(:,3);
GFP_Reelin_data = GFP(:,2);
GFP_Cal_data = GFP(:,4);

x = [1 2 3.5 4.5];
x1 = x-0.2;
x2 = x+0.2;
plot([x1(1);x2(1)],Reelin_GFP_data,'ok','MarkerSize',2)
hold on
plot([x1(2);x2(2)],Cal_GFP_data,'ok','MarkerSize',2)
hold on
plot([x1(3);x2(3)],GFP_Reelin_data,'ok','MarkerSize',2)
hold on
plot([x1(4);x2(4)],GFP_Cal_data,'ok','MarkerSize',2)
hold on
%%
subplot(2,1,2)
Reelin_ChR2_data = ChR2(:,1);
Cal_ChR2_data = ChR2(:,3);
ChR2_Reelin_data = ChR2(:,2);
ChR2_Cal_data = ChR2(:,4);

x = [1 2 3.5 4.5];
x1 = x-0.2;
x2 = x+0.2;
plot([x1(1);x2(1)],Reelin_ChR2_data,'ok','MarkerSize',2)
hold on
plot([x1(2);x2(2)],Cal_ChR2_data,'ok','MarkerSize',2)
hold on
plot([x1(3);x2(3)],ChR2_Reelin_data,'ok','MarkerSize',2)
hold on
plot([x1(4);x2(4)],ChR2_Cal_data,'ok','MarkerSize',2)