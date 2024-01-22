%% This script plot histology bar graphs
% goodF condition
clear;clc;
% load('raw_data.mat')

load('raw_data_individual.mat')
goodF.Reelin = mean(goodF_R);
goodN.Reelin = mean(goodN_R);
badF.Reelin = mean(badF_R);
badN.Reelin = mean(badN_R);

goodF.Cal = mean(goodF_C);
goodN.Cal = mean(goodN_C);
badF.Cal = mean(badF_C);
badN.Cal = mean(badN_C);

goodF.Gad = mean(goodF_G);
goodN.Gad = mean(goodN_G);
badF.Gad = mean(badF_G);
badN.Gad = mean(badN_G);

x = 1;
y1 = [goodF.Reelin,goodF.Cal,goodF.Gad];
if sum(y1)<1
    goodF.Unknown = 1-sum(y1);
    y = [goodF.Reelin,goodF.Cal,goodF.Gad,goodF.Unknown];
else
    y = [goodF.Reelin,goodF.Cal,goodF.Gad];
end

b1 = bar(x,y,'stacked','FaceColor','flat');
colormap summer;
color = summer(6);
b1(1).CData = color(1,:);
b1(2).CData = color(2,:);
b1(3).CData = color(3,:);

if sum(y1)<1
    b1(4).CData = color(4,:);
else
end
hold on
%% goodN condition
clear;clc;
load('raw_data.mat')
x = 2;
y1 = [goodN.Reelin,goodN.Cal,goodN.Gad];
if sum(y1)<1
    goodN.Unknown = 1-sum(y1);
    y = [goodN.Reelin,goodN.Cal,goodN.Gad,goodN.Unknown];
else
    y = [goodN.Reelin,goodN.Cal,goodN.Gad];
end

b1 = bar(x,y,'stacked','FaceColor','flat','HandleVisibility','off');
colormap summer;
color = summer(6);
b1(1).CData = color(1,:);
b1(2).CData = color(2,:);
b1(3).CData = color(3,:);

if sum(y1)<1
    b1(4).CData = color(4,:);
else
end

hold on
box off
ylim([0 1.1])
xticks([1 2])
xticklabels({'Familiar','Novel'})
yticks([0 0.25 0.5 0.75 1])
yticklabels({'0','0.25','0.5','0.75','1'})

legend('R+','C+','G+','U')
set(gcf, 'Position',  [100, 100, 300, 400]);

%% badF condition
clear;clc;
load('raw_data.mat')
figure;
x = 1;
y1 = [badF.Reelin,badF.Cal,badF.Gad];
if sum(y1)<1
    badF.Unknown = 1-sum(y1);
    y = [badF.Reelin,badF.Cal,badF.Gad,badF.Unknown];
else
    y = [badF.Reelin,badF.Cal,badF.Gad];
end

b1 = bar(x,y,'stacked','FaceColor','flat','HandleVisibility','off');
colormap spring;
color = spring(6);
b1(1).CData = color(1,:);
b1(2).CData = color(2,:);
b1(3).CData = color(3,:);

if sum(y1)<1
    b1(4).CData = color(4,:);
else
end
hold on
%% goodN condition
clear;clc;
load('raw_data.mat')
x = 2;
y1 = [badN.Reelin,badN.Cal,badN.Gad];
if sum(y1)<1
    badN.Unknown = 1-sum(y1);
    y = [badN.Reelin,badN.Cal,badN.Gad,badN.Unknown];
else
    y = [badN.Reelin,badN.Cal,badN.Gad];
end

b1 = bar(x,y,'stacked','FaceColor','flat');
colormap spring;
color = spring(6);
b1(1).CData = color(1,:);
b1(2).CData = color(2,:);
b1(3).CData = color(3,:);

if sum(y1)<1
    b1(4).CData = color(4,:);
else
end

hold on
box off
ylim([0 1.1])
xticks([1 2])
xticklabels({'Familiar','Novel'})
yticks([0 0.25 0.5 0.75 1])
yticklabels({'0','0.25','0.5','0.75','1'})

legend('R+','C+','G+','U')
set(gcf, 'Position',  [100, 100, 300, 400]);
