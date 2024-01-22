%% compare different cell types
% load the data
clear;clc;
goodF.R = [0,1,1,0.923100000000000,1,0.800000000000000,0.666700000000000,0.800000000000000,1,1,0.750000000000000];
goodN.R = [0.647100000000000,0.680000000000000,0.636400000000000,0.800000000000000,0.533300000000000,1,0.400000000000000,0.727300000000000,0.782600000000000,0.857100000000000,0.818200000000000,0.615400000000000,0.800000000000000,0.444400000000000,0.400000000000000,1];
badF.R = [0.800000000000000,0.666700000000000,0.500000000000000,1,0.200000000000000,0.500000000000000,0.800000000000000,1,0.727300000000000,0.750000000000000,0.888900000000000,0.750000000000000,0.571400000000000,0.636400000000000];
badN.R = [0.800000000000000,0.375000000000000,0.222200000000000,0.857100000000000,0.666700000000000,0.545500000000000,0.923100000000000,0.785700000000000,0.909100000000000,0.714300000000000,0.818200000000000];

goodF.C = [0.111100000000000,0.372500000000000,0.142900000000000,0,0,0];
goodN.C = [0,0,0,0,0.166700000000000];
badF.C = [0.0909000000000000,0,0.0526000000000000];
badN.C = [0,0.454500000000000,0.0667000000000000,0.0769000000000000];

goodF.G = [0.0909000000000000,0,0];
goodN.G = [0,0.0714000000000000,0];
badF.G = [0 0];
badN.G = [0 0];

%% compare good learners
envF = repmat('F',length([goodF.R goodF.C goodF.G]),1);
envN = repmat('N',length([goodN.R goodN.C goodN.G]),1);
envBF = repmat('F',length([badF.R badF.C badF.G]),1);
envBN = repmat('N',length([badN.R badN.C badN.G]),1);
% env = [envF; envN]; % two way anova
env = [envF; envN;envBF;envBN]; % three way anova

cellFR = repmat('R',length(goodF.R),1);
cellFC = repmat('C',length(goodF.C),1);
cellFG = repmat('G',length(goodF.G),1);
cellNR = repmat('R',length(goodN.R),1);
cellNC = repmat('C',length(goodN.C),1);
cellNG = repmat('G',length(goodN.G),1);

cellBFR = repmat('R',length(badF.R),1);
cellBFC = repmat('C',length(badF.C),1);
cellBFG = repmat('G',length(badF.G),1);
cellBNR = repmat('R',length(badN.R),1);
cellBNC = repmat('C',length(badN.C),1);
cellBNG = repmat('G',length(badN.G),1);
% cell = [cellFR;cellFC;cellFG;cellNR;cellNC;cellNG]; % two way anova
cell = [cellFR;cellFC;cellFG;cellNR;cellNC;cellNG;...
    cellBFR;cellBFC;cellBFG;cellBNR;cellBNC;cellBNG];% three way anova

goodmouse = repmat('g',length([goodF.R goodF.C goodF.G goodN.R goodN.C goodN.G] ),1);
badmouse = repmat('b',length([badF.R badF.C badF.G badN.R badN.C badN.G] ),1);

mouse = [goodmouse;badmouse];
% data = [goodF.R';goodF.C';goodF.G';goodN.R';goodN.C';goodN.G']; % two way anova
data = [goodF.R';goodF.C';goodF.G';goodN.R';goodN.C';goodN.G';...
    badF.R';badF.C';badF.G';badN.R';badN.C';badN.G']; % three way anova
% p = anovan(data,{env cell},'model','full','varnames',{'env','cell types'})
[p,tbl,stats] = anovan(data,{env cell mouse},'model','full','varnames',{'env','cell types','performance'}) ;% three way 
[results,~,~,gnames] = multcompare(stats,"Dimension",[1 2 3],"ctype",'bonferroni');

%% compare cell types
R = [goodF.R goodN.C badF.R badN.R];
C = [goodF.C goodN.C badF.C badN.C];
G = [goodF.G goodN.G badF.G badN.G];

meanR = nanmean(R);
semR = nanstd(R,1)/sqrt(length(R));
meanC = nanmean(C);
semC = nanstd(C,1)/sqrt(length(C));
meanG = nanmean(G);
semG = nanstd(G,1)/sqrt(length(G));

x1 = 1;
y1 = meanR;
x2 = 2;
y2 = meanC;
x3 = 3;
y3 = meanG;

b = bar(x1,y1,'FaceColor','flat','EdgeColor','flat');
hold on
e1 = semR;
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
e2 = semC;
er = errorbar(x2,y2,e2);
er.Color = [0 0 0];                            
er.LineStyle = 'none';
hold on
b = bar(x3,y3,'FaceColor','flat','EdgeColor','flat');
hold on
e3 = semG;
er = errorbar(x3,y3,e3);
er.Color = [0 0 0];                            
er.LineStyle = 'none';
b.CData = [1, 1, 1];
b.CData = [0.5,0.5,0.5];
b.CData = [1, 1, 1];
b.CData = [0.5,0.5,0.5];

xticks([1 2 3])
xticklabels({'Reelin+','Calbindin+','Gad67+'})

box off
hold on
ylabel('Predictive Licks')

set(gcf, 'Position',  [100, 100, 300, 400]);


[h1,p1] = ttest2(R,C);
[h2,p2] = ttest2(R,G);
[h3,p3] = ttest2(C,G);
data = [p1 p2 p3];
[isSignificant,adjusted_pvals,alpha]= bonferroni_holm(data,0.05)
