%% plotLocationConsistencyCombined

clear all; close all; clc


%% Plot day-by-day spatial correlation

folds = findSubF('RBRLoc.mat',1,[],0);

% useIdx = [2 3];
% type = 'bad';

useIdx = [1 4 5 6];
type = 'good';

maxDay = 10;

RBRLocComb = cell(1,maxDay);

for ii = 1:length(useIdx)
    load(folds{useIdx(ii)})
    
    for jj = 1:length(RBRLoc)
        RBRLocComb{jj} = [RBRLocComb{jj};RBRLoc{jj}];
    end
end

save(['combined\RBRLocComb_' type '.mat'],'RBRLocComb')



%%

load('D:\Mice\cueTemp.mat');

binWidth=5;
figure
for n=1:length(RBRLocComb)
    
    A = RBRLocComb{n};
    M=nanmean(A,1);
    E=nansem(A,1);
    
    subplot(3,4,n);
    errorbar([3:1:length(M)+2]*binWidth,M,E,'k')
    hold on
    plot([1:1:length(cueTemp)]*binWidth,cueTemp*max(M),'m');
    hold on
    r=366;
    line([r r],[0 max(M)],'Color','g','LineWidth',2);
    title(['day ',num2str(n)]);
    xlabel('Track position');
    ylabel('RBR Corr');
end

saveas(gcf,['combined\SpCorrDays_' type '.fig']);


%% Plot combined spatial correlation

load('D:\Mice\cueTemp.mat');
binWidth=5;
figure

c = 1;
d = 9;

A = cat(1,RBRLocComb{c:d});
M=nanmean(A,1);
E=nansem(A,1);

errorbar([3:1:length(M)+2]*binWidth,M,E,'k')
hold on
plot([1:1:length(cueTemp)]*binWidth,cueTemp*max(M),'m');
hold on
r=366;
line([r r],[0 max(M)],'Color','g','LineWidth',2);
title(['day' num2str(c) ' to day' num2str(d)]);
xlabel('Track position');
ylabel('RBR Corr');

saveas(gcf,['combined\SpCorrAll_' type '.fig']);


%% Plot day-by-day consistency

folds = findSubF('RBRAll.mat',1,[],0);

RBRComb = cell(1,maxDay);

for ii = 1:length(useIdx)
    load(folds{useIdx(ii)})
    
    for jj = 1:length(RBR)
        RBRComb{jj} = [RBRComb{jj};RBR{jj}];
    end
end

save(['combined\RBRComb_' type '.mat'],'RBRComb')


%%

M = [];
E = [];

for n=1:length(RBRComb)
    M(n) = nanmean(RBRComb{n});
    E(n) = nansem(RBRComb{n}');
end

figure;
errorbar(1:length(M),M,E,'k')

xlabel('Day');
ylabel('RBR Corr');
title('RBR Corr');

saveas(gcf,['combined\RBRDays_' type '.fig']);
