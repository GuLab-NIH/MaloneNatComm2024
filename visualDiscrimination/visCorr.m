%% VisCorr

clear all; close all; clc

load('dataLick.mat')
load('dataSlow.mat')
load('dataVis.mat')

% yLabs = {'d'' (day 1-10)','d'' (day 11-20)','d'' (day 1-20)'};
% ySets = {1:10,11:20,1:20};

yLabs = {'d-prime'};
ySets = {1:10};
nSets = length(ySets);

lickAvg = mean(dataLick,1);
slowAvg = mean(dataSlow,1);

visAvg = cell(1,nSets);

for ii = 1:nSets
    visAvg{ii} = mean(dataVis(ySets{ii},:),1);
end

figure; hold on

for ii = 1:2
    if ii==1
        curX = lickAvg;
        labX = 'Licking percentage';
    elseif ii==2
        curX = slowAvg;
        labX = 'Slow down percentile';
    end
    
    for jj = 1:nSets
        subplot(nSets,2,nSets*(ii-1)+jj); hold on
        
        scatter(curX,visAvg{jj})
        xlabel(labX)
        ylabel(yLabs{jj})
        
        [r,p] = corr(curX',visAvg{jj}');
        
%         title(['r=' num2str(r,2) ', p=' num2str(p,2)])
        axis('square')
        set(gca,'FontSize',18)
    end
end


sgtitle('Predictive Behavior vs. Visual Discrimination')




