%% visPlotSum

clear all; close all; clc

load('dataVis.mat')
useDays = 1:10;

cats = {[1 3:5 8 10:14],[2 6 7 9]};
catNames = {'Good Learners','Bad Learners'};
catCols = {[0 1 0],[1 0 1]};

yLab = 'd-prime';


%% Plot time course

figure; hold on

for ii = 1:2
    dataCur = dataVis(useDays,cats{ii});
    nMouse = size(dataCur,2);
    
    for jj = 1:nMouse
        plot(useDays,dataCur(:,jj),'-','LineWidth',1,...
            'Color',[catCols{ii} 0.2])
    end
end

for ii = 1:2
    dataCur = dataVis(useDays,cats{ii});
    nMouse = size(dataCur,2);
    
    dataAvg = mean(dataCur,2);
    dataStd = std(dataCur,0,2)/sqrt(nMouse);
    
    plot(useDays,dataAvg,'LineWidth',1.5,'Color',0.9*catCols{ii})
    errorbar(useDays,dataAvg,dataStd,'LineWidth',1.5,...
        'Color',0.9*catCols{ii},'LineStyle','none')
end

plot(useDays,ones(size(useDays)),'r','LineWidth',1)

% pVal = stat(dataVis(useDays,cats{1}),dataVis(useDays,cats{2}));

set(gca,'FontSize',18)
xlabel('Day')
ylabel(yLab)
xlim([min(useDays) max(useDays)])
xlim([0 11])
set(gca,'xtick',useDays)
axis('square')


%% Plot correlations

load('dataLick.mat')
load('dataSlow.mat')

lickAvg = mean(dataLick,1);
slowAvg = mean(dataSlow,1);
visAvg = mean(dataVis(useDays,:),1);

figure; hold on

for ii = 1:2
    subplot(1,2,ii); hold on
    
    if ii==1
        dataX = lickAvg;
        xLab = 'Licking percentage';
    elseif ii==2
        dataX = slowAvg;
        xLab = 'Slow down percentile';
    end
    
    for jj = 1:2
        curX = dataX(cats{jj});
        curY = visAvg(cats{jj});
        
        scatter(curX,curY,50,0.9*catCols{jj},'Filled')
        
        xlabel(xLab)
        ylabel(yLab)
    end
    
    [r,p] = corr(dataX',visAvg');
    
    disp([xLab ': r=' num2str(r,2) ', p=' num2str(p,2)])
    
%     title(['r=' num2str(r,2) ', p=' num2str(p,2)])
%     axis('square')
    set(gca,'FontSize',18)
    
end

