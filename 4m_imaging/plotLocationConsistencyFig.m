%% plotLocationConsistencyFig

clear all; close all; clc


%% Load RBRLoc by group

load('RBRLocComb_good.mat'); %mean of the consistency as a function of track location
locG = RBRLocComb;
meanRBR = cellfun(@(x) mean(x,1,'omitnan'),locG,'UniformOutput',false);
semRBR = cellfun(@(x) nansem(x,1),locG,'UniformOutput',false);
meanG = cat(1,meanRBR{:});
semG = cat(1,semRBR{:});

load('RBRLocComb_bad.mat');
locB = RBRLocComb;
meanRBR = cellfun(@(x) mean(x,1,'omitnan'),locB,'UniformOutput',false);
semRBR = cellfun(@(x) nansem(x,1),locB,'UniformOutput',false);
meanB = cat(1,meanRBR{:});
semB = cat(1,semRBR{:});


%% Plot day-by-day spatial correlation

load('D:\4m_imaging_analysis\cueTemp.mat');

binWidth=5;
day = 9;
colors = {'g','m'};

for n = day
    figure; hold on
    
    for ii = 1:2
        if ii==1
            useMean = meanG(n,:);
            useSEM = semG(n,:);
        else
            useMean= meanB(n,:);
            useSEM = semB(n,:);
        end
        
        errorbar((1:length(useMean))*binWidth,useMean,useSEM,colors{ii})
    end
    
    ylim([0 0.4])
    
    mx = max(ylim);
    plot((1:length(cueTemp))*binWidth,cueTemp*mx,'color',[0.5 0.5 0.5]);
    r=366;
    line([r r],[0 mx],'Color','r','LineWidth',1);
    title(['day ',num2str(n)]);
    xlabel('Track position');
    ylabel('RBR Corr');
    
end


%% Plot day-by-day in-cue/out-cue differences

day = 9;
colors = {'g','m'};

load('cueTemp.mat');
beforeReward = 67:73;%before reward
temp = cueTemp;
temp(beforeReward) = 2;
% tempUse = temp(3:78);
% tempUse = temp(1:76);
tempUse = temp(1:73);


outCue = find(tempUse==0);%out cue bins
inCueNoR = find(tempUse==1); %in cue but no reward
inRNoCue = find(tempUse==2); % in reward no cue

maxDays = 10;

R = cell(1,3);
M = zeros(maxDays,3);
E = zeros(maxDays,3);
p = zeros(maxDays,3);


for n=1:maxDays
    
    diffGBcur = nanmean(locG{n})-nanmean(locB{n});
    
    %use diffGB
    inCueNoRBinDiff = diffGBcur(inCueNoR);
    beforeRBinDiff = diffGBcur(inRNoCue);
    outCueBinDiff = diffGBcur(outCue);
    
    R{1}(:,n) = outCueBinDiff;
    R{2}(:,n) = inCueNoRBinDiff;
    R{3}(:,n) = beforeRBinDiff;
    
    M(n,1) = mean(beforeRBinDiff);
    M(n,2) = mean(inCueNoRBinDiff);
    M(n,3) = mean(outCueBinDiff);
    
    E(n,1) = nansem(beforeRBinDiff,2);
    E(n,2) = nansem(inCueNoRBinDiff,2);
    E(n,3) = nansem(outCueBinDiff,2);
    
    [~,p(n,1)] = ttest2(beforeRBinDiff,inCueNoRBinDiff);
    [~,p(n,2)] = ttest2(beforeRBinDiff,outCueBinDiff);
    [~,p(n,3)] = ttest2(inCueNoRBinDiff,outCueBinDiff);
    
    
    %%
    if n==day
        figure; hold on
        plot((1:length(useMean))*binWidth,diffGBcur,'k')
        
        
        title(['day' num2str(n)]);
        xlabel('Track position');
        ylabel('RBR Corr');
        ylim([-.05 0.3])
        
        mx = max(ylim);
        mn = min(ylim);
        
        plot((1:length(cueTemp))*binWidth,cueTemp*(mx-mn)+mn,'color',[0.5 0.5 0.5]);
        
        r=366;
        line([r r],[mn mx],'Color','r','LineWidth',1);
        
    end
    
end



%% Plot day-by-day bar graph

figure; hold on

b = bar(M,'grouped');

% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(M);

% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end

% Plot the errorbars
errorbar(x',M,E,'k','linestyle','none');

title('in-cue/out-cue differences');

% legend('Before Reward','In Cue','Out Cue')
saveas(gcf,'spCorrDifDays_bins');
ylim([0 0.4])


%% Combine all days

allRBR = cell(2,1);

allRBR{1} = cat(1,locG{:});
allRBR{2} = cat(1,locB{:});

meanAllRBR = cellfun(@(x) mean(x,1,'omitnan'),allRBR,'UniformOutput',false);
diffAll = meanAllRBR{1}-meanAllRBR{2};
    
%use diffGB
outCueBinDiffAll = diffAll(outCue);
inCueNoRBinDiffAll = diffAll(inCueNoR);
beforeRBinDiffAll = diffAll(inRNoCue);

M = [];
M(1) = mean(outCueBinDiffAll);
M(2) = mean(inCueNoRBinDiffAll);
M(3) = mean(beforeRBinDiffAll);

E = [];
E(1) = nansem(outCueBinDiffAll,2);
E(2) = nansem(inCueNoRBinDiffAll,2);
E(3) = nansem(beforeRBinDiffAll,2);
  
pAll = [];
[~,pAll(1)] = ttest2(inCueNoRBinDiffAll,outCueBinDiffAll);
[~,pAll(2)] = ttest2(beforeRBinDiffAll,outCueBinDiffAll);
[~,pAll(3)] = ttest2(beforeRBinDiffAll,inCueNoRBinDiffAll);
% [0.00421211158106985,0.000741715519512066,0.152959281017870]

figure; hold on

b = bar(M,'grouped');

% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(M);

% Get the x coordinate of the bars
x = b.XEndPoints;

% Plot the errorbars
errorbar(x',M,E,'k','linestyle','none');

title(['in-cue/out-cue differences All']);

