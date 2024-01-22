%%

clear; close all; clc

files = findSubF('behavior_NT.mat',1,[],0);

groups = {[1 4 5 6],[2 3]};
skipDays = {[0 4 4 0],[0 0]};

nG = length(groups);
nDays = 10;

combLick = cell(nG,1);
combSlow = cell(nG,1);
combLickMean = zeros(nG,nDays);
combSlowMean = zeros(nG,nDays);
combLickSEM = zeros(nG,nDays);
combSlowSEM = zeros(nG,nDays);

for ii = 1:nG
    nCur = length(groups{ii});
    
    tempLick = zeros(nCur,nDays);
    tempSlow = zeros(nCur,nDays);
    
    for jj = 1:nCur
        load(files{groups{ii}(jj)})
        
        skCur = skipDays{ii}(jj);
        if skCur~=0
            pLick(skCur) = [];
            pSlow(skCur) = [];
        end
        
        tempLick(jj,:) = pLick;
        tempSlow(jj,:) = pSlow;
    end
    
    combLick{ii} = tempLick*100;
    combSlow{ii} = tempSlow;
    combLickMean(ii,:) = nanmean(tempLick*100,1);
    combSlowMean(ii,:) = nanmean(tempSlow,1);
    combLickSEM(ii,:) = nansem(tempLick*100,1);
    combSlowSEM(ii,:) = nansem(tempSlow,1);
    
end


%% Plot combined behavior

figure; hold on

colors = {'g','m'};
ylabels = {'Predictive Licking (%)','Predictive Slowing (%)'};

ii = 2;

if ii==1
    useData = combLick;
    useMean = combLickMean;
    useSEM = combLickSEM;
else
    useData = combSlow;
    useMean = combSlowMean;
    useSEM = combSlowSEM;
end

for jj = 1:2
    errorbar(1:length(useMean),useMean(jj,:),useSEM(jj,:),colors{jj})
end

xlabel('Day');
ylabel(ylabels{ii});
ylim([30 90])
xlim([0 11])
title(ylabels{ii});

p = anovaRM2W(useData{1},useData{2});
[rG,pG] = corr((1:length(useMean))',useMean(1,:)');
[rB,pB] = corr((1:length(useMean))',useMean(2,:)');
disp([p pG pB])
disp([rG rB])