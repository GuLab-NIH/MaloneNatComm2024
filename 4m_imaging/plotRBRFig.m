%% Load RBR by group

clear; close all; clc

load('RBRComb_good.mat'); %mean of the consistency as a function of track location
RBRG = cat(2,RBRComb(:));
meanRBRG = cellfun(@(x) mean(x,1,'omitnan'),RBRComb);
semRBRG = cellfun(@(x) nansem(x,1),RBRComb);

load('RBRComb_bad.mat'); %mean of the consistency as a function of track location
RBRB = cat(2,RBRComb(:));
meanRBRB = cellfun(@(x) mean(x,1,'omitnan'),RBRComb);
semRBRB = cellfun(@(x) nansem(x,1),RBRComb);


%% Plot day-by-day consistency

figure; hold on

colors = {'g','m'};

for ii = 1:2
    if ii==1
        useMean = meanRBRG;
        useSEM = semRBRG;
    else
        useMean = meanRBRB;
        useSEM = semRBRB;
    end
    
    errorbar(1:length(useMean),useMean,useSEM,colors{ii})

end

xlabel('Day');
ylabel('RBR Corr');
title('RBR Corr');
ylim([0 0.6])
xlim([0 11])

saveas(gcf,'RBRDays_comb.fig');


%% Run anova

% create label matices
dayG = cell(size(RBRB));
dayB = cell(size(RBRB));

for ii = 1:length(dayG)
    RBRG{ii}(isnan(RBRG{ii})) = [];
    RBRB{ii}(isnan(RBRB{ii})) = [];
    
    dayG{ii} = ii*ones(size(RBRG{ii}));
    dayB{ii} = ii*ones(size(RBRB{ii}));
end

catG = cat(1,RBRG{:});
catB = cat(1,RBRB{:});

dayCatG = cat(1,dayG{:});
dayCatB = cat(1,dayB{:});
numG = 1*ones(size(catG));
numB = 2*ones(size(catB));

catAll = [catG;catB];
dayCatAll = [dayCatG;dayCatB];
numCat = [numG;numB];

[p,tbl,stats,terms] = anovan(catAll,{dayCatAll,numCat},'model','interaction','varnames',{'day','type'});


%% Calculate correlation stats

[rG,pG] = corr((1:10)',meanRBRG');
[rB,pB] = corr((1:10)',meanRBRB');

[~,pG2] = corr(dayCatG,catG);
[~,pB2] = corr(dayCatB,catB);

