%% dfof_byLoc
% Plots rolling average of dfof by track location for a batch of mice


clear all; close all; clc


%% Calculate spatial dfof rolling average for all FOV

base = 'D:\4m_imaging_analysis';
cd(base)

N=1;    %rolling average
type = 'dfofaveragesmooth';

% suff = '';
suff = '_sig';

svName = [type suff '_' num2str(N)];

folds = findSubF('pcaica',3,[],0);

load('D:\4m_imaging_analysis\cueTemp.mat');

for ff = 1:length(folds)
    cd(folds{ff})
    
    d = dir(['dfofaveragesmooth' suff '*_cells.mat']);
    load(d(1).name);
    
    if exist('dfofaveragesmooth_sig','var')
        dfofaveragesmooth = dfofaveragesmooth_sig;
    end
    
    dfofInfoLocation = struct();
    dfofInfoLocation.dfofMean = [];
    dfofInfoLocation.meanAll=[];
    dfofInfoLocation.semAll=[];
    
    B = movmean(dfofaveragesmooth,[0 N-1],1)';
    dfofInfoLocation.dfofMean = B(:,1:end-N+1);
    
    dfofInfoLocation.meanAll = nanmean(dfofInfoLocation.dfofMean,1);
    dfofInfoLocation.semAll = nansem(dfofInfoLocation.dfofMean,1);
    
    save(['dfofInfoLocation' suff '.mat'],'dfofInfoLocation');
    
%     figure; hold on
%     errorbar(1:length(dfofInfoLocation.meanAll),dfofInfoLocation.meanAll,dfofInfoLocation.semAll,'k')
%     plot(1:length(cueTemp),cueTemp*max(dfofInfoLocation.meanAll),'m');
%     
%     saveas(gcf,'dfofInfoLocation.fig');
%     close
    
    cd(base)
end


%% Calculate spatial dfof per mouse

base = 'D:\4m_imaging_analysis';
cd(base)

load('D:\4m_imaging_analysis\cueTemp.mat');

% Find all mice
d = dir('*22*');
mouseNames = {d.name}';

skipDays = {[],[],[],4,4,[]};


for mm = 1:length(mouseNames)
    %% Identify all fov per mouse per day
    
    cd(mouseNames{mm})
    
    % Find all days
    d = dir();
    isub = [d(:).isdir]; %# returns logical vector
    nameDays = {d(isub).name}';
    nameDays(ismember(nameDays,{'.','..'})) = [];
    
    % Find all locs
    allLocs = {};
    
    for ff = 1:length(nameDays)
        cd(nameDays{ff})
        
        d = dir('loc*');
        isub = [d(:).isdir]; %# returns logical vector
        subLocs = {d(isub).name}';
        
        for gg = 1:length(subLocs)
            curF = [pwd '\' subLocs{gg} '\pcaica'];
            if isfolder(curF)
                allLocs{ff,gg} = curF;
            end
        end
        
        cd([base '\' mouseNames{mm}])
    end

    
    %% Calculate day-by-day spatial correlation
    
    DFOFLoc = cell(1,size(allLocs,1));
    
    p=pwd;
    binWidth=5;
%     figure
    nn = 1;
    for n=1:size(allLocs,1)
        if ismember(n,skipDays{mm})
            continue
        end
        A=[];
        for m=1:size(allLocs,2)
            if isempty(allLocs{n,m})
                continue
            end
            cd(allLocs{n,m});
            load(['dfofInfoLocation' suff '.mat']);
            load('roiIdxUse.mat')
            B = dfofInfoLocation.dfofMean(roiIdxUse,:);
            
            A(end+1:end+size(B,1),:)=B;
        end
%         M=nanmean(A,1);
%         E=nansem(A,1);
%         
%         subplot(3,4,nn);
%         errorbar([3:1:length(M)+2]*binWidth,M,E,'k')
%         hold on
%         plot((1:length(cueTemp))*binWidth,cueTemp*max(M),'m');
%         hold on
%         r=366;
%         line([r r],[0 max(M)],'Color','g','LineWidth',2);
%         title(['day ',num2str(nn)]);
%         xlabel('Track position');
%         ylabel('DFOF');
        
        DFOFLoc{n} = A;
        
        nn = nn+1;
    end
    
    cd([base '\' mouseNames{mm}])

%     saveas(gcf,'SpDFOFDays.fig');
%     close;
        
    DFOFLoc(skipDays{mm}) = [];
    save('DFOFLoc.mat','DFOFLoc')
    
    cd(base)
    
end


%% Combine data learner group

folds = findSubF('DFOFLoc.mat',1,[],0);

maxDay = 10;

useIdx = {[1 4 5 6],[2 3]};
type = {'good','bad'};
DFOFLocComb = cell(2,maxDay);


for t = 1:2
    for ii = 1:length(useIdx{t})
        load(folds{useIdx{t}(ii)})
        
        for jj = 1:length(DFOFLoc)
            DFOFLocComb{t,jj} = [DFOFLocComb{t,jj};DFOFLoc{jj}];
        end
    end   
end

meanDFOF = cellfun(@(x) mean(x,1,'omitnan'),DFOFLocComb,'UniformOutput',false);
semDFOF = cellfun(@(x) nansem(x,1),DFOFLocComb,'UniformOutput',false);

% save('combined\DFOFLocComb.mat','DFOFLocComb','type')


%% Plot day-by-day dfof

load('D:\4m_imaging_analysis\cueTemp.mat');

binWidth=5;
day = 9;
colors = {'g','m'};

figure; hold on

for ii = 1:2
    useMean = meanDFOF{ii,day};
    useSEM = semDFOF{ii,day};
    
    errorbar((1:length(useMean))*binWidth,useMean,useSEM,colors{ii})
end

% ylim([0 0.15])

mx = max(ylim);
plot((1:length(cueTemp))*binWidth,cueTemp*mx,'color',[0.5 0.5 0.5]);
r=366;
line([r r],[0 mx],'Color','r','LineWidth',1);
title(['day ',num2str(day)]);
xlabel('Track position');
ylabel('DFOF');

% calculate difference
diffGBcur = meanDFOF{1,day}-meanDFOF{2,day};

figure; hold on
plot((1:length(useMean))*binWidth,diffGBcur,'k')

title(['day' num2str(day)]);
xlabel('Track position');
ylabel('RBR Corr');
% ylim([-.05 0.3])

mx = max(ylim);
mn = min(ylim);

plot((1:length(cueTemp))*binWidth,cueTemp*(mx-mn)+mn,'color',[0.5 0.5 0.5]);

r=366;
line([r r],[mn mx],'Color','r','LineWidth',1);
line([min(xlim) max(xlim)],[0 0],'Color','k','LineStyle',':','LineWidth',0.5)


%% Plot day-by-day bar graph

load('D:\4m_imaging_analysis\cueTemp.mat');

% process cue locations
beforeReward = 67:73;%before reward
temp = cueTemp;
temp(beforeReward) = 2;
tempUse = temp(1:73);

% define catergory bins
outCue = find(tempUse==0);%out cue bins
inCueNoR = find(tempUse==1); %in cue but no reward
inRNoCue = find(tempUse==2); % in reward no cue

nDays = 10;

R = cell(1,3);
M = zeros(nDays,3);
E = zeros(nDays,3);
p = zeros(nDays,3);

% claculate day-by-day data
for n=1:nDays
    
    diffGBcur = meanDFOF{1,n}-meanDFOF{2,n};
    
    %use diffGB
    outCueBinDiff = diffGBcur(outCue);
    inCueNoRBinDiff = diffGBcur(inCueNoR);
    beforeRBinDiff = diffGBcur(inRNoCue);
    
    R{1}(:,n) = outCueBinDiff;
    R{2}(:,n) = inCueNoRBinDiff;
    R{3}(:,n) = beforeRBinDiff;
    
    M(n,1) = mean(outCueBinDiff);
    M(n,2) = mean(inCueNoRBinDiff);
    M(n,3) = mean(beforeRBinDiff);
    
    E(n,1) = nansem(outCueBinDiff,2);
    E(n,2) = nansem(inCueNoRBinDiff,2);
    E(n,3) = nansem(beforeRBinDiff,2);
    
    [~,p(n,1)] = ttest2(inCueNoRBinDiff,outCueBinDiff);
    [~,p(n,2)] = ttest2(beforeRBinDiff,outCueBinDiff);
    [~,p(n,3)] = ttest2(beforeRBinDiff,inCueNoRBinDiff);
end

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

title([svName ': in-cue/out-cue differences']);

% legend('Before Reward','In Cue','Out Cue')
saveas(gcf,['dfofFigs\binsByDay_' svName]);
% ylim([0 0.4])


%% Combine all days

allDFOF = cell(2,1);
for ii = 1:2
    allDFOF{ii} = cat(1,DFOFLocComb{ii,:});
end

meanAllDFOF = cellfun(@(x) mean(x,1,'omitnan'),allDFOF,'UniformOutput',false);
semAllDFOF = cellfun(@(x) nansem(x,1),allDFOF,'UniformOutput',false);

diffAll = meanAllDFOF{1}-meanAllDFOF{2};
    
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
% [0.228244197848860,0.00614356649957464,0.00305677567595555]

figure; hold on

b = bar(M,'grouped');

% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(M);

% Get the x coordinate of the bars
x = b.XEndPoints;

% Plot the errorbars
errorbar(x',M,E,'k','linestyle','none');

title([svName ': in-cue/out-cue differences All']);
