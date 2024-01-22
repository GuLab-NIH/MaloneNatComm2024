%% plotLocationConsistency

clear all; clc


%% Find location folders

p = pwd;

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
        curF = [p '\' nameDays{ff} '\' subLocs{gg} '\pcaica'];
        if isfolder(curF)
            allLocs{ff,gg} = curF;
        end
    end
    
    cd(p)
end

% allLocs = allLocs(1:8,:);
skipDays = 4;
% skipDays = [];

%% Plot day-by-day spatial correlation

RBRLoc = cell(1,size(allLocs,1));

load('D:\Mice\cueTemp.mat');
p=pwd;
binWidth=5;
figure
nn = 1;
for n=1:size(allLocs,1)
    if ismember(n,skipDays)
        continue
    end
    A=[];
    for m=1:size(allLocs,2)
        if isempty(allLocs{n,m})
            continue
        end
        cd(allLocs{n,m});
        load('RunByRun_dfof\corrInfoLocation.mat');
        load('roiIdxUse.mat')
        B = corrInfoLocation.toOthersMean(roiIdxUse,:);
        
        A(end+1:end+size(B,1),:)=B;
    end
    M=nanmean(A,1);
    E=nansem(A,1);
    
    subplot(3,4,nn);
    errorbar([3:1:length(M)+2]*binWidth,M,E,'k')
    hold on
    plot([1:1:length(cueTemp)]*binWidth,cueTemp*max(M),'m');
    hold on
    r=366;
    line([r r],[0 max(M)],'Color','g','LineWidth',2);
    title(['day ',num2str(nn)]);
    xlabel('Track position');
    ylabel('RBR Corr');
    
    RBRLoc{n} = A;
    
    nn = nn+1;
end

cd(p)
saveas(gcf,'SpCorrDays.fig');


RBRLoc(skipDays) = [];
save('RBRLoc.mat','RBRLoc')


%% Plot day-by-day consistency

RBRLoc = cell(1,size(allLocs,1));

p=pwd;
binWidth=5;

M = [];
E = [];

for n=1:size(allLocs,1)
    A=[];
    for m=1:size(allLocs,2)
        if isempty(allLocs{n,m})
            continue
        end
        cd(allLocs{n,m});
        load('RunByRun_dfof\corrInfo.mat');
        load('roiIdxUse.mat')
        B = corrInfo.meantoOthers(roiIdxUse);
        
        A = [A,B];
    end
    M(n) = nanmean(A);
    E(n) = nansem(A);
    
    RBR{n} = A';
end

M(skipDays) = [];
E(skipDays) = [];

cd(p)
figure;
errorbar(1:length(M),M,E,'k')

xlabel('Day');
ylabel('RBR Corr');
title('RBR Corr');

saveas(gcf,'RBRDays.fig');

RBR(skipDays) = [];
RBRM = cellfun(@(x) mean(x,'omitnan'),RBR);
RBRSEM = cellfun(@(x) std(x,'omitnan')/sqrt(sum(~isnan(x))),RBR);

save('RBRAll.mat','RBR','RBRM','RBRSEM')


%% Plot combined spatial correlation

load('D:\Mice\cueTemp.mat');
p=pwd;
binWidth=5;
figure,

A = [];
c = 1;
d = 10;%size(allLocs,1);

for n=c:d
    for m=1:size(allLocs,2)
        if isempty(allLocs{n,m})
            continue
        end
        cd(allLocs{n,m});
        load('RunByRun_dfof\corrInfoLocation.mat');
        load('roiIdxUse.mat')
        B = corrInfoLocation.toOthersMean(roiIdxUse,:);
        
        A(end+1:end+size(B,1),:)=B;
    end
end

cd(p)

allC = A;
save('allC.mat','allC')

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

saveas(gcf,'SpCorrAll.fig');