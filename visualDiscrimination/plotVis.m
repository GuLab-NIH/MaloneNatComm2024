

%%

clear; close all; clc

buffSt = 1;
buffEnd = 1;
timeCue = 6;
rewCue = 4;
punCue = 5;

nRoll = 50;
plotOn = 1;

d = dir('*T*.txt');

file = d(20).name;

if isfile([file(1:end-3) 'mat'])
    load([file(1:end-3) 'mat'],'logData')
else
    logData = readLog(file,10);
    save([file(1:end-3) 'mat'],'logData')
end

params = logParamsVis(logData);

lickIdxSt = params.lickIdx;
tLickSt = params.tLick;
worldIdx = params.worldIdx;

t = params.t;
R = logData(:,8);
L = logData(:,9);
W = logData(:,10);

lickIdxAll = find(L>0.5);
tLickAll = t(lickIdxAll);

len = length(t);
nLick = length(lickIdxAll);


%% Process worlds

lickWorld = W(lickIdxAll);

if plotOn
    
    figure; hold on
    
    plot(t/60,W)
    
    starts = ismember(lickIdxAll,lickIdxSt);
    
    scatter(tLickAll(~starts)/60,lickWorld(~starts),25,'ko','Filled')
    scatter(tLickSt/60,lickWorld(starts),'ro','LineWidth',2)
end


%% Process trials

trials.N = size(worldIdx{1,2},1)-1;
trials.idx = zeros(trials.N,2);
trials.idx(:,1) = worldIdx{1,2}(1:end-1,1);
trials.idx(:,2) = worldIdx{1,2}(2:end,1)-1;

trials.cueIdx = zeros(trials.N,1);
trials.cue = zeros(trials.N,1);

trials.lickRng = zeros(trials.N,2);
trials.isLick = zeros(trials.N,1);

for i = 1:trials.N
    
    idx1 = trials.idx(i,1);
    idx2 = trials.idx(i,2);
    
    idxCue = find(W(idx1:idx2)~=1,1)+idx1-1;
    trials.cueIdx(i,1) = idxCue;
    trials.cue(i) = W(idxCue);
    
    
    idxL1 = find(t(idx1:idx2)>=t(idxCue)+buffSt,1)+idx1-1;
    idxL2 = find(t(idx1:idx2)<=t(idxCue)+timeCue-buffEnd,1,'last')+idx1-1;
    
    trials.lickRng(i,1) = idxL1;
    trials.lickRng(i,2) = idxL2;
    
    trials.isLick(i) = any(lickIdxAll>=idxL1 & lickIdxAll<=idxL2);
    
    if any(W(idx1:idx2)==punCue)
        trials.isLick(i) = 1;
    end
    
end

trials.type = ones(trials.N,1);
trials.type(trials.cue~=rewCue) = -1;


%% Calculate rolling stats

roll.idx = zeros(trials.N-nRoll+1,2);
roll.idx(:,1) = 1:trials.N-nRoll+1;
roll.idx(:,2) = nRoll:trials.N;
roll.x = mean(roll.idx,2);
roll.n = length(roll.x);

roll.HR = zeros(roll.n,1);
roll.FA = zeros(roll.n,1);
roll.dprime = zeros(roll.n,1);

for i = 1:roll.n
    RisLick = trials.isLick(roll.idx(i,1):roll.idx(i,2));
    Rtype = trials.type(roll.idx(i,1):roll.idx(i,2));
    
    roll.HR(i) = sum(Rtype==1 & RisLick==1)/sum(Rtype==1);
    roll.FA(i) = sum(Rtype==-1 & RisLick==1)/sum(Rtype==-1);
    
end

roll.HR(roll.HR==1) = 0.99;
roll.FA(roll.FA==1) = 0.99;
roll.HR(roll.HR==0) = 0.01;
roll.FA(roll.FA==0) = 0.01;

roll.dprime = norminv(roll.HR)-norminv(roll.FA);

if plotOn
    figure; hold on
    
    plot(roll.x,roll.HR,'-g')
    plot(roll.x,roll.FA,'-r')
    
    plot(roll.x,roll.dprime,'.b')
    
    legend('Hit Rate','False Alarms','d''')
end


%% Calculate stats

RC = sum(trials.type==1 & trials.isLick==1);
RI = sum(trials.type==1 & trials.isLick~=1);
RT = RC + RI;
HR = RC/RT;

if HR==1
    HR = 0.99;
elseif HR==0
    HR = 0.01;
end

EC = sum(trials.type==-1 & trials.isLick~=1);
EI = sum(trials.type==-1 & trials.isLick==1);
ET = EC + EI;
FA = EI/ET;

if FA==1
    FA = 0.99;
elseif FA==0
    FA = 0.01;
end

stats = struct();

stats.nTotal = trials.N;
stats.nRew = RT;
stats.perCorrect = (RC+EC)/(RT+ET);
stats.HR = HR;
stats.FA = FA;
stats.dprime = norminv(HR)-norminv(FA);

fprintf(['nRuns = '  num2str(stats.nTotal) '\n'])
fprintf(['d'' = '  num2str(stats.dprime,2) '\n\n'])
