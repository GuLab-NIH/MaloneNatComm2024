function [trials,roll,stats] = plotVisFun(file,st,trialsN,params,plotOn,rollON,lickTh)
%%

buffSt = params.buffSt;
buffEnd = params.buffEnd;
timeCue = params.timeCue;
rewCue = 4;
punCue = 5;

if nargin<7 || isempty(lickTh)
    lickTh = 0.5;
end

nRoll = 50;

fields = {'like','unlike'};
fCues = [3,2];

if isfile([file(1:end-3) 'mat'])
    load([file(1:end-3) 'mat'],'logData')
else
    logData = readLog(file,10);
    save([file(1:end-3) 'mat'],'logData')
end

params = logParamsVis(logData,lickTh);

lickIdxSt = params.lickIdx;
tLickSt = params.tLick;
worldIdx = params.worldIdx;

t = params.t;
R = logData(:,8);
L = logData(:,9);
W = logData(:,10);

lickIdxAll = find(L>lickTh);
tLickAll = t(lickIdxAll);


%% Process worlds

lickWorld = W(lickIdxAll);

if plotOn
    
    figure; hold on
    
    plot(t,W)
    
    starts = ismember(lickIdxAll,lickIdxSt);
    
    scatter(tLickAll(~starts),lickWorld(~starts),25,'ko','Filled')
    scatter(tLickSt,lickWorld(starts),'ro','LineWidth',2)
end


%% Process trials

trials = struct();

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

roll = struct();

if rollON
    
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
        roll.dprime(i) = norminv(roll.HR(i))-norminv(roll.FA(i));
    end
    
    if plotOn
        figure; hold on
        
        plot(roll.x,roll.HR,'-g')
        plot(roll.x,roll.FA,'-r')
        
        plot(roll.x,roll.dprime,'.b')
        
        legend('Hit Rate','False Alarms','d''')
    end
    
end


%% Limit trials

if trialsN==-1
    trialsN = trials.N;
end

% beginning
% st = 1;
% trialsN = 60;
trialsType = trials.type(st:st+trialsN-1);
trialsIsLick = trials.isLick(st:st+trialsN-1);
trialsCue = trials.cue(st:st+trialsN-1);

% end
% trials.N = 75;
% trials.type = trials.type(end-trials.N+1:end);
% trials.isLick = trials.isLick(end-trials.N+1:end);


%% Calculate stats for all cues

RC = sum(trialsType==1 & trialsIsLick==1);
RI = sum(trialsType==1 & trialsIsLick~=1);
RT = RC + RI;
HR = RC/RT;
if HR==1
    HR = 0.99;
elseif HR==0
    HR = 0.01;
end

EC = sum(trialsType==-1 & trialsIsLick~=1);
EI = sum(trialsType==-1 & trialsIsLick==1);
ET = EC + EI;

FA = EI/ET;
if FA==1
    FA = 0.99;
elseif FA==0
    FA = 0.01;
end

stats = struct();

stats.all.nTotal = trials.N;
stats.all.nRew = RT;
stats.all.perCorrect = (RC+EC)/(RT+ET);
stats.all.HR = HR;
stats.all.FA = FA;
stats.all.dprime = norminv(HR)-norminv(FA);

if plotOn
    fprintf(['nRuns = '  num2str(stats.all.nTotal) '\n'])
    fprintf(['d'' = '  num2str(stats.all.dprime,2) '\n\n'])
end


%% Calculate stats for specific cues

for i = 1:length(fields)
    
    EC = sum(trialsCue==fCues(i) & trialsIsLick~=1);
    EI = sum(trialsCue==fCues(i) & trialsIsLick==1);
    ET = EC + EI;
    
    FA = EI/ET;
    if FA==1
        FA = 0.99;
    elseif FA==0
        FA = 0.01;
    end
    
    stats.(fields{i}).nTotal = RT+ET;
    stats.(fields{i}).nRew = RT;
    stats.(fields{i}).perCorrect = (RC+EC)/(RT+ET);
    stats.(fields{i}).HR = HR;
    stats.(fields{i}).FA = FA;
    stats.(fields{i}).dprime = norminv(HR)-norminv(FA);
    
end


