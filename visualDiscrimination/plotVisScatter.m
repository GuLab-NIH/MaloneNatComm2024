%% plotVisScatter
% genterate scatter plot of lick locations by trial for each world.

% clear; close all; clc


%% Set parameters

day = 20;                            % training day

buffITI = 2;                        % amount of ITI included in plot
timeCue = 6;                        % total cue display time
buffSt = 1;                         % time buffer before reward zone
buffEnd = 1;                        % time buffer after reward zone
timeTot = buffITI + timeCue;        % total plot time

times = 0:0.005:timeTot;        % post-interpolation plot times

cueW = [4 2];                   % reward and false alarm cue worlds
cueTitle = {'Reward','False Alarm'};

% day 1 params
if day==1
    timeCue = 5;
    buffSt = 2.5;
    buffEnd = 0.5;
end

actWin = timeCue-buffSt-buffEnd;    % active window time


%% Load Data

% identfy valid log files
d = dir('*T*.txt');

% select file to plot (manually selected)
file = d(day).name;

% load or creat .mat file from log file
if isfile([file(1:end-3) 'mat'])
    load([file(1:end-3) 'mat'],'logData')
else
    logData = readLog(file,10);
    save([file(1:end-3) 'mat'],'logData')
end

% extract parameters from log file
params = logParamsVis(logData);

% extract required parameters
t = params.t;
L = logData(:,9);
W = logData(:,10);

worldIdx = params.worldIdx;
N = size(worldIdx{1,2},1)-1;

% binarize lick signal
L = L>=0.5;


%% Process trials

% initialize trials struct
trials = struct();

cues = zeros(N,1);

for i = 1:N
    
    % calculate trial indices
    idxSt = worldIdx{1,2}(i,1);
    idxEnd = worldIdx{1,2}(i+1,1)-1;
    idxCue = find(W(idxSt:idxEnd)~=1,1)+idxSt-1;
    idxPre = find(t>(t(idxCue)-buffITI),1);
    
    % svae trial indices
    trials(i).idx = [idxSt idxPre idxCue idxEnd];
    
    % calculate trial cue
    trials(i).cue = W(idxCue);
    cues(i) = W(idxCue);

    % calculate trial times
    tTrial = t(idxPre:idxEnd)-t(idxPre);
    
    % calculate trial lick signal
    licks = double(L(idxPre:idxEnd));

    % interpolate lick signal and re-binarize
    tLick = interp1(tTrial,licks,times)>=0.5;
    
    d = [false diff(tLick)==0];
    
    tLick(d & tLick==1) = 0;
    
    % save trial lick signal
    trials(i).tLick = tLick;
    
end


%% Plot licks by trial

figure

colors = [205,255,204; 255,204,203; 224,237,246]/255;

for i = 1:2
    cueCur = cueW(i);
    cueIdx = find(cues==cueCur);
    cueN = length(cueIdx);
        
    subplot(1,2,i); hold on
    
    pos = [0 0 timeCue cueN];
    rectangle('Position',pos,'FaceColor',colors(3,:),...
        'EdgeColor','none')
    
    pos = [buffSt 0 actWin cueN];
    rectangle('Position',pos,'FaceColor',colors(i,:),...
        'EdgeColor','none')
    
    for j = 1:cueN
        curTimes = times(trials(cueIdx(j)).tLick==1);
        curTrial = j*ones(size(curTimes));
        scatter(curTimes-buffITI,curTrial,36,'k.');
    end
    
    if i==1
        plot([timeCue-buffEnd timeCue-buffEnd],[0 cueN],'r','LineWidth',2)
    end
    
    axis('square')
    xlim([-buffITI,timeCue])
    ylim([1 cueN])
    
    title(cueTitle{i})
    
end
    







