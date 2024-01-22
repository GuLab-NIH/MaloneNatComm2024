%% visAccRead.m
% analyze visual accuity discrimination task. Calculates the percent of
% trials correct and d-prime statistics


%% Set paramters

wReward = 4;
wError = 2:3;

d = dir('*T*.txt');

logData = readLog(d(end).name,10);

params = logParamsVis(logData);

trials = struct();


%% Identify reward trial outcomes

trials.reward.correct = 0;
trials.reward.incorrect = 0;

for i = 1:length(wReward)
    wIDXS = params.worldIdx{wReward(i),2};
    
    for j = 1:size(wIDXS,1)
        wIdx = wIDXS(j,:);
        
        if any(params.lickIdx >= wIdx(1) & params.lickIdx <= wIdx(2))
            trials.reward.correct = trials.reward.correct+1;
        else
            trials.reward.incorrect = trials.reward.incorrect+1;
        end
    end
end

trials.reward.total = trials.reward.correct + trials.reward.incorrect;


%% Identify error trial outcomes

trials.error.correct = 0;
trials.error.incorrect = 0;

for i = 1:length(wError)
    wIDXS = params.worldIdx{wError(i),2};
    
    for j = 1:size(wIDXS,1)
        wIdx = wIDXS(j,:);
        
        if any(params.lickIdx >= wIdx(1) & params.lickIdx <= wIdx(2))
            trials.error.incorrect = trials.error.incorrect+1;
        else
            trials.error.correct = trials.error.correct+1;
        end
    end
end

trials.error.total = trials.error.correct + trials.error.incorrect;


%% Calculate stats

RC = trials.reward.correct;
RI = trials.reward.incorrect;
RT = trials.reward.total;
HR = RC/RT;

EC = trials.error.correct;
EI = trials.error.incorrect;
ET = trials.error.total;
FA = EI/ET;

stats.perCorrect = (RC+EC)/(RT+ET);
stats.HR = HR;
stats.FA = FA;
stats.dprime = norminv(HR)-norminv(FA);

