function [params] = logParamsVis(logData,lickTh)
%this code extract parameters from virmen logs on linear track
% Output:
%params.y:y position
%params.lickIdx=licking indices
%params.rewardIdx=the indices when the reward is delivered

%% Read data

a = logData;
t = (a(:,1)-a(1,1))*24*60*60;   % time
r = a(:,8);                     % reward
b = a(:,9);                     % licking
w = a(:,10);                    % world


%% Find reward locations

rewardIdx = find(r==1);         % reward indices
tReward = t(rewardIdx);         % time of reward


%% Find lick locations

if nargin<2 || isempty(lickTh)
    lickTh = 0.5;
end

b(b<lickTh)=0;
b(b>=lickTh)=1;%turning licking signal to 1 and 0

if any(b>0.5)
    n=contiguous(b,1);  % find contiguous 1
    nn=n{1,2};
    lickIdx=nn(:,1);    % licking indices
    tLick=t(lickIdx);   % time of licking
else
    lickIdx=nan;
    tLick=nan;
end


%% Find world numbers

worldIdx = contiguous(w);


%% save parameters

params.t = t;
params.rewardIdx = rewardIdx;
params.tReward = tReward;
params.lickIdx = lickIdx;
params.tLick = tLick;
params.worldIdx = worldIdx;


end
