%% plotVisTime

% clear; close all; clc


%% Initialize parameters

start = 1;
stop = -1;
stop = 98;
plotOn = 0;
rollOn = 0;
stDay = 1;

params = struct();

% set day 1 params
params(1).buffSt = 2.5;
params(1).buffEnd = 0.5;
params(1).timeCue = 5;

% set day 2 params
params(2).buffSt = 1;
params(2).buffEnd = 1;
params(2).timeCue = 6;

d = dir('*T*.txt');
nFile = length(stDay:length(d));

days = struct();
HR = struct();
FA = struct();
dPrime = struct();
trials = struct();


%% Read files

fields = {'all','like','unlike'};

mx = 4;

for i = 1:length(fields)
    days.(fields{i}) = zeros(nFile,1);
    HR.(fields{i}) = zeros(nFile,1);
    FA.(fields{i}) = zeros(nFile,1);
    dPrime.(fields{i}) = zeros(nFile,1);
    trials.(fields{i}) = zeros(nFile,1);
    
    for f = 1:nFile
        if f==1; p = 1; else; p = 2; end
        
        [~,~,stats] = plotVisFun(d(f+stDay-1).name,start,stop,...,
            params(p),rollOn,plotOn);
        
        days.(fields{i})(f) = stats.(fields{i}).nTotal;
        HR.(fields{i})(f) = stats.(fields{i}).HR;
        FA.(fields{i})(f) = stats.(fields{i}).FA;
        dPrime.(fields{i})(f) = stats.(fields{i}).dprime;
        trials.(fields{i})(f) = stats.(fields{i}).nTotal;
        
        if  stats.(fields{i}).dprime~=Inf
            mx = max(mx,ceil(stats.(fields{i}).dprime));
        end
        
    end
        
end


%% Save results

data = struct();

data.days = days;
data.HR = HR;
data.FA = FA;
data.dPrime = dPrime;
data.trials = trials;

save('data.mat','data');


%% Plot results

figure

for i = 1:length(fields)
    subplot(1,length(fields),i)
    hold on
    
    plot(HR.(fields{i}),'-g')
    plot(FA.(fields{i}),'-r')
    
    scatter(1:nFile,dPrime.(fields{i}),'b')
    
    title(['Visual Discrimination: ' fields{i}])
    
    legend('Hit Rate','False Alarms','d''','Location','northwest')
    ylim([0 mx])
    set(gca,'FontSize',18)
end


savefig('VisTime_sub.fig')
