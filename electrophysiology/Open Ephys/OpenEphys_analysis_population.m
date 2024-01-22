%% Population analysis for OpenEphys
% This analysis calculate the fold changes in firing rate from all the channels
% as a function of probe depth. Run this analysis after OpenEphys_analysis1
clear;clc;
firing.stim = nan(1,6);
firing.nostim = nan(1,6);

for recording = 1:6
    recording
    [fileName pathName]=uigetfile('*.mat','Select .mat file',...
        'D:\Data\optogenetics\Open Ephys\analized data\20211027\Analysis');

    file_location = [pathName fileName];
    load(file_location);
    
    firing.stim(recording) = nanmean(ephys.spike_stim_epoch,'all');
    firing.nostim(recording) = nanmean(ephys.spike_nostim_epoch,'all');
    
    % stimulus
    firing.spikeRateMean_2ms10Hz_stim{recording} = ephys.spikeRateMean_2ms10Hz_stim(:);
    firing.spikeRateMean_2ms20Hz_stim{recording} = ephys.spikeRateMean_2ms20Hz_stim(:);
    firing.spikeRateMean_2ms40Hz_stim{recording} = ephys.spikeRateMean_2ms40Hz_stim(:);
    
    firing.spikeRateMean_5ms10Hz_stim{recording} = ephys.spikeRateMean_5ms10Hz_stim(:);
    firing.spikeRateMean_5ms20Hz_stim{recording} = ephys.spikeRateMean_5ms20Hz_stim(:);
    firing.spikeRateMean_5ms40Hz_stim{recording} = ephys.spikeRateMean_5ms40Hz_stim(:);
    
    firing.spikeRateMean_10ms10Hz_stim{recording} = ephys.spikeRateMean_10ms10Hz_stim(:);
    firing.spikeRateMean_10ms20Hz_stim{recording} = ephys.spikeRateMean_10ms20Hz_stim(:);
    firing.spikeRateMean_10ms40Hz_stim{recording} = ephys.spikeRateMean_10ms40Hz_stim(:);
    
    firing.spikeRateMean_20ms10Hz_stim{recording} = ephys.spikeRateMean_20ms10Hz_stim(:);
    firing.spikeRateMean_20ms20Hz_stim{recording} = ephys.spikeRateMean_20ms20Hz_stim(:);
    firing.spikeRateMean_20ms40Hz_stim{recording} = ephys.spikeRateMean_20ms40Hz_stim(:);
    
    % no stimulus
    firing.spikeRateMean_2ms10Hz_nostim{recording} = ephys.spikeRateMean_2ms10Hz_nostim(:);
    firing.spikeRateMean_2ms20Hz_nostim{recording} = ephys.spikeRateMean_2ms20Hz_nostim(:);
    firing.spikeRateMean_2ms40Hz_nostim{recording} = ephys.spikeRateMean_2ms40Hz_nostim(:);
    
    firing.spikeRateMean_5ms10Hz_nostim{recording} = ephys.spikeRateMean_5ms10Hz_nostim(:);
    firing.spikeRateMean_5ms20Hz_nostim{recording} = ephys.spikeRateMean_5ms20Hz_nostim(:);
    firing.spikeRateMean_5ms40Hz_nostim{recording} = ephys.spikeRateMean_5ms40Hz_nostim(:);
    
    firing.spikeRateMean_10ms10Hz_nostim{recording} = ephys.spikeRateMean_10ms10Hz_nostim(:);
    firing.spikeRateMean_10ms20Hz_nostim{recording} = ephys.spikeRateMean_10ms20Hz_nostim(:);
    firing.spikeRateMean_10ms40Hz_nostim{recording} = ephys.spikeRateMean_10ms40Hz_nostim(:);
    
    firing.spikeRateMean_20ms10Hz_nostim{recording} = ephys.spikeRateMean_20ms10Hz_nostim(:);
    firing.spikeRateMean_20ms20Hz_nostim{recording} = ephys.spikeRateMean_20ms20Hz_nostim(:);
    firing.spikeRateMean_20ms40Hz_nostim{recording} = ephys.spikeRateMean_20ms40Hz_nostim(:);
    
end
%% plot the fold changes as a function of depth

firing.fold_changes = (firing.stim - firing.nostim)./firing.nostim;

fig = figure;
x = [7 6 5 4 3 2];
graph = barh(x,firing.fold_changes);

box off
yticks([2 3 4 5 6 7])
yticklabels({'2.1','2.0','1.9','1.8','1.7','1.6'})
ylabel({'Probe depth (mm)'})
xlabel({'Fold changes in number of spikes'})

graph.FaceColor = [.5 .5 .5];
graph.EdgeColor = [.5 .5 .5];

fig_name = strcat('20211027_fold_changes_removewavers');
saveas(fig,fig_name)

%% Plot overall rate changes, separating different stimuli conditions
for recording = 1:6 
    if recording ==1
        % stimulus
        spikeRateMean_2ms10Hz_stim = [firing.spikeRateMean_2ms10Hz_stim{recording}]';
        spikeRateMean_2ms20Hz_stim = [firing.spikeRateMean_2ms20Hz_stim{recording}]';
        spikeRateMean_2ms40Hz_stim = [firing.spikeRateMean_2ms40Hz_stim{recording}]';
        
        spikeRateMean_5ms10Hz_stim = [firing.spikeRateMean_5ms10Hz_stim{recording}]';
        spikeRateMean_5ms20Hz_stim = [firing.spikeRateMean_5ms20Hz_stim{recording}]';
        spikeRateMean_5ms40Hz_stim = [firing.spikeRateMean_5ms40Hz_stim{recording}]';
        
        spikeRateMean_10ms10Hz_stim = [firing.spikeRateMean_10ms10Hz_stim{recording}]';
        spikeRateMean_10ms20Hz_stim = [firing.spikeRateMean_10ms20Hz_stim{recording}]';
        spikeRateMean_10ms40Hz_stim = [firing.spikeRateMean_10ms40Hz_stim{recording}]';
        
        spikeRateMean_20ms10Hz_stim = [firing.spikeRateMean_20ms10Hz_stim{recording}]';
        spikeRateMean_20ms20Hz_stim = [firing.spikeRateMean_20ms20Hz_stim{recording}]';
        spikeRateMean_20ms40Hz_stim = [firing.spikeRateMean_20ms40Hz_stim{recording}]';
        
        % no stimulus
        spikeRateMean_2ms10Hz_nostim = [firing.spikeRateMean_2ms10Hz_nostim{recording}]';
        spikeRateMean_2ms20Hz_nostim = [firing.spikeRateMean_2ms20Hz_nostim{recording}]';
        spikeRateMean_2ms40Hz_nostim = [firing.spikeRateMean_2ms40Hz_nostim{recording}]';
        
        spikeRateMean_5ms10Hz_nostim = [firing.spikeRateMean_5ms10Hz_nostim{recording}]';
        spikeRateMean_5ms20Hz_nostim = [firing.spikeRateMean_5ms20Hz_nostim{recording}]';
        spikeRateMean_5ms40Hz_nostim = [firing.spikeRateMean_5ms40Hz_nostim{recording}]';
        
        spikeRateMean_10ms10Hz_nostim = [firing.spikeRateMean_10ms10Hz_nostim{recording}]';
        spikeRateMean_10ms20Hz_nostim = [firing.spikeRateMean_10ms20Hz_nostim{recording}]';
        spikeRateMean_10ms40Hz_nostim = [firing.spikeRateMean_10ms40Hz_nostim{recording}]';
        
        spikeRateMean_20ms10Hz_nostim = [firing.spikeRateMean_20ms10Hz_nostim{recording}]';
        spikeRateMean_20ms20Hz_nostim = [firing.spikeRateMean_20ms20Hz_nostim{recording}]';
        spikeRateMean_20ms40Hz_nostim = [firing.spikeRateMean_20ms40Hz_nostim{recording}]';
    else
        % stimulus
        spikeRateMean_2ms10Hz_stim = [spikeRateMean_2ms10Hz_stim [firing.spikeRateMean_2ms10Hz_stim{recording}]'];
        spikeRateMean_2ms20Hz_stim = [spikeRateMean_2ms20Hz_stim [firing.spikeRateMean_2ms20Hz_stim{recording}]'];
        spikeRateMean_2ms40Hz_stim = [spikeRateMean_2ms40Hz_stim [firing.spikeRateMean_2ms40Hz_stim{recording}]'];
        
        spikeRateMean_5ms10Hz_stim = [spikeRateMean_5ms10Hz_stim [firing.spikeRateMean_5ms10Hz_stim{recording}]'];
        spikeRateMean_5ms20Hz_stim = [spikeRateMean_5ms20Hz_stim [firing.spikeRateMean_5ms20Hz_stim{recording}]'];
        spikeRateMean_5ms40Hz_stim = [spikeRateMean_5ms40Hz_stim [firing.spikeRateMean_5ms40Hz_stim{recording}]'];
        
        spikeRateMean_10ms10Hz_stim = [spikeRateMean_10ms10Hz_stim [firing.spikeRateMean_10ms10Hz_stim{recording}]'];
        spikeRateMean_10ms20Hz_stim = [spikeRateMean_10ms20Hz_stim [firing.spikeRateMean_10ms20Hz_stim{recording}]'];
        spikeRateMean_10ms40Hz_stim = [spikeRateMean_10ms40Hz_stim [firing.spikeRateMean_10ms40Hz_stim{recording}]'];
        
        spikeRateMean_20ms10Hz_stim = [spikeRateMean_20ms10Hz_stim [firing.spikeRateMean_20ms10Hz_stim{recording}]'];
        spikeRateMean_20ms20Hz_stim = [spikeRateMean_20ms20Hz_stim [firing.spikeRateMean_20ms20Hz_stim{recording}]'];
        spikeRateMean_20ms40Hz_stim = [spikeRateMean_20ms40Hz_stim [firing.spikeRateMean_20ms40Hz_stim{recording}]'];
       
        % no stimulus
        spikeRateMean_2ms10Hz_nostim = [spikeRateMean_2ms10Hz_nostim [firing.spikeRateMean_2ms10Hz_nostim{recording}]'];
        spikeRateMean_2ms20Hz_nostim = [spikeRateMean_2ms20Hz_nostim [firing.spikeRateMean_2ms20Hz_nostim{recording}]'];
        spikeRateMean_2ms40Hz_nostim = [spikeRateMean_2ms40Hz_nostim [firing.spikeRateMean_2ms40Hz_nostim{recording}]'];
        
        spikeRateMean_5ms10Hz_nostim = [spikeRateMean_5ms10Hz_nostim [firing.spikeRateMean_5ms10Hz_nostim{recording}]'];
        spikeRateMean_5ms20Hz_nostim = [spikeRateMean_5ms20Hz_nostim [firing.spikeRateMean_5ms20Hz_nostim{recording}]'];
        spikeRateMean_5ms40Hz_nostim = [spikeRateMean_5ms40Hz_nostim [firing.spikeRateMean_5ms40Hz_nostim{recording}]'];
        
        spikeRateMean_10ms10Hz_nostim = [spikeRateMean_10ms10Hz_nostim [firing.spikeRateMean_10ms10Hz_nostim{recording}]'];
        spikeRateMean_10ms20Hz_nostim = [spikeRateMean_10ms20Hz_nostim [firing.spikeRateMean_10ms20Hz_nostim{recording}]'];
        spikeRateMean_10ms40Hz_nostim = [spikeRateMean_10ms40Hz_nostim [firing.spikeRateMean_10ms40Hz_nostim{recording}]'];
        
        spikeRateMean_20ms10Hz_nostim = [spikeRateMean_20ms10Hz_nostim [firing.spikeRateMean_20ms10Hz_nostim{recording}]'];
        spikeRateMean_20ms20Hz_nostim = [spikeRateMean_20ms20Hz_nostim [firing.spikeRateMean_20ms20Hz_nostim{recording}]'];
        spikeRateMean_20ms40Hz_nostim = [spikeRateMean_20ms40Hz_nostim [firing.spikeRateMean_20ms40Hz_nostim{recording}]'];
    end
end

%% Plot all recordings based on stimulus parameters
fig = figure;

% for 2 ms 10 Hz
subplot(3,4,1)
x = [1 2];
y = [mean(spikeRateMean_2ms10Hz_stim) mean(spikeRateMean_2ms10Hz_nostim)];
value = bar(x,y,'FaceColor',[.7 .7 .7],'EdgeColor','none');

box off
xticks([1 2])
xticklabels({'Stim','No Stim'})
ylabel({'10 Hz';'Firing rate (Hz)'})

hold on
xdata = ones(length(spikeRateMean_2ms10Hz_stim),2);
xdata = x.*xdata;
data1 = spikeRateMean_2ms10Hz_stim;
data2 = spikeRateMean_2ms10Hz_nostim;

for i=1:length(data1)  % iterate over number of bar objects
      plot(xdata(i,1),data1(i),'o','MarkerSize',5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k')
      hold on
      plot(xdata(i,2),data2(i),'o','MarkerSize',5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k')
      hold on
      plot([data1(i) data2(i)],'Color','k')
end
box off
xlim([0.5 2.5])
title('2 ms')

% for 2 ms 20 Hz
subplot(3,4,5)
x = [1 2];
y = [mean(spikeRateMean_2ms20Hz_stim) mean(spikeRateMean_2ms20Hz_nostim)];
value = bar(x,y,'FaceColor',[.7 .7 .7],'EdgeColor','none');

box off
xticks([1 2])
xticklabels({'Stim','No Stim'})
ylabel({'20 Hz';'Firing rate (Hz)'})

hold on
xdata = ones(length(spikeRateMean_2ms20Hz_stim),2);
xdata = x.*xdata;
data1 = spikeRateMean_2ms20Hz_stim;
data2 = spikeRateMean_2ms20Hz_nostim;

for i=1:length(data1)  % iterate over number of bar objects
      plot(xdata(i,1),data1(i),'o','MarkerSize',5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k')
      hold on
      plot(xdata(i,2),data2(i),'o','MarkerSize',5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k')
      hold on
      plot([data1(i) data2(i)],'Color','k')
end
box off
xlim([0.5 2.5])

% for 2 ms 40 Hz
subplot(3,4,9)
x = [1 2];
y = [mean(spikeRateMean_2ms40Hz_stim) mean(spikeRateMean_2ms40Hz_nostim)];
value = bar(x,y,'FaceColor',[.7 .7 .7],'EdgeColor','none');

box off
xticks([1 2])
xticklabels({'Stim','No Stim'})
ylabel({'40 Hz';'Firing rate (Hz)'})

hold on
xdata = ones(length(spikeRateMean_2ms40Hz_stim),2);
xdata = x.*xdata;
data1 = spikeRateMean_2ms40Hz_stim;
data2 = spikeRateMean_2ms40Hz_nostim;

for i=1:length(data1)  % iterate over number of bar objects
      plot(xdata(i,1),data1(i),'o','MarkerSize',5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k')
      hold on
      plot(xdata(i,2),data2(i),'o','MarkerSize',5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k')
      hold on
      plot([data1(i) data2(i)],'Color','k')
end
box off
xlim([0.5 2.5])

% for 5 ms 10 Hz
subplot(3,4,2)
x = [1 2];
y = [mean(spikeRateMean_5ms10Hz_stim) mean(spikeRateMean_5ms10Hz_nostim)];
value = bar(x,y,'FaceColor',[.7 .7 .7],'EdgeColor','none');

box off
xticks([1 2])
xticklabels({'Stim','No Stim'})
ylabel('Firing rate (Hz)')

hold on
xdata = ones(length(spikeRateMean_5ms10Hz_stim),2);
xdata = x.*xdata;
data1 = spikeRateMean_5ms10Hz_stim;
data2 = spikeRateMean_5ms10Hz_nostim;

for i=1:length(data1)  % iterate over number of bar objects
      plot(xdata(i,1),data1(i),'o','MarkerSize',5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k')
      hold on
      plot(xdata(i,2),data2(i),'o','MarkerSize',5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k')
      hold on
      plot([data1(i) data2(i)],'Color','k')
end
box off
xlim([0.5 2.5])
title('5 ms')

% for 5 ms 20 Hz
subplot(3,4,6)
x = [1 2];
y = [mean(spikeRateMean_5ms20Hz_stim) mean(spikeRateMean_5ms20Hz_nostim)];
value = bar(x,y,'FaceColor',[.7 .7 .7],'EdgeColor','none');

box off
xticks([1 2])
xticklabels({'Stim','No Stim'})
ylabel('Firing rate (Hz)')

hold on
xdata = ones(length(spikeRateMean_5ms20Hz_stim),2);
xdata = x.*xdata;
data1 = spikeRateMean_5ms20Hz_stim;
data2 = spikeRateMean_5ms20Hz_nostim;

for i=1:length(data1)  % iterate over number of bar objects
      plot(xdata(i,1),data1(i),'o','MarkerSize',5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k')
      hold on
      plot(xdata(i,2),data2(i),'o','MarkerSize',5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k')
      hold on
      plot([data1(i) data2(i)],'Color','k')
end
box off
xlim([0.5 2.5])

% for 5 ms 40 Hz
subplot(3,4,10)
x = [1 2];
y = [mean(spikeRateMean_5ms40Hz_stim) mean(spikeRateMean_5ms40Hz_nostim)];
value = bar(x,y,'FaceColor',[.7 .7 .7],'EdgeColor','none');

box off
xticks([1 2])
xticklabels({'Stim','No Stim'})
ylabel('Firing rate (Hz)')

hold on
xdata = ones(length(spikeRateMean_5ms40Hz_stim),2);
xdata = x.*xdata;
data1 = spikeRateMean_5ms40Hz_stim;
data2 = spikeRateMean_5ms40Hz_nostim;

for i=1:length(data1)  % iterate over number of bar objects
      plot(xdata(i,1),data1(i),'o','MarkerSize',5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k')
      hold on
      plot(xdata(i,2),data2(i),'o','MarkerSize',5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k')
      hold on
      plot([data1(i) data2(i)],'Color','k')
end
box off
xlim([0.5 2.5])


% for 10 ms 10 Hz
subplot(3,4,3)
x = [1 2];
y = [mean(spikeRateMean_10ms10Hz_stim) mean(spikeRateMean_10ms10Hz_nostim)];
value = bar(x,y,'FaceColor',[.7 .7 .7],'EdgeColor','none');

box off
xticks([1 2])
xticklabels({'Stim','No Stim'})
ylabel('Firing rate (Hz)')

hold on
xdata = ones(length(spikeRateMean_10ms10Hz_stim),2);
xdata = x.*xdata;
data1 = spikeRateMean_10ms10Hz_stim;
data2 = spikeRateMean_10ms10Hz_nostim;

for i=1:length(data1)  % iterate over number of bar objects
      plot(xdata(i,1),data1(i),'o','MarkerSize',5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k')
      hold on
      plot(xdata(i,2),data2(i),'o','MarkerSize',5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k')
      hold on
      plot([data1(i) data2(i)],'Color','k')
end
box off
xlim([0.5 2.5])
title('10 ms')

% for 10 ms 20 Hz
subplot(3,4,7)
x = [1 2];
y = [mean(spikeRateMean_10ms20Hz_stim) mean(spikeRateMean_10ms20Hz_nostim)];
value = bar(x,y,'FaceColor',[.7 .7 .7],'EdgeColor','none');

box off
xticks([1 2])
xticklabels({'Stim','No Stim'})
ylabel('Firing rate (Hz)')

hold on
xdata = ones(length(spikeRateMean_10ms20Hz_stim),2);
xdata = x.*xdata;
data1 = spikeRateMean_10ms20Hz_stim;
data2 = spikeRateMean_10ms20Hz_nostim;

for i=1:length(data1)  % iterate over number of bar objects
      plot(xdata(i,1),data1(i),'o','MarkerSize',5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k')
      hold on
      plot(xdata(i,2),data2(i),'o','MarkerSize',5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k')
      hold on
      plot([data1(i) data2(i)],'Color','k')
end
box off
xlim([0.5 2.5])


% for 10 ms 40 Hz
subplot(3,4,11)
x = [1 2];
y = [mean(spikeRateMean_10ms40Hz_stim) mean(spikeRateMean_10ms40Hz_nostim)];
value = bar(x,y,'FaceColor',[.7 .7 .7],'EdgeColor','none');

box off
xticks([1 2])
xticklabels({'Stim','No Stim'})
ylabel('Firing rate (Hz)')

hold on
xdata = ones(length(spikeRateMean_10ms40Hz_stim),2);
xdata = x.*xdata;
data1 = spikeRateMean_10ms40Hz_stim;
data2 = spikeRateMean_10ms40Hz_nostim;

for i=1:length(data1)  % iterate over number of bar objects
      plot(xdata(i,1),data1(i),'o','MarkerSize',5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k')
      hold on
      plot(xdata(i,2),data2(i),'o','MarkerSize',5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k')
      hold on
      plot([data1(i) data2(i)],'Color','k')
end
box off
xlim([0.5 2.5])

% for 20 ms 10 Hz
subplot(3,4,4)
x = [1 2];
y = [mean(spikeRateMean_20ms10Hz_stim) mean(spikeRateMean_20ms10Hz_nostim)];
value = bar(x,y,'FaceColor',[.7 .7 .7],'EdgeColor','none');

box off
xticks([1 2])
xticklabels({'Stim','No Stim'})
ylabel('Firing rate (Hz)')

hold on
xdata = ones(length(spikeRateMean_20ms10Hz_stim),2);
xdata = x.*xdata;
data1 = spikeRateMean_20ms10Hz_stim;
data2 = spikeRateMean_20ms10Hz_nostim;

for i=1:length(data1)  % iterate over number of bar objects
      plot(xdata(i,1),data1(i),'o','MarkerSize',5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k')
      hold on
      plot(xdata(i,2),data2(i),'o','MarkerSize',5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k')
      hold on
      plot([data1(i) data2(i)],'Color','k')
end
box off
xlim([0.5 2.5])
title('20 ms')

% for 20 ms 20 Hz
subplot(3,4,8)
x = [1 2];
y = [mean(spikeRateMean_20ms20Hz_stim) mean(spikeRateMean_20ms20Hz_nostim)];
value = bar(x,y,'FaceColor',[.7 .7 .7],'EdgeColor','none');

box off
xticks([1 2])
xticklabels({'Stim','No Stim'})
ylabel('Firing rate (Hz)')

hold on
xdata = ones(length(spikeRateMean_20ms20Hz_stim),2);
xdata = x.*xdata;
data1 = spikeRateMean_20ms20Hz_stim;
data2 = spikeRateMean_20ms20Hz_nostim;

for i=1:length(data1)  % iterate over number of bar objects
      plot(xdata(i,1),data1(i),'o','MarkerSize',5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k')
      hold on
      plot(xdata(i,2),data2(i),'o','MarkerSize',5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k')
      hold on
      plot([data1(i) data2(i)],'Color','k')
end
box off
xlim([0.5 2.5])

% for 20 ms 40 Hz
subplot(3,4,12)
x = [1 2];
y = [mean(spikeRateMean_20ms40Hz_stim) mean(spikeRateMean_20ms40Hz_nostim)];
value = bar(x,y,'FaceColor',[.7 .7 .7],'EdgeColor','none');

box off
xticks([1 2])
xticklabels({'Stim','No Stim'})
ylabel('Firing rate (Hz)')

hold on
xdata = ones(length(spikeRateMean_20ms40Hz_stim),2);
xdata = x.*xdata;
data1 = spikeRateMean_20ms40Hz_stim;
data2 = spikeRateMean_20ms40Hz_nostim;

for i=1:length(data1)  % iterate over number of bar objects
      plot(xdata(i,1),data1(i),'o','MarkerSize',5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k')
      hold on
      plot(xdata(i,2),data2(i),'o','MarkerSize',5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k')
      hold on
      plot([data1(i) data2(i)],'Color','k')
end
box off
xlim([0.5 2.5])


sgtitle('Firing Rate of MUA within 1s of stimulus vs no stimulus')
fig_name = strcat('Summary_firingRate_removewavers');
fig_name = strrep(fig_name,'.','');
saveas(fig,fig_name)

%% Pool all the stimuli together and plot
Stim = [spikeRateMean_2ms10Hz_stim spikeRateMean_2ms20Hz_stim spikeRateMean_2ms40Hz_stim ...
    spikeRateMean_5ms10Hz_stim spikeRateMean_5ms20Hz_stim spikeRateMean_5ms40Hz_stim ...
    spikeRateMean_10ms10Hz_stim spikeRateMean_10ms20Hz_stim spikeRateMean_10ms40Hz_stim ...
    spikeRateMean_20ms10Hz_stim spikeRateMean_20ms20Hz_stim spikeRateMean_20ms40Hz_stim];

NoStim = [spikeRateMean_2ms10Hz_nostim spikeRateMean_2ms20Hz_nostim spikeRateMean_2ms40Hz_nostim ...
    spikeRateMean_5ms10Hz_nostim spikeRateMean_5ms20Hz_nostim spikeRateMean_5ms40Hz_nostim ...
    spikeRateMean_10ms10Hz_nostim spikeRateMean_10ms20Hz_nostim spikeRateMean_10ms40Hz_nostim ...
    spikeRateMean_20ms10Hz_nostim spikeRateMean_20ms20Hz_nostim spikeRateMean_20ms40Hz_nostim];


fig = figure;
x = [1 2];
y = [mean(Stim) mean(NoStim)];
value = bar(x,y,'FaceColor',[.7 .7 .7],'EdgeColor','none');

hold on
errorbar(x,y,[std(Stim)/sqrt(length(Stim)) std(NoStim)/sqrt(length(NoStim))],'Color','k','CapSize',0)
box off
xticks([1 2])
xticklabels({'Stim','No Stim'})
ylabel('Firing rate (Hz)')

% hold on
% xdata = ones(length(Stim),2);
% xdata = x.*xdata;
% data1 = Stim;
% data2 = NoStim;
% 
% for i=1:length(data1)  % iterate over number of bar objects
%       plot(xdata(i,1),data1(i),'o','MarkerSize',5,...
%         'MarkerEdgeColor','k',...
%         'MarkerFaceColor','k')
%       hold on
%       plot(xdata(i,2),data2(i),'o','MarkerSize',5,...
%         'MarkerEdgeColor','k',...
%         'MarkerFaceColor','k')
%       hold on
%       plot([data1(i) data2(i)],'Color','k')
% end

box off
xlim([0.5 2.5])


sgtitle('Firing Rate of MUA within 1s of stimulus vs no stimulus')
fig_name = strcat('Summary_firingRate_all_stimulus_removeWave');
fig_name = strrep(fig_name,'.','');
saveas(fig,fig_name)
[h,p,ci,stats] = ttest(Stim, NoStim);
p

%% Scatter plot
fig = figure;
%plot(NoStim, Stim,'.b','MarkerSize',5)
plot(NoStim, Stim,'.','Color',[.5 .5 .5],'MarkerSize',5)

hold on
x = linspace(0,400);
y = linspace(0,400);
plot(x,y,'--k');
hold on

axis tight
xlim([0 400])
ylim([0 400])
box off
ylabel('MUA Firing Rate-Light On (Hz)')
xlabel('MUA Firing Rate-Light Off (Hz)')
set(gcf,'position',[10,10,650,550])

sgtitle('Firing Rate of MUA stimulus vs no stimulus')
fig_name = strcat('scatter_firingRate_all_stimulus');
fig_name = strrep(fig_name,'.','');
saveas(fig,fig_name)

%% Save data
prompt = {'Reulst name to save'};
dlgTitle = 'User-defined input';
nLines = 1;
def = {'20211027_tdT'};
answer = inputdlg(prompt, dlgTitle, nLines, def);

file_name = answer{1};
file_name = [file_name,'.','mat'];
save(file_name,'firing','p','-v7.3')
% clear;
% clc;

