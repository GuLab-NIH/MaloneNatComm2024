% load ephys data
load('2021110_AAV8_hsyn_ChR2_x6.mat')
stim_firing = [];
no_stim_firing = [];

for i = 1:7
    if i==1
        stim_firing = firing.spikeRateMean_2ms10Hz_stim{i}';
        no_stim_firing = firing.spikeRateMean_2ms10Hz_nostim{i}';
    else
        stim_firing = [stim_firing firing.spikeRateMean_2ms10Hz_stim{i}'];
        no_stim_firing = [no_stim_firing firing.spikeRateMean_2ms10Hz_nostim{i}'];
    end
end

mean_stim_firing = mean(stim_firing);
sem_stim_firing = std(stim_firing)/sqrt(length(stim_firing));
mean_no_firing = mean(no_stim_firing);
sem_no_firing = std(no_stim_firing)/sqrt(length(no_stim_firing));

%%
fig = figure;

data = [mean_no_firing mean_stim_firing];
error = [sem_no_firing sem_stim_firing];

x = 1:2;
bar(x,data);

hold on
er= errorbar(x,data,error);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

[h,p]=ttest(stim_firing,no_stim_firing)
title(['p=',num2str(p)])
box off
xticks([1 2])
xticklabels({'No Stim','Stim'})
ylabel('Firing rate (Hz)')
yticks([0 40 80])
yticklabels({'0','40','80'})