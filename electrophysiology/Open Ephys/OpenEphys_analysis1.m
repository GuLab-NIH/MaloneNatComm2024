%% Open ephys data analysis
% First run "load_open_ephys_binary.m" to extract data for each recording,
% then modify and run "Extract_and_save_20211021.m" to organize the data format 
% then load the data and stimulus timestmaps
clear;clc;
[fileName pathName]=uigetfile('*.mat','Select .mat file',...
    'Z:\SNM\labMembers\NT\Documentation\data\Optogenetics\ephys\DataFiles\20211110_AAV8_hsyn_ChR2_dilute6_UsedForPaper');

file_location = [pathName fileName];
load(file_location);
% Change the name to "data_file" and "stim_file"

%% Select channels with signals
% View all the channels in the Open Ephys GUI then decide the channels
% with signals
prompt = {'Which channels of interest?','Reference channels (0 for no reference)'};
dlgTitle = 'channels with signals';
nLines = 1;
def = {'6, 13, 17, 22, 27','0'};
answer = inputdlg(prompt, dlgTitle, nLines, def);

ephys.channel = str2num(answer{1});
ephys.ref = str2num(answer{2});


voltage = zeros(length(ephys.channel),length(data_file.Data));
ref_volt = zeros(length(ephys.ref),length(data_file.Data));

for channel_interest = 1:length(ephys.channel)
    voltage(channel_interest,:) = data_file.Data(channel_interest,:);
end

if ephys.ref>0
    for channel_ref = 1:length(ephys.ref)
        ref_volt(channel_ref,:) = data_file.Data(channel_ref,:);
    end
else
end

ephys.stim = stim_file.Timestamps;
ephys.time = data_file.Timestamps; 

% find the median from all the channels and subtract
ephys.median_noise = median(data_file.Data);
voltage = voltage - ephys.median_noise;
if ephys.ref>0
    ref_volt = ref_volt - ephys.median_noise;
end


%% Segment the pulse trains
PulseTrain_interval = [0 diff(ephys.stim)'];
PulseSeg = [1 find(PulseTrain_interval>5*30000)]; % the sampling rate is 30k/
% sec, so what it does is to find the index when the train interval is 
% larger than 5 sec, we set the train interval to be 10 sec long

% Get the pulse duration and frequency
number = 1;
first_stimON = PulseSeg; % the 1st LED ON of each pulse train
ephys.stimuli_info = zeros(length(PulseSeg),2); % the pulse duration and freque

for stim_num = 1:length(PulseSeg) % should be 36 stimuli
    if stim_num < length(PulseSeg)
        PulseTrain_time = [PulseSeg(stim_num):PulseSeg(stim_num+1)-1];% find the timestamps within one pulse train set
        Pulse_interval = diff(ephys.stim(PulseTrain_time));
        PulseDur = round(mean(Pulse_interval(1:2:end)))/30000; % the pulse duration in sec
        PulseFeq = length(PulseTrain_time)/2; % the pulse frequency in Hz
        
        ephys.stimON{stim_num} = ephys.stim(PulseTrain_time(1:2:end)); % the timestamps for the stim ON within one pulse train
        ephys.stimOFF{stim_num} = ephys.stim(PulseTrain_time(2:2:end)); % the timestamps for the stim OFF within one pulse train
        ephys.stimuli_info(number,:) = [PulseDur PulseFeq]; 
        number = number+1;
    else
        PulseTrain_time = [PulseSeg(stim_num):length(PulseTrain_interval)];% find the timestamps within the last pulse train set
        Pulse_interval = diff(ephys.stim(PulseTrain_time));
        PulseDur = round(mean(Pulse_interval(1:2:end)))/30000; % the pulse duration in sec
        PulseFeq = length(PulseTrain_time)/2; % the pulse frequency in Hz
        
        ephys.stimON{stim_num} = ephys.stim(PulseTrain_time(1:2:end)); % the timestamps for the stim ON within one pulse train
        ephys.stimOFF{stim_num} = ephys.stim(PulseTrain_time(2:2:end)); % the timestamps for the stim OFF within one pulse train
        ephys.stimuli_info(number,:) = [PulseDur PulseFeq]; 
    end
       
end

%% Segment the filtered signals from given channels based on stimuli

% filter the reference signals for removing electric noises
for channel_ref = 1:length(ephys.ref)
    ephys.ref_volt(channel_ref,:) = ref_volt(channel_ref,:);
    ephys.refSpike_filtered(channel_ref,:) = bandpass(ref_volt(channel_ref,:),[250 8000],30000); % 250-8000 Hz for spikes
end

% Filter the signals based on frequency
for channel_interest = 1:length(ephys.channel)
    ephys.raw(channel_interest,:) = voltage(channel_interest,:);
    ephys.Spike_filtered(channel_interest,:) = bandpass(voltage(channel_interest,:),[250 8000],30000); % 250-8000 Hz for spikes
end
% For LFP signals
Hd = LFP_bandpass_filter; % call the filter function
for channel_interest = 1:length(ephys.channel)
    ephys.LFP_filtered(channel_interest,:) = filter(Hd,voltage(channel_interest,:)); % 0.7-170 Hz for LFP
end
% Segment the recording from 0.5 sec before the stim ON, 8 sec after the
% stim OFF
ephys.start = zeros(1,length(PulseSeg)); % the timestamps 0.5 sec before each pulse train
ephys.stop = zeros(1,length(PulseSeg)); % the timestamps 8 sec after each pulse train

for stim_num = 1:length(PulseSeg) % should be 36 stimuli
    for channel_interest = 1:length(ephys.channel)
        ephys.start(stim_num) = ephys.stim(first_stimON(stim_num)) - 0.5*30000 - ephys.time(1); % 30000 per sec for sampling rate
        ephys.stop(stim_num) = ephys.start(stim_num)+ 1*30000 + 8*30000;% the start time plus 1 sec of
        % stimulus and 8 sec after stimulus off
        
        % segment the filtered data
        ephys.Spike_filtered_seg(channel_interest,stim_num,:) = ephys.Spike_filtered(channel_interest,ephys.start(stim_num): ephys.stop(stim_num))*0.195; %*0.195 convert to uV
        ephys.LFP_filtered_seg(channel_interest,stim_num,:) = ephys.LFP_filtered(channel_interest,ephys.start(stim_num): ephys.stop(stim_num))*0.195 ; %* convert to uV
        ephys.raw_seg(channel_interest,stim_num,:) = ephys.raw(channel_interest,ephys.start(stim_num): ephys.stop(stim_num))*0.195; %*0.195 convert to uV
    end
    
    for channel_ref = 1:length(ephys.ref)              
        % segment the filtered data for reference channels
        ephys.ref_Spike_filtered_seg(channel_ref,stim_num,:) = ephys.refSpike_filtered(channel_ref,ephys.start(stim_num): ephys.stop(stim_num))*0.195; %*0.195 convert to uV
        ephys.ref_raw_seg(channel_ref,stim_num,:) = ephys.ref_volt(channel_ref,ephys.start(stim_num): ephys.stop(stim_num))*0.195; %*0.195 convert to uV
    end
    
end

%% Set the threshold
prompt = {'Standard deviation threshold? (3)','Time bin (5 ms)','Threshold for electric noise (std = 9)'};
dlgTitle = 'User-defined input';
nLines = 1;
def = {'3','0.005','9'};
answer = inputdlg(prompt, dlgTitle, nLines, def);

threshold = str2num(answer{1});
bin = str2num(answer{2});
threshold_noise = str2num(answer{3});


% % calculate the threshold
% for channel_interest = 1:length(ephys.channel)
%     ephys.threshold(channel_interest) = threshold.*std(voltage(channel_interest,:)*0.195,0,'all');
% end

for stim_num = 1:length(PulseSeg) % should be 36 stimuli
    for channel_interest = 1:length(ephys.channel)
        % Calculate the std for each channel
        ephys.threshold(channel_interest,stim_num) = threshold.*std(ephys.Spike_filtered_seg(channel_interest,stim_num,:),0,3);
%         
        % z-score
        %z_trace(channel_interest,stim_num,:) = squeeze(zscore(ephys.Spike_filtered_seg(channel_interest,stim_num,:)))';
        
        % find spike that passes the threshold
        ephys.spike_binary(channel_interest,stim_num,:) = zeros(1,1,length(ephys.Spike_filtered_seg(channel_interest,stim_num,:)));
        trace(1,:) = -ephys.Spike_filtered_seg(channel_interest,stim_num,:);
        [pks,locs] = findpeaks(trace); % find local peaks of the reversed waveforms
        threshold_spike = find(trace(locs)>ephys.threshold(channel_interest,stim_num));        
        
        ephys.spike_binary(channel_interest,stim_num,locs(threshold_spike)) = 1;
        
        % Remove electric noises
        noise_spike = find(trace(locs)>threshold_noise.*std(ephys.Spike_filtered_seg(channel_interest,stim_num,:),0,3)); 
        ephys.spike_binary(channel_interest,stim_num,locs(noise_spike)) = 0;
    end
    
    % find spikes in reference channels that passes the threshold =>
    % not real spike but electrical noises
    for channel_ref = 1:length(ephys.ref)
        ephys.ref_volt_spike_binary(channel_ref,stim_num,:) = zeros(1,1,length(ephys.ref_Spike_filtered_seg(channel_ref,stim_num,:)));
        ref_trace(1,:) = -ephys.ref_Spike_filtered_seg(channel_ref,stim_num,:);
        [pks,locs] = findpeaks(ref_trace); % find local peaks of the reversed waveforms
        threshold_spike = find(ref_trace(locs)>ephys.threshold(channel_ref,stim_num));        
        
        ephys.ref_volt_spike_binary(channel_interest,stim_num,locs(threshold_spike)) = 1;   
    end
    
    % Remove the electric noises by comparing the ref channels
    for channel_interest = 1:length(ephys.channel)
        noise = find(ephys.spike_binary(channel_interest,stim_num,:).*ephys.ref_volt_spike_binary(1,stim_num,:));
        ephys.spike_binary(channel_interest,stim_num,noise) = 0;
    end
end 


%% Calculate the spike numbers within the 1 sec stimulus and without stimulus
prompt = {'Start time window for spontaneous spike (sec, relative to stimulus off)',...
    'Stot time window for spontaneous spike (sec, relative to stimulus off)'};
dlgTitle = 'User-defined input';
nLines = 1;
def = {'1','2'};
answer = inputdlg(prompt, dlgTitle, nLines, def);

no_stimulus_window = str2num(answer{2}) - str2num(answer{1});
for stim_num = 1:length(PulseSeg) % should be 36 stimuli
    ephys.no_stimulus_start(stim_num) = str2num(answer{1})*30000+0.5*30000+30000;
    % the default is 1 seconds after the stimulus if off, 30000 samples = 1 sec

    ephys.no_stimulus_stop(stim_num) = str2num(answer{2})*30000+0.5*30000+30000;
    % the default is 3 seconds after the stimulus if off, 30000 samples = 1 sec
    
    for channel_interest = 1:length(ephys.channel)
        ephys.spike_stim_epoch(channel_interest,stim_num) = sum(ephys.spike_binary(channel_interest,stim_num,0.5*30000:0.5*30000+30000));
%         x1 = double([ephys.stimON{stim_num}]-[ephys.stimON{stim_num}(1)])+3000;
%         x2 = double([ephys.stimOFF{stim_num}]-[ephys.stimON{stim_num}(1)])+3000+0.002*30000; 
%         
%         loop = length(x1);
%         for pulse = 1:loop
%             pulse_spike(1,loop) = sum(ephys.spike_binary(channel_interest,stim_num,x1(loop):x2(loop)));
%         end%         
        
        ephys.spike_nostim_epoch(channel_interest,stim_num) = sum(ephys.spike_binary(channel_interest,stim_num,ephys.no_stimulus_start(stim_num):ephys.no_stimulus_stop(stim_num)));
        ephys.spikeRate_nostim_epoch(channel_interest,stim_num) = ephys.spike_nostim_epoch(channel_interest,stim_num)/no_stimulus_window;
        total_spike = ephys.spike_stim_epoch(channel_interest,stim_num) + ephys.spike_nostim_epoch(channel_interest,stim_num);
        if total_spike <=5
            ephys.spike_stim_epoch(channel_interest,stim_num) = nan;
            ephys.spike_nostim_epoch(channel_interest,stim_num) = nan;
            ephys.spikeRate_nostim_epoch(channel_interest,stim_num) = nan;
        end
    end
end

ephys.spikeRate_2ms10Hz_stim = nan(length(ephys.channel),3);
ephys.spikeRate_2ms20Hz_stim = nan(length(ephys.channel),3);
ephys.spikeRate_2ms40Hz_stim = nan(length(ephys.channel),3);

ephys.spikeRate_5ms10Hz_stim = nan(length(ephys.channel),3);
ephys.spikeRate_5ms20Hz_stim = nan(length(ephys.channel),3);
ephys.spikeRate_5ms40Hz_stim = nan(length(ephys.channel),3);

ephys.spikeRate_10ms10Hz_stim = nan(length(ephys.channel),3);
ephys.spikeRate_10ms20Hz_stim = nan(length(ephys.channel),3);
ephys.spikeRate_10ms40Hz_stim = nan(length(ephys.channel),3);

ephys.spikeRate_20ms10Hz_stim = nan(length(ephys.channel),3);
ephys.spikeRate_20ms20Hz_stim = nan(length(ephys.channel),3);
ephys.spikeRate_20ms40Hz_stim = nan(length(ephys.channel),3);

ephys.spikeRate_2ms10Hz_nostim = nan(length(ephys.channel),3);
ephys.spikeRate_2ms20Hz_nostim = nan(length(ephys.channel),3);
ephys.spikeRate_2ms40Hz_nostim = nan(length(ephys.channel),3);

ephys.spikeRate_5ms10Hz_nostim = nan(length(ephys.channel),3);
ephys.spikeRate_5ms20Hz_nostim = nan(length(ephys.channel),3);
ephys.spikeRate_5ms40Hz_nostim = nan(length(ephys.channel),3);

ephys.spikeRate_10ms10Hz_nostim = nan(length(ephys.channel),3);
ephys.spikeRate_10ms20Hz_nostim = nan(length(ephys.channel),3);
ephys.spikeRate_10ms40Hz_nostim = nan(length(ephys.channel),3);

ephys.spikeRate_20ms10Hz_nostim = nan(length(ephys.channel),3);
ephys.spikeRate_20ms20Hz_nostim = nan(length(ephys.channel),3);
ephys.spikeRate_20ms40Hz_nostim = nan(length(ephys.channel),3);


% organize spike rate based on different pulse duration and frequency
% 2 ms 10 Hz
pulse_dur = ismember(ephys.stimuli_info(:,1),0.002)';
pulse_freq = ismember(ephys.stimuli_info(:,2),10)';
dur2ms_10Hz = pulse_dur.*pulse_freq; % for 2 ms and 10 Hz

for channel_interest = 1:length(ephys.channel)
    spike_locs = find(dur2ms_10Hz);
    ephys.spikeRate_2ms10Hz_stim(channel_interest,:) = ephys.spike_stim_epoch(channel_interest,spike_locs);
    ephys.spikeRate_2ms10Hz_nostim(channel_interest,:) = ephys.spikeRate_nostim_epoch(channel_interest,spike_locs);
    ephys.spikeRateMean_2ms10Hz_stim(channel_interest) = nanmean(ephys.spike_stim_epoch(channel_interest,spike_locs));
    ephys.spikeRateMean_2ms10Hz_nostim(channel_interest) = nanmean(ephys.spike_nostim_epoch(channel_interest,spike_locs));
end

% 2 ms 20 Hz
pulse_dur = ismember(ephys.stimuli_info(:,1),0.002)';
pulse_freq = ismember(ephys.stimuli_info(:,2),20)';
dur2ms_20Hz = pulse_dur.*pulse_freq; % for 2 ms and 20 Hz

for channel_interest = 1:length(ephys.channel)
    spike_locs = find(dur2ms_20Hz);
    ephys.spikeRate_2ms20Hz_stim(channel_interest,:) = ephys.spike_stim_epoch(channel_interest,spike_locs);
    ephys.spikeRate_2ms20Hz_nostim(channel_interest,:) = ephys.spikeRate_nostim_epoch(channel_interest,spike_locs);
    ephys.spikeRateMean_2ms20Hz_stim(channel_interest) = nanmean(ephys.spike_stim_epoch(channel_interest,spike_locs));
    ephys.spikeRateMean_2ms20Hz_nostim(channel_interest) = nanmean(ephys.spike_nostim_epoch(channel_interest,spike_locs));
end

% 2 ms 40 Hz
pulse_dur = ismember(ephys.stimuli_info(:,1),0.002)';
pulse_freq = ismember(ephys.stimuli_info(:,2),40)';
dur2ms_40Hz = pulse_dur.*pulse_freq; % for 2 ms and 40 Hz

for channel_interest = 1:length(ephys.channel)
    spike_locs = find(dur2ms_40Hz);
    ephys.spikeRate_2ms40Hz_stim(channel_interest,:) = ephys.spike_stim_epoch(channel_interest,spike_locs);
    ephys.spikeRate_2ms40Hz_nostim(channel_interest,:) = ephys.spikeRate_nostim_epoch(channel_interest,spike_locs);
    ephys.spikeRateMean_2ms40Hz_stim(channel_interest) = nanmean(ephys.spike_stim_epoch(channel_interest,spike_locs));
    ephys.spikeRateMean_2ms40Hz_nostim(channel_interest) = nanmean(ephys.spike_nostim_epoch(channel_interest,spike_locs));
end

% 5 ms 10 Hz
pulse_dur = ismember(ephys.stimuli_info(:,1),0.005)';
pulse_freq = ismember(ephys.stimuli_info(:,2),10)';
dur5ms_10Hz = pulse_dur.*pulse_freq; % for 5 ms and 10 Hz

for channel_interest = 1:length(ephys.channel)
    spike_locs = find(dur5ms_10Hz);
    ephys.spikeRate_5ms10Hz_stim(channel_interest,:) = ephys.spike_stim_epoch(channel_interest,spike_locs);
    ephys.spikeRate_5ms10Hz_nostim(channel_interest,:) = ephys.spikeRate_nostim_epoch(channel_interest,spike_locs);
    ephys.spikeRateMean_5ms10Hz_stim(channel_interest) = nanmean(ephys.spike_stim_epoch(channel_interest,spike_locs));
    ephys.spikeRateMean_5ms10Hz_nostim(channel_interest) = nanmean(ephys.spike_nostim_epoch(channel_interest,spike_locs));
end

% 5 ms 20 Hz
pulse_dur = ismember(ephys.stimuli_info(:,1),0.005)';
pulse_freq = ismember(ephys.stimuli_info(:,2),20)';
dur5ms_20Hz = pulse_dur.*pulse_freq; % for 5 ms and 20 Hz

for channel_interest = 1:length(ephys.channel)
    spike_locs = find(dur5ms_20Hz);
    ephys.spikeRate_5ms20Hz_stim(channel_interest,:) = ephys.spike_stim_epoch(channel_interest,spike_locs);
    ephys.spikeRate_5ms20Hz_nostim(channel_interest,:) = ephys.spikeRate_nostim_epoch(channel_interest,spike_locs);
    ephys.spikeRateMean_5ms20Hz_stim(channel_interest) = nanmean(ephys.spike_stim_epoch(channel_interest,spike_locs));
    ephys.spikeRateMean_5ms20Hz_nostim(channel_interest) = nanmean(ephys.spike_nostim_epoch(channel_interest,spike_locs));
end

% 5 ms 40 Hz
pulse_dur = ismember(ephys.stimuli_info(:,1),0.005)';
pulse_freq = ismember(ephys.stimuli_info(:,2),40)';
dur5ms_40Hz = pulse_dur.*pulse_freq; % for 5 ms and 40 Hz

for channel_interest = 1:length(ephys.channel)
    spike_locs = find(dur5ms_40Hz);
    ephys.spikeRate_5ms40Hz_stim(channel_interest,:) = ephys.spike_stim_epoch(channel_interest,spike_locs);
    ephys.spikeRate_5ms40Hz_nostim(channel_interest,:) = ephys.spikeRate_nostim_epoch(channel_interest,spike_locs);
    ephys.spikeRateMean_5ms40Hz_stim(channel_interest) = nanmean(ephys.spike_stim_epoch(channel_interest,spike_locs));
    ephys.spikeRateMean_5ms40Hz_nostim(channel_interest) = nanmean(ephys.spike_nostim_epoch(channel_interest,spike_locs));
end

% 10 ms 10 Hz
pulse_dur = ismember(ephys.stimuli_info(:,1),0.01)';
pulse_freq = ismember(ephys.stimuli_info(:,2),10)';
dur10ms_10Hz = pulse_dur.*pulse_freq; % for 10 ms and 10 Hz

for channel_interest = 1:length(ephys.channel)
    spike_locs = find(dur10ms_10Hz);
    ephys.spikeRate_10ms10Hz_stim(channel_interest,:) = ephys.spike_stim_epoch(channel_interest,spike_locs);
    ephys.spikeRate_10ms10Hz_nostim(channel_interest,:) = ephys.spikeRate_nostim_epoch(channel_interest,spike_locs);
    ephys.spikeRateMean_10ms10Hz_stim(channel_interest) = nanmean(ephys.spike_stim_epoch(channel_interest,spike_locs));
    ephys.spikeRateMean_10ms10Hz_nostim(channel_interest) = nanmean(ephys.spike_nostim_epoch(channel_interest,spike_locs));
end

% 10 ms 20 Hz
pulse_dur = ismember(ephys.stimuli_info(:,1),0.01)';
pulse_freq = ismember(ephys.stimuli_info(:,2),20)';
dur10ms_20Hz = pulse_dur.*pulse_freq; % for 10 ms and 20 Hz

for channel_interest = 1:length(ephys.channel)
    spike_locs = find(dur10ms_20Hz);
    ephys.spikeRate_10ms20Hz_stim(channel_interest,:) = ephys.spike_stim_epoch(channel_interest,spike_locs);
    ephys.spikeRate_10ms20Hz_nostim(channel_interest,:) = ephys.spikeRate_nostim_epoch(channel_interest,spike_locs);
    ephys.spikeRateMean_10ms20Hz_stim(channel_interest) = nanmean(ephys.spike_stim_epoch(channel_interest,spike_locs));
    ephys.spikeRateMean_10ms20Hz_nostim(channel_interest) = nanmean(ephys.spike_nostim_epoch(channel_interest,spike_locs));
end

% 10 ms 40 Hz
pulse_dur = ismember(ephys.stimuli_info(:,1),0.01)';
pulse_freq = ismember(ephys.stimuli_info(:,2),40)';
dur10ms_40Hz = pulse_dur.*pulse_freq; % for 10 ms and 40 Hz

for channel_interest = 1:length(ephys.channel)
    spike_locs = find(dur10ms_40Hz);
    ephys.spikeRate_10ms40Hz_stim(channel_interest,:) = ephys.spike_stim_epoch(channel_interest,spike_locs);
    ephys.spikeRate_10ms40Hz_nostim(channel_interest,:) = ephys.spikeRate_nostim_epoch(channel_interest,spike_locs);
    ephys.spikeRateMean_10ms40Hz_stim(channel_interest) = nanmean(ephys.spike_stim_epoch(channel_interest,spike_locs));
    ephys.spikeRateMean_10ms40Hz_nostim(channel_interest) = nanmean(ephys.spike_nostim_epoch(channel_interest,spike_locs));
end

% 20 ms 10 Hz
pulse_dur = ismember(ephys.stimuli_info(:,1),0.02)';
pulse_freq = ismember(ephys.stimuli_info(:,2),10)';
dur20ms_10Hz = pulse_dur.*pulse_freq; % for 20 ms and 10 Hz

for channel_interest = 1:length(ephys.channel)
    spike_locs = find(dur20ms_10Hz);
    ephys.spikeRate_20ms10Hz_stim(channel_interest,:) = ephys.spike_stim_epoch(channel_interest,spike_locs);
    ephys.spikeRate_20ms10Hz_nostim(channel_interest,:) = ephys.spikeRate_nostim_epoch(channel_interest,spike_locs);
    ephys.spikeRateMean_20ms10Hz_stim(channel_interest) = nanmean(ephys.spike_stim_epoch(channel_interest,spike_locs));
    ephys.spikeRateMean_20ms10Hz_nostim(channel_interest) = nanmean(ephys.spike_nostim_epoch(channel_interest,spike_locs));
end

% 20 ms 20 Hz
pulse_dur = ismember(ephys.stimuli_info(:,1),0.02)';
pulse_freq = ismember(ephys.stimuli_info(:,2),20)';
dur20ms_20Hz = pulse_dur.*pulse_freq; % for 20 ms and 10 Hz

for channel_interest = 1:length(ephys.channel)
    spike_locs = find(dur20ms_20Hz);
    ephys.spikeRate_20ms20Hz_stim(channel_interest,:) = ephys.spike_stim_epoch(channel_interest,spike_locs);
    ephys.spikeRate_20ms20Hz_nostim(channel_interest,:) = ephys.spikeRate_nostim_epoch(channel_interest,spike_locs);
    ephys.spikeRateMean_20ms20Hz_stim(channel_interest) = nanmean(ephys.spike_stim_epoch(channel_interest,spike_locs));
    ephys.spikeRateMean_20ms20Hz_nostim(channel_interest) = nanmean(ephys.spike_nostim_epoch(channel_interest,spike_locs));
end

% 20 ms 40 Hz
pulse_dur = ismember(ephys.stimuli_info(:,1),0.02)';
pulse_freq = ismember(ephys.stimuli_info(:,2),40)';
dur20ms_40Hz = pulse_dur.*pulse_freq; % for 20 ms and 40 Hz

for channel_interest = 1:length(ephys.channel)
    spike_locs = find(dur20ms_40Hz);
    ephys.spikeRate_20ms40Hz_stim(channel_interest,:) = ephys.spike_stim_epoch(channel_interest,spike_locs);
    ephys.spikeRate_20ms40Hz_nostim(channel_interest,:) = ephys.spikeRate_nostim_epoch(channel_interest,spike_locs);
    ephys.spikeRateMean_20ms40Hz_stim(channel_interest) = nanmean(ephys.spike_stim_epoch(channel_interest,spike_locs));
    ephys.spikeRateMean_20ms40Hz_nostim(channel_interest) = nanmean(ephys.spike_nostim_epoch(channel_interest,spike_locs));
end

%******************************************************************************
% Plot the firing rate in two conditions

background = [ephys.spikeRateMean_2ms10Hz_nostim; ephys.spikeRateMean_2ms20Hz_nostim;...
    ephys.spikeRateMean_2ms40Hz_nostim; ephys.spikeRateMean_5ms10Hz_nostim;...
    ephys.spikeRateMean_5ms20Hz_nostim; ephys.spikeRateMean_5ms40Hz_nostim;...
    ephys.spikeRateMean_10ms10Hz_nostim; ephys.spikeRateMean_10ms20Hz_nostim;...
    ephys.spikeRateMean_10ms40Hz_nostim; ephys.spikeRateMean_20ms10Hz_nostim;...
    ephys.spikeRateMean_20ms20Hz_nostim; ephys.spikeRateMean_20ms40Hz_nostim];
    
background_spikeRate = mean(background);    
    
fig = figure;

% for 2 ms 10 Hz
subplot(3,4,1)
x = [1 2];
y = [mean(ephys.spikeRateMean_2ms10Hz_stim) mean(background_spikeRate)];
value = bar(x,y,'FaceColor',[.7 .7 .7],'EdgeColor','none');

box off
xticks([1 2])
xticklabels({'Stim','No Stim'})
ylabel({'10 Hz';'Firing rate (Hz)'})

hold on
xdata = ones(length(ephys.channel),2);
xdata = x.*xdata;
data1 = ephys.spikeRateMean_2ms10Hz_stim;
data2 = background_spikeRate;

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
y = [mean(ephys.spikeRateMean_2ms20Hz_stim) mean(background_spikeRate)];
value = bar(x,y,'FaceColor',[.7 .7 .7],'EdgeColor','none');

box off
xticks([1 2])
xticklabels({'Stim','No Stim'})
ylabel({'20 Hz';'Firing rate (Hz)'})

hold on
xdata = ones(length(ephys.channel),2);
xdata = x.*xdata;
data1 = ephys.spikeRateMean_2ms20Hz_stim;
data2 = background_spikeRate;

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
y = [mean(ephys.spikeRateMean_2ms40Hz_stim) mean(background_spikeRate)];
value = bar(x,y,'FaceColor',[.7 .7 .7],'EdgeColor','none');

box off
xticks([1 2])
xticklabels({'Stim','No Stim'})
ylabel({'40 Hz';'Firing rate (Hz)'})

hold on
xdata = ones(length(ephys.channel),2);
xdata = x.*xdata;
data1 = ephys.spikeRateMean_2ms40Hz_stim;
data2 = background_spikeRate;

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
y = [mean(ephys.spikeRateMean_5ms10Hz_stim) mean(background_spikeRate)];
value = bar(x,y,'FaceColor',[.7 .7 .7],'EdgeColor','none');

box off
xticks([1 2])
xticklabels({'Stim','No Stim'})
ylabel('Firing rate (Hz)')

hold on
xdata = ones(length(ephys.channel),2);
xdata = x.*xdata;
data1 = ephys.spikeRateMean_5ms10Hz_stim;
data2 = background_spikeRate;

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
y = [mean(ephys.spikeRateMean_5ms20Hz_stim) mean(background_spikeRate)];
value = bar(x,y,'FaceColor',[.7 .7 .7],'EdgeColor','none');

box off
xticks([1 2])
xticklabels({'Stim','No Stim'})
ylabel('Firing rate (Hz)')

hold on
xdata = ones(length(ephys.channel),2);
xdata = x.*xdata;
data1 = ephys.spikeRateMean_5ms20Hz_stim;
data2 = background_spikeRate;

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
y = [mean(ephys.spikeRateMean_5ms40Hz_stim) mean(background_spikeRate)];
value = bar(x,y,'FaceColor',[.7 .7 .7],'EdgeColor','none');

box off
xticks([1 2])
xticklabels({'Stim','No Stim'})
ylabel('Firing rate (Hz)')

hold on
xdata = ones(length(ephys.channel),2);
xdata = x.*xdata;
data1 = ephys.spikeRateMean_5ms40Hz_stim;
data2 = background_spikeRate;

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
y = [mean(ephys.spikeRateMean_10ms10Hz_stim) mean(background_spikeRate)];
value = bar(x,y,'FaceColor',[.7 .7 .7],'EdgeColor','none');

box off
xticks([1 2])
xticklabels({'Stim','No Stim'})
ylabel('Firing rate (Hz)')

hold on
xdata = ones(length(ephys.channel),2);
xdata = x.*xdata;
data1 = ephys.spikeRateMean_10ms10Hz_stim;
data2 = background_spikeRate;

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
y = [mean(ephys.spikeRateMean_10ms20Hz_stim) mean(background_spikeRate)];
value = bar(x,y,'FaceColor',[.7 .7 .7],'EdgeColor','none');

box off
xticks([1 2])
xticklabels({'Stim','No Stim'})
ylabel('Firing rate (Hz)')

hold on
xdata = ones(length(ephys.channel),2);
xdata = x.*xdata;
data1 = ephys.spikeRateMean_10ms20Hz_stim;
data2 = background_spikeRate;

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
y = [mean(ephys.spikeRateMean_10ms40Hz_stim) mean(background_spikeRate)];
value = bar(x,y,'FaceColor',[.7 .7 .7],'EdgeColor','none');

box off
xticks([1 2])
xticklabels({'Stim','No Stim'})
ylabel('Firing rate (Hz)')

hold on
xdata = ones(length(ephys.channel),2);
xdata = x.*xdata;
data1 = ephys.spikeRateMean_10ms40Hz_stim;
data2 = background_spikeRate;

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
y = [mean(ephys.spikeRateMean_20ms10Hz_stim) mean(background_spikeRate)];
value = bar(x,y,'FaceColor',[.7 .7 .7],'EdgeColor','none');

box off
xticks([1 2])
xticklabels({'Stim','No Stim'})
ylabel('Firing rate (Hz)')

hold on
xdata = ones(length(ephys.channel),2);
xdata = x.*xdata;
data1 = ephys.spikeRateMean_20ms10Hz_stim;
data2 = background_spikeRate;

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
y = [mean(ephys.spikeRateMean_20ms20Hz_stim) mean(background_spikeRate)];
value = bar(x,y,'FaceColor',[.7 .7 .7],'EdgeColor','none');

box off
xticks([1 2])
xticklabels({'Stim','No Stim'})
ylabel('Firing rate (Hz)')

hold on
xdata = ones(length(ephys.channel),2);
xdata = x.*xdata;
data1 = ephys.spikeRateMean_20ms20Hz_stim;
data2 = background_spikeRate;

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
y = [mean(ephys.spikeRateMean_20ms40Hz_stim) mean(background_spikeRate)];
value = bar(x,y,'FaceColor',[.7 .7 .7],'EdgeColor','none');

box off
xticks([1 2])
xticklabels({'Stim','No Stim'})
ylabel('Firing rate (Hz)')

hold on
xdata = ones(length(ephys.channel),2);
xdata = x.*xdata;
data1 = ephys.spikeRateMean_20ms40Hz_stim;
data2 = background_spikeRate;

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
fig_name = strcat('Summary_firingRate');
fig_name = strrep(fig_name,'.','');
saveas(fig,fig_name)

%% figure plot the raw data, filtered voltage trace and PSTH (Optional)
bin_window = bin*30000; % for PSTH

for stim_num = 1:36%length(PulseSeg) % should be 36 stimuli
    for channel_interest = 1:length(ephys.channel)
        fig = figure;
        
        % **************for plotting raw and filtered trace************
        subplot(3,1,1)
        y_small = min(squeeze(ephys.Spike_filtered_seg(channel_interest, stim_num,12000:135000)));
        y_large = max(squeeze(ephys.Spike_filtered_seg(channel_interest, stim_num,12000:135000)));
        
        % plot stim
        x1 = double([ephys.stimON{stim_num}]-[ephys.stimON{stim_num}(1)])+3000;
        x2 = double([ephys.stimOFF{stim_num}]-[ephys.stimON{stim_num}(1)])+3000;
        
        for pulse_num = 1:length(x1)
            x = [x1(pulse_num) x1(pulse_num) x2(pulse_num) x2(pulse_num)];
            y = [y_small-500 y_large+500 y_large+500 y_small-500];
            h = fill(x,y,'r','EdgeColor','none'); % for antibiased
            h(1).FaceColor = [0,0.4,0.9];
            h(1).FaceAlpha = 0.1;
            hold on
        end
        
        ax = gca;
        yyaxis left
        plot(squeeze(ephys.Spike_filtered_seg(channel_interest,stim_num,12000:135000)),'Color',[1.00,0.52,0.69],'LineWidth',1);
        axis tight
        xl = xlim;
        ylabel({'Spike channel (uV)';'250-8000 Hz'},'Color',[1.00,0.52,0.695])
        ax.YAxis(1).Color = [1.00,0.52,0.695];
        hold on
        spike_time = find(squeeze(ephys.spike_binary(channel_interest,stim_num,12000:135000))==1);
        spike_time = spike_time';
        plot(spike_time,ones(1,length(spike_time))*y_large+3,'o','MarkerSize',2,...
            'MarkerEdgeColor','k','MarkerFaceColor','k');
        hold on
        xticks([0.1 1.1 2.1 3.1 4.1]*30000)
        xticklabels({'0','1','2','3','4'})
        xlabel('Time (s)')
        ylim([y_small-5 y_large+8])
        yyaxis right
        plot(squeeze(ephys.raw_seg(channel_interest,stim_num,12000:135000)),'Color','k','LineWidth',1);
        ylabel({'Raw voltage(uV)'},'Color','k')
        ax.YAxis(2).Color = [0,0,0];
        
        y_small = min(squeeze(ephys.raw_seg(channel_interest,stim_num,12000:135000)));
        y_large = max(squeeze(ephys.raw_seg(channel_interest,stim_num,12000:135000)));
        
        ylim([y_small-10 y_large+10])
        title('Raw and filtered traces');
        box off
        
        %*************plot filtered spike and LFP********************
        subplot(3,1,2)        
        y_small = min(squeeze(ephys.Spike_filtered_seg(channel_interest, stim_num,12000:135000)));
        y_large = max(squeeze(ephys.Spike_filtered_seg(channel_interest, stim_num,12000:135000)));
        
        % plot stim
        x1 = double([ephys.stimON{stim_num}]-[ephys.stimON{stim_num}(1)])+3000;
        x2 = double([ephys.stimOFF{stim_num}]-[ephys.stimON{stim_num}(1)])+3000;
        
        for pulse_num = 1:length(x1)
            x = [x1(pulse_num) x1(pulse_num) x2(pulse_num) x2(pulse_num)];
            y = [y_small-500 y_large+500 y_large+500 y_small-500];
            h = fill(x,y,'r','EdgeColor','none'); % for antibiased
            h(1).FaceColor = [0,0.4,0.9];
            h(1).FaceAlpha = 0.1;
            hold on
        end
        
        ax = gca;
        yyaxis left
        plot(squeeze(ephys.Spike_filtered_seg(channel_interest,stim_num,12000:135000)),'Color',[1.00,0.52,0.69],'LineWidth',1);
        axis tight
        xl = xlim;
        ylabel({'Spike channel (uV)';'250-8000 Hz'},'Color',[1.00,0.52,0.695])
        ax.YAxis(1).Color = [1.00,0.52,0.695];
        hold on
        spike_time = find(squeeze(ephys.spike_binary(channel_interest,stim_num,12000:135000))==1);
        spike_time = spike_time';
        plot(spike_time,ones(1,length(spike_time))*y_large+5,'o','MarkerSize',2,...
            'MarkerEdgeColor','k','MarkerFaceColor','k');
        hold on
        yline(-ephys.threshold(channel_interest,stim_num),'--b');
        hold on
        xticks([0.1 1.1 2.1 3.1 4.1]*30000)
        xticklabels({'0','1','2','3','4'})
        xlabel('Time (s)')
        ylim([y_small-5 y_large+5])
        yyaxis right
        plot(squeeze(ephys.LFP_filtered_seg(channel_interest,stim_num,12000:135000)),'Color','k','LineWidth',1);
        ylabel({'LFP channel (uV)';'0.7-170 Hz'},'Color','k')
        ax.YAxis(2).Color = [0,0,0];
        
        y_small = min(squeeze(ephys.LFP_filtered_seg(channel_interest,stim_num,12000:135000)));
        y_large = max(squeeze(ephys.LFP_filtered_seg(channel_interest,stim_num,12000:135000)));
        
        ylim([y_small-10 y_large+10])
        title('Spike LFP');
        box off
        
        %*************plot PSTH*******************************************
        subplot(3,1,3)
        y_small = 0;
        y_large = max(movmean(squeeze(ephys.spike_binary(channel_interest,stim_num,12000:135000)),bin_window)*150)+5;
        
        % plot stim
        x1 = double([ephys.stimON{stim_num}]-[ephys.stimON{stim_num}(1)])+3000;
        x2 = double([ephys.stimOFF{stim_num}]-[ephys.stimON{stim_num}(1)])+3000;
        
        for pulse_num = 1:length(x1)
            x = [x1(pulse_num) x1(pulse_num) x2(pulse_num) x2(pulse_num)];
            y = [y_small y_large y_large y_small];
            h = fill(x,y,'r','EdgeColor','none'); % for antibiased
            h(1).FaceColor = [0,0.4,0.9];
            h(1).FaceAlpha = 0.1;
            hold on
        end
        
        ax = gca;
        plot(movmean(squeeze(ephys.spike_binary(channel_interest,stim_num,12000:135000)),bin_window)*150,'Color',[1.00,0.52,0.69],'LineWidth',1);
        
        axis tight
        xl = xlim;
        ylabel({'Number of spikes'},'Color',[1.00,0.52,0.695])
        ax.YAxis(1).Color = [1.00,0.52,0.695];
        hold on
        xticks([0.1 1.1 2.1 3.1 4.1]*30000)
        xticklabels({'0','1','2','3','4'})
        xlabel('Time (s)')
        ylim([y_small y_large])        
        
        title('PSTH');
        box off
        
        info = ['Channel', ' ', num2str(ephys.channel(channel_interest)),', ','Pulse', ' ',num2str(ephys.stimuli_info(stim_num,1)),' ','sec,', ' ', num2str(ephys.stimuli_info(stim_num,2)),' ','Hz'];
        sgtitle(info);
        box off
        
        set(gcf,'position',[10,10,550,950])

        fig_name = strcat('channel',num2str(ephys.channel(channel_interest)),'Pulse',num2str(ephys.stimuli_info(stim_num,1)),'s', num2str(ephys.stimuli_info(stim_num,2)),'Hz','stim',num2str(stim_num));
        fig_name = strrep(fig_name,'.','');
        saveas(fig,fig_name)
        close;
    end
end


%% figure plot the PSTH with user-defined bin (Optional)
% bin_window = bin*30000;
% 
% for stim_num = 1:12%length(PulseSeg) % should be 36 stimuli
%     for channel_interest = 1:length(ephys.channel)
%         fig = figure;
%         info = ['PSTH',' ','Channel', ' ', num2str(ephys.channel(channel_interest)),', ','Pulse', ' ',num2str(ephys.stimuli_info(stim_num,1)),' ','sec,', ' ', num2str(ephys.stimuli_info(stim_num,2)),' ','Hz'];
%         
%         y_small = 0;
%         y_large = max(movmean(squeeze(ephys.spike_binary(channel_interest,stim_num,12000:135000)),bin_window)*150)+5;
%         
%         % plot stim
%         x1 = double([ephys.stimON{stim_num}]-[ephys.stimON{stim_num}(1)])+3000;
%         x2 = double([ephys.stimOFF{stim_num}]-[ephys.stimON{stim_num}(1)])+3000;
%         
%         for pulse_num = 1:length(x1)
%             x = [x1(pulse_num) x1(pulse_num) x2(pulse_num) x2(pulse_num)];
%             y = [y_small y_large y_large y_small];
%             h = fill(x,y,'r','EdgeColor','none'); % for antibiased
%             h(1).FaceColor = [0,0.4,0.9];
%             h(1).FaceAlpha = 0.1;
%             hold on
%         end
%         
%         ax = gca;
%         plot(movmean(squeeze(ephys.spike_binary(channel_interest,stim_num,12000:135000)),bin_window)*150,'Color',[1.00,0.52,0.69],'LineWidth',1);
%         
%         axis tight
%         xl = xlim;
%         ylabel({'Number of spikes'},'Color',[1.00,0.52,0.695])
%         ax.YAxis(1).Color = [1.00,0.52,0.695];
%         hold on
%         xticks([0.1 1.1 2.1 3.1 4.1]*30000)
%         xticklabels({'0','1','2','3','4'})
%         xlabel('Time (s)')
%         ylim([y_small y_large])        
%         
%         title(info);
%         box off
%         
%         fig_name = strcat('PSTH',num2str(ephys.channel(channel_interest)),'Pulse',num2str(ephys.stimuli_info(stim_num,1)),'s', num2str(ephys.stimuli_info(stim_num,2)),'Hz');
%         fig_name = strrep(fig_name,'.','');
%         saveas(fig,fig_name)
%         close;
%     end
% end
% 
% 
% %% figure plot raw data
% for stim_num = 1:12%length(PulseSeg) % should be 36 stimuli
%     for channel_interest = 1:length(ephys.channel)
%         fig = figure;
%         info = ['Raw trace',' ','Channel', ' ', num2str(ephys.channel(channel_interest)),', ','Pulse', ' ',num2str(ephys.stimuli_info(stim_num,1)),' ','sec,', ' ', num2str(ephys.stimuli_info(stim_num,2)),' ','Hz'];
%         
%         y_small = min(squeeze(ephys.Spike_filtered_seg(channel_interest, stim_num,12000:end)));
%         y_large = max(squeeze(ephys.Spike_filtered_seg(channel_interest, stim_num,12000:end)));
%         
%         % plot stim
%         x1 = double([ephys.stimON{stim_num}]-[ephys.stimON{stim_num}(1)])+3000;
%         x2 = double([ephys.stimOFF{stim_num}]-[ephys.stimON{stim_num}(1)])+3000;
%         
%         for pulse_num = 1:length(x1)
%             x = [x1(pulse_num) x1(pulse_num) x2(pulse_num) x2(pulse_num)];
%             y = [y_small-500 y_large+500 y_large+500 y_small-500];
%             h = fill(x,y,'r','EdgeColor','none'); % for antibiased
%             h(1).FaceColor = [0,0.4,0.9];
%             h(1).FaceAlpha = 0.1;
%             hold on
%         end
%         
%         ax = gca;
%         yyaxis left
%         plot(squeeze(ephys.Spike_filtered_seg(channel_interest,stim_num,12000:end)),'Color',[1.00,0.52,0.69],'LineWidth',1);
%         axis tight
%         xl = xlim;
%         ylabel({'Spike channel (uV)';'250-8000 Hz'},'Color',[1.00,0.52,0.695])
%         ax.YAxis(1).Color = [1.00,0.52,0.695];
%         hold on
%         spike_time = find(squeeze(ephys.spike_binary(channel_interest,stim_num,12000:end))==1);
%         spike_time = spike_time';
%         plot(spike_time,ones(1,length(spike_time))*y_large+3,'o','MarkerSize',2,...
%             'MarkerEdgeColor','k','MarkerFaceColor','k');
%         hold on
%         xticks([0.1 0.6 1.1]*30000)
%         xticklabels({'0','0.5','1'})
%         xlabel('Time (s)')
%         ylim([y_small-5 y_large+8])
%         yyaxis right
%         plot(squeeze(ephys.raw_seg(channel_interest,stim_num,12000:end)),'Color','k','LineWidth',1);
%         ylabel({'LFP channel (uV)';'0.7-170 Hz'},'Color','k')
%         ax.YAxis(2).Color = [0,0,0];
%         
%         y_small = min(squeeze(ephys.raw_seg(channel_interest,stim_num,12000:end)));
%         y_large = max(squeeze(ephys.raw_seg(channel_interest,stim_num,12000:end)));
%         
%         ylim([y_small-10 y_large+10])
%         title(info);
%         box off
%         
%         fig_name = strcat('Raw',num2str(ephys.channel(channel_interest)),'Pulse',num2str(ephys.stimuli_info(stim_num,1)),'s', num2str(ephys.stimuli_info(stim_num,2)),'Hz');
%         fig_name = strrep(fig_name,'.','');
%         saveas(fig,fig_name)
%         close;
%     end
% end

%% SETTING FILE- AND PATHNAME FOR SAVING RESULTS
% [saveFile savePath] = uiputfile('*.mat','save results as',...
%     'D:\Data\optogenetics\Open Ephys\analized data\20211021\analysis');

prompt = {'Reulst name to save'};
dlgTitle = 'User-defined input';
nLines = 1;
def = {'20211027_Record3'};
answer = inputdlg(prompt, dlgTitle, nLines, def);

file_name = answer{1};
file_name = [file_name,'.','mat'];
save(file_name,'ephys','-v7.3')
clear;
clc;