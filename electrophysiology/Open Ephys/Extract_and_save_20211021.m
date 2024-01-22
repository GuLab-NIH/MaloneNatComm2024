%% To extract the voltage from 32 channels and timestamps for stimuli for experiment 20211021

Record1_data_20211021 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\max_ephys1_2021-10-21_16-34-10\Record Node 110\experiment1\recording1\structure.oebin','continuous',1);
Record1_stim_20211021 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\max_ephys1_2021-10-21_16-34-10\Record Node 110\experiment1\recording1\structure.oebin','events',1);

Record2_data_20211021 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\max_ephys2_2021-10-21_17-13-44\Record Node 110\experiment1\recording1\structure.oebin','continuous',1);
Record2_stim_20211021 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\max_ephys2_2021-10-21_17-13-44\Record Node 110\experiment1\recording1\structure.oebin','events',1);

Record3_data_20211021 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\max_ephys3_2021-10-21_17-32-54\Record Node 110\experiment1\recording1\structure.oebin','continuous',1);
Record3_stim_20211021 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\max_ephys3_2021-10-21_17-32-54\Record Node 110\experiment1\recording1\structure.oebin','events',1);

Record4_data_20211021 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\max_ephys4_2021-10-21_17-51-54\Record Node 110\experiment1\recording1\structure.oebin','continuous',1);
Record4_stim_20211021 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\max_ephys4_2021-10-21_17-51-54\Record Node 110\experiment1\recording1\structure.oebin','events',1);

Record5_data_20211021 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\max_ephys5_2021-10-21_18-10-38\Record Node 110\experiment1\recording1\structure.oebin','continuous',1);
Record5_stim_20211021 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\max_ephys5_2021-10-21_18-10-38\Record Node 110\experiment1\recording1\structure.oebin','events',1);

Record6_data_20211021 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\max_ephys6_2021-10-21_18-29-20\Record Node 110\experiment1\recording1\structure.oebin','continuous',1);
Record6_stim_20211021 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\max_ephys6_2021-10-21_18-29-20\Record Node 110\experiment1\recording1\structure.oebin','events',1);

Record7_data_20211021 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\max_ephys7_2021-10-21_18-47-40\Record Node 110\experiment1\recording1\structure.oebin','continuous',1);
Record7_stim_20211021 = load_open_ephys_binary('D:\Data\optogenetics\Open Ephys\DataFiles\max_ephys7_2021-10-21_18-47-40\Record Node 110\experiment1\recording1\structure.oebin','events',1);

%% To save the voltage from 32 channels and timestamps for stimuli for experiment 20211021
save('Record1_20211021','Record1_data_20211021','Record1_stim_20211021','-v7.3');
save('Record2_20211021','Record2_data_20211021','Record2_stim_20211021','-v7.3');
save('Record3_20211021','Record3_data_20211021','Record3_stim_20211021','-v7.3');
save('Record4_20211021','Record4_data_20211021','Record4_stim_20211021','-v7.3');
save('Record5_20211021','Record5_data_20211021','Record5_stim_20211021','-v7.3');
save('Record6_20211021','Record6_data_20211021','Record6_stim_20211021','-v7.3');
save('Record7_20211021','Record7_data_20211021','Record7_stim_20211021','-v7.3');
