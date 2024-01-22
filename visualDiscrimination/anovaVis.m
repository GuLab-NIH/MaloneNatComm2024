%%

load('dataVis.mat')

useDays = 1:10;
cats = {[1 3:5 8 10:14],[2 6 7 9]};
catNames = {'Good Learners','Bad Learners'};

data1 = dataVis(useDays,cats{1})';
data2 = dataVis(useDays,cats{2})';

[pA,pMC] = anovaRM2W(data1,data2);

% pMC = pMC';

