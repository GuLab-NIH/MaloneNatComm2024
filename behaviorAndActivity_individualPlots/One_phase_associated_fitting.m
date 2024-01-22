%% curve fitting for neural and behavioral activity
% example
fake_data = [1 2 3 4 5 6 7 8 9];

raw_data = fake_data;
x = 1:1:length(fake_data);
y = raw_data;
[param,stat]=sigm_fit(x,y,[],[],[]);
corrIndv95 = (max(stat.ypred)-min(stat.ypred))*0.95+min(stat.ypred);
% find the closesd day
func = @(xval)param(1)+(param(2)-param(1))./(1+10.^((param(3)-xval)*param(4)))-corrIndv95;
fitting95_day.matrixCorrGood = fzero(func,1);