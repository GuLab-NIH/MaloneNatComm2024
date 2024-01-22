%% Sigmoid Calculation
%
% Y=Bottom + (Top-Bottom)/(1+(IC50/X)^HillSlope)
%

clear

% manually input parameters 
params = [];


%% Calculate y values

bottom = params(1,:);
top = params(2,:);
IC50 = params(3,:);
slope = params(4,:);

N = size(params,2);

x = 1:0.1:20;
y = zeros(N,length(x));

for i = 1:N
    y(i,:) = bottom(i) + (top(i)-bottom(i))./(1+(IC50(i)./x).^slope(i));
end


%% Find learning threshold day

thresh = 2;

thDay = zeros(N,1);

for i = 1:N
    thIdx = find(y(i,:)>=thresh,1);
    thDay(i) = x(thIdx);
end


%% Plot results

figure; hold on

plot(x,y)



