%%

clear; close all; clc

fig = openfig('channel21Pulse0002s10Hzstim7_simplified.fig');


%%

scaleX = 30000;

dataObjs = axObjs.Children;

Y = dataObjs(1).YData';
Y = Y(round(min(xlim)):round(max(xlim)));

X = dataObjs(1).XData';
X = X(round(min(xlim)):round(max(xlim)))/scaleX;

outData = [X Y];

stimData= zeros(10,2);
del = [];
for ii = 1:10
    stimData(ii,:) = dataObjs(ii+1).XData([1 3]);
    if ~inrange(stimData(ii,1),xlim)
        del(end+1) = ii;
    end
end

stimData(del,:) = [];
stimData = flipud(stimData)/scaleX;




