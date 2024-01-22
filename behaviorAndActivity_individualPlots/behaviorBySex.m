%mice are in this order:

clear; close all; clc

% % Slowing
% ylab = 'Predictive Slowing (%)';
% A = [78.15801024
% 62.52792071
% 60.00126751
% 79.71455351
% 85.03011564
% 74.70156268
% 72.00241561
% ];
% B = [81.85264945
% 74.96015417
% 75.06798356
% 68.97577327
% 56.49964455
% 65.69598949
% 71.69336514
% 73.16132358
% ];

% Licking
ylab = 'Predictive licking (%)';
A = [5.281352134
23.98993002
16.994132
45.28930982
56.10216657
58.51162399
36.06570096
];
B = [31.82425214
22.06346899
80.46918768
55.17327634
3.430833247
10.73593074
27.0095018
61.34037169
];

meanA = nanmean(A(:,1));
meanB = nanmean(B(:,1));

semA = nansem(A,1);
semB = nansem(B,1);

figure; hold on

x1 = [1 2];
y1 = [meanA meanB];
e1 = [semA semB];

h = bar({'Female','Male'},y1,'FaceColor','none','EdgeColor',[0 0 0 ]);
er = errorbar(x1,y1,e1,'Color',[0 0 0],'LineStyle','none');
h.CData = [0.5,0.5,0.5];


% simple version
scatter(repmat(h.XEndPoints(1), length(A), 1),A,30,'MarkerFaceColor','k',...
    'MarkerEdgeColor','none','XJitter','randn','XJitterWidth',.1)
scatter(repmat(h.XEndPoints(2), length(B), 1),B,30,'MarkerFaceColor','k',...
    'MarkerEdgeColor','none','XJitter','randn','XJitterWidth',.1)
ylabel(ylab)
ylim([0 100])
