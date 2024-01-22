%% visAll.m

clear all; close all; clc

good = 1;
bad = 2;

folders = cell(1,2);
folders{good} = {'210206','210208','210209','210413'};
folders{bad} = {'210207','210519-1','210519-2'};

dPrime = cell(1,2);
nTrial  = cell(1,2);

for i = good:bad
    
    for f = 1:length(folders{i})
        load([folders{i}{f} '\data_full.mat'])
        
        dPrime{i}(f,1:length(data.dPrime.all)) = data.dPrime.all;
        nTrial{i}(f,1:length(data.trials.all)) = data.trials.all;
    end
    
%     for f = 1:length(folders{i})
%         load([folders{i}{f} '\data.mat'])
%         
%         lenMax = size(dPrime{i},2);
%         lenCur = size(data.dPrime.all,1);
%         
%         if lenMax>lenCur
%             dPrime{i}(f,lenCur+1:lenMax) = mean(maxk(dPrime{i}(f,:),4));
%         end
%     end
    
    dPrimeAvg(i,:) = mean(dPrime{i},1);
    dPrimeErr(i,:) = std(dPrime{i},1)/sqrt(size(dPrime{i},1));
    
end

dataAll = struct();
dataAll.dPrime = dPrime;
dataAll.dPrimeAvg = dPrimeAvg;
dataAll.dPrimeErr = dPrimeErr;
dataAll.nTrial = nTrial;

save('data.mat','dataAll')


%% Plot combined data

close all

colors = {'g','m'};
labels = {'good learners','bad learners'};
nDays = size(dPrimeAvg(:,:),2);

FIG = figure; hold on

for i = good:bad
    plot(dPrimeAvg(i,:),colors{i},'LineWidth',2,'DisplayName',labels{i});
    
%     for f = 1:length(folders{i})
%         scatter(1:nDays,dPrime{i}(f,:),colors{i},'HandleVisibility','off');
%     end
    
    errorbar(dPrimeAvg(i,:),dPrimeErr(i,:),'k','linestyle','none',...
        'LineWidth',1,'HandleVisibility','off')
end

plot([0 size(dPrimeAvg,2)],[1 1],'r','LineWidth',1,...
    'HandleVisibility','off')

title('Visual discrimination task learning')
xlabel('Training Day')
ylabel('d''')

set(get(gca,'YLabel'),'Rotation',0)
set(gca,'fontsize',22')

axis('square')
ylim([-0.5 4.5])
legend('Location','northwest')
legend boxoff

hold off

savefig(FIG,'visTrainingAll.fig')


%% Plot correlation data

allTrial = [];
allDPrime = [];

for i = 1:2
    allTrial = [allTrial; nTrial{i}(:)];
    allDPrime = [allDPrime; dPrime{i}(:)];
end

idx = find(allDPrime>2);
allDPrime = allDPrime(idx);
allTrial = allTrial(idx);

figure; hold on
scatter(allTrial,allDPrime)

[R,p] = corrcoef(allTrial,allDPrime);

disp(['Correlation: ' num2str(R(2,1))])
disp(['p-value: ' num2str(p(2,1))])


%%

saveas(FIG,'visTrainingAll.tif')

