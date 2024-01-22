%%

clear all; close all; clc

d = dir('data_*');
nFile = length(d);

PR = [0.78,0.977];

outCol = [1 2 3 4];
outType = {'S-S','S-O','O-O','O-S'};

%%

dists = cell(nFile,4);

for f = 1:nFile
    load(d(f).name)
    
    if f==1 && strcmp(pwd,'D:\Novelty_Analysis\clusterMan\good')
        pr = PR(1);
    else
        pr = PR(2);
    end
    
    %%
        
    centLocs = allCen;
    popA = cellsShr;
    popB = cellsOth;
    szA = length(popA);
    szB = length(popB);
    
    szMx = max(szA,szB);
    
    SS = pdist2(centLocs(popA,:),centLocs(popA,:),'euclidean','Smallest',szMx)*pr;
    SO = pdist2(centLocs(popA,:),centLocs(popB,:),'euclidean','Smallest',szMx)*pr;
    OO = pdist2(centLocs(popB,:),centLocs(popB,:),'euclidean','Smallest',szMx)*pr;
    OS = pdist2(centLocs(popB,:),centLocs(popA,:),'euclidean','Smallest',szMx)*pr;
    
    for ss = 1:szMx
        if ss<=size(SS,1)-1
            dists{f,outCol(1)}(:,ss) = mean(SS(2:ss+1,:),1,'omitnan')';
        end
        
        if ss<=size(SO,1)
            dists{f,outCol(2)}(:,ss) = mean(SO(1:ss,:),1,'omitnan')';
        end
        
        if ss<=size(OO,1)-1
            dists{f,outCol(3)}(:,ss) = mean(OO(2:ss+1,:),1,'omitnan')';
        end
        
        if ss<=size(OS,1)
            dists{f,outCol(4)}(:,ss) = mean(OS(1:ss,:),1,'omitnan')';
        end
    end
    
end

save('distsAll.mat','dists')


%% Set parameters

load('distsAll.mat')

close all

outType = {'S-S','S-O','O-O','O-S'};
limX = 200;
smOn = 0;

useIdx = 1:9;%1:length(dists);
distsUse = dists(useIdx,:);

figure
useT = erase(pwd,'D:\Novelty_Analysis\clusterMan\');
ttl = [useT ' learners'];
sgtitle('Cluster Anlaysis')


%% Plot data by FOV

distsAll = cell(size(distsUse));

% calculate raw means by FOV
for ii = 1:size(distsUse,1)
    curTTL = ['datatset ' num2str(ii)];
    curDist = distsUse(ii,:);

    for jj = 1:size(distsUse,2)
        distsAll{ii,jj} = mean(distsUse{ii,jj},1,'omitnan');
    end
end

% calculate number of elements by FOV
nA = cellfun(@size,distsAll,'UniformOutput',false)';
nA = vertcat(nA{:});
nN = max(nA(:,2));

% set values beyond element number
for ii = 1:size(distsUse,1)
    for jj = 1:size(distsUse,2)
        if length(distsAll{ii,jj})<nN
            if smOn
                distsAll{ii,jj}(end+1:nN) = distsAll{ii,jj}(end);
            else
                distsAll{ii,jj}(end+1:nN) = NaN;
            end
        end
    end
end

% combine data by distance outType
distsCom = cell(1,size(distsUse,2));
for jj = 1:size(distsUse,2)
    distsCom{jj} = vertcat(distsAll{:,jj});
end

% plot by FOV
subplot(1,2,1); hold on;
plotPDISTns(distsCom,outType,[ttl ' by FOV'],1)
xlim([1 limX])
% legend(outType)


%% Plot data by cell

DISTSF = cell(1,size(distsUse,2));

for ii = 1:size(distsUse,2)
    % calculate raw means by cell
    for jj = 1:size(distsUse,1)
        DISTSF{ii}(end+1:end+size(distsUse{jj,ii},1),1:size(distsUse{jj,ii},2)) = distsUse{jj,ii};
    end
    
    % set values beyond element number
    for jj = 1:size(DISTSF{ii})
        % set post element value
        if smOn
            rep = DISTSF{ii}(jj,find(DISTSF{ii}(jj,:)~=0,1,'last'));
        else
            rep = NaN;
        end
        
        % implement post element value
        DISTSF{ii}(jj,DISTSF{ii}(jj,:)==0) = rep;
        
        if size(DISTSF{ii},2)<nN
            DISTSF{ii}(jj,end+1:nN) = rep;
        end
    end
end

% plot by cell
subplot(1,2,2); hold on;
plotPDISTns(DISTSF,outType,[ttl ' by cell'])
xlim([1 limX])


%% Save figures

savefig(['clusters_' useT '_' num2str(limX) '.fig'])


