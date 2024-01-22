%%

% clear; close all; clc

codeT = 'Fx';
d = dir(['data' codeT '_*']);
nFile = length(d);

PR = [0.78,0.977];

globOn = 1;
nType = 3;
nSubsets = 1000;

%%

dists = cell(nFile,nType);

for f = 1:nFile
    load(d(f).name)
    
    if f==1 && strcmp(pwd,'D:\Novelty_Analysis\clusterMan\good')
        pr = PR(1);
    else
        pr = PR(2);
    end
    
    %%
    
    if globOn
        load('minGlob.mat')
        useSz = minGlob;
    else
        useSz = min(size(allCen.common,1),size(allCen.unique,1));
    end
        
    cenCom = allCen.common;
    cenUnq = allCen.unique;
    
     distsCur = cell(1,4);
     
     for jj = 1:nSubsets
         popCom = randsample(size(cenCom,1),useSz);
         popUnq = randsample(size(cenUnq,1),useSz);
         
         useCom = cenCom(popCom,:);
         useUnq = cenUnq(popUnq,:);

         distsCom = pdist2(useCom,useCom,'euclidean','Smallest',useSz)*pr;
         distsUnq = pdist2(useUnq,useUnq,'euclidean','Smallest',useSz)*pr;
         distsCU = pdist2(useCom,useUnq,'euclidean','Smallest',useSz)*pr;
         distsUC = pdist2(useUnq,useCom,'euclidean','Smallest',useSz)*pr;
         
         for ss = 1:useSz-1
             distsCur{1}(:,ss,jj) = mean(distsCom(2:ss+1,:),1,'omitnan')';
             distsCur{2}(:,ss,jj) = mean(distsUnq(2:ss+1,:),1,'omitnan')';
             
             distsCur{3}(:,ss,jj) = mean(distsCU(1:ss,:),1,'omitnan')';
             distsCur{4}(:,ss,jj) = mean(distsUC(1:ss,:),1,'omitnan')';
         end
     end
     
     for kk = 1:3
         if kk<=2
             meanDists = mean(distsCur{kk},3);
         else
             meanCU = mean(distsCur{3},3);
             meanUC = mean(distsCur{4},3);
             meanDists = cat(1,meanCU,meanUC);
         end
         dists{f,kk} = meanDists;
     end
end

% save('distsAll.mat','dists')


%% Set parameters

% load('distsAll.mat')

% figure use distsFOVFx_type_global
%%
close all

outType = {'Common','Other','Cross'};
limX = 10;
smOn = 0;

useIdx = 1:size(dists,1);
% useIdx = randsample(size(dists,1),9);
distsUse = dists(useIdx,:);

figure; hold on
ylim([30 120])

useT = erase(pwd,'Z:\SNM\labMembers\TM\Novelty_Analysis\clusterMan\');
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
% subplot(1,2,1); hold on;
% plotPDISTall(distsCom,outType,[ttl ' by FOV'],1);
pVals = plotPDISTallFigure(distsCom,outType,[ttl ' by FOV'],1);

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
plotPDISTall(DISTSF,outType,[ttl ' by cell'])
xlim([1 limX])


%% Save figures

if globOn
    sType = 'global';
else
    sType = 'local';
end

% savefig(['clusters' codeT '_' useT '_' sType '.fig'])
% save(['distsFOV' codeT '_' useT '_' sType '.mat'],'dists','distsCom')

