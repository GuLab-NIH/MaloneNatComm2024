%%

good = 1;
bad = 0;

use = good;

if use
    cd('Z:\labMembers\YG\data\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells')
    outFolder = 'D:\Novelty_Analysis\clusterMan\good';
else
    cd('Z:\labMembers\YG\data\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\commonUniqueCells')
    outFolder = 'D:\Novelty_Analysis\clusterMan\bad';
end

rmN = 'E:';
preN = 'Z:\labMembers\YG\data';

load('foldersAllN.mat')

p = pwd;

%%

for ii = st:length(foldersAllN)
    
    newF = [preN erase(foldersAllN{ii},rmN)];
    
    cd(newF)
    
    %% Identify common cells
    
    load('commonCells.mat')
    
    if size(commonCells,1)==1
        disp('Error\: insufficient cells')
        return
    else
        anchor = commonCells([1,end],:);
    end
    
    
    %% Identify ID groups
    
    
    
    
    %
    
    load('minDistManualAct.mat')
    
    cellsAll = minDistManualAct(:,1);
    cellN = length(cellsAll);
    cellsActive = cellsAll(minDistManualAct(:,2)<=thrSize);
    idxActive = minDistManualAct(cellsActive,3);
    
    load('allROIsManual.mat')
    allCen = zeros(cellN,2);
    for jj = 1:cellN
        cenCur = regionprops(allROIsManual(:,:,jj),'Centroid');
        allCen(jj,:) = cenCur.Centroid;
    end
    
    fid=fopen('ROISizes.m');
    lineNum = 82;
    manFile = textscan(fid,'%s',1,'delimiter','\n','headerlines',lineNum-1);
    fclose(fid);
    
    disp(manFile{1}{1})
    cd('D:\Novelty_Analysis\Alignments')
    
    
    %% Manually move to correct shared folder
    
    load('novel.mat')
    idx = input('Select FOV: ');
    idxShared = novel(idx).share;
    
    cellsShr = cellsActive(ismember(idxActive,idxShared));
    cellsOth = setdiff(cellsAll,cellsShr);
    
    save([outFolder '\data_' num2str(ii)],'cellsShr','cellsOth','allCen')
    
    
    %%
    
    st = st+1;
    cd(p)
end

