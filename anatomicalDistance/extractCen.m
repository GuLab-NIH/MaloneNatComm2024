%%

clear; close all; clc

good = 1;
bad = 0;

use = bad;

if use
    cd('Z:\SNM\labMembers\YG\data\learningAnalysis\summaryManyMice_includingOldEnvGoodMoreMice\commonUniqueCells')
    outFolder = 'Z:\SNM\labMembers\TM\Novelty_Analysis\clusterMan\good';
else
    cd('Z:\SNM\labMembers\YG\data\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\commonUniqueCells')
    outFolder = 'Z:\SNM\labMembers\TM\Novelty_Analysis\clusterMan\bad';
end

rmN = 'E:';
preN = 'Z:\SNM\labMembers\YG\data';

ancNew = [0 0;1 0];

load('foldersAllN.mat')

p = pwd;

%%

for ii = 1:length(foldersAllN)
    
    newF1 = [preN erase(foldersAllN{ii},rmN)];
    cd(newF1)
    
    %% Identify common cells
    
    load('commonCells.mat')
    
    if size(commonCells,1)==1
        disp('Error\: insufficient cells')
        return
    else
        anchor = commonCells([1,end],:);
    end
    
    
    %% Identify ID groups
    
    load('trackIDsAll.mat')
    
    IDs.common = trackIDsAll{1};
    
    IDs.unique = trackIDsAll(5:15);
        
    for jj = 1:size(commonCells,2)
        IDs.unique{jj}(ismember(IDs.unique{jj},commonCells(:,jj))) = [];
    end
    
    
    %% Process all FOVs
    
    load('useFolders.mat')
    
    allCen = struct();
    allCen.common = [];
    allCen.unique = [];
    
    for jj = 1:length(useFolders)
        
        newF2 = [preN erase(useFolders{jj},rmN)];
        cd(newF2)
        
        %% Calculate centroids
        
        load('allROIs.mat')
        curCen = zeros(size(roi,3),2);
        for kk = 1:size(roi,3)
            cenCur = regionprops(roi(:,:,kk),'Centroid');
            if isempty(cenCur)
                curCen(kk,:) = [NaN NaN];
            else
                curCen(kk,:) = cenCur.Centroid;
            end
        end
        
        if jj==1
            allCen.common = curCen(IDs.common,:);
            
            ancGlob = curCen(anchor(:,1),:);
            tformGlob = fitgeotrans(ancGlob,ancNew,'NonreflectiveSimilarity');
        elseif jj==2
            tformFilt = fitgeotrans(ancGlob,curCen(anchor(:,2),:),...
                'NonreflectiveSimilarity');
        end
        
        
        %% Transform centroids
        
        ancCur = curCen(anchor(:,jj),:);
        
        tform = fitgeotrans(ancCur,ancNew,'NonreflectiveSimilarity');
        
        X = transformPointsForward(tform,curCen);
        
        
        %% Update unique centroids
        
        curUnq = IDs.unique{jj};
        
        allCen.unique(end+1:end+length(curUnq),1:2) = X(curUnq,:);
    
        
    end
    
    
    %%
    
    unqTempD0 = transformPointsInverse(tformGlob,allCen.unique);
    
    limX = [0;0;512;512];
    limY = [0;512;512;0];
    
    % filter day 0
    inD0 = inpolygon(unqTempD0(:,1),unqTempD0(:,2),limX,limY);
%     inD0 = true(size(unqTempD0(:,1)));
    
    % filter day 1
    unqTempD1 = transformPointsForward(tformFilt,unqTempD0);
    inD1 = inpolygon(unqTempD1(:,1),unqTempD1(:,2),limX,limY);
%     inD1 = true(size(unqTempD0(:,1)));
    
    allCen.unique = unqTempD0(inD0&inD1,:);
    
    save([outFolder '\dataFx_' num2str(ii)],'allCen')
    
    cd(newF1)
    
end

cd(p)

