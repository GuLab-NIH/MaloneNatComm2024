%use aligned common cells to register all FOVS
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\foldersAll.mat');
% foldersAll{3}='E:\ID20210206\findCells\loc2\subsetFOVsWithOldNo11\analyzeCells';
%above: fix one mismatch between fodlersAll and foldersAllN
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\commonUniqueCells\foldersAllN.mat');
p=pwd;
%there are total 23 FOVs in good performers. align individual FOVs to the
%first day
rotation={};%each cell is one FOV, withn that, each cell is the rotation of the roi to the first fov on day 0.
translation={};%same. each cell is one fov, withn that, each cell is the translation of the roi to the first fov on day 0.
newCOMs={};%new COMs according to the coms on day0 (old env)

for n=1:length(foldersAll);%last three FOVs are the same FOV
    disp(n)
    cd(foldersAll{n});
    load('commonCells.mat');
    load('useFolders.mat');
    rotation{n}={};
    translation{n}={};
    newCOMs{n}={};
    %use day 0 as reference
    
    cd(useFolders{1});
    load('COMs.mat');%first column is x (column), second column is y (row)
    oCOMs=COMs;
    i=commonCells(:,1); %use the first FOV (day 1) as a reference.
    iCOMs=oCOMs(i,:);      
    
    for m=1:length(useFolders);
        disp(m)
        cd(useFolders{m});
        ii=commonCells(:,m);
        load('COMs.mat');
        imCOMs=COMs(ii,:);
        [U, r, lrms] = Kabsch(imCOMs', iCOMs');
        rotation{n}{m}=U;
        translation{n}{m}=r;
        %convert rois
        newCOMs{n}{m}=[];
        
        for i=1:size(COMs,1);
             a=COMs(i,:)';
             b=U*a+r;
             newCOMs{n}{m}(i,:)=b';
        end            
    end
end
cd(p);

save('newCOMs.mat','newCOMs');
save('rotation.mat','rotation');
save('translation.mat','translation');

%% get the COMs of all cells

fovI=[1 2 1 2 1 2 3 4 5 6 7 8 9 10 11];%the FOV should be used
newCOMsTrackCells={};
for n=1:length(foldersAllN);
    cd(foldersAllN{n});
    load('trackIDsAll.mat');
    newCOMsTrackCells{n}={};
    for m=1:length(trackIDsAll);
        newCOMsTrackCells{n}{m}=newCOMs{n}{fovI(m)}(trackIDsAll{m},:);
    end
end
cd(p);
 save('newCOMsTrackCells.mat','newCOMsTrackCells');   
 %these are new coms of tracked cells
 %% see how they overlap with the manual cells
 load('newCOMsTrackCells.mat');
 newCOMsManualMatch={};
load('E:\learningAnalysis\summaryManyMice_includingOldEnvBadMoreMice\foldersAll.mat');
p=pwd;
for n=1:length(foldersAll);
    disp(n)
     newCOMsManualMatch{n}={};
     
    cd(foldersAll{n});
    load('pyrSte/allCells/comManual.mat');
    load('pyrSte/allCells/params.mat');
    for m=1:length(newCOMsTrackCells{n});
         newCOMsManualMatch{n}{m}=[];
         if ~isempty(newCOMsTrackCells{n}{m})
        for i=1:size(comManual,1);
           newCOMsManualMatch{n}{m}(i,1)=i;%first number is the manual cell ID
           newCOMsManualMatch{n}{m}(i,2)=params.longAxisum(i);%second number is the manual cell param
           allDis=[];
           for ii=1:size(newCOMsTrackCells{n}{m},1);
               allDis(ii)=pdist([comManual(i,:);newCOMsTrackCells{n}{m}(ii,:)]);
           end
           
           [newCOMsManualMatch{n}{m}(i,3),newCOMsManualMatch{n}{m}(i,4)]=min(allDis);%SECOND number is the manual cell and track cell distance, third number is the order in the track ID
        
        %change to um
        if n==1;
           newCOMsManualMatch{n}{m}(i,3)=newCOMsManualMatch{n}{m}(i,3)*0.7813;
        else
            newCOMsManualMatch{n}{m}(i,3)=newCOMsManualMatch{n}{m}(i,3)*0.9766;
        end
        end
         end
    end
end
cd(p);
        
%REMOVE cells that are too far
thresh=10;%distance smaller than 10um is the same cell
for n=1:length(newCOMsManualMatch);
    for m=1:length(newCOMsManualMatch{n});
        if ~isempty(newCOMsManualMatch{n}{m});
        d=newCOMsManualMatch{n}{m}(:,3);
        i=find(d<10);
        newCOMsManualMatch{n}{m}=newCOMsManualMatch{n}{m}(i,:);
        end
    end
end

 save('newCOMsManualMatch.mat','newCOMsManualMatch');   
%% determine ste pyra
load('E:\learningAnalysis\allStePyr\AddMoreMice\thresh.mat');
mustD=1;
steThresh=thresh+mustD;
pyrThresh=thresh-mustD;

stePyrAllCells={};

for n=1:length(newCOMsManualMatch);
    stePyrAllCells{n}={};
    for m=1:length(newCOMsManualMatch{n});
        stePyrAllCells{n}{m}=[];
       
        if ~isempty(newCOMsManualMatch{n}{m})
             yn=[];
        d=newCOMsManualMatch{n}{m}(:,2);
        for i=1:length(d);
            if d(i)>steThresh;
                yn(i)=1;
            elseif d(i)<pyrThresh;
                yn(i)=0;
            else
                yn(i)=nan;
            end
        end
        stePyrAllCells{n}{m}=yn(~isnan(yn));
        else
           stePyrAllCells{n}{m}=[];
        end
    end
end
     save('stePyrAllCells.mat','stePyrAllCells');       
stePer=[];
pyrPer=[];

for n=1:length(stePyrAllCells);
    for m=1:length(stePyrAllCells{n});
        i=stePyrAllCells{n}{m};
        if ~isempty(i);
        stePer(n,m)=length(find(i))/length(i);
         pyrPer(n,m)=1-length(find(i))/length(i);
        else
          stePer(n,m)=nan;
          pyrPer(n,m)=nan;
        end
    end
end

figure,
errorbar([1:1:15],nanmean(stePer,1),nansem(stePer,1));
hold on
errorbar([1:1:15],nanmean(pyrPer,1),nansem(pyrPer,1));
title('fractionOfStePyrAllTrackCells');
legend('ste','pyr','Location','Northwest')

saveas(gcf,'fractionOfStePyrAllTrackCells.fig');