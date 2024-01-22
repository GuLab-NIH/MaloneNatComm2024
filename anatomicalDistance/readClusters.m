%%

cType = {'dataR_*.mat','dataFx_*.mat','data_*.mat'};
cType = {'dataFx_*.mat'};

cDraw = {'g.','m.','k.'};
cSz = [36 100 25];
nC = length(cType);

d = cell(1,nC);
for ii = 1:nC
    d{ii} = dir(cType{ii});
end

for ff = [3 9 11 22 23]%1:length(d{1})
    F = figure; hold on
    title(['FOV ' num2str(ff)])
    
    for ii = 1:nC
        load(d{ii}(ff).name)
        
        if ii==1
            scatter(allCen.common(:,1),allCen.common(:,2),'k.')
        end
        
        scatter(allCen.unique(:,1),allCen.unique(:,2),cSz(ii),cDraw{ii})
    end
    
    axis('square')
    axis([0 512 0 512])
    axis('off')
    
%     pause
%     close
end

% scatter(unqTempD0(:,1),unqTempD0(:,2),'bo')
% 
% scatter(unqTempD0(inD0,1),unqTempD0(inD0,2),'k.')
% 
% scatter(unqTempD0(inD1,1),unqTempD0(inD1,2),'k.')
% 
% scatter(unqTempD0(inD0&inD1,1),unqTempD0(inD0&inD1,2),'b.')




