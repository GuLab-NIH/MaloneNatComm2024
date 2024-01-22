%% Set parameters

clear; close all; clc

useT = erase(pwd,'Z:\labMembers\TM\Novelty_Analysis\clusterMan\');
anaT = 'norm2other';
ttl = ['Fx-global_' useT '_' anaT];

useFile = ['distsFOVFx_' useT '_global.mat'];
load(useFile)

outType = {'Common','Other','Cross'};

xSet = 1:38;
xN = length(xSet);

XLIM = [min(xSet)-0.5 max(xSet)+0.5];

normIdx = 2;

nSubsets = 1;
subSz = 9;
pGlobal = zeros(nSubsets,xN);

%%

for kk = 1:nSubsets
    useIdx = 1:size(dists,1);
%     useIdx = randsample(size(dists,1),subSz);

    
    %% Plot data by FOV
    
    if normIdx==0
        normF = ones(1,xN);
        YLIM = [30 120];
    else
        normF = mean(distsCom{normIdx}(useIdx,xSet),1);
        YLIM = [0.7 1.1];
    end
    
    useCom = cell(1,length(outType));
    
    for jj = 1:length(outType)
        useCom{jj} = distsCom{jj}(useIdx,xSet)./normF;
    end
    
    if normIdx==2
        useCom = useCom(1:2);
        outType = outType(1:2);
    end
    
    % plot by FOV
    figure; hold on
    % plotPDISTall(distsCom,outType,[ttl ' by FOV'],1);
    pVals = plotPDISTallFigure(useCom,outType,[ttl ' by FOV'],xSet);
    
    xlim(XLIM)
%     set(gca,'XLabels',xSet)
    ylim(YLIM)
    
    pGlobal(kk,:) = pVals{1};
    
    %% Save figures
    
    if nSubsets==1
        ii = 1;
        while true
            cName = [ttl '_' num2str(ii) '.fig'];
            if isfile(cName)
                ii = ii+1;
            else
%                 savefig(cName)
                break
            end
        end
    else
        close
    end
    
end


%% Test significance percentage

pThresh = 0.05;

pSig = pGlobal*size(pGlobal,2)>pThresh;
% pSig = pGlobal>pThresh;

pSigPer = mean(pSig,1)*100;

