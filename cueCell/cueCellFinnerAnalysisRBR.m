p=pwd;
cd ..\
load('foldersAll.mat');
cd(p);
cueRRBR={};%each cell is a day
cueLRBR={};
for n=1:11;
    cueRRBR{n}={};
    cueLRBR{n}={};
end
for n=1:length(foldersAll);
    disp(n)
    cd(foldersAll{n});
        load('reIdentifyCueCellsUniThreshAddingNewMice\indices.mat');
        load('reIdentifyCueCellsUniThreshAddingNewMice\useFolders.mat');
        load('reIdentifyCueCellsUniThreshAddingNewMice\commonCells.mat');
        
        iCueR=indices.iCueR;
        if ~isempty(iCueR);
        i=commonCells(iCueR,:);
               for ii=1:length(useFolders);
                     disp(ii)
                cd(useFolders{ii});
                iii=i(:,ii);
                load('RunByRun_sig\dfofMInterpM_sig.mat');
                cueRRBR{ii}(end+1:end+length(iii))=dfofMInterpM_sig(iii);
               end
        end
        
         iCueL=indices.iCueL;
        if ~isempty(iCueL);
        i=commonCells(iCueL,:);
               for ii=1:length(useFolders);
                     disp(ii)
                cd(useFolders{ii});
                iii=i(:,ii);
                load('RunByRun_sig\dfofMInterpM_sig.mat');
                cueLRBR{ii}(end+1:end+length(iii))=dfofMInterpM_sig(iii);
               end
        end
end
      
cd(p);
save('cueRRBR.mat','cueRRBR');
save('cueLRBR.mat','cueLRBR');
% for n=1:length(cueRRBR);
%     for m=1:length(cueRRBR{n});
%         A=cueRRBR{n}{m};
%         f=floor(size(A,1)/2);
%         cueRRBR{n}{m}=A(1:f*2-1,:)+A(2:f*2,:);
%     end
% end

%% lag and amplitudes
lagRRBR={};
peakRRBR={};
lagbins=60;
load('E:\ID20201118\20210103\1stloc1\pcaica\cueAnalysisNew_sig\tempR.mat');
w=20;%moving window 20 bins
cueRRBRNorm={};

for n=1:length(cueRRBR);
        lagRRBR{n}={};
    peakRRBR{n}={};
    cueRRBRNorm{n}={};
    for m=1:length(cueRRBR{n});
        A=cueRRBR{n}{m};
        peakRRBR{n}{m}=[];
        lagRRBR{n}{m}=[];
          cueRRBRNorm{n}{m}=[];
        for i=1:size(A,1);
            a=A(i,:);
           [~,lagRRBR{n}{m}(i,1),peakBins] = get_cuescoreYGModified_lessLagsNew_ExportPeakBins(a,tempR,lagbins); 
           %normalize AA;
           
            d=[];
        for iii=1:200-w+1;
            d(iii)=std(a(iii:iii+w-1));
        end
           [~,ii]=min(d);
            md=mean(a(ii:ii+w-1));
        a=a-md;
        a=a/max(a);
        cueRRBRNorm{n}{m}(i,:)=a;
        
        for ii=1:length(peakBins);
            p=a(peakBins{ii});
            if ~isempty(p);
            peakRRBR{n}{m}(i,ii)=max(p);
            else
             peakRRBR{n}{m}(i,ii)=nan;
            end
        end
        end
    end
end

save('cueRRBRNorm.mat','cueRRBRNorm');
save('lagRRBR.mat','lagRRBR');
save('peakRRBR.mat','peakRRBR');

%%

lagLRBR={};
peakLRBR={};
lagbins=60;
load('E:\ID20201118\20210103\1stloc1\pcaica\cueAnalysisNew_sig\tempL.mat');
w=20;%moving window 20 bins
cueLRBRNorm={};

for n=1:length(cueLRBR);
        lagLRBR{n}={};
    peakLRBR{n}={};
    cueLRBRNorm{n}={};
    for m=1:length(cueLRBR{n});
        A=cueLRBR{n}{m};
        peakLRBR{n}{m}=[];
        lagLRBR{n}{m}=[];
          cueLRBRNorm{n}{m}=[];
        for i=1:size(A,1);
            a=A(i,:);
           [~,lagLRBR{n}{m}(i,1),peakBins] = get_cuescoreYGModified_lessLagsNew_ExportPeakBins(a,tempL,lagbins); 
           %normalize AA;
           
            d=[];
        for iii=1:200-w+1;
            d(iii)=std(a(iii:iii+w-1));
        end
           [~,ii]=min(d);
            md=mean(a(ii:ii+w-1));
        a=a-md;
        a=a/max(a);
        cueLRBRNorm{n}{m}(i,:)=a;
        
        for ii=1:length(peakBins);
            p=a(peakBins{ii});
            if ~isempty(p);
            peakLRBR{n}{m}(i,ii)=max(p);
            else
             peakLRBR{n}{m}(i,ii)=nan;
            end
        end
        end
    end
end

save('cueLRBRNorm.mat','cueLRBRNorm');
save('lagLRBR.mat','lagLRBR');
save('peakLRBR.mat','peakLRBR');
%%
lagRRBRSTD=[];

for n=1:length(lagRRBR);
    for m=1:length(lagRRBR{n});
       lagRRBRSTD(m,n)=nanstd(lagRRBR{n}{m});
    end
end
save('lagRRBRSTD.mat','lagRRBRSTD');
lagLRBRSTD=[];

for n=1:length(lagLRBR);
    for m=1:length(lagLRBR{n});
       lagLRBRSTD(m,n)=nanstd(lagLRBR{n}{m});
    end
end
save('lagLRBRSTD.mat','lagLRBRSTD');

lagRLRBRSTD=[lagRRBRSTD;lagLRBRSTD];
save('lagRLRBRSTD.mat','lagRLRBRSTD');


peakRRBRAllCues={};

for n=1:length(peakRRBR)
    A=peakRRBR{n};
    B={};
    for m=1:length(A);
        B{m}=nanstd(A{m},1);
    end
    peakRRBRAllCues{n}=cell2mat(B)';
end
peakRRBRAllCues=cell2mat(peakRRBRAllCues);
save('peakRRBRAllCues.mat','peakRRBRAllCues');

peakLRBRAllCues={};

for n=1:length(peakLRBR)
    A=peakLRBR{n};
    B={};
    for m=1:length(A);
        B{m}=nanstd(A{m},1);
    end
    peakLRBRAllCues{n}=cell2mat(B)';
end
peakLRBRAllCues=cell2mat(peakLRBRAllCues);
save('peakLRBRAllCues.mat','peakLRBRAllCues');

peakRLRBRAllCues=[peakRRBRAllCues;peakLRBRAllCues];
save('peakRLRBRAllCues.mat','peakRLRBRAllCues');


figure,
subplot(121)%remove old env
errorbar([1:1:10],nanmean(lagRLRBRSTD(:,2:end),1),nansem(lagRLRBRSTD(:,2:end),1))
 [r,p]=corr([1:1:10]',nanmean(lagRLRBRSTD(:,2:end),1)','tail','right');        
title(['lag std p=',num2str(p)]);       
subplot(122)%remove old env
errorbar([1:1:10],nanmean(peakRLRBRAllCues(:,2:end),1),nansem(peakRLRBRAllCues(:,2:end),1))
[r,p]=corr([1:1:10]',nanmean(peakRLRBRAllCues(:,2:end),1)','tail','right'); 
title(['peak std p=',num2str(p)]);       

saveas(gcf,'lagPeakVariationRBR.fig');

%% without normalizing run by run

%% lag and amplitudes
lagRRBRNonorm={};
peakRRBRNonorm={};
lagbins=60;
load('E:\ID20201118\20210103\1stloc1\pcaica\cueAnalysisNew_sig\tempR.mat');
w=20;%moving window 20 bins

for n=1:length(cueRRBR);
        lagRRBRNonorm{n}={};
    peakRRBRNonorm{n}={};
   
    for m=1:length(cueRRBR{n});
        A=cueRRBR{n}{m};
        peakRRBRNonorm{n}{m}=[];
        lagRRBRNonorm{n}{m}=[];
         
        for i=1:size(A,1);
            a=A(i,:);
           [~,lagRRBRNonorm{n}{m}(i,1),peakBins] = get_cuescoreYGModified_lessLagsNew_ExportPeakBins(a,tempR,lagbins); 
           %normalize AA;
           
             d=[];
        for iii=1:200-w+1;
            d(iii)=std(a(iii:iii+w-1));
        end
           [~,ii]=min(d);
            md=mean(a(ii:ii+w-1));
        a=a-md;
        
        for ii=1:length(peakBins);
            p=a(peakBins{ii});
            if ~isempty(p);
            peakRRBRNonorm{n}{m}(i,ii)=max(p);
            else
             peakRRBRNonorm{n}{m}(i,ii)=nan;
            end
        end
        end
    end
end

save('lagRRBRNonorm.mat','lagRRBRNonorm');
save('peakRRBRNonorm.mat','peakRRBRNonorm');

%%

lagLRBRNonorm={};
peakLRBRNonorm={};
lagbins=60;
load('E:\ID20201118\20210103\1stloc1\pcaica\cueAnalysisNew_sig\tempL.mat');
w=20;%moving window 20 bins

for n=1:length(cueLRBR);
        lagLRBRNonorm{n}={};
    peakLRBRNonorm{n}={};
   
    for m=1:length(cueLRBR{n});
        A=cueLRBR{n}{m};
        peakLRBRNonorm{n}{m}=[];
        lagLRBRNonorm{n}{m}=[];
       
        for i=1:size(A,1);
            a=A(i,:);
           [~,lagLRBRNonorm{n}{m}(i,1),peakBins] = get_cuescoreYGModified_lessLagsNew_ExportPeakBins(a,tempL,lagbins); 
           %normalize AA;
       
           
            d=[];
        for iii=1:200-w+1;
            d(iii)=std(a(iii:iii+w-1));
        end
           [~,ii]=min(d);
            md=mean(a(ii:ii+w-1));
        a=a-md;
      
           
        for ii=1:length(peakBins);
            p=a(peakBins{ii});
            if ~isempty(p);
            peakLRBRNonorm{n}{m}(i,ii)=max(p);
            else
             peakLRBRNonorm{n}{m}(i,ii)=nan;
            end
        end
        end
    end
end

save('lagLRBRNonorm.mat','lagLRBRNonorm');
save('peakLRBRNonorm.mat','peakLRBRNonorm');
%%
lagRRBRSTDNonorm=[];

for n=1:length(lagRRBRNonorm);
    for m=1:length(lagRRBRNonorm{n});
       lagRRBRSTDNonorm(m,n)=nanstd(lagRRBRNonorm{n}{m});
    end
end
save('lagRRBRSTDNonorm.mat','lagRRBRSTDNonorm');
lagLRBRSTDNonorm=[];

for n=1:length(lagLRBRNonorm);
    for m=1:length(lagLRBRNonorm{n});
       lagLRBRSTDNonorm(m,n)=nanstd(lagLRBRNonorm{n}{m});
    end
end
save('lagLRBRSTDNonorm.mat','lagLRBRSTDNonorm');

lagRLRBRSTDNonorm=[lagRRBRSTDNonorm;lagLRBRSTDNonorm];
save('lagRLRBRSTDNonorm.mat','lagRLRBRSTDNonorm');


peakRRBRAllCuesNonorm={};

for n=1:length(peakRRBRNonorm)
    A=peakRRBRNonorm{n};
    B={};
    for m=1:length(A);
        B{m}=nanstd(A{m},1);
    end
    peakRRBRAllCuesNonorm{n}=cell2mat(B)';
end
peakRRBRAllCuesNonorm=cell2mat(peakRRBRAllCuesNonorm);
save('peakRRBRAllCuesNonorm.mat','peakRRBRAllCuesNonorm');

peakLRBRAllCuesNonorm={};

for n=1:length(peakLRBRNonorm)
    A=peakLRBRNonorm{n};
    B={};
    for m=1:length(A);
        B{m}=nanstd(A{m},1);
    end
    peakLRBRAllCuesNonorm{n}=cell2mat(B)';
end
peakLRBRAllCuesNonorm=cell2mat(peakLRBRAllCuesNonorm);
save('peakLRBRAllCuesNonorm.mat','peakLRBRAllCuesNonorm');

peakRLRBRAllCuesNonorm=[peakRRBRAllCuesNonorm;peakLRBRAllCuesNonorm];
save('peakRLRBRAllCuesNonorm.mat','peakRLRBRAllCuesNonorm');

%calculate correlation coefficient per cue using mean of cue
peakRRBRAllCuesNonormCV={};

for n=1:length(peakRRBRNonorm)
    A=peakRRBRNonorm{n};
    B={};
    for m=1:length(A);
        B{m}=nanstd(A{m},1)./nanmean(A{m},1);
    end
    peakRRBRAllCuesNonormCV{n}=cell2mat(B)';
end
peakRRBRAllCuesNonormCV=cell2mat(peakRRBRAllCuesNonormCV);
save('peakRRBRAllCuesNonormCV.mat','peakRRBRAllCuesNonormCV');

peakLRBRAllCuesNonormCV={};

for n=1:length(peakLRBRNonorm)
    A=peakLRBRNonorm{n};
    B={};
    for m=1:length(A);
        B{m}=nanstd(A{m},1)./nanmean(A{m},1);
    end
    peakLRBRAllCuesNonormCV{n}=cell2mat(B)';
end
peakLRBRAllCuesNonormCV=cell2mat(peakLRBRAllCuesNonormCV);
save('peakLRBRAllCuesNonormCV.mat','peakLRBRAllCuesNonormCV');

peakRLRBRAllCuesNonormCV=[peakRRBRAllCuesNonormCV;peakLRBRAllCuesNonormCV];
save('peakRLRBRAllCuesNonormCV.mat','peakRLRBRAllCuesNonormCV');


%% mean ddfof based on dfofaveragesmooh

p=pwd;
cueR={};%each cell is a day
cueL={};
for n=1:11;
    cueR{n}=[];
    cueL{n}=[];
end
for n=1:length(foldersAll);
    disp(n)
    cd(foldersAll{n});
        load('reIdentifyCueCellsUniThreshAddingNewMice\indices.mat');
        load('reIdentifyCueCellsUniThreshAddingNewMice\useFolders.mat');
        load('reIdentifyCueCellsUniThreshAddingNewMice\commonCells.mat');
        
        iCueR=indices.iCueR;
        if ~isempty(iCueR);
        i=commonCells(iCueR,:);
               for ii=1:length(useFolders);
                     disp(ii)
                cd(useFolders{ii});
                iii=i(:,ii);
                d=dir('dfofaveragesmooth_*');
                for k=1:length(d);
                    load(d(k).name);
                end
                cueR{ii}(end+1:end+length(iii),:)=dfofaveragesmooth_sig(:,iii)';
               end
        end
        
         iCueL=indices.iCueL;
        if ~isempty(iCueL);
        i=commonCells(iCueL,:);
               for ii=1:length(useFolders);
                     disp(ii)
                cd(useFolders{ii});
                iii=i(:,ii);
                 d=dir('dfofaveragesmooth_*');
                for k=1:length(d);
                    load(d(k).name);
                end
                cueL{ii}(end+1:end+length(iii),:)=dfofaveragesmooth_sig(:,iii)';
               end
        end
end
cd(p)
cueRMeanF=[];
cueLMeanF=[];

for n=1:length(cueR);
    cueRMeanF(:,n)=mean(cueR{n},2);
     cueLMeanF(:,n)=mean(cueL{n},2);
end

cueRLMeanF=[cueRMeanF;cueLMeanF];

peakRLRBRAllCellsNonorm=[];
for n=1:size(cueRLMeanF,1);
    peakRLRBRAllCellsNonorm(n,:)=mean(peakRLRBRAllCuesNonorm(4*(n-1)+1:4*n,:),1);
end
cueRLCV=peakRLRBRAllCellsNonorm./cueRLMeanF;
save('cueRLCV.mat','cueRLCV');


%% %mean dfof based on run by run
cueLMeanFRBR=[];
for n=1:length(cueLRBR);
    for m=1:length(cueLRBR{n});
        cueLMeanFRBR(m,n)=mean(mean(cueLRBR{n}{m}));
    end
end

cueRMeanFRBR=[];
for n=1:length(cueRRBR);
    for m=1:length(cueRRBR{n});
        cueRMeanFRBR(m,n)=mean(mean(cueRRBR{n}{m}));
    end
end

cueRLMeanFRBR=[cueRMeanFRBR;cueLMeanFRBR];
%caculate coefficient of variation

cueRLCVRBR=peakRLRBRAllCellsNonorm./cueRLMeanFRBR;

%%
figure,
subplot(421)%remove old env
errorbar([1:1:10],nanmean(lagRLRBRSTDNonorm(:,2:end),1),nansem(lagRLRBRSTDNonorm(:,2:end),1))
 [r,p]=corr([1:1:10]',nanmean(lagRLRBRSTDNonorm(:,2:end),1)','tail','left');        
title(['lag std p=',num2str(p)]);       
subplot(422)%remove old env
errorbar([1:1:10],nanmean(peakRLRBRAllCuesNonorm(:,2:end),1),nansem(peakRLRBRAllCuesNonorm(:,2:end),1))
[r,p]=corr([1:1:10]',nanmean(peakRLRBRAllCuesNonorm(:,2:end),1)','tail','left');
title(['peak std p=',num2str(p)]);       

sig=[];
for n=1:9;
    sig(1,n)=ttest(peakRLRBRAllCuesNonorm(:,2),peakRLRBRAllCuesNonorm(:,n+2));
end

subplot(424)%remove old env
errorbar([1:1:10],nanmean(peakRLRBRAllCuesNonormCV(:,2:end),1),nansem(peakRLRBRAllCuesNonormCV(:,2:end),1))
[r,p]=corr([1:1:10]',nanmean(peakRLRBRAllCuesNonormCV(:,2:end),1)','tail','left');
title(['peak CV p=',num2str(p)]);       

sig=[];
for n=1:9;
    sig(1,n)=ttest(peakRLRBRAllCuesNonorm(:,2),peakRLRBRAllCuesNonorm(:,n+2));
end

subplot(425);
errorbar([1:1:10],nanmean(cueRLMeanF(:,2:end),1),nansem(cueRLMeanF(:,2:end),1))
[r,p]=corr([1:1:10]',nanmean(cueRLMeanF(:,2:end),1)','tail','left');

title(['MEAN dfof p=',num2str(p)]);       


subplot(426);
errorbar([1:1:10],nanmean(cueRLMeanFRBR(:,2:end),1),nansem(cueRLMeanFRBR(:,2:end),1))
[r,p]=corr([1:1:10]',nanmean(cueRLMeanFRBR(:,2:end),1)','tail','left');
title(['MEAN dfof based on rur p=',num2str(p)]);     

subplot(427)%remove old env
errorbar([1:1:10],nanmean(cueRLCV(:,2:end),1),nansem(cueRLCV(:,2:end),1))
[r,p]=corr([1:1:10]',nanmean(cueRLCV(:,2:end),1)','tail','left');
title(['CV meanF p=',num2str(p)]);       

subplot(428)%remove old env
errorbar([1:1:10],nanmean(cueRLCVRBR(:,2:end),1),nansem(cueRLCVRBR(:,2:end),1))
[r,p]=corr([1:1:10]',nanmean(cueRLCVRBR(:,2:end),1)','tail','left');
title(['CV meanFRBR p=',num2str(p)]);  

saveas(gcf,'lagPeakVariationRBRNonorm.fig');
%%

%run by run mean activity variation with time: no trend

A=[];

for n=1:11;
    A{n}=[cueRRBR{n} cueLRBR{n}];
end

B=[];

for n=1:11;
    for m=1:23;
        C=A{n}{m};
        mc=mean(C,2);
        B(m,n)=std(mc)/mean(mc);
%   B(m,n)=std(mc);
    end
end

figure,errorbar([1:1:10],mean(B(:,2:end),1),nansem(B(:,2:end),1))

%%
%run by run mean activity variation and activity correlation with
%individual runs: : no trend

B=[];
for n=1:23;
    B{n}=[];
    for m=2:11;
        C=A{m}{n};
        B{n}(end+1:end+size(C,1),:)=C;
    end
end
    

C=nan(164,23);
D=nan(164,23);
for n=1:23;
    AA=B{n};
    a=mean(AA,2);
    C(1:length(a)-1,n)=diff(a);
    for m=1:length(a)-1;
        D(m,n)=corr(AA(m,:)',AA(m+1,:)');
    end
end

figure,
subplot(121)
errorbar([1:1:164],nanmean(C,2),nansem(C,2))
title('mean f diff');
subplot(122)
errorbar([1:1:164],nanmean(D,2),nansem(D,2))
title('correlation to next');

