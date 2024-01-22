load('foldersAllIndvFOV.mat');
fieldD={};%day by day distribution
fieldDA={};%day by day distribution and amplitude
p=pwd;

  for n=1:11;    
      fieldD{n}=[];
    fieldDA{n}=[];
      for m=1:length(foldersAllIndvFOV);
      
        cd(foldersAllIndvFOV{m}{n});
        load('PValueClassifier_KY2_6_sig\allCells.mat');
        da=allCells.dfofaveragesmoothFields';
        d=da;
        d(d>0)=1;
        fieldD{n}(end+1:end+size(d,1),:)=d;
        fieldDA{n}(end+1:end+size(da,1),:)=da;
    end
end
cd(p);   
%%
figure
load('E:\ID20201118\20201224\loc1\pcaica\cueAnalysisNew_sig\tempRL.mat');
oldTempRL=tempRL;
load('E:\ID20201118\20210103\1stloc1\pcaica\cueAnalysisNew_sig\tempRL');

for n=1:length(fieldD);
    subplot(4, 6, n);
    imagesc(fieldD{n});
    subplot(4, 6, 12+n);
    A=sum(fieldD{n},1)/size(fieldD{n},1);
    plot([1:1:200],A);
    hold on
    if n==1;
    plot([1:1:200],oldTempRL*max(A));
        hold on
    line([10 10],[0 max(A)],'Color','green');
    hold on
    line([110 110],[0 max(A)],'Color','green');
    else 
       plot([1:1:200],tempRL*max(A));  
    hold on
    line([178 178],[0 max(A)],'Color','green');
     end
    xlim([1 200]);
    ylim([min(A) max(A)*1.05]);
end
tightfig;
saveas(gcf,'fieldDistri_allCells.fig');

    
figure
for n=1:length(fieldDA);
      subplot(4, 6, n);
    imagesc(fieldDA{n});
    subplot(4, 6, 12+n);
    A=sum(fieldDA{n},1)/size(fieldDA{n},1);
    plot([1:1:200],A);
    hold on
 if n==1;
    plot([1:1:200],oldTempRL*max(A));
        hold on
    line([10 10],[0 max(A)],'Color','green');
    hold on
    line([110 110],[0 max(A)],'Color','green');
    else 
       plot([1:1:200],tempRL*max(A));  
    hold on
    line([178 178],[0 max(A)],'Color','green');
     end
      xlim([1 200]);
    ylim([min(A) max(A)*1.05]);
end   
tightfig;
saveas(gcf,'fieldAmp_allCells.fig');

save('fieldAllCells.mat','fieldD','fieldDA');

%% plot field distri and amp of common cells
load('fieldDistri_zeroOne.mat');
load('fieldDistriAmp.mat');

for n=1:length(fieldDistri_zeroOne);
    subplot(4, 6, n);
    imagesc(fieldDistri_zeroOne{n});
    subplot(4, 6, 12+n);
    A=sum(fieldDistri_zeroOne{n},1)/size(fieldDistri_zeroOne{n},1);
    plot([1:1:200],A);
    hold on
    if n==1;
    plot([1:1:200],oldTempRL*max(A));
        hold on
    line([10 10],[0 max(A)],'Color','green');
    hold on
    line([110 110],[0 max(A)],'Color','green');
    else 
       plot([1:1:200],tempRL*max(A));  
    hold on
    line([178 178],[0 max(A)],'Color','green');
     end
    xlim([1 200]);
    ylim([min(A) max(A)*1.05]);
end
tightfig;

saveas(gcf,'fieldDistri_commonCells.fig');

    
figure
for n=1:length(fieldDistriAmp);
      subplot(4, 6, n);
    imagesc(fieldDistriAmp{n});
    subplot(4, 6, 12+n);
    A=sum(fieldDistriAmp{n},1)/size(fieldDistriAmp{n},1);
    plot([1:1:200],A);
    hold on
 if n==1;
    plot([1:1:200],oldTempRL*max(A));
        hold on
    line([10 10],[0 max(A)],'Color','green');
    hold on
    line([110 110],[0 max(A)],'Color','green');
    else 
       plot([1:1:200],tempRL*max(A));  
    hold on
    line([178 178],[0 max(A)],'Color','green');
     end
      xlim([1 200]);
    ylim([min(A) max(A)*1.05]);
end   
tightfig;
saveas(gcf,'fieldAmp_commonCells.fig');

