%in the "summary.m", the field was actually number of in field bins, here
%let's calculate the real fields
p=pwd;
load('foldersAll.mat');
load('foldersAllIndvFOV.mat');
load('commonCellsAllIdx.mat');

nFieldReal=[];

for m=1:length(foldersAllIndvFOV);
    disp(m)
    N=[];
    for n=1:length(foldersAllIndvFOV{m});
        disp(n)
        nn=[];
         cd(foldersAllIndvFOV{m}{n});
  B=commonCellsAllIdx{m}(:,n);
   load('PValueClassifier_KY2_6_sig\allCellsCorrected.mat');
   for i=1:length(B);
       if ~isnan(allCellsCorrected.fieldCenters{B(i)});
       nn(i,1)=length(allCellsCorrected.fieldCenters{B(i)});
       else
           nn(i,1)=nan;
       end
   end
   N(:,n)=nn;
    end
    
    nFieldReal(end+1:end+size(N,1),:)=N;
end
cd(p);
       save('nFieldReal.mat','nFieldReal');
figure,plot(nanmean(nFieldReal,1))
title('nFields');
saveas(gcf,'nFieldsReal.fig');