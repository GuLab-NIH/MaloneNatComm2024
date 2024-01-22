% name='dfofaveragesmooth_*';
% fn= dir(name);
% for n=1:length(fn);
% load(fn(n).name);
% end
load('D:\Mice\cueTemp.mat');

cd('RunByRun_dfof');  
load('dfofMInterpM.mat');
load('corrInfo.mat');
clear useRun
useRun=find(corrInfo.noNaN(:,1));
if isempty(useRun)
    useRun=find(corrInfo.noNaN(:,2));
end

if length(useRun)>15;
useRun=useRun(1:15);%only run the first 15
end

N=5;%rolling average every three bins
corrInfoLocation=[];
corrInfoLocation.toOthers={};
corrInfoLocation.toOthersMean=[];
corrInfoLocation.meanAll=[];
corrInfoLocation.semAll=[];

for n=1:length(dfofMInterpM);
    disp(n)
    corrInfoLocation.toOthers{n}=[];
    A=dfofMInterpM{n}(useRun,:);
    for m=1:size(A,2)-(N-1);
        nrow=0;
        for i=1:size(A,1)-1;
            for ii=i+1:size(A,1);
            c=corr(A(i,m:m+N-1)',A(ii,m:m+N-1)');
            nrow=nrow+1;
            corrInfoLocation.toOthers{n}(nrow,m)=c;
            end
        end
    end
     corrInfoLocation.toOthersMean(n,:)=nanmean(corrInfoLocation.toOthers{n},1);
end

corrInfoLocation.meanAll=nanmean(corrInfoLocation.toOthersMean,1);
corrInfoLocation.semAll=nansem(corrInfoLocation.toOthersMean,1);

save('corrInfoLocation.mat','corrInfoLocation');

figure,errorbar([1:1:length(corrInfoLocation.meanAll)],corrInfoLocation.meanAll,corrInfoLocation.semAll,'k')
hold on
plot([1:1:length(cueTemp)],cueTemp*max(corrInfoLocation.meanAll),'m');
% hold on
% line([890/5 890/5],[0 max(corrInfoLocation.meanAll)]);

saveas(gcf,'corrInfoLocation.fig');
close