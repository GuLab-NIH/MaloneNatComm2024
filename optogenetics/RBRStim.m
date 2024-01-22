%make the excel data (20230724-RBR-new) into matlab file: RBR=[];
A=RBR(:,31:48);
A(A==0)=nan;
RBR(:,31:48)=A;
save('RBR.mat','RBR')

m=nanmean(RBR,1);
s=nansem(RBR,1);
figure
errorbar([1:1:43],m(1:43),s(1:43),'k')
hold on
plot([11:1:20],m(11:20),'b.','MarkerSize',15);
saveas(gcf,'RBR.fig');

%% same thing for CS (20230724-RBR-3x5cm-2ms-ChR2-G-CS): RBRCS

A=RBRCS(:,31:52);
A(A==0)=nan;
RBRCS(:,31:52)=A;
save('RBRCS.mat','RBRCS')

m=nanmean(RBRCS,1);
s=nansem(RBRCS,1);
figure
errorbar([1:1:43],m(1:43),s(1:43),'k')
hold on
plot([11:1:20],m(11:20),'b.','MarkerSize',15);
saveas(gcf,'RBRCS.fig');
%% PLOT TOGETHER
m=nanmean(RBRCS,1);
s=nansem(RBRCS,1);
figure
subplot(121)
errorbar([1:1:43],m(1:43),s(1:43),'k')
hold on
plot([11:1:20],m(11:20),'b.','MarkerSize',15);
ylim([50 95])
m=nanmean(RBR,1);
s=nansem(RBR,1);
subplot(122)
errorbar([1:1:43],m(1:43),s(1:43),'k')
hold on
plot([11:1:20],m(11:20),'b.','MarkerSize',15);
ylim([50 95])

saveas(gcf,'RBRRSCS.fig')

%% plot control CHR2
%copy 20230726-RBR-3x5cm-2ms-GFP-G-CS to RBRCSGFP=[];
%copy 20230726-RBR-3x5cm-2ms-GFP-G-RS to RBRRSGFP=[];
A=RBRCSGFP(:,31:46);
A(A==0)=nan;
RBRCSGFP(:,31:46)=A;
save('RBRCSGFP.mat','RBRCSGFP')

A=RBRRSGFP(:,31:43);
A(A==0)=nan;
RBRRSGFP(:,31:43)=A;
save('RBRRSGFP.mat','RBRRSGFP')

%% 
load('RBRCS.mat');
load('RBR.mat');
load('RBRCSGFP.mat')
load('RBRRSGFP.mat')

m=nanmean(RBRCS,1);
s=nansem(RBRCS,1);

N=40;

figure
subplot(221)
errorbar([1:1:N],m(1:N),s(1:N),'k')
hold on
plot([11:1:20],m(11:20),'b.','MarkerSize',15);
ylim([50 95])
title('chr2 CS')
xlim([0 N])

m=nanmean(RBR,1);
s=nansem(RBR,1);
subplot(222)
errorbar([1:1:N],m(1:N),s(1:N),'k')
hold on
plot([11:1:20],m(11:20),'b.','MarkerSize',15);
ylim([50 95])
title('chr2 RS')
xlim([0 N])

N=40;
m=nanmean(RBRCSGFP,1);
s=nansem(RBRCSGFP,1);
subplot(223)
errorbar([1:1:N],m(1:N),s(1:N),'k')
hold on
plot([11:1:20],m(11:20),'b.','MarkerSize',15);
ylim([50 95])
title('gfp CS')
xlim([0 N])

m=nanmean(RBRRSGFP,1);
s=nansem(RBRRSGFP,1);
subplot(224)
errorbar([1:1:N],m(1:N),s(1:N),'k')
hold on
plot([11:1:20],m(11:20),'b.','MarkerSize',15);
ylim([50 95])
title('gfp RS')
xlim([0 N])

saveas(gcf,'RBRRSCSCHR2GFP.fig')

%
%% 
load('RBRCS.mat');
load('RBR.mat');
load('RBRCSGFP.mat')
load('RBRRSGFP.mat')

%PLOT TOGETHER

N=40;

figure
subplot(221)
A=RBR(:,1:N);
for n=1:size(A,1);
    ma=mean(A(n,[1:10]));
    A(n,:)=A(n,:)/ma;
end

B=RBRCS(:,1:N);
for n=1:size(B,1);
    ma=mean(B(n,[1:10]));
    B(n,:)=B(n,:)/ma;
end
semshade(B,0.2,'r',[1:40])
hold on
semshade(A,0.2,'b',[1:40])
ylim([0.6 1.3])
title('chr2 CS RS')
xlim([0 N])
p1=[];
for n=1:N;
    [r,p1(n)]= ttest2(A(:,n),B(:,n));
end

for n=1:N;
    if p1(n)<=0.05;
        hold on
        plot(n,1.2,'k*')
    end
end

subplot(222)
A=RBRRSGFP(:,1:N);
for n=1:size(A,1);
    ma=mean(A(n,[1:10]));
    A(n,:)=A(n,:)/ma;
end

B=RBRCSGFP(:,1:N);
for n=1:size(B,1);
    ma=mean(B(n,[1:10]));
    B(n,:)=B(n,:)/ma;
end
semshade(B,0.2,'r',[1:40])
hold on
semshade(A,0.2,'b',[1:40])
ylim([0.6 1.3])
title('gfp CS RS')
xlim([0 N])
p2=[];
for n=1:N;
    [r,p2(n)]= ttest2(A(:,n),B(:,n));
end
for n=1:N;
    if p2(n)<=0.05;
        hold on
        plot(n,1.2,'k*')
    end
end



subplot(223)
A=RBRCS(:,1:N);
for n=1:size(A,1);
    ma=mean(A(n,[1:10]));
    A(n,:)=A(n,:)/ma;
end

B=RBRCSGFP(:,1:N);
for n=1:size(B,1);
    ma=mean(B(n,[1:10]));
    B(n,:)=B(n,:)/ma;
end
semshade(B,0.2,'r',[1:40])
hold on
semshade(A,0.2,'b',[1:40])
ylim([0.6 1.3])
title('CS CHR2 GFP')
xlim([0 N])
p3=[];
for n=1:N;
    [r,p3(n)]= ttest2(A(:,n),B(:,n));
end
for n=1:N;
    if p3(n)<=0.05;
        hold on
        plot(n,1.2,'k*')
    end
end



subplot(224)
A=RBR(:,1:N);
for n=1:size(A,1);
    ma=mean(A(n,[1:10]));
    A(n,:)=A(n,:)/ma;
end

B=RBRRSGFP(:,1:N);
for n=1:size(B,1);
    ma=mean(B(n,[1:10]));
    B(n,:)=B(n,:)/ma;
end
semshade(B,0.2,'r',[1:40])
hold on
semshade(A,0.2,'b',[1:40])
ylim([0.6 1.3])
title('RS CHR2 GFP')
xlim([0 N])
p4=[];
for n=1:N;
    [r,p4(n)]= ttest2(A(:,n),B(:,n));
end
for n=1:N;
    if p4(n)<=0.05;
        hold on
        plot(n,1.2,'k*')
    end
end

saveas(gcf,'RBR_together.fig')

