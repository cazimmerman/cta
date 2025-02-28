function fig_ED10(data_path)

disp('Generating panels for Extended Data Figure 10...')
%% Fig ED10a

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 10a','FontWeight','bold')
load([data_path,'\Neuropixels\CGRP stim CFA\cgrp-stim-retrieval-summary.mat'],'data')

hold on
axis square
A = []; B = [];
for i = unique(data.meta.session)
    A(i) = sum(data.stats.significant & data.stats.preference>0 & data.meta.session==i)./sum(data.meta.session==i)*100;
    B(i) = sum(data.stats_retrieval.significant & data.stats_retrieval.preference>0 & data.meta.session==i)./sum(data.meta.session==i)*100;
end
for i = 1:length(A)
    plot([1 2],[A(i),B(i)],'color',[.85 .85 .85],'linewidth',1)
end
xx = [NaN mean(A)-std(A)./sqrt(length(A)) mean(A) mean(A)+std(A)./sqrt(length(A)) NaN];
plot([-0.125 0.125]+1,[xx(3) xx(3)],'color','k','LineWidth',1)
plot([1 1],[xx(2) xx(4)],'color','k','LineWidth',1)
xx = [NaN mean(B)-std(B)./sqrt(length(B)) mean(B) mean(B)+std(B)./sqrt(length(B)) NaN];
plot([-0.125 0.125]+2,[xx(3) xx(3)],'color','k','LineWidth',1)
plot([1 1]*2,[xx(2) xx(4)],'color','k','LineWidth',1)
ylabel('Flavor-preferring neurons (%)'); ylim([0 60]); yticks(0:10:60);
xlim([.5 2.5]); xticks([1 2]); xticklabels({'Conditioning','Retrieval'});
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

[p,~,stats]=signrank(A,B,'method','exact');
StatsTbl = table({'ED 10a'},{'Conditioning vs. Retrieval'},{'Wilcoxon signed-rank'},{'N/A'},{[length(A)]},stats.signedrank,p, ...
    'VariableNames',{'Figure panel','Group','Statistical test','Multiple comparisons','Sample size','Test statistic','P-value'});
%% Fig ED10b

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 10b','FontWeight','bold')
load([data_path,'\Neuropixels\CGRP stim CFA\cgrp-stim-retrieval-summary.mat'],'data')

cmap = struct;
idx1 = find(data.stats.significant & data.stats.preference>0);
idx2 = find(data.stats.significant & data.stats.preference<0);
idx = sort([idx1,idx2]);
PCAdata.novel = data.psth.smooth.novel(idx,501:end) - mean([data.psth.smooth.novel(idx,501:600),data.psth.smooth.water(idx,501:600)],2);
PCAdata.water = data.psth.smooth.water(idx,501:end) - mean([data.psth.smooth.novel(idx,501:600),data.psth.smooth.water(idx,501:600)],2);
PCAdata.cgrp  = data.psth.smooth.cgrp(idx,101:600)  - mean(data.psth.smooth.cgrp(idx,101:200),2);
PCAdata.novel_retrieval = data.psth.smooth.novel_retrieval(idx,501:end) - mean([data.psth.smooth.novel(idx,501:600),data.psth.smooth.water(idx,501:600)],2);
PCAdata.water_retrieval = data.psth.smooth.water_retrieval(idx,501:end) - mean([data.psth.smooth.novel(idx,501:600),data.psth.smooth.water(idx,501:600)],2);
t=max(abs([PCAdata.novel,PCAdata.water]),[],2);
PCAdata.novel = PCAdata.novel./t;
PCAdata.water = PCAdata.water./t;
PCAdata.novel_retrieval = PCAdata.novel_retrieval./t;
PCAdata.water_retrieval = PCAdata.water_retrieval./t;
t=max(abs([PCAdata.cgrp]),[],2);
PCAdata.cgrp = PCAdata.cgrp./t;
cmap.novel = cbrewer('seq','Reds',size(PCAdata.novel,2)*1.5,'spline'); cmap.novel = cmap.novel(size(PCAdata.novel,2)/2+1:size(PCAdata.novel,2)*1.5,:);
cmap.water = cbrewer('seq','Blues',size(PCAdata.water,2)*1.5,'spline'); cmap.water = cmap.water(size(PCAdata.water,2)/2+1:size(PCAdata.water,2)*1.5,:);
cmap.cgrp = cbrewer('seq','Greens',size(PCAdata.cgrp,2)*1.5,'spline'); cmap.cgrp = cmap.cgrp(size(PCAdata.cgrp,2)/2+1:size(PCAdata.cgrp,2)*1.5,:); cmap.cgrp(cmap.cgrp<0) = 0;

x = [PCAdata.novel(:,501:1000) PCAdata.water(:,501:1000)]; x = x-mean(x);
[coeff,score,latent,tsquared,explained,mu] = pca(x,'Centered',false);

PCA_CFA = struct;
ax1 = subplot(1,2,1);
axis square
hold on
pc1 = []; pc2 = [];
for i = 1:size(PCAdata.water,2)
    x = PCAdata.water(:,i); x = x-mean(x);
    y = x'*score; pc1(i) = y(1); pc2(i) = y(2);
    scatter(pc1(i),pc2(i),100,cmap.water(i,:),'filled')
    if i>1
        plot([pc1(i-1) pc1(i)],[pc2(i-1) pc2(i)],'Color',cmap.water(i,:),'LineWidth',2)
    end
end
pc1 = []; pc2 = [];
for i = 1:size(PCAdata.novel,2)
    x = PCAdata.novel(:,i); x = x-mean(x);
    y = x'*score; pc1(i) = y(1); pc2(i) = y(2);
    scatter(pc1(i),pc2(i),100,cmap.novel(i,:),'filled')
    if i>1
        plot([pc1(i-1) pc1(i)],[pc2(i-1) pc2(i)],'Color',cmap.novel(i,:),'LineWidth',2)
    end
end
PCA_CFA.novel = pc2;
pc1 = []; pc2 = [];
for i = 1:size(PCAdata.cgrp,2)
    x = PCAdata.cgrp(:,i); x = x-mean(x);
    y = x'*score; pc1(i) = y(1); pc2(i) = y(2);
    scatter(pc1(i),pc2(i),100,cmap.cgrp(i,:),'filled')
    if i>1
        plot([pc1(i-1) pc1(i)],[pc2(i-1) pc2(i)],'Color',cmap.cgrp(i,:),'LineWidth',2)
    end
end
x1 = xlim; y1 = ylim;
xticks([]); yticks([]);
xlabel('PC1'); ylabel('PC2');
set(gca,'FontSize',12,'LineWidth',1)
title('CGRP stim conditioning day','FontWeight','normal')
hold off

ax2 = subplot(1,2,2);
axis square
hold on
pc1 = []; pc2 = [];
for i = 1:size(PCAdata.water_retrieval,2)
    x = PCAdata.water_retrieval(:,i); x = x-mean(x);
    y = x'*score; pc1(i) = y(1); pc2(i) = y(2);
    scatter(pc1(i),pc2(i),100,cmap.water(i,:),'filled')
    if i>1
        plot([pc1(i-1) pc1(i)],[pc2(i-1) pc2(i)],'Color',cmap.water(i,:),'LineWidth',2)
    end
end
pc1 = []; pc2 = [];
for i = 1:size(PCAdata.novel_retrieval,2)
    x = PCAdata.novel_retrieval(:,i); x = x-mean(x);
    y = x'*score; pc1(i) = y(1); pc2(i) = y(2);
    scatter(pc1(i),pc2(i),100,cmap.novel(i,:),'filled')
    if i>1
        plot([pc1(i-1) pc1(i)],[pc2(i-1) pc2(i)],'Color',cmap.novel(i,:),'LineWidth',2)
    end
end
PCA_CFA.retrieval = pc2;
y2 = ylim; x2 = xlim;
xticks([]); yticks([]);
xlabel('PC1'); ylabel('PC2');
set(gca,'FontSize',12,'LineWidth',1)
title('Retrieval day','FontWeight','normal')
hold off

x = [min([x1,x2]) max([x1 x2])];
y = [min([y1,y2]) max([y1 y2])];
set(ax1,'xlim',x,'ylim',y)
set(ax2,'xlim',x,'ylim',y)
%% Fig ED10c

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 10c','FontWeight','bold')
load([data_path,'\Neuropixels\CGRP-CEA stim CFA\cgrp-cea-stim-retrieval-summary.mat'],'data')

idx = find(data.stats.significant & data.stats.preference>0);
[~,idx2] = sort(data.stats.response.cgrp_period(idx),'descend'); idx = idx(idx2);
A1 = data.psth.smooth.novel(idx,:);
A2 = data.psth.smooth.water(idx,:);
A4 = data.psth.smooth.novel_retrieval(idx,:);
A5 = data.psth.smooth.water_retrieval(idx,:);
A6 = data.psth.cgrp_period(idx,:);

idx = find(data.stats.significant & data.stats.preference<0);
[~,idx2] = sort(data.stats.response.cgrp_period(idx),'descend'); idx = idx(idx2);
B1 = data.psth.smooth.novel(idx,:);
B2 = data.psth.smooth.water(idx,:);
B4 = data.psth.smooth.novel_retrieval(idx,:);
B5 = data.psth.smooth.water_retrieval(idx,:);
B6 = data.psth.cgrp_period(idx,:);

idx = find(~data.stats.significant);
[~,idx2] = sort(data.stats.response.cgrp_period(idx),'descend'); idx = idx(idx2);
C1 = data.psth.smooth.novel(idx,:);
C2 = data.psth.smooth.water(idx,:);
C4 = data.psth.smooth.novel_retrieval(idx,:);
C5 = data.psth.smooth.water_retrieval(idx,:);
C6 = data.psth.cgrp_period(idx,:);

pl1 = [A1(:,501:2000);nan(15,1500);B1(:,501:2000);nan(15,1500);C1(:,501:2000)];
pl2 = [A2(:,501:2000);nan(15,1500);B2(:,501:2000);nan(15,1500);C2(:,501:2000)];
pl = [pl1,nan(size(pl1,1),75),pl2];
cmap = flipud(cbrewer('div','RdBu',1000,'spline')); cmap(cmap<0) = 0;

subplot(1,3,1)
hold on
heatmap(flipud(pl),[],[],[],'Colormap',cmap,'ColorLevels',1000,'MaxColorValue',.5,'MinColorValue',-.5,'NaNColor',[1 1 1]);
nomod = size(C1,1); nopref = size(B1,1); flavorpref = size(A1,1);
yticks([nomod./2 nomod+nopref./2+15 nomod+nopref+flavorpref./2+30]+0.5)
ytickangle(90)
yticklabels({'Non-selective','Water','Flavor'})
xticks([750 750+1500+75])
xticklabels({'Flavor','Water'})
title('CGRP^{CEA} stim conditioning day')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0, 0],'TickDir','out')
hold off

pl6 = [A6;nan(15,75);B6;nan(15,75);C6];
subplot(1,3,2)
hold on
heatmap(flipud(pl6),[],[],[],'Colormap',cmap,'ColorLevels',1000,'MaxColorValue',.5,'MinColorValue',-.5,'NaNColor',[1 1 1]);
xticks([0:15:75]+.5)
title('CGRP^{CEA} stim conditioning day')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0, 0],'TickDir','out')
hold off

pl1 = [A4(:,501:2000);nan(15,1500);B4(:,501:2000);nan(15,1500);C4(:,501:2000)];
pl2 = [A5(:,501:2000);nan(15,1500);B5(:,501:2000);nan(15,1500);C5(:,501:2000)];
pl = [pl1,nan(size(pl1,1),75),pl2];
subplot(1,3,3)
hold on
heatmap(flipud(pl),[],[],[],'Colormap',cmap,'ColorLevels',1000,'MaxColorValue',.5,'MinColorValue',-.5,'NaNColor',[1 1 1]);
xticks([750 750+1500+75])
xticklabels({'Flavor','Water'})
title('Retrieval day')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0, 0],'TickDir','out')
hold off
%% Fig ED10d

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 10d','FontWeight','bold')
load([data_path,'\Neuropixels\CGRP-CEA stim CFA\cgrp-cea-stim-retrieval-summary.mat'],'data')

hold on
axis square
A = []; B = [];
for i = unique(data.meta.session)
    A(i) = sum(data.stats.significant & data.stats.preference>0 & data.meta.session==i)./sum(data.meta.session==i)*100;
    B(i) = sum(data.stats_retrieval.significant & data.stats_retrieval.preference>0 & data.meta.session==i)./sum(data.meta.session==i)*100;
end
for i = 1:length(A)
    plot([1 2],[A(i),B(i)],'color',[.85 .85 .85],'linewidth',1)
end
xx = [NaN mean(A)-std(A)./sqrt(length(A)) mean(A) mean(A)+std(A)./sqrt(length(A)) NaN];
plot([-0.125 0.125]+1,[xx(3) xx(3)],'color','k','LineWidth',1)
plot([1 1],[xx(2) xx(4)],'color','k','LineWidth',1)
xx = [NaN mean(B)-std(B)./sqrt(length(B)) mean(B) mean(B)+std(B)./sqrt(length(B)) NaN];
plot([-0.125 0.125]+2,[xx(3) xx(3)],'color','k','LineWidth',1)
plot([1 1]*2,[xx(2) xx(4)],'color','k','LineWidth',1)
ylabel('Flavor-preferring neurons (%)'); ylim([0 60]); yticks(0:10:60);
xlim([.5 2.5]); xticks([1 2]); xticklabels({'Conditioning','Retrieval'});
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

[p,~,stats]=signrank(A,B,'method','exact');
StatsTbl(end+1,:) = table({'ED 10d'},{'Conditioning vs. Retrieval'},{'Wilcoxon signed-rank'},{'N/A'},{[length(A)]},stats.signedrank,p);
%% Fig ED10e

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 10e','FontWeight','bold')
load([data_path,'\Neuropixels\CGRP-CEA stim CFA\cgrp-cea-stim-retrieval-summary.mat'],'data')

subplot(1,3,1)
hold on
axis square
times = -5:.01:10;
times = times(1:end-1) + mean(diff(times))/2;
idx = find(data.stats.significant & data.stats.preference>0);
[~,idx2] = sort(data.stats.response.cgrp_period(idx),'descend'); idx = idx(idx2(1:floor(length(idx2)*0.1)));
A1 = data.psth.smooth.novel(idx,501:end);
A2 = data.psth.smooth.novel_retrieval(idx,501:end);
fill([times fliplr(times)],[mean(A1)+std(A1)/sqrt(size(A1,1)) fliplr(mean(A1)-std(A1)/sqrt(size(A1,1)))],[.85 .85 .85],'LineStyle','none')
a = plot(times,mean(A1),'k','LineWidth',1);
fill([times fliplr(times)],[mean(A2)+std(A2)/sqrt(size(A2,1)) fliplr(mean(A2)-std(A2)/sqrt(size(A2,1)))],[252 216 213]/255,'LineStyle','none')
b = plot(times,mean(A2),'color',[229 45 38]/255,'LineWidth',1);
ylabel('Spiking (σ)'); ylim([-.2 .8]); yticks(-.2:.2:.8);
xlabel('Time (s)'); xlim([-5 10]); xticks(-5:5:10)
title('High CGRP^{CEA} response')
legend([a,b],{'CGRP^{CEA} stim conditioning day','Retrieval day'})
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

subplot(1,3,2)
hold on
axis square
times = -5:.01:10;
times = times(1:end-1) + mean(diff(times))/2;
idx = find(data.stats.significant & data.stats.preference>0);
[~,idx2] = sort(data.stats.response.cgrp_period(idx),'descend'); idx = idx(idx2(floor(length(idx2)*0.1)+1:end));
A1 = data.psth.smooth.novel(idx,501:end);
A2 = data.psth.smooth.novel_retrieval(idx,501:end);
fill([times fliplr(times)],[mean(A1)+std(A1)/sqrt(size(A1,1)) fliplr(mean(A1)-std(A1)/sqrt(size(A1,1)))],[.85 .85 .85],'LineStyle','none')
a = plot(times,mean(A1),'k','LineWidth',1);
fill([times fliplr(times)],[mean(A2)+std(A2)/sqrt(size(A2,1)) fliplr(mean(A2)-std(A2)/sqrt(size(A2,1)))],[252 216 213]/255,'LineStyle','none')
b = plot(times,mean(A2),'color',[229 45 38]/255,'LineWidth',1);
ylabel('Spiking (σ)'); ylim([-.2 .8]); yticks(-.2:.2:.8);
xlabel('Time (s)'); xlim([-5 10]); xticks(-5:5:10)
title('Low CGRP^{CEA} response')
legend([a,b],{'CGRP^{CEA} stim conditioning day','Retrieval day'})
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

subplot(1,3,3)
hold on
axis square
for i = 0:.1:2.9
    plot([i i],[.980 .995],'Color',[54 161 86]/255,'LineWidth',2)
end
times = -1:.01:4;
times = times(1:end-1) + mean(diff(times))/2;
idx = find(data.stats.significant & data.stats.preference>0);
[~,idx2] = sort(data.stats.response.cgrp_period(idx),'descend'); idx = idx(idx2(floor(length(idx2)*0.1)+1:end));
A1 = data.psth.raw.cgrp(idx,101:600);
fill([times fliplr(times)],[mean(A1)+std(A1)/sqrt(size(A1,1)) fliplr(mean(A1)-std(A1)/sqrt(size(A1,1)))],[.85 .85 .85],'LineStyle','none')
b = plot(times,mean(A1),'k','LineWidth',1);
idx = find(data.stats.significant & data.stats.preference>0);
[~,idx2] = sort(data.stats.response.cgrp_period(idx),'descend'); idx = idx(idx2(1:floor(length(idx2)*0.1)));
A1 = data.psth.raw.cgrp(idx,101:600);
fill([times fliplr(times)],[mean(A1)+std(A1)/sqrt(size(A1,1)) fliplr(mean(A1)-std(A1)/sqrt(size(A1,1)))],[252 216 213]/255,'LineStyle','none')
a = plot(times,mean(A1),'color',[229 45 38]/255,'LineWidth',1);
ylim([-.2 1])
yticks(-.2:.2:1)
legend([a,b],{'High CGRP^{CEA} response','Low CGRP^{CEA} response'})
xticks(-1:1:4)
xlim([-1 4])
xlabel('Time (s)')
ylabel('Spiking (σ)')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off
%% Fig ED10f

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 10f','FontWeight','bold')
load([data_path,'\Neuropixels\LiCl CFA\licl-control-retrieval-summary.mat'],'data')

subplot(1,2,1)
hold on
axis square
times = -5:.01:10;
times = times(1:end-1) + mean(diff(times))/2;
idx = find(data.stats.significant & data.stats.preference>0);
A1 = data.psth.smooth.novel(idx,501:end);
A2 = data.psth.smooth.novel_retrieval(idx,501:end);
fill([times fliplr(times)],[mean(A1)+std(A1)/sqrt(size(A1,1)) fliplr(mean(A1)-std(A1)/sqrt(size(A1,1)))],[.85 .85 .85],'LineStyle','none')
a = plot(times,mean(A1),'k','LineWidth',1);
fill([times fliplr(times)],[mean(A2)+std(A2)/sqrt(size(A2,1)) fliplr(mean(A2)-std(A2)/sqrt(size(A2,1)))],[252 216 213]/255,'LineStyle','none')
b = plot(times,mean(A2),'color',[229 45 38]/255,'LineWidth',1);
ylabel('Spiking (σ)'); ylim([-.1 .3]); yticks(-.1:.1:.3);
xlabel('Time (s)'); xlim([-5 10]); xticks(-5:5:10)
legend([a,b],{'LiCl conditioning day','Retrieval day'})
title(['Control mice',char(10),'All Flavor-preferring'])
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

subplot(1,2,2)
axis square
hold on
idx = find(data.stats.significant & data.stats.preference>0);
A = data.stats.response.novel(idx);
B = data.stats_retrieval.response.novel(idx);
simpleboxplot(1,A,'k')
simpleboxplot(2,B,[229 45 38]/255)
ylabel('Flavor response (σ)'); ylim([-.1 .3]); yticks(ylim)
xticks([1,2])
xticklabels({'LiCl conditioning day','Retrieval day'})
xlim([0.25 2.75])
title(['Control mice',char(10),'All Flavor-preferring'])
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

[p,~,stat] = signrank(A,B,'method','exact');
StatsTbl(end+1,:) = table({'ED 10f'},{'Conditioning vs. Retrieval'},{'Wilcoxon signed-rank'},{'N/A'},{[length(A)]},stat.signedrank,p);
%% Fig ED10g

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 10g','FontWeight','bold')
load([data_path,'\Neuropixels\LiCl CFA\licl-tacasp3-retrieval-summary.mat'],'data')

subplot(1,2,1)
hold on
axis square
times = -5:.01:10;
times = times(1:end-1) + mean(diff(times))/2;
idx = find(data.stats.significant & data.stats.preference>0);
A1 = data.psth.smooth.novel(idx,501:end);
A2 = data.psth.smooth.novel_retrieval(idx,501:end);
fill([times fliplr(times)],[mean(A1)+std(A1)/sqrt(size(A1,1)) fliplr(mean(A1)-std(A1)/sqrt(size(A1,1)))],[.85 .85 .85],'LineStyle','none')
a = plot(times,mean(A1),'k','LineWidth',1);
fill([times fliplr(times)],[mean(A2)+std(A2)/sqrt(size(A2,1)) fliplr(mean(A2)-std(A2)/sqrt(size(A2,1)))],[252 216 213]/255,'LineStyle','none')
b = plot(times,mean(A2),'color',[229 45 38]/255,'LineWidth',1);
ylabel('Spiking (σ)'); ylim([-.1 .3]); yticks(-.1:.1:.3);
xlabel('Time (s)'); xlim([-5 10]); xticks(-5:5:10)
legend([a,b],{'LiCl conditioning day','Retrieval day'})
title(['CGRP neuron ablation',char(10),'All Flavor-preferring'])
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

subplot(1,2,2)
axis square
hold on
idx = find(data.stats.significant & data.stats.preference>0);
A = data.stats.response.novel(idx);
B = data.stats_retrieval.response.novel(idx);
simpleboxplot(1,A,'k')
simpleboxplot(2,B,[229 45 38]/255)
ylabel('Flavor response (σ)'); ylim([-.1 .3]); yticks(ylim)
xticks([1,2])
xticklabels({'LiCl conditioning day','Retrieval day'})
xlim([0.25 2.75])
title(['CGRP neuron ablation',char(10),'All Flavor-preferring'])
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

[p,~,stat] = signrank(A,B,'method','exact');
StatsTbl(end+1,:) = table({'ED 10g'},{'Conditioning vs. Retrieval'},{'Wilcoxon signed-rank'},{'N/A'},{[length(A)]},stat.signedrank,p);
%% Fig ED10h

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 10h','FontWeight','bold')
load([data_path,'\Neuropixels\Familiarization\familiarization-summary.mat'],'data')

idx = find(data.stats.significant & data.stats.preference>0);
[~,idx2] = sort(data.stats.response.novel(idx),'descend'); idx = idx(idx2);
A1 = data.psth.smooth.novel(idx,:);
A2 = data.psth.smooth.water(idx,:);
A4 = data.psth.smooth.novel_retrieval(idx,:);
A5 = data.psth.smooth.water_retrieval(idx,:);

idx = find(data.stats.significant & data.stats.preference<0);
[~,idx2] = sort(data.stats.response.water(idx),'descend'); idx = idx(idx2);
B1 = data.psth.smooth.novel(idx,:);
B2 = data.psth.smooth.water(idx,:);
B4 = data.psth.smooth.novel_retrieval(idx,:);
B5 = data.psth.smooth.water_retrieval(idx,:);

idx = find(~data.stats.significant);
[~,idx2] = sort(mean([data.stats.response.novel(idx);data.stats.response.water(idx)]),'descend'); idx = idx(idx2);
C1 = data.psth.smooth.novel(idx,:);
C2 = data.psth.smooth.water(idx,:);
C4 = data.psth.smooth.novel_retrieval(idx,:);
C5 = data.psth.smooth.water_retrieval(idx,:);

pl1 = [A1(:,501:2000);nan(15,1500);B1(:,501:2000);nan(15,1500);C1(:,501:2000)];
pl2 = [A2(:,501:2000);nan(15,1500);B2(:,501:2000);nan(15,1500);C2(:,501:2000)];
pl = [pl1,nan(size(pl1,1),75),pl2];
cmap = flipud(cbrewer('div','RdBu',1000,'spline')); cmap(cmap<0) = 0;

subplot(1,2,1)
hold on
heatmap(flipud(pl),[],[],[],'Colormap',cmap,'ColorLevels',1000,'MaxColorValue',.5,'MinColorValue',-.5,'NaNColor',[1 1 1]);
nomod = size(C1,1); nopref = size(B1,1); flavorpref = size(A1,1);
yticks([nomod./2 nomod+nopref./2+15 nomod+nopref+flavorpref./2+30]+0.5)
ytickangle(90)
yticklabels({'Non-selective','Water','Flavor'})
xticks([750 750+1500+75])
xticklabels({'Flavor','Water'})
title('Novel day')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0, 0],'TickDir','out')
hold off

pl1 = [A4(:,501:2000);nan(15,1500);B4(:,501:2000);nan(15,1500);C4(:,501:2000)];
pl2 = [A5(:,501:2000);nan(15,1500);B5(:,501:2000);nan(15,1500);C5(:,501:2000)];
pl = [pl1,nan(size(pl1,1),75),pl2];
subplot(1,2,2)
hold on
heatmap(flipud(pl),[],[],[],'Colormap',cmap,'ColorLevels',1000,'MaxColorValue',.5,'MinColorValue',-.5,'NaNColor',[1 1 1]);
xticks([750 750+1500+75])
xticklabels({'Flavor','Water'})
title('Familiar day')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0, 0],'TickDir','out')
hold off
%% Fig ED10i

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 10i','FontWeight','bold')
load([data_path,'\Neuropixels\Familiarization\familiarization-summary.mat'],'data')

hold on
axis square
A = []; B = [];
for i = unique(data.meta.session)
    A(i) = sum(data.stats.significant & data.stats.preference>0 & data.meta.session==i)./sum(data.meta.session==i)*100;
    B(i) = sum(data.stats_retrieval.significant & data.stats_retrieval.preference>0 & data.meta.session==i)./sum(data.meta.session==i)*100;
end
for i = 1:length(A)
    plot([1 2],[A(i),B(i)],'color',[.85 .85 .85],'linewidth',1)
end
xx = [NaN mean(A)-std(A)./sqrt(length(A)) mean(A) mean(A)+std(A)./sqrt(length(A)) NaN];
plot([-0.125 0.125]+1,[xx(3) xx(3)],'color','k','LineWidth',1)
plot([1 1],[xx(2) xx(4)],'color','k','LineWidth',1)
xx = [NaN mean(B)-std(B)./sqrt(length(B)) mean(B) mean(B)+std(B)./sqrt(length(B)) NaN];
plot([-0.125 0.125]+2,[xx(3) xx(3)],'color','k','LineWidth',1)
plot([1 1]*2,[xx(2) xx(4)],'color','k','LineWidth',1)
ylabel('Flavor-preferring neurons (%)'); ylim([0 60]); yticks(0:10:60);
xlim([.5 2.5]); xticks([1 2]); xticklabels({'Novel','Familiar'});
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

[p,~,stats]=signrank(A,B,'method','exact');
StatsTbl(end+1,:) = table({'ED 10i'},{'Novel vs. Familiar'},{'Wilcoxon signed-rank'},{'N/A'},{[length(A)]},stats.signedrank,p);
%% Fig ED10j

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 10j','FontWeight','bold')
load([data_path,'\Neuropixels\Familiarization\familiarization-summary.mat'],'data')

subplot(1,2,1)
hold on
axis square
times = -5:.01:10;
times = times(1:end-1) + mean(diff(times))/2;
idx = find(data.stats.significant & data.stats.preference<0);
A1 = data.psth.smooth.novel(idx,501:end);
A2 = data.psth.smooth.novel_retrieval(idx,501:end);
fill([times fliplr(times)],[mean(A1)+std(A1)/sqrt(size(A1,1)) fliplr(mean(A1)-std(A1)/sqrt(size(A1,1)))],[.85 .85 .85],'LineStyle','none')
a = plot(times,mean(A1),'k','LineWidth',1);
fill([times fliplr(times)],[mean(A2)+std(A2)/sqrt(size(A2,1)) fliplr(mean(A2)-std(A2)/sqrt(size(A2,1)))],[216 231 243]/255,'LineStyle','none')
b = plot(times,mean(A2),'color',[55 136 193]/255,'LineWidth',1);
ylim([-.1 .1])
yticks(-.1:.05:.1)
xticks(-5:5:10)
xlim([-5 10])
xlabel('Time (s)')
ylabel('Spiking (σ)')
title('All Water-preferring')
legend([a,b],{'Novel day','Familiar day'})
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

subplot(1,2,2)
axis square
hold on
idx = find(data.stats.significant & data.stats.preference<0);
A = data.stats.response.novel(idx);
B = data.stats_retrieval.response.novel(idx);
simpleboxplot(1,A,'k')
simpleboxplot(2,B,[55 136 193]/255)
ylabel('Flavor response (σ)'); ylim([-.2 .2]); yticks(ylim)
xticks([1,2])
xticklabels({'Novel day','Familiar day'})
xlim([0.25 2.75])
title('All Water-preferring')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

[p,~,stat] = signrank(A,B,'method','exact');
StatsTbl(end+1,:) = table({'ED 10j'},{'Novel vs. Familiar'},{'Wilcoxon signed-rank'},{'N/A'},{[length(A)]},stat.signedrank,p);
%% Fig ED10k

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 10k','FontWeight','bold')
load([data_path,'\Neuropixels\Familiarization\familiarization-summary.mat'],'data')

cmap = struct;
idx1 = find(data.stats.significant & data.stats.preference>0);
idx2 = find(data.stats.significant & data.stats.preference<0);
idx = sort([idx1,idx2]);
PCAdata.novel = data.psth.smooth.novel(idx,501:end) - mean([data.psth.smooth.novel(idx,501:600),data.psth.smooth.water(idx,501:600)],2);
PCAdata.water = data.psth.smooth.water(idx,501:end) - mean([data.psth.smooth.novel(idx,501:600),data.psth.smooth.water(idx,501:600)],2);
PCAdata.novel_retrieval = data.psth.smooth.novel_retrieval(idx,501:end) - mean([data.psth.smooth.novel(idx,501:600),data.psth.smooth.water(idx,501:600)],2);
PCAdata.water_retrieval = data.psth.smooth.water_retrieval(idx,501:end) - mean([data.psth.smooth.novel(idx,501:600),data.psth.smooth.water(idx,501:600)],2);
t=max(abs([PCAdata.novel,PCAdata.water]),[],2);
PCAdata.novel = PCAdata.novel./t;
PCAdata.water = PCAdata.water./t;
PCAdata.novel_retrieval = PCAdata.novel_retrieval./t;
PCAdata.water_retrieval = PCAdata.water_retrieval./t;
cmap.novel = cbrewer('seq','Reds',size(PCAdata.novel,2)*1.5,'spline'); cmap.novel = cmap.novel(size(PCAdata.novel,2)/2+1:size(PCAdata.novel,2)*1.5,:);
cmap.water = cbrewer('seq','Blues',size(PCAdata.water,2)*1.5,'spline'); cmap.water = cmap.water(size(PCAdata.water,2)/2+1:size(PCAdata.water,2)*1.5,:);

x = [PCAdata.novel(:,501:1000) PCAdata.water(:,501:1000)]; x = x-mean(x);
[coeff,score,latent,tsquared,explained,mu] = pca(x,'Centered',false);

PCA_FAM = struct;
ax1 = subplot(1,2,1);
axis square
hold on
pc1 = []; pc2 = [];
for i = 1:size(PCAdata.water,2)
    x = PCAdata.water(:,i); x = x-mean(x);
    y = x'*score; pc1(i) = y(1); pc2(i) = y(2);
    scatter(pc1(i),pc2(i),100,cmap.water(i,:),'filled')
    if i>1
        plot([pc1(i-1) pc1(i)],[pc2(i-1) pc2(i)],'Color',cmap.water(i,:),'LineWidth',2)
    end
end
pc1 = []; pc2 = [];
for i = 1:size(PCAdata.novel,2)
    x = PCAdata.novel(:,i); x = x-mean(x);
    y = x'*score; pc1(i) = y(1); pc2(i) = y(2);
    scatter(pc1(i),pc2(i),100,cmap.novel(i,:),'filled')
    if i>1
        plot([pc1(i-1) pc1(i)],[pc2(i-1) pc2(i)],'Color',cmap.novel(i,:),'LineWidth',2)
    end
end
PCA_FAM.novel = pc2;
x1 = xlim; y1 = ylim;
xticks([]); yticks([]);
xlabel('PC1'); ylabel('PC2');
set(gca,'FontSize',12,'LineWidth',1)
title('Novel day','FontWeight','normal')
hold off

ax2 = subplot(1,2,2);
axis square
hold on
pc1 = []; pc2 = [];
for i = 1:size(PCAdata.water_retrieval,2)
    x = PCAdata.water_retrieval(:,i); x = x-mean(x);
    y = x'*score; pc1(i) = y(1); pc2(i) = y(2);
    scatter(pc1(i),pc2(i),100,cmap.water(i,:),'filled')
    if i>1
        plot([pc1(i-1) pc1(i)],[pc2(i-1) pc2(i)],'Color',cmap.water(i,:),'LineWidth',2)
    end
end
pc1 = []; pc2 = [];
for i = 1:size(PCAdata.novel_retrieval,2)
    x = PCAdata.novel_retrieval(:,i); x = x-mean(x);
    y = x'*score; pc1(i) = y(1); pc2(i) = y(2);
    scatter(pc1(i),pc2(i),100,cmap.novel(i,:),'filled')
    if i>1
        plot([pc1(i-1) pc1(i)],[pc2(i-1) pc2(i)],'Color',cmap.novel(i,:),'LineWidth',2)
    end
end
PCA_FAM.familiar = pc2;
y2 = ylim; x2 = xlim;
xticks([]); yticks([]);
xlabel('PC1'); ylabel('PC2');
set(gca,'FontSize',12,'LineWidth',1)
title('Familiar day','FontWeight','normal')
hold off

x = [min([x1,x2]) max([x1 x2])];
y = [min([y1,y2]) max([y1 y2])];
set(ax1,'xlim',x,'ylim',y)
set(ax2,'xlim',x,'ylim',y)
%% Fig ED10l

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 10l','FontWeight','bold')

cmap.novel = cbrewer('seq','Greys',size(PCAdata.novel,2)*1.5,'spline'); cmap.novel = cmap.novel(size(PCAdata.novel,2)/2+1:size(PCAdata.novel,2)*1.5,:);
cmap.retrieval = cbrewer('seq','Reds',size(PCAdata.novel,2)*1.5,'spline'); cmap.retrieval = cmap.retrieval(size(PCAdata.novel,2)/2+1:size(PCAdata.novel,2)*1.5,:);
cmap.familiar = cbrewer('seq','Blues',size(PCAdata.novel,2)*1.5,'spline'); cmap.familiar = cmap.familiar(size(PCAdata.novel,2)/2+1:size(PCAdata.novel,2)*1.5,:);

subplot(1,2,1)
axis square
hold on
times = -5:.01:10;
times = times(1:end-1) + median(diff(times))/2;
for i = 1:length(times)
    scatter(times(i),PCA_CFA.novel(i),100,cmap.novel(i,:),'filled')
end
for i = 1:length(times)
    scatter(times(i),PCA_CFA.retrieval(i),100,cmap.retrieval(i,:),'filled')
end
xlabel('Time (s)'); xlim([-5 10]); xticks(-5:5:10);
ylabel('PC2'); yticks([]);
title('Novel → CFA Retrieval')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

subplot(1,2,2)
axis square
hold on
times = -5:.01:10;
times = times(1:end-1) + median(diff(times))/2;
for i = 1:length(times)
    scatter(times(i),PCA_FAM.novel(i),100,cmap.novel(i,:),'filled')
end
for i = 1:length(times)
    scatter(times(i),PCA_FAM.familiar(i),100,cmap.familiar(i,:),'filled')
end
xlabel('Time (s)'); xlim([-5 10]); xticks(-5:5:10);
ylabel('PC2'); yticks([]);
title('Novel → Familiar')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

drawnow
%% Stats table

StatsTbl{cellfun(@(x) isequal(x,'Pearson correlation'),StatsTbl{:,3}),6} = round(StatsTbl{cellfun(@(x) isequal(x,'Pearson correlation'),StatsTbl{:,3}),6},3);
StatsTbl{cellfun(@(x) ~isequal(x,'Pearson correlation'),StatsTbl{:,3}),6} = round(StatsTbl{cellfun(@(x) ~isequal(x,'Pearson correlation'),StatsTbl{:,3}),6},2);
StatsTbl{:,7} = round(StatsTbl{:,7},2,'significant');
StatsTbl{StatsTbl{:,7}>0.05,8} = {'NS'};
StatsTbl{StatsTbl{:,7}<=0.05,8} = {'*'};
StatsTbl{StatsTbl{:,7}<=0.01,8} = {'**'};
StatsTbl{StatsTbl{:,7}<=0.001,8} = {'***'};
StatsTbl{StatsTbl{:,7}<=0.0001,8} = {'****'};
StatsTbl = renamevars(StatsTbl,'Var8','Significant');
disp(StatsTbl); disp(' ');