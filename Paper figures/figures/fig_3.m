function fig_3(data_path)

disp('Generating panels for Figure 3...')
%% Fig 3c

openfig([data_path,'\Neuropixels\CGRP stim CFA\cgrp-stim-trajectories.fig']);
set(gcf,'Position', get(0, 'Screensize'))
%% Fig 3d

figure('Position', get(0, 'Screensize'))
sgtitle('Figure 3d','FontWeight','bold')
load([data_path,'\Neuropixels\CGRP stim CFA\cgrp-stim-conditioning-summary.mat'],'data')

idx = find(data.stats.significant & data.stats.preference>0);
[~,idx2] = sort(data.stats.response.novel(idx),'descend'); idx = idx(idx2);
A1 = data.psth.smooth.novel(idx,:);
A2 = data.psth.smooth.water(idx,:);
A5 = data.psth.cgrp_period(idx,:);

idx = find(data.stats.significant & data.stats.preference<0);
[~,idx2] = sort(data.stats.response.water(idx),'descend'); idx = idx(idx2);
B1 = data.psth.smooth.novel(idx,:);
B2 = data.psth.smooth.water(idx,:);
B5 = data.psth.cgrp_period(idx,:);

idx = find(~data.stats.significant);
[~,idx2] = sort(mean([data.stats.response.novel(idx);data.stats.response.water(idx)]),'descend'); idx = idx(idx2);
C1 = data.psth.smooth.novel(idx,:);
C2 = data.psth.smooth.water(idx,:);
C5 = data.psth.cgrp_period(idx,:);

pl1 = [A1(:,501:2000);nan(15,1500);B1(:,501:2000);nan(15,1500);C1(:,501:2000)];
pl2 = [A2(:,501:2000);nan(15,1500);B2(:,501:2000);nan(15,1500);C2(:,501:2000)];
pl = [pl1,nan(size(pl1,1),75),pl2];
cmap = flipud(cbrewer('div','RdBu',1000,'spline')); cmap(cmap<0) = 0;

subplot(1,2,1)
hold on
heatmap(flipud(pl),[],[],[],'Colormap',cmap,'ColorLevels',1000,'MaxColorValue',.5,'MinColorValue',-.5,'NaNColor',[1 1 1]);
nomod = size(C1,1); nopref = size(B1,1); familiarpref = size(A1,1);
yticks([nomod./2 nomod+nopref./2+15 nomod+nopref+familiarpref./2+30]+0.5)
ytickangle(90)
yticklabels({['Non-selective'],['Water'],['Flavor']})
xticks([750 750+1500+75])
xticklabels({'Novel flavor','Water'})
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0, 0],'TickDir','out')
hold off

pl2 = [A5;nan(15,75);B5;nan(15,75);C5];
subplot(1,2,2)
hold on
heatmap(flipud(pl2),[],[],[],'Colormap',cmap,'ColorLevels',1000,'MaxColorValue',.5,'MinColorValue',-.5,'NaNColor',[1 1 1]);
ytickangle(90)
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off
%% Fig 3e

figure('Position', get(0, 'Screensize'))
sgtitle('Figure 3e','FontWeight','bold')
load([data_path,'\Neuropixels\CGRP stim CFA\cgrp-stim-conditioning-summary.mat'],'data')

subplot(1,2,1)
axis square
hold on
A = [data.psth.drinking_period(:,6:end),data.psth.cgrp_period(:,31:end)];
t = [0:1:89]+.5;
plot([0 0],[-.1 .3],'k','LineWidth',1)
idx = find(~data.stats.significant);
fill([t fliplr(t)],[nanmean(A(idx,:))+nanstd(A(idx,:))/sqrt(size(A(idx,:),1)) fliplr(nanmean(A(idx,:))-nanstd(A(idx,:))/sqrt(size(A(idx,:),1)))],[.85 .85 .85],'LineStyle','none');
a = plot(t,nanmean(A(idx,:)),'k','LineWidth',1);
idx = find(data.stats.significant & data.stats.preference<0);
fill([t fliplr(t)],[nanmean(A(idx,:))+nanstd(A(idx,:))/sqrt(size(A(idx,:),1)) fliplr(nanmean(A(idx,:))-nanstd(A(idx,:))/sqrt(size(A(idx,:),1)))],[216 231 243]/255,'LineStyle','none');
b = plot(t,nanmean(A(idx,:)),'color',[55 136 193]/255,'LineWidth',1);
idx = find(data.stats.significant & data.stats.preference>0);
fill([t fliplr(t)],[nanmean(A(idx,:))+nanstd(A(idx,:))/sqrt(size(A(idx,:),1)) fliplr(nanmean(A(idx,:))-nanstd(A(idx,:))/sqrt(size(A(idx,:),1)))],[252 216 213]/255,'LineStyle','none');
c = plot(t,nanmean(A(idx,:)),'color',[229 45 38]/255,'LineWidth',1);
legend([c,b,a],{'Flavor-pref','Water-pref','Non-selective'})
ylim([-.05 .15])
yticks(-.05:.05:.15)
xticks(0:15:90)
xlim([0 90])
xlabel('Time (min)')
ylabel('Spiking (σ)')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

subplot(1,2,2)
X.Novel = data.stats.response.cgrp_period(find(data.stats.significant & data.stats.preference>0));
X.Water = data.stats.response.cgrp_period(find(data.stats.significant & data.stats.preference<0));
X.Neither = data.stats.response.cgrp_period(find(~data.stats.significant));
axis square
hold on
simpleboxplot(1,X.Novel,[229 45 38]/255)
simpleboxplot(2,X.Water,[55 136 193]/255)
simpleboxplot(3,X.Neither,[0 0 0])
ylim([-.1 .2])
yticks([-.1 .2])
xticks(1:3)
xlim([.25 3.75])
ylabel('CGRP response (σ)')
xticklabels({'Flavor-pref','Water-pref','Non-selective'})
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

stat = []; p = [];
[p(1),~,s] = ranksum(X.Novel,X.Water,'method','approximate'); stat(1) = s.zval;
[p(2),~,s] = ranksum(X.Novel,X.Neither,'method','approximate'); stat(2) = s.zval;
[p(3),~,s] = ranksum(X.Water,X.Neither,'method','approximate'); stat(3) = s.zval;
p = multicmp(p,'up',0.05);
StatsTbl = table({'3e'},{'Flavor-pref vs. Water-pref'},{'Wilcoxon rank-sum (approx.)'},{'3 pairs of neuron groups'},{[length(X.Novel) length(X.Water)]},stat(1),p(1), ...
    'VariableNames',{'Figure panel','Group','Statistical test','Multiple comparisons','Sample size','Test statistic','P-value'});
StatsTbl(end+1,:) = table({'3e'},{'Flavor-pref vs. Non-selective'},{'Wilcoxon rank-sum (approx.)'},{'3 pairs of neuron groups'},{[length(X.Novel) length(X.Neither)]},stat(2),p(2));
StatsTbl(end+1,:) = table({'3e'},{'Water-pref vs. Non-selective'},{'Wilcoxon rank-sum (approx.)'},{'3 pairs of neuron groups'},{[length(X.Water) length(X.Neither)]},stat(3),p(3));
%% Fig 3f

figure('Position', get(0, 'Screensize'))
sgtitle('Figure 3f','FontWeight','bold')
load([data_path,'\Neuropixels\CGRP stim CFA\cgrp-stim-conditioning-summary.mat'],'data')

idx = find(data.stats.significant & data.stats.preference>0);
[~,idx2] = sort(data.stats.response.novel(idx),'descend'); idx = idx(idx2);
A3 = data.psth.raw.cgrp(idx,:);

idx = find(data.stats.significant & data.stats.preference<0);
[~,idx2] = sort(data.stats.response.water(idx),'descend'); idx = idx(idx2);
B3 = data.psth.raw.cgrp(idx,:);

idx = find(~data.stats.significant);
[~,idx2] = sort(mean([data.stats.response.novel(idx);data.stats.response.water(idx)]),'descend'); idx = idx(idx2);
C3 = data.psth.raw.cgrp(idx,:);

pl3 = [A3(:,101:600);nan(15,500);B3(:,101:600);nan(15,500);C3(:,101:600)];
cmap = flipud(cbrewer('div','RdBu',1000,'spline')); cmap(cmap<0) = 0;

subplot(1,3,1)
axis square
hold on
heatmap(flipud(pl3),[],[],[],'Colormap',cmap,'ColorLevels',1000,'MaxColorValue',.5,'MinColorValue',-.5,'NaNColor',[1 1 1]);
nomod = size(C3,1); nopref = size(B3,1); familiarpref = size(A3,1);
yticks([nomod./2 nomod+nopref./2+15 nomod+nopref+familiarpref./2+30]+0.5)
ytickangle(90)
yticklabels({['Non-selective'],['Water'],['Flavor']})
xticks([750 750+1500+75])
xticklabels({'Novel flavor','Water'})
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0, 0],'TickDir','out')
hold off

subplot(1,3,2)
axis square
hold on
t = [-1:.01:3.99]+.005;
idx = find(~data.stats.significant);
fill([t fliplr(t)],[nanmean(data.psth.raw.cgrp(idx,101:600))+nanstd(data.psth.raw.cgrp(idx,101:600))/sqrt(size(data.psth.raw.cgrp(idx,101:600),1)) fliplr(nanmean(data.psth.raw.cgrp(idx,101:600))-nanstd(data.psth.raw.cgrp(idx,101:600))/sqrt(size(data.psth.raw.cgrp(idx,101:600),1)))],[.9 .9 .9],'LineStyle','none');
a = plot(t,nanmean(data.psth.raw.cgrp(idx,101:600)),'k','LineWidth',1);
idx = find(data.stats.significant & data.stats.preference<0);
fill([t fliplr(t)],[nanmean(data.psth.raw.cgrp(idx,101:600))+nanstd(data.psth.raw.cgrp(idx,101:600))/sqrt(size(data.psth.raw.cgrp(idx,101:600),1)) fliplr(nanmean(data.psth.raw.cgrp(idx,101:600))-nanstd(data.psth.raw.cgrp(idx,101:600))/sqrt(size(data.psth.raw.cgrp(idx,101:600),1)))],[216 231 243]/255,'LineStyle','none');
b = plot(t,nanmean(data.psth.raw.cgrp(idx,101:600)),'color',[55 136 193]/255,'LineWidth',1);
idx = find(data.stats.significant & data.stats.preference>0);
fill([t fliplr(t)],[nanmean(data.psth.raw.cgrp(idx,101:600))+nanstd(data.psth.raw.cgrp(idx,101:600))/sqrt(size(data.psth.raw.cgrp(idx,101:600),1)) fliplr(nanmean(data.psth.raw.cgrp(idx,101:600))-nanstd(data.psth.raw.cgrp(idx,101:600))/sqrt(size(data.psth.raw.cgrp(idx,101:600),1)))],[252 216 213]/255,'LineStyle','none');
c = plot(t,nanmean(data.psth.raw.cgrp(idx,101:600)),'color',[229 45 38]/255,'LineWidth',1);
for i = 0:.1:2.9
    plot([i i],[.290 .295],'Color',[54 161 86]/255,'LineWidth',2)
end
ylim([-.1 .3])
yticks(-.1:.1:.3)
xticks(-1:1:4)
xlim([-1 4])
xlabel('Time (s)')
ylabel('Spiking (σ)')
legend([c,b,a],{'Flavor-pref','Water-pref','Non-selective'})
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

subplot(1,3,3)
X.Novel = data.stats.response.cgrp(find(data.stats.significant & data.stats.preference>0));
X.Water = data.stats.response.cgrp(find(data.stats.significant & data.stats.preference<0));
X.Neither = data.stats.response.cgrp(find(~data.stats.significant));
axis square
hold on
simpleboxplot(1,X.Novel,[229 45 38]/255)
simpleboxplot(2,X.Water,[55 136 193]/255)
simpleboxplot(3,X.Neither,[0 0 0])
ylim([-.1 .2])
yticks([-.1 .2])
xticks(1:3)
xlim([.25 3.75])
ylabel('CGRP stim bout response (σ)')
xticklabels({'Flavor-pref','Water-pref','Non-selective'})
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

stat = []; p = [];
[p(1),~,s] = ranksum(X.Novel,X.Water,'method','approximate'); stat(1) = s.zval;
[p(2),~,s] = ranksum(X.Novel,X.Neither,'method','approximate'); stat(2) = s.zval;
[p(3),~,s] = ranksum(X.Water,X.Neither,'method','approximate'); stat(3) = s.zval;
p = multicmp(p,'up',0.05);
StatsTbl(end+1,:) = table({'3f'},{'Flavor-pref vs. Water-pref'},{'Wilcoxon rank-sum (approx.)'},{'3 pairs of neuron groups'},{[length(X.Novel) length(X.Water)]},stat(1),p(1));
StatsTbl(end+1,:) = table({'3f'},{'Flavor-pref vs. Non-selective'},{'Wilcoxon rank-sum (approx.)'},{'3 pairs of neuron groups'},{[length(X.Novel) length(X.Neither)]},stat(2),p(2));
StatsTbl(end+1,:) = table({'3f'},{'Water-pref vs. Non-selective'},{'Wilcoxon rank-sum (approx.)'},{'3 pairs of neuron groups'},{[length(X.Water) length(X.Neither)]},stat(3),p(3));
%% Fig 3j

figure('Position', get(0, 'Screensize'))
sgtitle('Figure 3j','FontWeight','bold')
load([data_path,'\Neuropixels\CGRP stim CFA\cgrp-stim-conditioning-summary.mat'],'data')

cmap = struct;
idx1 = find(data.stats.significant & data.stats.preference>0);
idx2 = find(data.stats.significant & data.stats.preference<0);
idx = sort([idx1,idx2]);
PCAdata.novel = data.psth.smooth.novel(idx,501:end) - mean([data.psth.smooth.novel(idx,501:600),data.psth.smooth.water(idx,501:600)],2);
PCAdata.water = data.psth.smooth.water(idx,501:end) - mean([data.psth.smooth.novel(idx,501:600),data.psth.smooth.water(idx,501:600)],2);
PCAdata.cgrp  = data.psth.smooth.cgrp(idx,101:600)  - mean(data.psth.smooth.cgrp(idx,101:200),2);
t=max(abs([PCAdata.novel,PCAdata.water]),[],2);
PCAdata.novel = PCAdata.novel./t;
PCAdata.water = PCAdata.water./t;
t=max(abs([PCAdata.cgrp]),[],2);
PCAdata.cgrp = PCAdata.cgrp./t;
cmap.novel = cbrewer('seq','Reds',size(PCAdata.novel,2)*1.5,'spline'); cmap.novel = cmap.novel(size(PCAdata.novel,2)/2+1:size(PCAdata.novel,2)*1.5,:);
cmap.water = cbrewer('seq','Blues',size(PCAdata.water,2)*1.5,'spline'); cmap.water = cmap.water(size(PCAdata.water,2)/2+1:size(PCAdata.water,2)*1.5,:);
cmap.cgrp = cbrewer('seq','Greens',size(PCAdata.cgrp,2)*1.5,'spline'); cmap.cgrp = cmap.cgrp(size(PCAdata.cgrp,2)/2+1:size(PCAdata.cgrp,2)*1.5,:); cmap.cgrp(cmap.cgrp<0) = 0;
x = [PCAdata.novel(:,501:1000) PCAdata.water(:,501:1000)]; x = x-mean(x);
[coeff,score,latent,tsquared,explained,mu] = pca(x,'Centered',false);
axis square
hold on
plot(1:10,cumsum(explained(1:10)),'k','linewidth',1)
scatter(1:10,cumsum(explained(1:10)),100,'k','MarkerFacecolor','k')
xlim([0 10]); xticks(0:5:10); xlabel('Principal components');
ylim([0 100]); yticks(0:50:100); ylabel('Var. explained (%)');
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off
%% Fig 3k

figure('Position', get(0, 'Screensize'))
sgtitle('Figure 3k','FontWeight','bold')
water = struct; novel = struct; cgrp = struct;

subplot(1,2,1)
axis square
hold on
pc1 = []; pc2 = [];
for i = 1:size(PCAdata.water,2)
    x = PCAdata.water(:,i); x = x-mean(x);
    y = x'*score; pc1(i) = y(1); pc2(i) = y(2);
    scatter(pc1(i),pc2(i),100,cmap.water(i,:),'filled')
    if i>1
        plot([pc1(i-1) pc1(i)],[pc2(i-1) pc2(i)],'Color',cmap.water(i,:),'LineWidth',1)
    end
end
water.pc1 = pc1;
water.pc2 = pc2;
for i = 1:size(PCAdata.novel,2)
    x = PCAdata.novel(:,i); x = x-mean(x);
    y = x'*score; pc1(i) = y(1); pc2(i) = y(2);
    scatter(pc1(i),pc2(i),100,cmap.novel(i,:),'filled')
    if i>1
        plot([pc1(i-1) pc1(i)],[pc2(i-1) pc2(i)],'Color',cmap.novel(i,:),'LineWidth',1)
    end
end
novel.pc1 = pc1;
novel.pc2 = pc2;
for i = 1:size(PCAdata.cgrp,2)
    x = PCAdata.cgrp(:,i); x = x-mean(x);
    y = x'*score; pc1(i) = y(1); pc2(i) = y(2);
    scatter(pc1(i),pc2(i),100,cmap.cgrp(i,:),'filled')
    if i>1
        plot([pc1(i-1) pc1(i)],[pc2(i-1) pc2(i)],'Color',cmap.cgrp(i,:),'LineWidth',1)
    end
end
cgrp.pc1 = pc1;
cgrp.pc2 = pc2;
xticks([])
yticks([])
xlabel('PC1')
ylabel('PC2')
set(gca,'FontSize',12,'LineWidth',1)
hold off

subplot(2,2,2)
hold on
times = -5:.01:10;
times = times(1:end-1) + median(diff(times))/2;
for i = 1:length(times)
    scatter(times(i),water.pc1(i),100,cmap.water(i,:),'filled')
end
for i = 1:length(times)
    scatter(times(i),novel.pc1(i),100,cmap.novel(i,:),'filled')
end
times = -1:.01:4;
times = times(1:end-1) + median(diff(times))/2;
for i = 1:length(times)
    scatter(times(i),cgrp.pc1(i),100,cmap.cgrp(i,:),'filled')
end
y = ylim;
ylim(y)
yticks([])
xticks(-5:5:15)
xlabel('Time (s)')
ylabel('PC1')
xlim([-5 10])
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

subplot(2,2,4)
hold on
times = -5:.01:10;
times = times(1:end-1) + median(diff(times))/2;
for i = 1:length(times)
    scatter(times(i),water.pc2(i),100,cmap.water(i,:),'filled')
end
for i = 1:length(times)
    scatter(times(i),novel.pc2(i),100,cmap.novel(i,:),'filled')
end
times = -1:.01:4;
times = times(1:end-1) + median(diff(times))/2;
for i = 1:length(times)
    scatter(times(i),cgrp.pc2(i),100,cmap.cgrp(i,:),'filled')
end
y = ylim;
ylim(y)
yticks([])
xticks(-5:5:15)
xlabel('Time (s)')
ylabel('PC2')
xlim([-5 10])
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off
%% Fig 3l

figure('Position', get(0, 'Screensize'))
sgtitle('Figure 3l','FontWeight','bold')

example_mice = [3,4,5,7];
for k = 1:3
    idx1 = find(data.stats.significant & data.stats.preference>0 & data.meta.session==example_mice(k));
    idx2 = find(data.stats.significant & data.stats.preference<0 & data.meta.session==example_mice(k));
    idx = sort([idx1,idx2]);
    PCAdata.novel = data.psth.smooth.novel(idx,501:end) - mean([data.psth.smooth.novel(idx,501:600),data.psth.smooth.water(idx,501:600)],2);
    PCAdata.water = data.psth.smooth.water(idx,501:end) - mean([data.psth.smooth.novel(idx,501:600),data.psth.smooth.water(idx,501:600)],2);
    PCAdata.cgrp  = data.psth.smooth.cgrp(idx,101:600)  - mean(data.psth.smooth.cgrp(idx,101:200),2);
    t=max(abs([PCAdata.novel,PCAdata.water]),[],2);
    PCAdata.novel = PCAdata.novel./t;
    PCAdata.water = PCAdata.water./t;
    t=max(abs([PCAdata.cgrp]),[],2);
    PCAdata.cgrp = PCAdata.cgrp./t;
    cmap = struct;
    cmap.novel = cbrewer('seq','Reds',size(PCAdata.novel,2)*1.5,'spline'); cmap.novel = cmap.novel(size(PCAdata.novel,2)/2+1:size(PCAdata.novel,2)*1.5,:);
    cmap.water = cbrewer('seq','Blues',size(PCAdata.water,2)*1.5,'spline'); cmap.water = cmap.water(size(PCAdata.water,2)/2+1:size(PCAdata.water,2)*1.5,:);
    cmap.cgrp = cbrewer('seq','Greens',size(PCAdata.cgrp,2)*1.5,'spline'); cmap.cgrp = cmap.cgrp(size(PCAdata.cgrp,2)/2+1:size(PCAdata.cgrp,2)*1.5,:); cmap.cgrp(cmap.cgrp<0) = 0;
    x = [PCAdata.novel(:,501:1000) PCAdata.water(:,501:1000)]; x = x-mean(x);
    [coeff,score,latent,tsquared,explained,mu] = pca(x,'Centered',false);
    subplot(2,2,k)
    axis square
    hold on
    pc1 = []; pc2 = []; pc3 = [];
    for i = 1:size(PCAdata.water,2)
        x = PCAdata.water(:,i); x = x-mean(x);
        y = x'*score; pc1(i) = y(1); pc2(i) = y(2);  pc3(i) = y(3);
        scatter(pc1(i),pc2(i),100,cmap.water(i,:),'filled')
        if i>1
            plot([pc1(i-1) pc1(i)],[pc2(i-1) pc2(i)],'Color',cmap.water(i,:),'LineWidth',1)
        end
    end
    for i = 1:size(PCAdata.novel,2)
        x = PCAdata.novel(:,i); x = x-mean(x);
        y = x'*score; pc1(i) = y(1); pc2(i) = y(2);  pc3(i) = y(3);
        scatter(pc1(i),pc2(i),100,cmap.novel(i,:),'filled')
        if i>1
            plot([pc1(i-1) pc1(i)],[pc2(i-1) pc2(i)],'Color',cmap.novel(i,:),'LineWidth',1)
        end
    end
    for i = 1:size(PCAdata.cgrp,2)
        x = PCAdata.cgrp(:,i); x = x-mean(x);
        y = x'*score; pc1(i) = y(1); pc2(i) = y(2);  pc3(i) = y(3);
        scatter(pc1(i),pc2(i),100,cmap.cgrp(i,:),'filled')
        if i>1
            plot([pc1(i-1) pc1(i)],[pc2(i-1) pc2(i)],'Color',cmap.cgrp(i,:),'LineWidth',1)
        end
    end
    %     if k < 3
    %         set(gca,'YDir','reverse')
    %     end
    xticks([])
    yticks([])
    xlabel('PC1')
    ylabel('PC2')
    set(gca,'FontSize',12,'LineWidth',1)
    title(['Animal ',num2str(k)],'FontWeight','normal')
    hold off
end

idx1 = find(data.stats.significant & data.stats.preference>0 & data.meta.session==example_mice(4));
idx2 = find(data.stats.significant & data.stats.preference<0 & data.meta.session==example_mice(4));
idx = sort([idx1,idx2]);
PCAdata.novel = data.psth.smooth.novel(idx,501:end) - mean([data.psth.smooth.novel(idx,501:600),data.psth.smooth.water(idx,501:600)],2);
PCAdata.water = data.psth.smooth.water(idx,501:end) - mean([data.psth.smooth.novel(idx,501:600),data.psth.smooth.water(idx,501:600)],2);
PCAdata.cgrp  = data.psth.smooth.cgrp(idx,101:600)  - mean(data.psth.smooth.cgrp(idx,101:200),2);
t=max(abs([PCAdata.novel,PCAdata.water]),[],2);
PCAdata.novel = PCAdata.novel./t;
PCAdata.water = PCAdata.water./t;
t=max(abs([PCAdata.cgrp]),[],2);
PCAdata.cgrp = PCAdata.cgrp./t;
cmap.novel = cbrewer('seq','Reds',size(PCAdata.novel,2)*1.5,'spline'); cmap.novel = cmap.novel(size(PCAdata.novel,2)/2+1:size(PCAdata.novel,2)*1.5,:);
cmap.water = cbrewer('seq','Blues',size(PCAdata.water,2)*1.5,'spline'); cmap.water = cmap.water(size(PCAdata.water,2)/2+1:size(PCAdata.water,2)*1.5,:);
cmap.cgrp = cbrewer('seq','Greens',size(PCAdata.cgrp,2)*1.5,'spline'); cmap.cgrp = cmap.cgrp(size(PCAdata.cgrp,2)/2+1:size(PCAdata.cgrp,2)*1.5,:); cmap.cgrp(cmap.cgrp<0) = 0;
x = [PCAdata.novel(:,501:1000) PCAdata.water(:,501:1000)]; x = x-mean(x);
[coeff,score,latent,tsquared,explained,mu] = pca(x,'Centered',false);
ax4 = subplot(2,2,4);
axis square
hold on
pc1 = []; pc2 = []; pc3 = [];
for i = 1:size(PCAdata.water,2)
    x = PCAdata.water(:,i); x = x-mean(x);
    y = x'*score; pc1(i) = y(1); pc2(i) = y(2);  pc3(i) = y(3);
    scatter3(pc1(i),pc2(i),pc3(i),64,cmap.water(i,:),'filled')
    if i>1
        plot3([pc1(i-1) pc1(i)],[pc2(i-1) pc2(i)],[pc3(i-1) pc3(i)],'Color',cmap.water(i,:),'LineWidth',2)
    end
end
for i = 1:size(PCAdata.novel,2)
    x = PCAdata.novel(:,i); x = x-mean(x);
    y = x'*score; pc1(i) = y(1); pc2(i) = y(2);  pc3(i) = y(3);
    scatter3(pc1(i),pc2(i),pc3(i),64,cmap.novel(i,:),'filled')
    if i>1
        plot3([pc1(i-1) pc1(i)],[pc2(i-1) pc2(i)],[pc3(i-1) pc3(i)],'Color',cmap.novel(i,:),'LineWidth',2)
    end
end
for i = 1:size(PCAdata.cgrp,2)
    x = PCAdata.cgrp(:,i); x = x-mean(x);
    y = x'*score; pc1(i) = y(1); pc2(i) = y(2);  pc3(i) = y(3);
    scatter3(pc1(i),pc2(i),pc3(i),64,cmap.cgrp(i,:),'filled')
    if i>1
        plot3([pc1(i-1) pc1(i)],[pc2(i-1) pc2(i)],[pc3(i-1) pc3(i)],'Color',cmap.cgrp(i,:),'LineWidth',2)
    end
end
xticks([])
yticks([])
zticks([])
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
title('Animal 4')
set(gca,'FontSize',12,'LineWidth',1)
view(ax4,-15,-45)
hold off
%% Fig 3n

figure('Position', get(0, 'Screensize'))
sgtitle('Figure 3n','FontWeight','bold')
load([data_path,'\Neuropixels\LiCl CFA\licl-control-conditioning-summary.mat'],'data')

subplot(1,2,1)
axis square
hold on
A = [data.psth.drinking_period(:,6:end),data.psth.licl_period(:,31:end)];
t = [0:1:89]+.5;
plot([0 0],[-.1 .3],'k','LineWidth',1)
idx = find(~data.stats.significant);
fill([t fliplr(t)],[nanmean(A(idx,:))+nanstd(A(idx,:))/sqrt(size(A(idx,:),1)) fliplr(nanmean(A(idx,:))-nanstd(A(idx,:))/sqrt(size(A(idx,:),1)))],[.85 .85 .85],'LineStyle','none');
a = plot(t,nanmean(A(idx,:)),'k','LineWidth',1);
idx = find(data.stats.significant & data.stats.preference<0);
fill([t fliplr(t)],[nanmean(A(idx,:))+nanstd(A(idx,:))/sqrt(size(A(idx,:),1)) fliplr(nanmean(A(idx,:))-nanstd(A(idx,:))/sqrt(size(A(idx,:),1)))],[216 231 243]/255,'LineStyle','none');
b = plot(t,nanmean(A(idx,:)),'color',[55 136 193]/255,'LineWidth',1);
idx = find(data.stats.significant & data.stats.preference>0);
fill([t fliplr(t)],[nanmean(A(idx,:))+nanstd(A(idx,:))/sqrt(size(A(idx,:),1)) fliplr(nanmean(A(idx,:))-nanstd(A(idx,:))/sqrt(size(A(idx,:),1)))],[252 216 213]/255,'LineStyle','none');
c = plot(t,nanmean(A(idx,:)),'color',[229 45 38]/255,'LineWidth',1);
legend([c,b,a],{'Flavor-pref','Water-pref','Non-selective'})
ylim([-.1 .3])
yticks(-.1:.1:.3)
xticks(0:15:90)
xlim([0 90])
xlabel('Time (min)')
ylabel('Spiking (σ)')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

subplot(1,2,2)
X.Novel = data.stats.response.licl(find(data.stats.significant & data.stats.preference>0));
X.Water = data.stats.response.licl(find(data.stats.significant & data.stats.preference<0));
X.Neither = data.stats.response.licl(find(~data.stats.significant));
axis square
hold on
simpleboxplot(1,X.Novel,[229 45 38]/255)
simpleboxplot(2,X.Water,[55 136 193]/255)
simpleboxplot(3,X.Neither,[0 0 0])
ylim([-.2 .5])
yticks([-.2 .5])
xticks(1:3)
xlim([.25 3.75])
ylabel('LiCl response (σ)')
xticklabels({'Flavor-pref','Water-pref','Non-selective'})
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

stat = []; p = [];
[p(1),~,s] = ranksum(X.Novel,X.Water,'method','approximate'); stat(1) = s.zval;
[p(2),~,s] = ranksum(X.Novel,X.Neither,'method','approximate'); stat(2) = s.zval;
[p(3),~,s] = ranksum(X.Water,X.Neither,'method','approximate'); stat(3) = s.zval;
p = multicmp(p,'up',0.05);
StatsTbl(end+1,:) = table({'3n'},{'Flavor-pref vs. Water-pref'},{'Wilcoxon rank-sum (approx.)'},{'3 pairs of neuron groups'},{[length(X.Novel) length(X.Water)]},stat(1),p(1));
StatsTbl(end+1,:) = table({'3n'},{'Flavor-pref vs. Non-selective'},{'Wilcoxon rank-sum (approx.)'},{'3 pairs of neuron groups'},{[length(X.Novel) length(X.Neither)]},stat(2),p(2));
StatsTbl(end+1,:) = table({'3n'},{'Water-pref vs. Non-selective'},{'Wilcoxon rank-sum (approx.)'},{'3 pairs of neuron groups'},{[length(X.Water) length(X.Neither)]},stat(3),p(3));
%% Fig 3q

figure('Position', get(0, 'Screensize'))
sgtitle('Figure 3q','FontWeight','bold')
load([data_path,'\Neuropixels\LiCl CFA\licl-tacasp3-conditioning-summary.mat'],'data')

subplot(1,2,1)
axis square
hold on
A = [data.psth.drinking_period(:,6:end),data.psth.licl_period(:,31:end)];
t = [0:1:89]+.5;
plot([0 0],[-.1 .3],'k','LineWidth',1)
idx = find(~data.stats.significant);
fill([t fliplr(t)],[nanmean(A(idx,:))+nanstd(A(idx,:))/sqrt(size(A(idx,:),1)) fliplr(nanmean(A(idx,:))-nanstd(A(idx,:))/sqrt(size(A(idx,:),1)))],[.85 .85 .85],'LineStyle','none');
a = plot(t,nanmean(A(idx,:)),'k','LineWidth',1);
idx = find(data.stats.significant & data.stats.preference<0);
fill([t fliplr(t)],[nanmean(A(idx,:))+nanstd(A(idx,:))/sqrt(size(A(idx,:),1)) fliplr(nanmean(A(idx,:))-nanstd(A(idx,:))/sqrt(size(A(idx,:),1)))],[216 231 243]/255,'LineStyle','none');
b = plot(t,nanmean(A(idx,:)),'color',[55 136 193]/255,'LineWidth',1);
idx = find(data.stats.significant & data.stats.preference>0);
fill([t fliplr(t)],[nanmean(A(idx,:))+nanstd(A(idx,:))/sqrt(size(A(idx,:),1)) fliplr(nanmean(A(idx,:))-nanstd(A(idx,:))/sqrt(size(A(idx,:),1)))],[252 216 213]/255,'LineStyle','none');
c = plot(t,nanmean(A(idx,:)),'color',[229 45 38]/255,'LineWidth',1);
legend([c,b,a],{'Flavor-pref','Water-pref','Non-selective'})
ylim([-.1 .3])
yticks(-.1:.1:.3)
xticks(0:15:90)
xlim([0 90])
xlabel('Time (min)')
ylabel('Spiking (σ)')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

subplot(1,2,2)
X.Novel = data.stats.response.licl(find(data.stats.significant & data.stats.preference>0));
X.Water = data.stats.response.licl(find(data.stats.significant & data.stats.preference<0));
X.Neither = data.stats.response.licl(find(~data.stats.significant));
axis square
hold on
simpleboxplot(1,X.Novel,[229 45 38]/255)
simpleboxplot(2,X.Water,[55 136 193]/255)
simpleboxplot(3,X.Neither,[0 0 0])
ylim([-.2 .5])
yticks([-.2 .5])
xticks(1:3)
xlim([.25 3.75])
ylabel('LiCl response (σ)')
xticklabels({'Flavor-pref','Water-pref','Non-selective'})
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

stat = []; p = [];
[p(1),~,s] = ranksum(X.Novel,X.Water,'method','approximate'); stat(1) = s.zval;
[p(2),~,s] = ranksum(X.Novel,X.Neither,'method','approximate'); stat(2) = s.zval;
[p(3),~,s] = ranksum(X.Water,X.Neither,'method','approximate'); stat(3) = s.zval;
p = multicmp(p,'up',0.05);
StatsTbl(end+1,:) = table({'3q'},{'Flavor-pref vs. Water-pref'},{'Wilcoxon rank-sum (approx.)'},{'3 pairs of neuron groups'},{[length(X.Novel) length(X.Water)]},stat(1),p(1));
StatsTbl(end+1,:) = table({'3q'},{'Flavor-pref vs. Non-selective'},{'Wilcoxon rank-sum (approx.)'},{'3 pairs of neuron groups'},{[length(X.Novel) length(X.Neither)]},stat(2),p(2));
StatsTbl(end+1,:) = table({'3q'},{'Water-pref vs. Non-selective'},{'Wilcoxon rank-sum (approx.)'},{'3 pairs of neuron groups'},{[length(X.Water) length(X.Neither)]},stat(3),p(3));

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