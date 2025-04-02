function fig_ED9(data_path)

disp('Generating panels for Extended Data Figure 9...')
%% Fig ED9a

openfig([data_path,'/Neuropixels/CGRP-CEA stim CFA/cgrp-cea-stim-trajectories.fig']);
set(gcf,'Position', get(0, 'Screensize'))
%% Fig ED9b

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 9b','FontWeight','bold')
load([data_path,'/Neuropixels/CGRP-CEA stim CFA/cgrp-cea-stim-conditioning-summary.mat'],'data')

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
yticklabels({['Nonselective'],['Water'],['Flavour']})
xticks([750 750+1500+75])
xticklabels({'Novel flavour','Water'})
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0, 0],'TickDir','out')
hold off

pl2 = [A5;nan(15,75);B5;nan(15,75);C5];
subplot(1,2,2)
hold on
heatmap(flipud(pl2),[],[],[],'Colormap',cmap,'ColorLevels',1000,'MaxColorValue',.5,'MinColorValue',-.5,'NaNColor',[1 1 1]);
ytickangle(90)
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off
%% Fig ED9c

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 9c','FontWeight','bold')
load([data_path,'/Neuropixels/CGRP-CEA stim CFA/cgrp-cea-stim-conditioning-summary.mat'],'data')

subplot(1,2,1)
axis square
hold on
A = [data.psth.drinking_period(:,6:end),data.psth.cgrp_period(:,31:end)];
t = [0:1:89]+.5;
plot([0 0],[-.1 .3],'k','LineWidth',1)
idx = find(~data.stats.significant);
fill([t fliplr(t)],[mean(A(idx,:))+std(A(idx,:))/sqrt(size(A(idx,:),1)) fliplr(mean(A(idx,:))-std(A(idx,:))/sqrt(size(A(idx,:),1)))],[.85 .85 .85],'LineStyle','none');
a = plot(t,mean(A(idx,:)),'k','LineWidth',1);
idx = find(data.stats.significant & data.stats.preference<0);
fill([t fliplr(t)],[mean(A(idx,:))+std(A(idx,:))/sqrt(size(A(idx,:),1)) fliplr(mean(A(idx,:))-std(A(idx,:))/sqrt(size(A(idx,:),1)))],[216 231 243]/255,'LineStyle','none');
b = plot(t,mean(A(idx,:)),'color',[55 136 193]/255,'LineWidth',1);
idx = find(data.stats.significant & data.stats.preference>0);
fill([t fliplr(t)],[mean(A(idx,:))+std(A(idx,:))/sqrt(size(A(idx,:),1)) fliplr(mean(A(idx,:))-std(A(idx,:))/sqrt(size(A(idx,:),1)))],[252 216 213]/255,'LineStyle','none');
c = plot(t,mean(A(idx,:)),'color',[229 45 38]/255,'LineWidth',1);
legend([c,b,a],{'Flavour-pref','Water-pref','Nonselective'})
ylim([-.05 .2])
yticks(-.05:.05:.2)
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
ylim([-.1 .3])
yticks([-.1 .3])
xticks(1:3)
xlim([.25 3.75])
ylabel('CGRP^{CEA} response (σ)')
xticklabels({'Flavour-pref','Water-pref','Nonselective'})
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

stat = []; p = [];
[p(1),~,s] = ranksum(X.Novel,X.Water,'method','approximate'); stat(1) = s.zval;
[p(2),~,s] = ranksum(X.Novel,X.Neither,'method','approximate'); stat(2) = s.zval;
[p(3),~,s] = ranksum(X.Water,X.Neither,'method','approximate'); stat(3) = s.zval;
p = multicmp(p,'up',0.05);
StatsTbl = table({'ED 9c'},{'Flavour-pref vs. Water-pref'},{'Wilcoxon rank-sum (approx.)'},{'3 pairs of neuron groups'},{[length(X.Novel) length(X.Water)]},stat(1),p(1), ...
    'VariableNames',{'Figure panel','Group','Statistical test','Multiple comparisons','Sample size','Test statistic','P-value'});
StatsTbl(end+1,:) = table({'ED 9c'},{'Flavour-pref vs. Nonselective'},{'Wilcoxon rank-sum (approx.)'},{'3 pairs of neuron groups'},{[length(X.Novel) length(X.Neither)]},stat(2),p(2));
StatsTbl(end+1,:) = table({'ED 9c'},{'Water-pref vs. Nonselective'},{'Wilcoxon rank-sum (approx.)'},{'3 pairs of neuron groups'},{[length(X.Water) length(X.Neither)]},stat(3),p(3));
%% Fig ED9d

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 9d','FontWeight','bold')
load([data_path,'/Neuropixels/CGRP-CEA stim CFA/cgrp-cea-stim-conditioning-summary.mat'],'data')

idx = find(data.stats.significant & data.stats.preference>0);
[~,idx2] = sort(data.stats.response.novel(idx),'descend'); idx = idx(idx2);
A3 = data.psth.raw.cgrp(idx,:);

idx = find(data.stats.significant & data.stats.preference<0);
[~,idx2] = sort(data.stats.response.water(idx),'descend'); idx = idx(idx2);
B3 = data.psth.raw.cgrp(idx,:);

idx = find(~data.stats.significant);
[~,idx2] = sort(mean([data.stats.response.novel(idx);data.stats.response.water(idx)]),'descend'); idx = idx(idx2);
C3 = data.psth.raw.cgrp(idx,:);

pl3 = [A3(:,101:600);nan(18,500);B3(:,101:600);nan(18,500);C3(:,101:600)];
cmap = flipud(cbrewer('div','RdBu',1000,'spline')); cmap(cmap<0) = 0;

subplot(1,2,1)
axis square
hold on
t = [-1:.01:3.99]+.005;
idx = find(~data.stats.significant);
fill([t fliplr(t)],[mean(data.psth.raw.cgrp(idx,101:600))+std(data.psth.raw.cgrp(idx,101:600))/sqrt(size(data.psth.raw.cgrp(idx,101:600),1)) fliplr(mean(data.psth.raw.cgrp(idx,101:600))-std(data.psth.raw.cgrp(idx,101:600))/sqrt(size(data.psth.raw.cgrp(idx,101:600),1)))],[.9 .9 .9],'LineStyle','none');
a = plot(t,mean(data.psth.raw.cgrp(idx,101:600)),'k','LineWidth',1);
idx = find(data.stats.significant & data.stats.preference<0);
fill([t fliplr(t)],[mean(data.psth.raw.cgrp(idx,101:600))+std(data.psth.raw.cgrp(idx,101:600))/sqrt(size(data.psth.raw.cgrp(idx,101:600),1)) fliplr(mean(data.psth.raw.cgrp(idx,101:600))-std(data.psth.raw.cgrp(idx,101:600))/sqrt(size(data.psth.raw.cgrp(idx,101:600),1)))],[216 231 243]/255,'LineStyle','none');
b = plot(t,mean(data.psth.raw.cgrp(idx,101:600)),'color',[55 136 193]/255,'LineWidth',1);
idx = find(data.stats.significant & data.stats.preference>0);
fill([t fliplr(t)],[mean(data.psth.raw.cgrp(idx,101:600))+std(data.psth.raw.cgrp(idx,101:600))/sqrt(size(data.psth.raw.cgrp(idx,101:600),1)) fliplr(mean(data.psth.raw.cgrp(idx,101:600))-std(data.psth.raw.cgrp(idx,101:600))/sqrt(size(data.psth.raw.cgrp(idx,101:600),1)))],[252 216 213]/255,'LineStyle','none');
c = plot(t,mean(data.psth.raw.cgrp(idx,101:600)),'color',[229 45 38]/255,'LineWidth',1);
for i = 0:.1:2.9
    plot([i i],[.290 .295],'Color',[54 161 86]/255,'LineWidth',2)
end
ylim([-.1 .3])
yticks(-.1:.1:.3)
xticks(-1:1:4)
xlim([-1 4])
xlabel('Time (s)')
ylabel('Spiking (σ)')
legend([c,b,a],{'Flavour-pref','Water-pref','Nonselective'})
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

subplot(1,2,2)
X.Novel = data.stats.response.cgrp(find(data.stats.significant & data.stats.preference>0));
X.Water = data.stats.response.cgrp(find(data.stats.significant & data.stats.preference<0));
X.Neither = data.stats.response.cgrp(find(~data.stats.significant));
axis square
hold on
simpleboxplot(1,X.Novel,[229 45 38]/255)
simpleboxplot(2,X.Water,[55 136 193]/255)
simpleboxplot(3,X.Neither,[0 0 0])
ylim([-.1 .3])
yticks([-.1 .3])
xticks(1:3)
xlim([.25 3.75])
ylabel('CGRP^{CEA} stim bout response (σ)')
xticklabels({'Flavour-pref','Water-pref','Nonselective'})
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

stat = []; p = [];
[p(1),~,s] = ranksum(X.Novel,X.Water,'method','approximate'); stat(1) = s.zval;
[p(2),~,s] = ranksum(X.Novel,X.Neither,'method','approximate'); stat(2) = s.zval;
[p(3),~,s] = ranksum(X.Water,X.Neither,'method','approximate'); stat(3) = s.zval;
p = multicmp(p,'up',0.05);
StatsTbl(end+1,:) = table({'ED 9d'},{'Flavour-pref vs. Water-pref'},{'Wilcoxon rank-sum (approx.)'},{'3 pairs of neuron groups'},{[length(X.Novel) length(X.Water)]},stat(1),p(1));
StatsTbl(end+1,:) = table({'ED 9d'},{'Flavour-pref vs. Nonselective'},{'Wilcoxon rank-sum (approx.)'},{'3 pairs of neuron groups'},{[length(X.Novel) length(X.Neither)]},stat(2),p(2));
StatsTbl(end+1,:) = table({'ED 9d'},{'Water-pref vs. Nonselective'},{'Wilcoxon rank-sum (approx.)'},{'3 pairs of neuron groups'},{[length(X.Water) length(X.Neither)]},stat(3),p(3));
%% Fig ED9e

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 9e','FontWeight','bold')
load([data_path,'/Neuropixels/CGRP-CEA stim CFA/cgrp-cea-stim-conditioning-summary.mat'],'data')

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
water = struct; novel = struct; cgrp = struct;

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
%% Fig ED9g

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 9g','FontWeight','bold')

Control = load([data_path,'/Neuropixels/LiCl CFA/licl-control-retrieval-behavior.mat'],'data');
taCasp3 = load([data_path,'/Neuropixels/LiCl CFA/licl-tacasp3-retrieval-behavior.mat'],'data');

axis square
hold on
fill([Control.data.time fliplr(Control.data.time)],[mean(Control.data.flavor)+std(Control.data.flavor)/sqrt(size(Control.data.flavor,1)) fliplr(mean(Control.data.flavor)-std(Control.data.flavor)/sqrt(size(Control.data.flavor,1)))]*0.02,[.85 .85 .85],'LineStyle','none');
a = plot(Control.data.time,mean(Control.data.flavor)*0.02,'Color','k','LineWidth',1);
fill([taCasp3.data.time fliplr(taCasp3.data.time)],[mean(taCasp3.data.flavor)+std(taCasp3.data.flavor)/sqrt(size(taCasp3.data.flavor,1)) fliplr(mean(taCasp3.data.flavor)-std(taCasp3.data.flavor)/sqrt(size(taCasp3.data.flavor,1)))]*0.02,[252 216 213]/255,'LineStyle','none');
b = plot(taCasp3.data.time,mean(taCasp3.data.flavor)*0.02,'Color',[229 45 38]/255,'LineWidth',1);
xlabel('Time (min)'); xlim([0 30]); xticks(0:10:30);
ylabel('Flavour acceptance (ml)'); ylim([0 0.6]); yticks(0:.2:.6);
legend([a,b],{'Control','taCasp3'})
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

[p,~,stats] = ranksum(taCasp3.data.flavor(:,end),Control.data.flavor(:,end),'method','exact');
StatsTbl(end+1,:) = table({'ED 9g'},{'taCasp3 vs. Control'},{'Wilcoxon rank-sum'},{'N/A'},{[length(taCasp3.data.flavor(:,end)) length(Control.data.flavor(:,end))]},stats.ranksum,p(1));
%% Fig ED9i

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 9i','FontWeight','bold')
load([data_path,'/Neuropixels/CGRP stim + LiCl/cgrp-stim-licl-summary.mat'],'data')

cmap = flipud(cbrewer('div','RdBu',1000,'spline')); cmap(cmap<0) = 0;

idx1 = find(ismember(data.stats.cgrp_cluster,[2,4]));
idx2 = find(ismember(data.stats.cgrp_cluster,[1,3]));
subplot(1,2,1)
hold on
pl = [data.psth.raw.cgrp(idx1,101:400);nan(10,300);data.psth.raw.cgrp(idx2,101:400)];
heatmap(flipud(pl),[],[],[],'Colormap',cmap,'ColorLevels',1000,'MaxColorValue',1,'MinColorValue',-1,'NaNColor',[1 1 1]);
yticks([length(idx2)./2 length(idx2)+length(idx1)./2+10]+0.5)
ytickangle(90)
yticklabels({'Other','CGRP-activated'})
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0, 0],'TickDir','out')
hold off

subplot(1,2,2)
hold on
pl2 = [data.psth.licl_period(idx1,:);nan(10,20);data.psth.licl_period(idx2,:)];
heatmap(flipud(pl2),[],[],[],'Colormap',cmap,'ColorLevels',1000,'MaxColorValue',1,'MinColorValue',-1,'NaNColor',[1 1 1]);
xlim([0 20]+.5)
ytickangle(90)
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0, 0],'TickDir','out')
hold off
%% Fig ED9j

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 9j','FontWeight','bold')
load([data_path,'/Neuropixels/CGRP stim + LiCl/cgrp-stim-licl-summary.mat'],'data')

subplot(1,2,1)
axis square
hold on
t = [-5:1:14]+.5;
idx = find(ismember(data.stats.cgrp_cluster,[1,3]));
fill([t fliplr(t)],[mean(data.psth.licl_period(idx,:))+std(data.psth.licl_period(idx,:))/sqrt(size(data.psth.licl_period(idx,:),1)) fliplr(mean(data.psth.licl_period(idx,:))-std(data.psth.licl_period(idx,:))/sqrt(size(data.psth.licl_period(idx,:),1)))],[.85 .85 .85],'LineStyle','none');
a = plot(t,mean(data.psth.licl_period(idx,:)),'k','LineWidth',1);
idx = find(ismember(data.stats.cgrp_cluster,[2,4]));
fill([t fliplr(t)],[mean(data.psth.licl_period(idx,:))+std(data.psth.licl_period(idx,:))/sqrt(size(data.psth.licl_period(idx,:),1)) fliplr(mean(data.psth.licl_period(idx,:))-std(data.psth.licl_period(idx,:))/sqrt(size(data.psth.licl_period(idx,:),1)))],[205 231 213]/255,'LineStyle','none');
b = plot(t,mean(data.psth.licl_period(idx,:)),'color',[54 161 86]/255,'LineWidth',1);
xlabel('Time (min)'); xlim([-5 15]); xticks(-30:5:15);
ylabel('Spiking (σ)'); ylim([-.05 .1]); yticks(-.1:.05:.3);
legend([b,a],{'CGRP-activated','Other'})
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

subplot(1,2,2)
hold on
axis square
simpleboxplot(1,data.stats.response.licl(ismember(data.stats.cgrp_cluster,[2,4])),[54 161 86]/255)
simpleboxplot(2,data.stats.response.licl(ismember(data.stats.cgrp_cluster,[1,3])),'k')
ylim([-.1 .2])
yticks([-.1 .2])
xlim([.25 2.75])
xticks(1:3)
xticklabels({'CGRP-activated','Other'})
ylabel('LiCl response (σ)')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

[p,~,stats] = ranksum(data.stats.response.licl(ismember(data.stats.cgrp_cluster,[2,4])),data.stats.response.licl(ismember(data.stats.cgrp_cluster,[1,3])),'method','approximate');
StatsTbl(end+1,:) = table({'ED 9j'},{'CGRP-activated vs. Other'},{'Wilcoxon rank-sum (approx.)'},{'N/A'},{[length(data.stats.response.licl(ismember(data.stats.cgrp_cluster,[2,4]))) length(data.stats.response.licl(ismember(data.stats.cgrp_cluster,[1,3])))]},stats.zval,p(1));

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