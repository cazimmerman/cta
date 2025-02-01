%% setup
addpath(genpath('Z:\Chris\matlab\cz\neuropixels-utils\'));
addpath(genpath('Z:\Chris\matlab\heatmap\'));
clearvars; close all; clc
min_amp = 20;
unittype = {'good','mua'};
regions = 'CEA'; % Options: 'CEA','notCEA','amygdala','all'

CGRP_cutoff = 0.2;
%% load data

library = struct;
library.CEA = {'CEAc','CEAl','CEAm','CEAv','CEAast'};%,'BMAa','COAa','IA','BLAp','PAA'};
library.amygdala = {'CEAc','CEAl','CEAm','CEAv','CEAast','BMAa','BMAp','BLAa','BLAp','BLAv','COAa','COAp','IA','MEA','LA','PAA','AAA','PA'};
library.all = cell(0,0);

data_in = cell(0,0);
load(['data/calca_302-CgrpLicl.mat'])
data_in{end+1} = data;
library.all = unique([library.all,data_in{end}.meta.location]);

load(['data/calca_304-CgrpLicl.mat'])
data_in{end+1} = data;
library.all = unique([library.all,data_in{end}.meta.location]);

load(['data/calca_305-CgrpLicl.mat'])
data_in{end+1} = data;
library.all = unique([library.all,data_in{end}.meta.location]);

load(['data/calca_306-CgrpLicl.mat'])
data_in{end+1} = data;
library.all = unique([library.all,data_in{end}.meta.location]);
library.notCEA = setdiff(library.all,library.CEA);
%% organize data
data = struct;

data.psth.full.cgrp = [];
data.psth.raw.cgrp = [];
data.psth.smooth.cgrp = [];
data.psth.cgrp_period = [];
data.psth.licl_period = [];

data.stats.response.cgrp_train = [];
data.stats.response.cgrp_period = [];
data.stats.response.licl = [];
data.stats.cgrp_cluster = [];

data.meta.session = [];
data.meta.cluster_id = [];
data.meta.mu_cgrp = [];
data.meta.sigma_cgrp = [];
data.meta.mu_baseline = [];
data.meta.sigma_baseline = [];

load('Z:/Chris/matlab/cz/neuropixels-cgrp-acute/GMModel.mat')

for i = 1:length(data_in)
    for j = 1:length(data_in{i}.stats.response.cgrp)
        if data_in{i}.meta.amplitude(j)>=min_amp && ismember(data_in{i}.meta.label{j},unittype) ...
                && ismember(data_in{i}.meta.location{j},library.(regions)) && data_in{i}.meta.sigma_baseline(j) ~= 0
            
            data.stats.response.cgrp_train(end+1)  = mean(mean(data_in{i}.psth.raw.cgrp(j,201:300))) - mean(mean(data_in{i}.psth.raw.cgrp(j,101:200)));
            data.stats.response.cgrp_period(end+1) = nanmean(data_in{i}.psth.cgrp_period(j,6:15))  - nanmean(data_in{i}.psth.cgrp_period(j,1:4));
            data.stats.response.licl(end+1)        = nanmean(data_in{i}.psth.licl_period(j,11:20)) - nanmean(data_in{i}.psth.licl_period(j,1:4));
            
            data.psth.full.cgrp(end+1,:,:) = data_in{i}.psth.full.cgrp(j,:,:);
            data.psth.raw.cgrp(end+1,:)    = data_in{i}.psth.raw.cgrp(j,:) - mean(data_in{i}.psth.raw.cgrp(j,101:200));
            data.psth.smooth.cgrp(end+1,:) = data_in{i}.psth.smooth.cgrp(j,:) - mean(data_in{i}.psth.smooth.cgrp(j,101:200));
            
            data.psth.cgrp_period(end+1,:) = data_in{i}.psth.cgrp_period(j,:) - mean(data_in{i}.psth.cgrp_period(j,1:4));
            data.psth.licl_period(end+1,:) = data_in{i}.psth.licl_period(j,:) - mean(data_in{i}.psth.licl_period(j,1:4));
            
            data.stats.cgrp_cluster(end+1) = cluster(GMModel,data.psth.raw.cgrp(end,201:300));
            
            data.meta.session(end+1) = i;
            data.meta.cluster_id(end+1) = data_in{i}.meta.cluster_id(j);
            data.meta.mu_cgrp(end+1) = data_in{i}.meta.mu_cgrp(j);
            data.meta.sigma_cgrp(end+1) = data_in{i}.meta.sigma_cgrp(j);
            data.meta.mu_baseline(end+1) = data_in{i}.meta.mu_baseline(j);
            data.meta.sigma_baseline(end+1) = data_in{i}.meta.sigma_baseline(j);
            
        end
    end
end

data.psth.cgrp_period(isinf(data.psth.cgrp_period)) = NaN;
data.psth.licl_period(isinf(data.psth.licl_period)) = NaN;
save(['data/cgrp-licl-',regions,'-',num2str(min_amp),'uV-',unittype{end},'.mat'],'data')
%% heatmap
close all
[~,idx] = sort(data.stats.response.cgrp_train,'descend');

idx1 = find(ismember(data.stats.cgrp_cluster,[2,4]));
size(idx1)
%[~,t] = sort(data.stats.response.cgrp_train(idx1),'descend');
%idx1 = idx1(t);
idx2 = find(ismember(data.stats.cgrp_cluster,[1,3]));
size(idx2)
%[~,t] = sort(data.stats.response.cgrp_train(idx2),'descend');
%idx2 = idx2(t);

pl = [data.psth.raw.cgrp(idx1,101:400);nan(10,300);data.psth.raw.cgrp(idx2,101:400)];
cmap = flipud(cbrewer('div','RdBu',1000,'spline')); cmap(cmap<0) = 0;
Clims = [-2 2];

figure('Position', get(0, 'Screensize'))

subplot(1,3,1)
hold on
heatmap(flipud(pl),[],[],[],'Colormap',cmap,'ColorLevels',1000,'MaxColorValue',1,'MinColorValue',-1,'NaNColor',[1 1 1]);
plot([101 101]+0.5,[0 size(pl,1)]+0.5,'k-','LineWidth',1)
plot([201 201]+0.5,[0 size(pl,1)]+0.5,'k-','LineWidth',1)
xticks([0:100:400]+.5)
ylabel([num2str(size(pl,1)),' Units'])
xticklabels({'-1','0','1','2'})
title({'CGRP Trains'})
set(gca,'FontSize',16,'LineWidth',1)
xlabel('Time (sec)')
hold off

pl2 = data.psth.cgrp_period(idx,:);
pl2 = [data.psth.cgrp_period(idx1,:);nan(10,20);data.psth.cgrp_period(idx2,:)];
subplot(1,3,2)
hold on
heatmap(flipud(pl2),[],[],[],'Colormap',cmap,'ColorLevels',1000,'MaxColorValue',1,'MinColorValue',-1,'NaNColor',[1 1 1]);
plot([5 5]+0.5,[0 size(pl2,1)]+0.5,'k-','LineWidth',1)
plot([15 15]+0.5,[0 size(pl2,1)]+0.5,'k-','LineWidth',1)
xticks([0:5:75]+.5)
xticklabels({'-5','0','5','10','15'})
ytickangle(90)
xtickangle(0)
xlabel('Time (min)')
set(gca,'FontSize',16,'LineWidth',1)
title('CGRP Period')
hold off

pl2 = data.psth.licl_period(idx,:);
pl2 = [data.psth.licl_period(idx1,:);nan(10,35);data.psth.licl_period(idx2,:)];
subplot(1,3,3)
hold on
heatmap(flipud(pl2),[],[],[],'Colormap',cmap,'ColorLevels',1000,'MaxColorValue',1,'MinColorValue',-1,'NaNColor',[1 1 1]);
plot([5 5]+0.5,[0 size(pl2,1)]+0.5,'k-','LineWidth',1)
xticks([0:5:75]+.5)
xlim([0 20]+.5)
xticklabels({'-5','0','5','10','15','20','25','30'})
ytickangle(90)
xtickangle(0)
xlabel('Time (min)')
set(gca,'FontSize',16,'LineWidth',1)
h = colorbar('FontSize',16);
ylabel(h,'Spiking (σ)');
title('LiCl Period')
h.Ticks = -1:.5:1;
h.TickLength = 0;
h.Box = 'on';
hold off

saveas(gcf,['plots-png/CGRP-LiCl/CGRP-LiCl-plot-1-',regions,'-',num2str(min_amp),'uV-',unittype{end}],'png')
set(gcf,'renderer','Painters')
saveas(gcf,['plots-eps/CGRP-LiCl/CGRP-LiCl-plot-1-',regions,'-',num2str(min_amp),'uV-',unittype{end}],'epsc')
%% summary plot - full experiment
figure('Position', get(0, 'Screensize'))

[~,idx] = sort(data.stats.response.cgrp_train,'descend'); idx1 = idx(1:floor(length(idx)*CGRP_cutoff));
%[~,idx] = sort(data.stats.response.cgrp_train,'ascend'); idx2 = idx(1:floor(length(idx)*CGRP_cutoff));

idx1 = find(ismember(data.stats.cgrp_cluster,[2,4]));

idx2 = find(ismember(data.stats.cgrp_cluster,[1,3]));
idx3 = setdiff(1:length(data.stats.response.cgrp_train),unique([idx1,idx2]));

subplot(1,3,1)
hold on
axis square
for i = 0:.1:0.9
    plot([i i],[.980 .995],'Color',[0 .5 1],'LineWidth',2)
end
times = -1:.01:4;
times = times(1:end-1) + mean(diff(times))/2;
idx = idx2;
A1 = data.psth.raw.cgrp(idx,101:600);
fill([times fliplr(times)],[mean(A1)+std(A1)/sqrt(size(A1,1)) fliplr(mean(A1)-std(A1)/sqrt(size(A1,1)))],'m','LineStyle','none','FaceAlpha',0.15)
c = plot(times,mean(A1),'m','LineWidth',1);
idx = idx3;
A1 = data.psth.raw.cgrp(idx,101:600);
fill([times fliplr(times)],[mean(A1)+std(A1)/sqrt(size(A1,1)) fliplr(mean(A1)-std(A1)/sqrt(size(A1,1)))],'k','LineStyle','none','FaceAlpha',0.15)
b = plot(times,mean(A1),'k','LineWidth',1);
idx = idx1;
A1 = data.psth.raw.cgrp(idx,101:600);
fill([times fliplr(times)],[mean(A1)+std(A1)/sqrt(size(A1,1)) fliplr(mean(A1)-std(A1)/sqrt(size(A1,1)))],'g','LineStyle','none','FaceAlpha',0.15)
a = plot(times,mean(A1),'g','LineWidth',1);
ylim([-.2 1])
yticks(-.2:.2:1)
legend([a,b,c],{'CGRP stim 1/2','CGRP stim 3','CGRP stim 4'},'box','off')
xticks(-1:1:2)
xlim([-1 2])
xlabel('Time (s)')
ylabel('Spiking (σ)')
title('CGRP Train Responses')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')

subplot(1,3,2)
axis square
hold on
t = [-5:1:14]+.5;
plot([0 0],[-.1 .3],'k','LineWidth',1)
plot([10 10],[-.1 .3],'k','LineWidth',1)

idx = idx3;
fill([t fliplr(t)],[nanmean(data.psth.cgrp_period(idx,:))+nanstd(data.psth.cgrp_period(idx,:))/sqrt(size(data.psth.cgrp_period(idx,:),1)) fliplr(nanmean(data.psth.cgrp_period(idx,:))-nanstd(data.psth.cgrp_period(idx,:))/sqrt(size(data.psth.cgrp_period(idx,:),1)))],[.9 .9 .9],'LineStyle','none');
a = plot(t,nanmean(data.psth.cgrp_period(idx,:)),'k','LineWidth',1);
%scatter(t,nanmean(data.psth.cgrp_period(idx,:)),64,'k','filled')

idx = idx2;
fill([t fliplr(t)],[nanmean(data.psth.cgrp_period(idx,:))+nanstd(data.psth.cgrp_period(idx,:))/sqrt(size(data.psth.cgrp_period(idx,:),1)) fliplr(nanmean(data.psth.cgrp_period(idx,:))-nanstd(data.psth.cgrp_period(idx,:))/sqrt(size(data.psth.cgrp_period(idx,:),1)))],[1 .9 1],'LineStyle','none');
b = plot(t,nanmean(data.psth.cgrp_period(idx,:)),'m','LineWidth',1);
%scatter(t,nanmean(data.psth.cgrp_period(idx,:)),64,'b','filled')

idx = idx1;
fill([t fliplr(t)],[nanmean(data.psth.cgrp_period(idx,:))+nanstd(data.psth.cgrp_period(idx,:))/sqrt(size(data.psth.cgrp_period(idx,:),1)) fliplr(nanmean(data.psth.cgrp_period(idx,:))-nanstd(data.psth.cgrp_period(idx,:))/sqrt(size(data.psth.cgrp_period(idx,:),1)))],[.9 1 .9],'LineStyle','none');
c = plot(t,nanmean(data.psth.cgrp_period(idx,:)),'g','LineWidth',1);
%scatter(t,nanmean(data.psth.cgrp_period(idx,:)),64,'r','filled')

plot([-4.5 -.5],[.275 .275],'k','LineWidth',10)
text(-2.5,.28,'Baseline','FontSize',16,'vert','bottom','horiz','center')
plot([.5 9.5],[.275 .275],'k','LineWidth',10)
text(5,.28,'CGRP','FontSize',16,'vert','bottom','horiz','center')
plot([10.5 14.5],[.275 .275],'k','LineWidth',10)
text(12.5,.28,'Delay','FontSize',16,'vert','bottom','horiz','center')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
ylim([-.1 .3])
yticks(-.1:.1:.3)
xticks(-20:5:30)
xlim([-5 15])
xlabel('Time (min)')
ylabel('Spiking (σ)')
title('CGRP Period Responses')
hold off

subplot(1,3,3)
axis square
hold on
t = [-5:1:29]+.5;
plot([0 0],[-.1 .3],'k','LineWidth',1)

idx = idx3;
fill([t fliplr(t)],[nanmean(data.psth.licl_period(idx,:))+nanstd(data.psth.licl_period(idx,:))/sqrt(size(data.psth.licl_period(idx,:),1)) fliplr(nanmean(data.psth.licl_period(idx,:))-nanstd(data.psth.licl_period(idx,:))/sqrt(size(data.psth.licl_period(idx,:),1)))],[.9 .9 .9],'LineStyle','none');
a = plot(t,nanmean(data.psth.licl_period(idx,:)),'k','LineWidth',1);
%scatter(t,nanmean(data.psth.licl_period(idx,:)),64,'k','filled')

idx = idx2;
fill([t fliplr(t)],[nanmean(data.psth.licl_period(idx,:))+nanstd(data.psth.licl_period(idx,:))/sqrt(size(data.psth.licl_period(idx,:),1)) fliplr(nanmean(data.psth.licl_period(idx,:))-nanstd(data.psth.licl_period(idx,:))/sqrt(size(data.psth.licl_period(idx,:),1)))],[1 .9 1],'LineStyle','none');
b = plot(t,nanmean(data.psth.licl_period(idx,:)),'m','LineWidth',1);
%scatter(t,nanmean(data.psth.licl_period(idx,:)),64,'b','filled')

idx = idx1;
fill([t fliplr(t)],[nanmean(data.psth.licl_period(idx,:))+nanstd(data.psth.licl_period(idx,:))/sqrt(size(data.psth.licl_period(idx,:),1)) fliplr(nanmean(data.psth.licl_period(idx,:))-nanstd(data.psth.licl_period(idx,:))/sqrt(size(data.psth.licl_period(idx,:),1)))],[.9 1 .9],'LineStyle','none');
c = plot(t,nanmean(data.psth.licl_period(idx,:)),'g','LineWidth',1);
%scatter(t,nanmean(data.psth.licl_period(idx,:)),64,'r','filled')

plot([-4.5 -.5],[.275 .275],'k','LineWidth',10)
text(-2.5,.28,'Delay','FontSize',16,'vert','bottom','horiz','center')
plot([.5 29.5],[.275 .275],'k','LineWidth',10)
text(15,.28,'LiCl','FontSize',16,'vert','bottom','horiz','center')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
ylim([-.05 .1])
yticks(-.1:.05:.3)
xticks(-30:5:15)
xlim([-5 15])
xlabel('Time (min)')
ylabel('Spiking (σ)')
title('LiCl Period Responses')
hold off

saveas(gcf,['plots-png/CGRP-LiCl/CGRP-LiCl-plot-2-',regions,'-',num2str(min_amp),'uV-',unittype{end}],'png')
set(gcf,'renderer','Painters')
saveas(gcf,['plots-eps/CGRP-LiCl/CGRP-LiCl-plot-2-',regions,'-',num2str(min_amp),'uV-',unittype{end}],'epsc')
%% scatter plot

figure('Position', get(0, 'Screensize'))

X.Set1 = [0 -1 -2];
X.Set3 = data.stats.response.licl(ismember(data.stats.cgrp_cluster,[1,3]));
X.Set4 = data.stats.response.licl(ismember(data.stats.cgrp_cluster,[2,4]));

[p,~,stats] = ranksum(X.Set4,X.Set3,'method','approximate')

subplot(1,3,1)
hold on
axis square
% plot([.25 3.75],[0 0],'k','LineWidth',1)
% errorbar([3:-1:1],[nanmean(X.Set1),nanmean(X.Set3),nanmean(X.Set4)],[nanstd(X.Set1)/sqrt(length(X.Set1)),nanstd(X.Set3)/sqrt(length(X.Set3)),nanstd(X.Set4)/sqrt(length(X.Set4))],'k','LineStyle','none','LineWidth',1,'CapSize',50)
% bar(3,nanmean(X.Set1),'FaceColor',[.75 .5 .75],'LineStyle','none');
% bar(2,nanmean(X.Set3),'FaceColor',[.74 .75 .75],'LineStyle','none');
% bar(1,nanmean(X.Set4),'FaceColor',[0 .75 0],'LineStyle','none');
simpleboxplot(1,X.Set4,[0 .75 0])
simpleboxplot(2,X.Set3,'k')
ylim([-.1 .2])
yticks(-.1:.1:.2)
xlim([.25 2.75])
xticks(1:3)
xticklabels({'1/2','3/4'})
xlabel('CGRP stim response type')
ylabel('LiCl Period Average (σ)')
text(0.1,0.95,['Rank-sum test: ',num2str(p,3)],'Units','Normalized','FontSize',12,'VerticalAlignment','top')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

idx = ~isnan(data.stats.response.cgrp_train) & ~isnan(data.stats.response.licl);
A = data.stats.response.cgrp_train(idx);
B = data.stats.response.licl(idx);

subplot(1,3,2)
[r,p]=corr(A',B');
hold on
axis square
a = scatter(A,B,64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([-.5 1.5])
ylim([-.5 1.5])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[y_fit] = polyval(p2,x,S);
fitresult = fit(A',B','poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[.9 .9 .9],'LineStyle','none'),
plot(x,y_fit,'k','LineWidth',1)
scatter(A,B,64,'k','filled','MarkerEdgeColor','w')
xlim([x(1) x(end)])
xlim([x(1) x(end)])
text(0.05,0.95,['r = ',num2str(r,3),char(10),'p = ',num2str(p,2)],'Units','Normalized','FontSize',12,'VerticalAlignment','top')
xlabel('CGRP Train Response (σ)')
ylabel('LiCl Period Average (σ)')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

A = data.stats.response.cgrp_period(idx);
subplot(1,3,3)
[r,p]=corr(A',B');
hold on
axis square
a = scatter(A,B,64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([-.5 1.5])
ylim([-.5 1.5])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[y_fit] = polyval(p2,x,S);
fitresult = fit(A',B','poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[.9 .9 .9],'LineStyle','none'),
plot(x,y_fit,'k','LineWidth',1)
scatter(A,B,64,'k','filled','MarkerEdgeColor','w')
xlim([x(1) x(end)])
xlim([x(1) x(end)])
text(0.05,0.95,['r = ',num2str(r,3),char(10),'p = ',num2str(p,2)],'Units','Normalized','FontSize',12,'VerticalAlignment','top')
xlabel('CGRP Period Average (σ)')
ylabel('LiCl Period Average (σ)')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

saveas(gcf,['plots-png/CGRP-LiCl/CGRP-LiCl-plot-3-',regions,'-',num2str(min_amp),'uV-',unittype{end}],'png')
set(gcf,'renderer','Painters')
saveas(gcf,['plots-eps/CGRP-LiCl/CGRP-LiCl-plot-3-',regions,'-',num2str(min_amp),'uV-',unittype{end}],'epsc')
%% save table

X.Activated = sort(data.stats.response.licl(ismember(data.stats.cgrp_cluster,[2,4])),'descend');
X.Other = sort(data.stats.response.licl(ismember(data.stats.cgrp_cluster,[1,3])),'descend');
group = cell(0,0);
out = [];
for i = 1:length(X.Activated)
    out(end+1) = X.Activated(i);
    group{end+1} = 'CGRP-activated';
end
for i = 1:length(X.Other)
    out(end+1) = X.Other(i);
    group{end+1} = 'Other';
end
tbl = table(group',out','VariableNames',{'Group','LiCl response'});
writetable(tbl,'Z:\Chris\matlab\cz\cta-source-data\Fig-ED9j.csv')