%% setup
addpath(genpath('Z:\Chris\matlab\cz\neuropixels-utils\'));
addpath(genpath('Z:\Chris\matlab\heatmap\'));
clearvars; close all; clc
min_amp = 20;
unittype = {'good','mua'};
regions = 'CEA'; % Options: 'CEA','notCEA','amygdala','all'
FDR = 0.05;

window = [0 10]; % Max: [-5 10]
%% load data

library = struct;
library.CEA = {'CEAc','CEAl','CEAm','CEAv','CEAast'};
library.amygdala = {'CEAc','CEAl','CEAm','CEAv','CEAast','BMAa','BMAp','BLAa','BLAp','BLAv','COAa','COAp','IA','MEA','LA','PAA','AAA','PA'};
library.all = cell(0,0);

data_in = cell(0,0);
load(['data/calca_461-pairing.mat'])
data_in{end+1} = data;
library.all = unique([library.all,data_in{end}.meta.location]);

load(['data/calca_462-pairing.mat'])
data_in{end+1} = data;
library.all = unique([library.all,data_in{end}.meta.location]);

load(['data/calca_463-pairing.mat'])
data_in{end+1} = data;
library.all = unique([library.all,data_in{end}.meta.location]);

load(['data/calca_913-pairing.mat'])
data_in{end+1} = data;
library.all = unique([library.all,data_in{end}.meta.location]);

library.notCEA = setdiff(library.all,library.CEA);
library.notAmygdala = setdiff(library.all,library.amygdala);
%% organize data
data = struct;

data.psth.full.novel = [];
data.psth.full.water = [];

data.psth.raw.novel = [];
data.psth.raw.water = [];

data.psth.smooth.novel = [];
data.psth.smooth.water = [];

data.psth.drinking_period = [];
data.psth.licl_period = [];

data.stats.Pvalue = [];
data.stats.FDR = FDR;
data.stats.significant = [];
data.stats.preference = [];
data.stats.response.novel = [];
data.stats.response.water = [];
data.stats.response.licl = [];

data.meta.session = [];
data.meta.cluster_id = [];
data.meta.label = cell(0,0);
data.meta.location = cell(0,0);
data.meta.amplitude = [];
data.meta.firing_rate = [];
data.meta.isi_viol = [];
data.meta.mu_drinking = [];
data.meta.sigma_drinking = [];
data.meta.mu_delay = [];
data.meta.sigma_delay = [];

idx = window*100+1000; idx(1) = idx(1)+1; idx = idx(1):idx(2);

for i = 1:length(data_in)
    for j = 1:length(data_in{i}.stats.Pvalue)
        if data_in{i}.meta.amplitude(j)>=min_amp && ismember(data_in{i}.meta.label{j},unittype) ...
           && ismember(data_in{i}.meta.location{j},library.(regions))% data_in{i}.meta.isi_viol(j)<=0.1
            
            data.stats.Pvalue(end+1) = ranksum(mean(squeeze(data_in{i}.psth.full.novel(j,:,idx)),2),mean(squeeze(data_in{i}.psth.full.water(j,:,idx)),2));
            data.stats.response.novel(end+1) = mean(data_in{i}.psth.raw.novel(j,idx));% - mean(data_in{i}.psth.raw.novel(j,501:600));
            data.stats.response.water(end+1) = mean(data_in{i}.psth.raw.water(j,idx));% - mean(data_in{i}.psth.raw.water(j,501:600));
            data.stats.response.licl(end+1) = mean(data_in{i}.psth.licl_period(j,36:46));
            data.stats.preference(end+1) = data.stats.response.novel(end) - data.stats.response.water(end);
            % data.stats.Pvalue(end+1) = data_in{i}.stats.Pvalue(j);
            % data.stats.preference(end+1) = data_in{i}.stats.preference(j);
            
            data.psth.full.novel(end+1,:,:) = data_in{i}.psth.full.novel(j,:,:);
            data.psth.full.water(end+1,:,:) = data_in{i}.psth.full.water(j,:,:);
            
            data.psth.raw.novel(end+1,:) = data_in{i}.psth.raw.novel(j,:);
            data.psth.raw.water(end+1,:) = data_in{i}.psth.raw.water(j,:);
            
            data.psth.smooth.novel(end+1,:) =  data_in{i}.psth.smooth.novel(j,:);
            data.psth.smooth.water(end+1,:) = data_in{i}.psth.smooth.water(j,:);
            
            data.psth.drinking_period(end+1,:) = data_in{i}.psth.drinking_period(j,:);
            data.psth.licl_period(end+1,:) = data_in{i}.psth.licl_period(j,:);
            
            data.meta.session(end+1) = i;
            data.meta.cluster_id(end+1) = data_in{i}.meta.cluster_id(j);
            data.meta.label{end+1} = data_in{i}.meta.label{j};
            data.meta.location{end+1} = data_in{i}.meta.location{j};
            data.meta.amplitude(end+1) = data_in{i}.meta.amplitude(j);
            data.meta.firing_rate(end+1) = data_in{i}.meta.firing_rate(j);
            data.meta.isi_viol(end+1) = data_in{i}.meta.isi_viol(j);
            data.meta.mu_drinking(end+1) = data_in{i}.meta.mu(j);
            data.meta.sigma_drinking(end+1) = data_in{i}.meta.sigma(j);
            data.meta.mu_delay(end+1) = data_in{i}.meta.mu_delay(j);
            data.meta.sigma_delay(end+1) = data_in{i}.meta.sigma_delay(j);
            
        end
    end
end
data.stats.significant = fdr_bky(data.stats.Pvalue,data.stats.FDR);
disp(['Number of significant neurons: ',num2str(sum(data.stats.significant)),'/',num2str(length(data.stats.significant)),' (',num2str(sum(data.stats.significant)/length(data.stats.significant)*100,'%0.01f'),'%)'])
save(['data/LiCl-pairing-GFP-',regions,'-',num2str(min_amp),'uV-',num2str(window(2)),'sec-',unittype{end},'.mat'],'data','-v7.3')
%% heatmap
close all
idx = find(data.stats.significant & data.stats.preference>0);
[~,idx2] = sort(data.stats.response.novel(idx),'descend'); idx = idx(idx2);
A1 = data.psth.smooth.novel(idx,:);
A2 = data.psth.smooth.water(idx,:);
A4 = data.psth.drinking_period(idx,:);
A5 = data.psth.licl_period(idx,:);

idx = find(data.stats.significant & data.stats.preference<0);
[~,idx2] = sort(data.stats.response.water(idx),'descend'); idx = idx(idx2);
B1 = data.psth.smooth.novel(idx,:);
B2 = data.psth.smooth.water(idx,:);
B4 = data.psth.drinking_period(idx,:);
B5 = data.psth.licl_period(idx,:);

idx = find(~data.stats.significant);
[~,idx2] = sort(mean([data.stats.response.novel(idx);data.stats.response.water(idx)]),'descend'); idx = idx(idx2);
C1 = data.psth.smooth.novel(idx,:);
C2 = data.psth.smooth.water(idx,:);
C4 = data.psth.drinking_period(idx,:);
C5 = data.psth.licl_period(idx,:);

pl1 = [A1(:,501:2000);nan(10,1500);B1(:,501:2000);nan(10,1500);C1(:,501:2000)];
pl2 = [A2(:,501:2000);nan(10,1500);B2(:,501:2000);nan(10,1500);C2(:,501:2000)];
pl = [pl1,nan(size(pl1,1),75),pl2];
cmap = flipud(cbrewer('div','RdBu',1000,'spline')); cmap(cmap<0) = 0;
Clims = [-2 2];

figure('Position', get(0, 'Screensize'))

subplot(2,2,1)
hold on
heatmap(flipud(pl),[],[],[],'Colormap',cmap,'ColorLevels',1000,'MaxColorValue',.5,'MinColorValue',-.5,'NaNColor',[1 1 1]);
plot([501 501]+0.5,[0 size(pl,1)]+0.5,'k-','LineWidth',1)
plot([501 501]+1501+75+0.5,[0 size(pl,1)]+0.5,'k-','LineWidth',1)
nomod = size(C1,1); nopref = size(B1,1); familiarpref = size(A1,1);
yticks([nomod./2 nomod+nopref./2+10 nomod+nopref+familiarpref./2+20]+0.5)
ytickangle(0)
yticklabels({['Non-selective (',num2str(nomod),')'],['Water (',num2str(nopref),')'],['Novel (',num2str(familiarpref),')']})
xticks([750 750+1500+75 250+3000+150])
xticklabels({'Novel','Water','CGRP'})
set(gca,'FontSize',16,'LineWidth',1)
h = colorbar('FontSize',16);
ylabel(h,'Spiking (σ)');
h.Ticks = -1:.5:1;
h.TickLength = 0;
h.Box = 'on';
hold off

pl2 = [A5;nan(10,75);B5;nan(10,75);C5];
subplot(2,2,2)
hold on
heatmap(flipud(pl2),[],[],[],'Colormap',cmap,'ColorLevels',1000,'MaxColorValue',.5,'MinColorValue',-.5,'NaNColor',[1 1 1]);
plot([30 30]+0.5,[0 size(pl2,1)]+0.5,'k-','LineWidth',1)
xticks([0:15:75]+.5)
xticklabels({'-30','-15','0','15','30','45'})
ytickangle(90)
xtickangle(0)
xlabel('Time (min)')
set(gca,'FontSize',16,'LineWidth',1)
hold off

subplot(2,4,5)

% %axis square
% hold on
% t = [-20:1:29]+.5;
% plot([0 0],[-.1 .3],'k','LineWidth',1)
% idx = find(~data.stats.significant);
% fill([t fliplr(t)],[nanmean(data.psth.drinking_period(idx,:))+nanstd(data.psth.drinking_period(idx,:))/sqrt(size(data.psth.drinking_period(idx,:),1)) fliplr(nanmean(data.psth.drinking_period(idx,:))-nanstd(data.psth.drinking_period(idx,:))/sqrt(size(data.psth.drinking_period(idx,:),1)))],[.9 .9 .9],'LineStyle','none');
% a = plot(t,nanmean(data.psth.drinking_period(idx,:)),'k','LineWidth',1);
% idx = find(data.stats.significant & data.stats.preference<0);
% fill([t fliplr(t)],[nanmean(data.psth.drinking_period(idx,:))+nanstd(data.psth.drinking_period(idx,:))/sqrt(size(data.psth.drinking_period(idx,:),1)) fliplr(nanmean(data.psth.drinking_period(idx,:))-nanstd(data.psth.drinking_period(idx,:))/sqrt(size(data.psth.drinking_period(idx,:),1)))],[.9 .9 1],'LineStyle','none');
% b = plot(t,nanmean(data.psth.drinking_period(idx,:)),'b','LineWidth',1);
% idx = find(data.stats.significant & data.stats.preference>0);
% fill([t fliplr(t)],[nanmean(data.psth.drinking_period(idx,:))+nanstd(data.psth.drinking_period(idx,:))/sqrt(size(data.psth.drinking_period(idx,:),1)) fliplr(nanmean(data.psth.drinking_period(idx,:))-nanstd(data.psth.drinking_period(idx,:))/sqrt(size(data.psth.drinking_period(idx,:),1)))],[1 .9 .9],'LineStyle','none');
% c = plot(t,nanmean(data.psth.drinking_period(idx,:)),'r','LineWidth',1);
% plot([-19 -1],[.265 .265],'k','LineWidth',10)
% text(-10,.275,'Drinking','FontSize',16,'vert','bottom','horiz','center')
% plot([1 29],[.265 .265],'k','LineWidth',10)
% text(15,.275,'Delay','FontSize',16,'vert','bottom','horiz','center')
% legend([c,b,a],{'Novel','Water','Non-selective'},'box','off','location','southwest')
% set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
% ylim([-.1 .3])
% yticks(-.1:.1:.3)
% xticks(-20:10:30)
% xlim([-20 30])
% xlabel('Time (min)')
% ylabel('Spiking (σ)')
% hold off
idx = find(data.stats.significant & data.stats.preference>0);
X.Novel = data.stats.response.licl(idx);
idx = find(data.stats.significant & data.stats.preference<0);
X.Water = data.stats.response.licl(idx);
idx = find(~data.stats.significant);
X.Neither = data.stats.response.licl(idx);
hold on
%axis square
% plot([.25 3.75],[0 0],'k','LineWidth',1)
% errorbar([1:3],[nanmean(X.Novel),nanmean(X.Water),nanmean(X.Neither)],[nanstd(X.Novel)/sqrt(length(X.Novel)),nanstd(X.Water)/sqrt(length(X.Water)),nanstd(X.Neither)/sqrt(length(X.Neither))],'k','LineStyle','none','LineWidth',1,'CapSize',50)
% bar(1,nanmean(X.Novel),'FaceColor',[1 0 0],'LineWidth',1);
% bar(2,nanmean(X.Water),'FaceColor',[.5 .5 1],'LineWidth',1);
% bar(3,nanmean(X.Neither),'FaceColor',[.75 .75 .75],'LineWidth',1);
simpleboxplot(1,X.Novel,[1 0 0])
simpleboxplot(2,X.Water,[.5 .5 1])
simpleboxplot(3,X.Neither,[.75 .75 .75])
ylim([-.2 .5])
yticks(-.2:.1:.5)
xticks(1:3)
xlim([.25 3.75])
[p1] = ranksum(X.Novel,X.Water);
[p2] = ranksum(X.Novel,X.Neither);
[p3] = ranksum(X.Water,X.Neither);
p = multicmp([p1 p2 p3],'up',0.05);
text(0.35,0.95,['Novel vs. Water: ',num2str(p(1),2),char(10),...
    'Novel vs. Non-selective:: ',num2str(p(2),2),char(10),...
    'Water vs. Non-selective: ',num2str(p(3),2)],'Units','Normalized','FontSize',12,'VerticalAlignment','top')
ylabel('LiCl Period Average (σ)')
xticklabels({'Novel','Water','Non-selective'})
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

subplot(2,4,6)
A = data.stats.response.licl;
B = data.stats.response.novel-data.stats.response.water;
[r,p]=corr(A',B');
hold on
a = scatter(A,B,64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B',1);
xlim([-1 3])
ylim([-1 1])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[y_fit] = polyval(p2,x,S);
fitresult = fit(A',B','poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[.9 .9 .9],'LineStyle','none'),
plot(x,y_fit,'k','LineWidth',1)
idx = find(~data.stats.significant);
scatter(A(idx),B(idx),128,'k','filled','MarkerEdgeColor','w')
idx = find(data.stats.significant & data.stats.preference<0);
scatter(A(idx),B(idx),128,'b','filled','MarkerEdgeColor','w')
idx = find(data.stats.significant & data.stats.preference>0);
scatter(A(idx),B(idx),128,'r','filled','MarkerEdgeColor','w')
xlim([x(1) x(end)])
ylim([-.5 .5])
yticks(-.5:.25:.5)
xlabel('LiCl Period Average (σ)')
ylabel('Novel–Water Response (σ)')
text(0.05,0.95,['r = ',num2str(r,3),char(10),'p = ',num2str(p,2)],'Units','Normalized','FontSize',12,'VerticalAlignment','top')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')

subplot(2,2,4)
%axis square
hold on
t = [-30:1:44]+.5;
plot([0 0],[-.1 .3],'k','LineWidth',1)
idx = find(~data.stats.significant);
fill([t fliplr(t)],[nanmean(data.psth.licl_period(idx,:))+nanstd(data.psth.licl_period(idx,:))/sqrt(size(data.psth.licl_period(idx,:),1)) fliplr(nanmean(data.psth.licl_period(idx,:))-nanstd(data.psth.licl_period(idx,:))/sqrt(size(data.psth.licl_period(idx,:),1)))],[.9 .9 .9],'LineStyle','none');
plot(t,nanmean(data.psth.licl_period(idx,:)),'k','LineWidth',1)
idx = find(data.stats.significant & data.stats.preference<0);
fill([t fliplr(t)],[nanmean(data.psth.licl_period(idx,:))+nanstd(data.psth.licl_period(idx,:))/sqrt(size(data.psth.licl_period(idx,:),1)) fliplr(nanmean(data.psth.licl_period(idx,:))-nanstd(data.psth.licl_period(idx,:))/sqrt(size(data.psth.licl_period(idx,:),1)))],[.9 .9 1],'LineStyle','none');
plot(t,nanmean(data.psth.licl_period(idx,:)),'b','LineWidth',1)
idx = find(data.stats.significant & data.stats.preference>0);
fill([t fliplr(t)],[nanmean(data.psth.licl_period(idx,:))+nanstd(data.psth.licl_period(idx,:))/sqrt(size(data.psth.licl_period(idx,:),1)) fliplr(nanmean(data.psth.licl_period(idx,:))-nanstd(data.psth.licl_period(idx,:))/sqrt(size(data.psth.licl_period(idx,:),1)))],[1 .9 .9],'LineStyle','none');
plot(t,nanmean(data.psth.licl_period(idx,:)),'r','LineWidth',1)
plot([-29 -1],[.265 .265],'k','LineWidth',10)
text(-15,.275,'Delay','FontSize',16,'vert','bottom','horiz','center')
plot([1 44],[.265 .265],'k','LineWidth',10)
text(22.5,.275,'LiCl','FontSize',16,'vert','bottom','horiz','center')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
ylim([-.1 .3])
yticks(-.1:.1:.3)
xticks(-30:15:45)
xlim([-30 45])
xlabel('Time (min)')
ylabel('Spiking (σ)')
hold off

saveas(gcf,['plots-png/LiCl-pairing-GFP/LiCl-pairing-GFP-plot-1-',regions,'-',num2str(min_amp),'uV-',num2str(window(2)),'sec-',unittype{end}],'png')
set(gcf,'renderer','Painters')
saveas(gcf,['plots-eps/LiCl-pairing-GFP/LiCl-pairing-GFP-plot-1-',regions,'-',num2str(min_amp),'uV-',num2str(window(2)),'sec-',unittype{end}],'epsc')
return
%% save table
idx = find(data.stats.significant & data.stats.preference>0);
X.Novel = sort(data.stats.response.licl(idx),'descend');
idx = find(data.stats.significant & data.stats.preference<0);
X.Water = sort(data.stats.response.licl(idx),'descend');
idx = find(~data.stats.significant);
X.Neither = sort(data.stats.response.licl(idx),'descend');
group = cell(0,0);
out = [];
for i = 1:length(X.Novel)
    out(end+1) = X.Novel(i);
    group{end+1} = 'Flavor-pref';
end
for i = 1:length(X.Water)
    out(end+1) = X.Water(i);
    group{end+1} = 'Water-pref';
end
for i = 1:length(X.Neither)
    out(end+1) = X.Neither(i);
    group{end+1} = 'Non-selective';
end
tbl = table(group',out','VariableNames',{'Group','LiCl response'});
writetable(tbl,'Z:\Chris\matlab\cz\cta-source-data\Fig-3n.csv')
%%
figure('Position', get(0, 'Screensize'))

subplot(1,2,1)
idx = find(data.stats.significant & data.stats.preference>0);
X.Novel = data.stats.response.licl(idx);
idx = find(data.stats.significant & data.stats.preference<0);
X.Water = data.stats.response.licl(idx);
idx = find(~data.stats.significant);
X.Neither = data.stats.response.licl(idx);
hold on
axis square
plot([.25 3.75],[0 0],'k','LineWidth',1)
errorbar([1:3],[nanmean(X.Novel),nanmean(X.Water),nanmean(X.Neither)],[nanstd(X.Novel)/sqrt(length(X.Novel)),nanstd(X.Water)/sqrt(length(X.Water)),nanstd(X.Neither)/sqrt(length(X.Neither))],'k','LineStyle','none','LineWidth',1,'CapSize',50)
bar(1,nanmean(X.Novel),'FaceColor',[1 0 0],'LineWidth',1);
bar(2,nanmean(X.Water),'FaceColor',[.5 .5 1],'LineWidth',1);
bar(3,nanmean(X.Neither),'FaceColor',[.75 .75 .75],'LineWidth',1);
ylim([-.1 .2])
yticks(-.1:.1:.2)
xticks(1:3)
xlim([.25 3.75])
[p1,~,stats1] = ranksum(X.Novel,X.Water,'method','approximate');
[p2,~,stats2] = ranksum(X.Novel,X.Neither,'method','approximate');
[p3,~,stats3] = ranksum(X.Water,X.Neither,'method','approximate');
p = multicmp([p1 p2 p3],'up',0.05);
text(0.65,0.95,['Novel vs. Water: ',num2str(p(1),2),char(10),...
    'Novel vs. Non-selective:: ',num2str(p(2),2),char(10),...
    'Water vs. Non-selective: ',num2str(p(3),2)],'Units','Normalized','FontSize',12,'VerticalAlignment','top')
ylabel('LiCl Period Average (σ)')
xticklabels({'Novel','Water','Non-selective'})
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

A = [data.psth.drinking_period(:,6:end),data.psth.licl_period(:,31:end)];
subplot(1,2,2)
hold on
axis square
t = [0:1:89]+.5;
plot([0 0],[-.1 .3],'k','LineWidth',1)
idx = find(~data.stats.significant);
fill([t fliplr(t)],[nanmean(A(idx,:))+nanstd(A(idx,:))/sqrt(size(A(idx,:),1)) fliplr(nanmean(A(idx,:))-nanstd(A(idx,:))/sqrt(size(A(idx,:),1)))],[.9 .9 .9],'LineStyle','none');
a = plot(t,nanmean(A(idx,:)),'k','LineWidth',1);
idx = find(data.stats.significant & data.stats.preference<0);
fill([t fliplr(t)],[nanmean(A(idx,:))+nanstd(A(idx,:))/sqrt(size(A(idx,:),1)) fliplr(nanmean(A(idx,:))-nanstd(A(idx,:))/sqrt(size(A(idx,:),1)))],[.9 .9 1],'LineStyle','none');
b = plot(t,nanmean(A(idx,:)),'b','LineWidth',1);
idx = find(data.stats.significant & data.stats.preference>0);
fill([t fliplr(t)],[nanmean(A(idx,:))+nanstd(A(idx,:))/sqrt(size(A(idx,:),1)) fliplr(nanmean(A(idx,:))-nanstd(A(idx,:))/sqrt(size(A(idx,:),1)))],[1 .9 .9],'LineStyle','none');
c = plot(t,nanmean(A(idx,:)),'r','LineWidth',1);
plot([1 14],[.275 .275],'k','LineWidth',10)
text(7.5,.28,'Drinking','FontSize',16,'vert','bottom','horiz','center')
plot([16 44],[.275 .275],'k','LineWidth',10)
text(30,.28,'Delay','FontSize',16,'vert','bottom','horiz','center')
plot([46 89],[.275 .275],'k','LineWidth',10)
text(67.5,.28,'LiCl','FontSize',16,'vert','bottom','horiz','center')
legend([c,b,a],{'Novel','Water','Non-selective'},'box','off','location','southwest')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
ylim([-.1 .3])
yticks(-.1:.1:.3)
xticks(0:15:90)
xlim([0 90])
xlabel('Time (min)')
ylabel('Spiking (σ)')
hold off
saveas(gcf,['plots-png/LiCl-pairing-GFP/LiCl-pairing-GFP-plot-2-',regions,'-',num2str(min_amp),'uV-',num2str(window(2)),'sec-',unittype{end}],'png')
set(gcf,'renderer','Painters')
saveas(gcf,['plots-eps/LiCl-pairing-GFP/LiCl-pairing-GFP-plot-2-',regions,'-',num2str(min_amp),'uV-',num2str(window(2)),'sec-',unittype{end}],'epsc')
%% summary plot - drinking responses
figure('Position', get(0, 'Screensize'))

subplot(1,3,1)
axis square
hold on
plot([0 0],[-.1 .4],'k','LineWidth',1)
times = -5:.01:10;
times = times(1:end-1) + mean(diff(times))/2;
idx = find(data.stats.significant & data.stats.preference>0);
A1 = data.psth.smooth.novel(idx,501:end);
B1 = data.psth.smooth.water(idx,501:end);
fill([times fliplr(times)],[mean(B1)+std(B1)/sqrt(size(B1,1)) fliplr(mean(B1)-std(B1)/sqrt(size(B1,1)))],'k','LineStyle','none','FaceAlpha',0.1)
fill([times fliplr(times)],[mean(A1)+std(A1)/sqrt(size(A1,1)) fliplr(mean(A1)-std(A1)/sqrt(size(A1,1)))],'r','LineStyle','none','FaceAlpha',0.1)
c = plot(times,mean(B1),'k','LineWidth',1);
a = plot(times,mean(A1),'r','LineWidth',1);
ylim([-.1 .4])
yticks(-.1:.1:.5)
xticks(-5:5:10)
xlim([-5 10])
xlabel('Time (s)')
ylabel('Spiking (σ)')
title('Novel-preferring Units')
legend([a,c],{'Novel','Water'},'box','off')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

subplot(1,3,2)
axis square
hold on
plot([0 0],[-.1 .4],'k','LineWidth',1)
times = -5:.01:10;
times = times(1:end-1) + mean(diff(times))/2;
idx = find(data.stats.significant & data.stats.preference<0);
A1 = data.psth.smooth.novel(idx,501:end);
B1 = data.psth.smooth.water(idx,501:end);
fill([times fliplr(times)],[mean(B1)+std(B1)/sqrt(size(B1,1)) fliplr(mean(B1)-std(B1)/sqrt(size(B1,1)))],'k','LineStyle','none','FaceAlpha',0.1)
fill([times fliplr(times)],[mean(A1)+std(A1)/sqrt(size(A1,1)) fliplr(mean(A1)-std(A1)/sqrt(size(A1,1)))],'r','LineStyle','none','FaceAlpha',0.1)
c = plot(times,mean(B1),'k','LineWidth',1);
a = plot(times,mean(A1),'r','LineWidth',1);
ylim([-.1 .4])
yticks(-.1:.1:.5)
xticks(-5:5:10)
xlim([-5 10])
xlabel('Time (s)')
ylabel('Spiking (σ)')
title('Water-preferring Units')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

subplot(1,3,3)
axis square
hold on
plot([0 0],[-.1 .4],'k','LineWidth',1)
times = -5:.01:10;
times = times(1:end-1) + mean(diff(times))/2;
idx = find(~data.stats.significant);
A1 = data.psth.smooth.novel(idx,501:end);
B1 = data.psth.smooth.water(idx,501:end);
fill([times fliplr(times)],[mean(B1)+std(B1)/sqrt(size(B1,1)) fliplr(mean(B1)-std(B1)/sqrt(size(B1,1)))],'k','LineStyle','none','FaceAlpha',0.1)
fill([times fliplr(times)],[mean(A1)+std(A1)/sqrt(size(A1,1)) fliplr(mean(A1)-std(A1)/sqrt(size(A1,1)))],'r','LineStyle','none','FaceAlpha',0.1)
c = plot(times,mean(B1),'k','LineWidth',1);
a = plot(times,mean(A1),'r','LineWidth',1);
ylim([-.1 .4])
yticks(-.1:.1:.5)
xticks(-5:5:10)
xlim([-5 10])
xlabel('Time (s)')
ylabel('Spiking (σ)')
title('Non-selective Units')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

saveas(gcf,['plots-png/LiCl-pairing-GFP/LiCl-pairing-GFP-plot-3-',regions,'-',num2str(min_amp),'uV-',num2str(window(2)),'sec-',unittype{end}],'png')
set(gcf,'renderer','Painters')
saveas(gcf,['plots-eps/LiCl-pairing-GFP/LiCl-pairing-GFP-plot-3-',regions,'-',num2str(min_amp),'uV-',num2str(window(2)),'sec-',unittype{end}],'epsc')
%% summary plot - reward responses

figure('Position', get(0, 'Screensize'))

idx = find(data.stats.significant & data.stats.preference>0);
X.Novel = data.stats.response.novel(idx);
idx = find(data.stats.significant & data.stats.preference<0);
X.Water = data.stats.response.novel(idx);
idx = find(~data.stats.significant);
X.Neither = data.stats.response.novel(idx);
subplot(2,2,1)
hold on
axis square
plot([.25 3.75],[0 0],'k','LineWidth',1)
errorbar([1:3],[mean(X.Novel),nanmean(X.Water),nanmean(X.Neither)],[nanstd(X.Novel)/sqrt(length(X.Novel)),nanstd(X.Water)/sqrt(length(X.Water)),nanstd(X.Neither)/sqrt(length(X.Neither))],'k','LineStyle','none','LineWidth',1,'CapSize',25)
bar(1,nanmean(X.Novel),'FaceColor',[1 0 0],'LineWidth',1);
bar(2,nanmean(X.Water),'FaceColor',[.5 .5 1],'LineWidth',1);
bar(3,nanmean(X.Neither),'FaceColor',[.75 .75 .75],'LineWidth',1);
ylim([-.2 .2])
yticks(-.2:.1:.2)
xticks(1:3)
xlim([.25 3.75])
[p1] = ranksum(X.Novel,X.Water);
[p2] = ranksum(X.Novel,X.Neither);
[p3] = ranksum(X.Water,X.Neither);
p = multicmp([p1 p2 p3],'up',0.05);
text(0.05,0.95,['Novel vs. Water: ',num2str(p(1),2),char(10),...
    'Novel vs. Non-selective:: ',num2str(p(2),2),char(10),...
    'Water vs. Non-selective: ',num2str(p(3),2)],'Units','Normalized','FontSize',12,'VerticalAlignment','top')
ylabel('Novel Response (σ)')
title('Novel Response')
xticklabels({'Novel','Water','Non-selective'})
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

idx = find(data.stats.significant & data.stats.preference>0);
X.Novel = data.stats.response.water(idx);
idx = find(data.stats.significant & data.stats.preference<0);
X.Water = data.stats.response.water(idx);
idx = find(~data.stats.significant);
X.Neither = data.stats.response.water(idx);
subplot(2,2,2)
hold on
axis square
plot([.25 3.75],[0 0],'k','LineWidth',1)
errorbar([1:3],[mean(X.Novel),nanmean(X.Water),nanmean(X.Neither)],[nanstd(X.Novel)/sqrt(length(X.Novel)),nanstd(X.Water)/sqrt(length(X.Water)),nanstd(X.Neither)/sqrt(length(X.Neither))],'k','LineStyle','none','LineWidth',1,'CapSize',25)
bar(1,nanmean(X.Novel),'FaceColor',[1 0 0],'LineWidth',1);
bar(2,nanmean(X.Water),'FaceColor',[.5 .5 1],'LineWidth',1);
bar(3,nanmean(X.Neither),'FaceColor',[.75 .75 .75],'LineWidth',1);
ylim([-.2 .2])
yticks(-.2:.1:.2)
xticks(1:3)
xlim([.25 3.75])
[p1] = ranksum(X.Novel,X.Water);
[p2] = ranksum(X.Novel,X.Neither);
[p3] = ranksum(X.Water,X.Neither);
p = multicmp([p1 p2 p3],'up',0.05);
text(0.05,0.95,['Novel vs. Water: ',num2str(p(1),2),char(10),...
    'Novel vs. Non-selective:: ',num2str(p(2),2),char(10),...
    'Water vs. Non-selective: ',num2str(p(3),2)],'Units','Normalized','FontSize',12,'VerticalAlignment','top')
ylabel('Water Response (σ)')
title('Water Response')
xticklabels({'Novel','Water','Non-selective'})
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

idx = find(data.stats.significant & data.stats.preference>0);
X.Novel = data.stats.preference(idx);
idx = find(data.stats.significant & data.stats.preference<0);
X.Water = data.stats.preference(idx);
idx = find(~data.stats.significant);
X.Neither = data.stats.preference(idx);
subplot(2,2,3)
hold on
axis square
plot([.25 3.75],[0 0],'k','LineWidth',1)
errorbar([1:3],[mean(X.Novel),nanmean(X.Water),nanmean(X.Neither)],[nanstd(X.Novel)/sqrt(length(X.Novel)),nanstd(X.Water)/sqrt(length(X.Water)),nanstd(X.Neither)/sqrt(length(X.Neither))],'k','LineStyle','none','LineWidth',1,'CapSize',25)
bar(1,nanmean(X.Novel),'FaceColor',[1 0 0],'LineWidth',1);
bar(2,nanmean(X.Water),'FaceColor',[.5 .5 1],'LineWidth',1);
bar(3,nanmean(X.Neither),'FaceColor',[.75 .75 .75],'LineWidth',1);
ylim([-.2 .2])
yticks(-.2:.1:.2)
xticks(1:3)
xlim([.25 3.75])
[p1] = ranksum(X.Novel,X.Water);
[p2] = ranksum(X.Novel,X.Neither);
[p3] = ranksum(X.Water,X.Neither);
p = multicmp([p1 p2 p3],'up',0.05);
text(0.05,0.95,['Novel vs. Water: ',num2str(p(1),2),char(10),...
    'Novel vs. Non-selective:: ',num2str(p(2),2),char(10),...
    'Water vs. Non-selective: ',num2str(p(3),2)],'Units','Normalized','FontSize',12,'VerticalAlignment','top')
ylabel('Novel–Water Response (σ)')
title('Novel–Water Response')
xticklabels({'Novel','Water','Non-selective'})
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

idx = find(data.stats.significant & data.stats.preference>0);
X.Novel = (length(idx))/length(data.stats.significant)*100;
idx = find(data.stats.significant & data.stats.preference<0);
X.Water = (length(idx))/length(data.stats.significant)*100;
idx = find(~data.stats.significant);
X.Neither = (length(idx))/length(data.stats.significant)*100;
subplot(2,2,4)
hold on
axis square
a=pie([X.Novel,X.Water,X.Neither]);
a(1).FaceColor = [1 0 0];
a(3).FaceColor = [.5 .5 1];
a(5).FaceColor = [.7 .7 .7];
title('Significant Units')
xticks([])
yticks([])
axis off
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

saveas(gcf,['plots-png/LiCl-pairing-GFP/LiCl-pairing-GFP-plot-4-',regions,'-',num2str(min_amp),'uV-',num2str(window(2)),'sec-',unittype{end}],'png')
set(gcf,'renderer','Painters')
saveas(gcf,['plots-eps/LiCl-pairing-GFP/LiCl-pairing-GFP-plot-4-',regions,'-',num2str(min_amp),'uV-',num2str(window(2)),'sec-',unittype{end}],'epsc')
%% scatter plot - LiCl response vs. response change
figure('Position', get(0, 'Screensize'))

A = data.stats.response.licl;
in = ~isnan(A);
B = data.stats.response.novel;
subplot(1,3,1)
[r,p]=corr(A(in)',B(in)');
hold on
axis square
a = scatter(A,B,64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A(in)',B(in)',1);
xlim([-1 3])
ylim([-1 1])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[y_fit] = polyval(p2,x,S);
fitresult = fit(A(in)',B(in)','poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[.9 .9 .9],'LineStyle','none'),
plot(x,y_fit,'k','LineWidth',1)
idx = find(~data.stats.significant);
scatter(A(idx),B(idx),64,'k','filled','MarkerEdgeColor','w')
idx = find(data.stats.significant & data.stats.preference<0);
scatter(A(idx),B(idx),64,'b','filled','MarkerEdgeColor','w')
idx = find(data.stats.significant & data.stats.preference>0);
scatter(A(idx),B(idx),64,'r','filled','MarkerEdgeColor','w')
xlim([x(1) x(end)])
xlabel('LiCl Period Average (σ)')
ylabel('Novel Response (σ)')
text(0.05,0.95,['r = ',num2str(r,3),char(10),'p = ',num2str(p,2)],'Units','Normalized','FontSize',12,'VerticalAlignment','top')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')

A = data.stats.response.licl;
B = data.stats.response.water;
subplot(1,3,2)
[r,p]=corr(A(in)',B(in)');
hold on
axis square
a = scatter(A,B,64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A(in)',B(in)',1);
xlim([-1 3])
ylim([-1 1])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[y_fit] = polyval(p2,x,S);
fitresult = fit(A(in)',B(in)','poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[.9 .9 .9],'LineStyle','none'),
plot(x,y_fit,'k','LineWidth',1)
idx = find(~data.stats.significant);
scatter(A(idx),B(idx),64,'k','filled','MarkerEdgeColor','w')
idx = find(data.stats.significant & data.stats.preference<0);
scatter(A(idx),B(idx),64,'b','filled','MarkerEdgeColor','w')
idx = find(data.stats.significant & data.stats.preference>0);
scatter(A(idx),B(idx),64,'r','filled','MarkerEdgeColor','w')
xlim([x(1) x(end)])
xlabel('LiCl Period Average (σ)')
ylabel('Water Response (σ)')
text(0.05,0.95,['r = ',num2str(r,3),char(10),'p = ',num2str(p,2)],'Units','Normalized','FontSize',12,'VerticalAlignment','top')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')

A = data.stats.response.licl;
B = data.stats.response.novel-data.stats.response.water;
subplot(1,3,3)
[r,p]=corr(A(in)',B(in)');
hold on
axis square
a = scatter(A,B,64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A(in)',B(in)',1);
xlim([-1 3])
ylim([-1 1])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[y_fit] = polyval(p2,x,S);
fitresult = fit(A(in)',B(in)','poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[.9 .9 .9],'LineStyle','none'),
plot(x,y_fit,'k','LineWidth',1)
idx = find(~data.stats.significant);
scatter(A(idx),B(idx),64,'k','filled','MarkerEdgeColor','w')
idx = find(data.stats.significant & data.stats.preference<0);
scatter(A(idx),B(idx),64,'b','filled','MarkerEdgeColor','w')
idx = find(data.stats.significant & data.stats.preference>0);
scatter(A(idx),B(idx),64,'r','filled','MarkerEdgeColor','w')
xlim([x(1) x(end)])
xlabel('LiCl Period Average (σ)')
ylabel('Novel–Water Response (σ)')
text(0.05,0.95,['r = ',num2str(r,3),char(10),'p = ',num2str(p,2)],'Units','Normalized','FontSize',12,'VerticalAlignment','top')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')

saveas(gcf,['plots-png/LiCl-pairing-GFP/LiCl-pairing-GFP-plot-5-',regions,'-',num2str(min_amp),'uV-',num2str(window(2)),'sec-',unittype{end}],'png')
set(gcf,'renderer','Painters')
saveas(gcf,['plots-eps/LiCl-pairing-GFP/LiCl-pairing-GFP-plot-5-',regions,'-',num2str(min_amp),'uV-',num2str(window(2)),'sec-',unittype{end}],'epsc')