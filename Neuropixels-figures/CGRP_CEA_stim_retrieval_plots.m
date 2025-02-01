%% setup
addpath(genpath('Z:\Chris\matlab\cz\neuropixels-utils\'));
addpath(genpath('Z:\Chris\matlab\heatmap\'));
clearvars; close all; clc
min_amp = 20;
unittype = {'good','mua'};
regions = 'CEA'; % Options: 'CEA','notCEA','amygdala','all'
FDR = 0.05;
CGRP_cutoff = 0.1;
window = [0 10]; % Max: [-5 10]
%% load data

library = struct;
library.CEA = {'CEAc','CEAl','CEAm','CEAv','CEAast'};
library.amygdala = {'CEAc','CEAl','CEAm','CEAv','CEAast','BMAa','BMAp','BLAa','BLAp','BLAv','COAa','COAp','IA','MEA','LA','PAA','AAA','PA'};
library.all = cell(0,0);

data_in = cell(0,0);
load(['data/calca_908-retrieval.mat'])
data_in{end+1} = data;
library.all = unique([library.all,data_in{end}.meta.location]);

load(['data/calca_911-retrieval.mat'])
data_in{end+1} = data;
library.all = unique([library.all,data_in{end}.meta.location]);

load(['data/calca_912-retrieval.mat'])
data_in{end+1} = data;
library.all = unique([library.all,data_in{end}.meta.location]);

load(['data/calca_924-retrieval.mat']')
data_in{end+1} = data;
library.all = unique([library.all,data_in{end}.meta.location]);

load(['data/calca_927-retrieval.mat'])
data_in{end+1} = data;
library.all = unique([library.all,data_in{end}.meta.location]);

load(['data/calca_929-retrieval.mat'])
data_in{end+1} = data;
library.all = unique([library.all,data_in{end}.meta.location]);

load(['data/calca_930-retrieval.mat'])
data_in{end+1} = data;
library.all = unique([library.all,data_in{end}.meta.location]);

load(['data/calca_933-retrieval.mat']')
data_in{end+1} = data;
library.all = unique([library.all,data_in{end}.meta.location]);

library.all = unique([library.all,data_in{end}.meta.location]);
library.notCEA = setdiff(library.all,library.CEA);
%% organize data

data = struct;

data.psth.full.novel = [];
data.psth.full.water = [];
data.psth.full.cgrp = cell(0,0);
data.psth.full.novel_retrieval = [];
data.psth.full.water_retrieval = [];

data.psth.raw.novel = [];
data.psth.raw.water = [];
data.psth.raw.cgrp = [];
data.psth.raw.novel_retrieval = [];
data.psth.raw.water_retrieval = [];

data.psth.smooth.novel = [];
data.psth.smooth.water = [];
data.psth.smooth.cgrp = [];
data.psth.smooth.novel_retrieval = [];
data.psth.smooth.water_retrieval = [];

data.psth.cgrp_period = [];

data.stats.Pvalue = [];
data.stats.FDR = FDR;
data.stats.significant = [];
data.stats.preference = [];
data.stats.response.novel = [];
data.stats.response.water = [];
data.stats.response.cgrp = [];
data.stats.response.cgrp_period = [];
data.stats.preference_subtracted = [];
data.stats.response.novel_subtracted = [];
data.stats.response.water_subtracted = [];

data.stats_retrieval.Pvalue = [];
data.stats_retrieval.FDR = FDR;
data.stats_retrieval.significant = [];
data.stats_retrieval.preference = [];
data.stats_retrieval.response.novel = [];
data.stats_retrieval.response.water = [];
data.stats_retrieval.preference_subtracted = [];
data.stats_retrieval.response.novel_subtracted = [];
data.stats_retrieval.response.water_subtracted = [];

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

FLAV = [];
WATR = [];

for i = 1:length(data_in)
    data.psth.full.cgrp{i} = [];
    for j = 1:length(data_in{i}.stats.Pvalue)
        if data_in{i}.meta.amplitude(j)>=min_amp && ismember(data_in{i}.meta.label{j},unittype) ...
           && ismember(data_in{i}.meta.location{j},library.(regions))% data_in{i}.meta.isi_viol(j)<=0.1
            
            data.stats.Pvalue(end+1) = ranksum(mean(squeeze(data_in{i}.psth.full.novel(j,:,idx)),2),mean(squeeze(data_in{i}.psth.full.water(j,:,idx)),2));
            data.stats.response.novel(end+1) = mean(data_in{i}.psth.raw.novel(j,idx));
            data.stats.response.water(end+1) = mean(data_in{i}.psth.raw.water(j,idx));
            data.stats.response.cgrp(end+1) = mean(data_in{i}.psth.raw.cgrp(j,201:500)) - mean(data_in{i}.psth.raw.cgrp(j,101:200));
            data.stats.response.cgrp_period(end+1) = mean(data_in{i}.psth.cgrp_period(j,31:end));
            data.stats.preference(end+1) = data.stats.response.novel(end) - data.stats.response.water(end);
            data.stats.response.novel_subtracted(end+1) = mean(data_in{i}.psth.raw.novel(j,idx)) - mean(data_in{i}.psth.raw.novel(j,01:500));
            data.stats.response.water_subtracted(end+1) = mean(data_in{i}.psth.raw.water(j,idx)) - mean(data_in{i}.psth.raw.water(j,01:500));
            data.stats.preference_subtracted(end+1) = data.stats.response.novel_subtracted(end) - data.stats.response.water_subtracted(end);
            % data.stats.Pvalue(end+1) = data_in{i}.stats.Pvalue(j);
            % data.stats.preference(end+1) = data_in{i}.stats.preference(j);
            
            %%%%%
            FLAV(end+1,:) = mean(data_in{i}.psth.full.novel(j,:,idx),3);
            WATR(end+1,:) = mean(data_in{i}.psth.full.water(j,:,idx),3);
            %%%%%

            data.stats_retrieval.Pvalue(end+1) = ranksum(mean(squeeze(data_in{i}.psth.full.novel_retrieval(j,:,idx)),2),mean(squeeze(data_in{i}.psth.full.water_retrieval(j,:,idx)),2));
            data.stats_retrieval.response.novel(end+1) = mean(data_in{i}.psth.raw.novel_retrieval(j,idx));
            data.stats_retrieval.response.water(end+1) = mean(data_in{i}.psth.raw.water_retrieval(j,idx));
            data.stats_retrieval.preference(end+1) = data.stats_retrieval.response.novel(end) - data.stats_retrieval.response.water(end);
            data.stats_retrieval.response.novel_subtracted(end+1) = mean(data_in{i}.psth.raw.novel_retrieval(j,idx)) - mean(data_in{i}.psth.raw.novel_retrieval(j,01:500));
            data.stats_retrieval.response.water_subtracted(end+1) = mean(data_in{i}.psth.raw.water_retrieval(j,idx)) - mean(data_in{i}.psth.raw.water_retrieval(j,01:500));
            data.stats_retrieval.preference_subtracted(end+1) = data.stats_retrieval.response.novel_subtracted(end) - data.stats_retrieval.response.water_subtracted(end);
            % data.stats_retrieval.Pvalue(end+1) = data_in{i}.stats_retrieval.Pvalue(j);
            % data.stats_retrieval.preference(end+1) = data_in{i}.stats_retrieval.preference(j);
            
            data.psth.full.novel(end+1,:,:) = data_in{i}.psth.full.novel(j,:,:);
            data.psth.full.water(end+1,:,:) = data_in{i}.psth.full.water(j,:,:);
            data.psth.full.cgrp{i}(end+1,:,:) = data_in{i}.psth.full.cgrp(j,:,:);
            data.psth.full.novel_retrieval(end+1,:,:) = data_in{i}.psth.full.novel_retrieval(j,:,:);
            data.psth.full.water_retrieval(end+1,:,:) = data_in{i}.psth.full.water_retrieval(j,:,:);

            data.psth.raw.novel(end+1,:) = data_in{i}.psth.raw.novel(j,:);
            data.psth.raw.water(end+1,:) = data_in{i}.psth.raw.water(j,:);
            data.psth.raw.cgrp(end+1,:) = data_in{i}.psth.raw.cgrp(j,:) - mean(data_in{i}.psth.raw.cgrp(j,101:200));
            data.psth.raw.novel_retrieval(end+1,:) = data_in{i}.psth.raw.novel_retrieval(j,:);
            data.psth.raw.water_retrieval(end+1,:) = data_in{i}.psth.raw.water_retrieval(j,:);
            
            data.psth.smooth.novel(end+1,:) =  data_in{i}.psth.smooth.novel(j,:);
            data.psth.smooth.water(end+1,:) = data_in{i}.psth.smooth.water(j,:);
            data.psth.smooth.cgrp(end+1,:) = data_in{i}.psth.smooth.cgrp(j,:) - mean(data_in{i}.psth.smooth.cgrp(j,101:200));
            data.psth.smooth.novel_retrieval(end+1,:) = data_in{i}.psth.smooth.novel_retrieval(j,:);
            data.psth.smooth.water_retrieval(end+1,:) = data_in{i}.psth.smooth.water_retrieval(j,:);
            
            data.psth.cgrp_period(end+1,:) = data_in{i}.psth.cgrp_period(j,:);

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
data.stats_retrieval.significant = fdr_bky(data.stats_retrieval.Pvalue,data.stats_retrieval.FDR);
disp(['Number of significant neurons: ',num2str(sum(data.stats.significant)),'/',num2str(length(data.stats.significant)),' (',num2str(sum(data.stats.significant)/length(data.stats.significant)*100,'%0.01f'),'%)'])
disp(['Number of significant neurons (retrieval): ',num2str(sum(data.stats_retrieval.significant)),'/',num2str(length(data.stats.significant)),' (',num2str(sum(data.stats_retrieval.significant)/length(data.stats_retrieval.significant)*100,'%0.01f'),'%)'])
save(['data/cgrp-chrmine-retrieval-',regions,'-',num2str(min_amp),'uV-',num2str(window(2)),'sec-',num2str(FDR*100),'pctFDR-',unittype{end},'.mat'],'data','-v7.3')
%% plot results

idx = find(data.stats.significant & data.stats.preference>0);
[~,idx2] = sort(data.stats.response.cgrp_period(idx),'descend'); idx = idx(idx2);
A1 = data.psth.smooth.novel(idx,:);
A2 = data.psth.smooth.water(idx,:);
A3 = data.psth.smooth.cgrp(idx,:);
A4 = data.psth.smooth.novel_retrieval(idx,:);
A5 = data.psth.smooth.water_retrieval(idx,:);
A6 = data.psth.cgrp_period(idx,:);

idx = find(data.stats.significant & data.stats.preference<0);
[~,idx2] = sort(data.stats.response.cgrp_period(idx),'descend'); idx = idx(idx2);
B1 = data.psth.smooth.novel(idx,:);
B2 = data.psth.smooth.water(idx,:);
B3 = data.psth.smooth.cgrp(idx,:);
B4 = data.psth.smooth.novel_retrieval(idx,:);
B5 = data.psth.smooth.water_retrieval(idx,:);
B6 = data.psth.cgrp_period(idx,:);

idx = find(~data.stats.significant);
[~,idx2] = sort(data.stats.response.cgrp_period(idx),'descend'); idx = idx(idx2);
C1 = data.psth.smooth.novel(idx,:);
C2 = data.psth.smooth.water(idx,:);
C3 = data.psth.smooth.cgrp(idx,:);
C4 = data.psth.smooth.novel_retrieval(idx,:);
C5 = data.psth.smooth.water_retrieval(idx,:);
C6 = data.psth.cgrp_period(idx,:);

pl1 = [A1(:,501:2000);nan(18,1500);B1(:,501:2000);nan(18,1500);C1(:,501:2000)];
pl2 = [A2(:,501:2000);nan(18,1500);B2(:,501:2000);nan(18,1500);C2(:,501:2000)];
pl3 = [A3(:,101:600);nan(18,500);B3(:,101:600);nan(18,500);C3(:,101:600)];
pl = [pl1,nan(size(pl1,1),75),pl2,nan(size(pl1,1),75),pl3];
cmap = flipud(cbrewer('div','RdBu',1000,'spline')); cmap(cmap<0) = 0;
Clims = [-2 2];

figure('Position', get(0, 'Screensize'))

subplot(2,3,1)
hold on
heatmap(flipud(pl),[],[],[],'Colormap',cmap,'ColorLevels',1000,'MaxColorValue',.5,'MinColorValue',-.5,'NaNColor',[1 1 1]);
plot([501 501]+0.5,[0 size(pl,1)]+0.5,'k-','LineWidth',1)
plot([501 501]+1501+75+0.5,[0 size(pl,1)]+0.5,'k-','LineWidth',1)
plot([101 101]+3002+150+0.5,[0 size(pl,1)]+0.5,'k-','LineWidth',1)
plot([401 401]+3002+150+0.5,[0 size(pl,1)]+0.5,'k-','LineWidth',1)
nomod = size(C1,1); nopref = size(B1,1); familiarpref = size(A1,1);
yticks([nomod./2 nomod+nopref./2+18 nomod+nopref+familiarpref./2+18*2]+0.5)
ytickangle(90)
yticklabels({['Non-selective (',num2str(nomod),')'],['Water (',num2str(nopref),')'],['Novel (',num2str(familiarpref),')']})
xticks([750 750+1500+75 250+3000+150])
xticklabels({'Novel','Water','CGRP'})
title('Conditioning day')
set(gca,'FontSize',16,'LineWidth',1)
hold off

pl6 = [A6;nan(18,75);B6;nan(18,75);C6];
subplot(2,3,2)
hold on
heatmap(flipud(pl6),[],[],[],'Colormap',cmap,'ColorLevels',1000,'MaxColorValue',.5,'MinColorValue',-.5,'NaNColor',[1 1 1]);
plot([30 30]+0.5,[0 size(pl2,1)]+0.5,'k-','LineWidth',1)
xticks([0:15:75]+.5)
xticklabels({'-30','-15','0','15','30','45'})
ytickangle(90)
xtickangle(0)
xlabel('Time (min)')
title('CGRP stimulation')
set(gca,'FontSize',16,'LineWidth',1)
hold off

pl1 = [A4(:,501:2000);nan(18,1500);B4(:,501:2000);nan(18,1500);C4(:,501:2000)];
pl2 = [A5(:,501:2000);nan(18,1500);B5(:,501:2000);nan(18,1500);C5(:,501:2000)];
pl = [pl1,nan(size(pl1,1),75),pl2];
subplot(2,3,3)
hold on
heatmap(flipud(pl),[],[],[],'Colormap',cmap,'ColorLevels',1000,'MaxColorValue',.5,'MinColorValue',-.5,'NaNColor',[1 1 1]);
plot([501 501]+0.5,[0 size(pl,1)]+0.5,'k-','LineWidth',1)
plot([501 501]+1501+75+0.5,[0 size(pl,1)]+0.5,'k-','LineWidth',1)
xticks([750 750+1500+75])
xticklabels({'Novel','Water'})
title('Retrieval day')
set(gca,'FontSize',16,'LineWidth',1)
hold off

idx = find(data.stats.significant & data.stats.preference>0);
[~,idx2] = sort(data.stats.response.cgrp_period(idx),'descend'); idx = idx(idx2);
A1 = data.psth.smooth.novel_retrieval(idx,:)-data.psth.smooth.novel(idx,:);
A2 = data.psth.smooth.water_retrieval(idx,:)-data.psth.smooth.water(idx,:);
A3 = data.psth.smooth.cgrp(idx,:);
A4 = data.psth.smooth.novel_retrieval(idx,:);
A5 = data.psth.smooth.water_retrieval(idx,:);

idx = find(data.stats.significant & data.stats.preference<0);
[~,idx2] = sort(data.stats.response.cgrp_period(idx),'descend'); idx = idx(idx2);
B1 = data.psth.smooth.novel_retrieval(idx,:)-data.psth.smooth.novel(idx,:);
B2 = data.psth.smooth.water_retrieval(idx,:)-data.psth.smooth.water(idx,:);
B3 = data.psth.smooth.cgrp(idx,:);
B4 = data.psth.smooth.novel_retrieval(idx,:);
B5 = data.psth.smooth.water_retrieval(idx,:);

idx = find(~data.stats.significant);
[~,idx2] = sort(data.stats.response.cgrp_period(idx),'descend'); idx = idx(idx2);
C1 = data.psth.smooth.novel_retrieval(idx,:)-data.psth.smooth.novel(idx,:);
C2 = data.psth.smooth.water_retrieval(idx,:)-data.psth.smooth.water(idx,:);
C3 = data.psth.smooth.cgrp(idx,:);
C4 = data.psth.smooth.novel_retrieval(idx,:);
C5 = data.psth.smooth.water_retrieval(idx,:);

pl1 = [A1(:,501:2000);nan(18,1500);B1(:,501:2000);nan(18,1500);C1(:,501:2000)];
pl2 = [A2(:,501:2000);nan(18,1500);B2(:,501:2000);nan(18,1500);C2(:,501:2000)];
pl3 = [A3(:,101:600);nan(18,500);B3(:,101:600);nan(18,500);C3(:,101:600)];
pl = [pl1,nan(size(pl1,1),75),pl2];
cmap = flipud(cbrewer('div','RdBu',1000,'spline')); cmap(cmap<0) = 0;
Clims = [-2 2];

A = [];
B = [];
for i = unique(data.meta.session)
    A(i) = sum(data.stats.significant & data.stats.preference>0 & data.meta.session==i)./sum(data.meta.session==i)*100;
    B(i) = sum(data.stats_retrieval.significant & data.stats_retrieval.preference>0 & data.meta.session==i)./sum(data.meta.session==i)*100;
end
tbl = table(A'/100',B'/100,'VariableNames',{'Conditioning','Retrieval'});
writetable(tbl,'Z:\Chris\matlab\cz\cta-source-data\Fig-ED10d.csv')

subplot(2,6,7)
hold on
axis square
for i = 1:length(A)
    plot([1 2],[A(i),B(i)],'color',[.85 .85 .85],'linewidth',1)
end
scatter([1 2],[mean(A) mean(B)],64,'k','filled')
plot([1 1],[mean(A)+std(A)./length(A) mean(A)-std(A)./length(A)],'k','linewidth',2)
plot([2 2],[mean(B)+std(B)./length(B) mean(B)-std(B)./length(B)],'k','linewidth',2)
ylim([0 50])
xlim([.5 2.5])
xticks([1 2])
[p,~,stats]=signrank(A,B,'method','exact')
xticklabels({'Conditioning','Retrieval'})
ylabel('Flavor-preferring neurons (%)')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')

idx = find(data.stats.significant & data.stats.preference>0);
[~,idx2] = sort(data.stats.response.cgrp_period(idx),'descend'); idx_high = idx(idx2(1:floor(length(idx2)*CGRP_cutoff)));
[~,idx2] = sort(data.stats.response.cgrp_period(idx),'descend'); idx_low = idx(idx2(floor(length(idx2)*CGRP_cutoff)+1:end));

A = data.stats.response.novel(idx);
B = data.stats_retrieval.response.novel(idx);
p1 = signrank(A,B);
p2 = ranksum(A,B);
[~,p3] = kstest([A,B]);
disp([p1 p2 p3])
%return


subplot(2,6,8)
hold on
axis square
plot([0 0],[-1 1],'k','LineWidth',1)
times = -5:.01:10;
times = times(1:end-1) + mean(diff(times))/2;
idx = find(data.stats.significant & data.stats.preference>0);
A1 = data.psth.smooth.novel(idx,501:end);
A2 = data.psth.smooth.novel_retrieval(idx,501:end);
B1 = data.psth.smooth.water(idx,501:end);
B2 = data.psth.smooth.water_retrieval(idx,501:end);
fill([times fliplr(times)],[mean(A1)+std(A1)/sqrt(size(A1,1)) fliplr(mean(A1)-std(A1)/sqrt(size(A1,1)))],'r','LineStyle','none','FaceAlpha',0.1)
fill([times fliplr(times)],[mean(A2)+std(A2)/sqrt(size(A2,1)) fliplr(mean(A2)-std(A2)/sqrt(size(A2,1)))],'r','LineStyle','none','FaceAlpha',0.1)
a = plot(times,mean(A1),'r','LineWidth',1);
b = plot(times,mean(A2),'r--','LineWidth',1);
ylim([-.2 .8])
yticks(-.2:.2:.8)
xticks(-5:5:10)
xlim([-5 10])
xlabel('Time (s)')
ylabel('Spiking (σ)')
title('All Novel-preferring')
legend([a,b],{'Conditioning day','Retrieval day'})
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')

subplot(2,6,9)
hold on
axis square
plot([0 0],[-1 1],'k','LineWidth',1)
times = -5:.01:10;
times = times(1:end-1) + mean(diff(times))/2;
idx = find(data.stats.significant & data.stats.preference>0);
[~,idx2] = sort(data.stats.response.cgrp_period(idx),'descend'); idx = idx(idx2(1:floor(length(idx2)*CGRP_cutoff)));
A1 = data.psth.smooth.novel(idx,501:end);
A2 = data.psth.smooth.novel_retrieval(idx,501:end);
B1 = data.psth.smooth.water(idx,501:end);
B2 = data.psth.smooth.water_retrieval(idx,501:end);
fill([times fliplr(times)],[mean(A1)+std(A1)/sqrt(size(A1,1)) fliplr(mean(A1)-std(A1)/sqrt(size(A1,1)))],'r','LineStyle','none','FaceAlpha',0.1)
fill([times fliplr(times)],[mean(A2)+std(A2)/sqrt(size(A2,1)) fliplr(mean(A2)-std(A2)/sqrt(size(A2,1)))],'r','LineStyle','none','FaceAlpha',0.1)
a = plot(times,mean(A1),'r','LineWidth',1);
b = plot(times,mean(A2),'r--','LineWidth',1);
ylim([-.2 .8])
yticks(-.2:.2:.8)
xticks(-5:5:10)
xlim([-5 10])
xlabel('Time (s)')
ylabel('Spiking (σ)')
title('High CGRP response')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')

subplot(2,6,10)
hold on
axis square
plot([0 0],[-1 1],'k','LineWidth',1)
times = -5:.01:10;
times = times(1:end-1) + mean(diff(times))/2;
idx = find(data.stats.significant & data.stats.preference>0);
[~,idx2] = sort(data.stats.response.cgrp_period(idx),'descend'); idx = idx(idx2(floor(length(idx2)*CGRP_cutoff)+1:end));
A1 = data.psth.smooth.novel(idx,501:end);
A2 = data.psth.smooth.novel_retrieval(idx,501:end);
B1 = data.psth.smooth.water(idx,501:end);
B2 = data.psth.smooth.water_retrieval(idx,501:end);
fill([times fliplr(times)],[mean(A1)+std(A1)/sqrt(size(A1,1)) fliplr(mean(A1)-std(A1)/sqrt(size(A1,1)))],'r','LineStyle','none','FaceAlpha',0.1)
fill([times fliplr(times)],[mean(A2)+std(A2)/sqrt(size(A2,1)) fliplr(mean(A2)-std(A2)/sqrt(size(A2,1)))],'r','LineStyle','none','FaceAlpha',0.1)
a = plot(times,mean(A1),'r','LineWidth',1);
b = plot(times,mean(A2),'r--','LineWidth',1);
ylim([-.2 .8])
yticks(-.2:.2:.8)
xticks(-5:5:10)
xlim([-5 10])
xlabel('Time (s)')
ylabel('Spiking (σ)')
title('Low CGRP response')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')

idx = find(data.stats.significant & data.stats.preference>0);
A = data.stats.response.cgrp_period(idx);
ranksum(data.stats.response.novel(idx),data.stats_retrieval.response.novel(idx));
B = data.stats_retrieval.response.novel_subtracted(idx)-data.stats.response.novel_subtracted(idx);
subplot(2,6,11)
[r,p]=corr(A',B');
hold on
axis square
a = scatter(A,B,64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([-.5 2])
ylim([-.5 1.5])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[y_fit] = polyval(p2,x,S);
fitresult = fit(A',B','poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[1 .9 .9],'LineStyle','none'),
plot(x,y_fit,'r','LineWidth',1)
scatter(A,B,64,'r','filled','MarkerEdgeColor','w')
xlim([x(1) x(end)])
xlabel('CGRP response (σ)') 
ylabel('ΔNovel flavor (σ)')
title('All Novel-preferring')
text(0.05,0.95,['r = ',num2str(r,3),char(10),'p = ',num2str(p,2)],'Units','Normalized','FontSize',12,'VerticalAlignment','top')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

B = data.stats_retrieval.preference_subtracted(idx)-data.stats.preference_subtracted(idx);
subplot(2,6,12)
[r,p]=corr(A',B');
hold on
axis square
a = scatter(A,B,64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([-.5 2])
ylim([-.5 1.5])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[y_fit] = polyval(p2,x,S);
fitresult = fit(A',B','poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[1 .9 .9],'LineStyle','none'),
plot(x,y_fit,'r','LineWidth',1)
scatter(A,B,64,'r','filled','MarkerEdgeColor','w')
xlim([x(1) x(end)])
xlabel('CGRP response (σ)') 
ylabel('ΔSelectivity (σ)')
title('All Novel-preferring')
text(0.05,0.95,['r = ',num2str(r,3),char(10),'p = ',num2str(p,2)],'Units','Normalized','FontSize',12,'VerticalAlignment','top')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

saveas(gcf,['plots-png/CGRP-chrmine-retrieval/CGRP-chrmine-retrieval-',regions,'-',num2str(min_amp),'uV-',num2str(window(2)),'sec-',num2str(FDR*100),'pctFDR-',unittype{end}],'png')
set(gcf,'renderer','Painters')
saveas(gcf,['plots-eps/CGRP-chrmine-retrieval/CGRP-chrmine-retrieval-',regions,'-',num2str(min_amp),'uV-',num2str(window(2)),'sec-',num2str(FDR*100),'pctFDR-',unittype{end}],'epsc')
%% save table
idx = find(data.stats.significant & data.stats.preference>0);
A = data.stats.response.cgrp_period(idx);
B = data.stats_retrieval.response.novel_subtracted(idx)-data.stats.response.novel_subtracted(idx);
C = data.stats_retrieval.preference_subtracted(idx)-data.stats.preference_subtracted(idx);
X.Novel = [A; B; C];

idx = find(data.stats.significant & data.stats.preference<0);
A = data.stats.response.cgrp_period(idx);
B = data.stats_retrieval.response.novel_subtracted(idx)-data.stats.response.novel_subtracted(idx);
C = data.stats_retrieval.preference_subtracted(idx)-data.stats.preference_subtracted(idx);
X.Water = [A; B; C];

idx = find(~data.stats.significant);
A = data.stats.response.cgrp_period(idx);
B = data.stats_retrieval.response.novel_subtracted(idx)-data.stats.response.novel_subtracted(idx);
C = data.stats_retrieval.preference_subtracted(idx)-data.stats.preference_subtracted(idx);
X.Neither = [A; B; C];

group = cell(0,0);
out = [];
for i = 1:length(X.Novel)
    out(end+1,:) = X.Novel(:,i);
    group{end+1} = 'Flavor-pref';
end
% for i = 1:length(X.Water)
%     out(end+1,:) = X.Water(:,i);
%     group{end+1} = 'Water-pref';
% end
% for i = 1:length(X.Neither)
%     out(end+1,:) = X.Neither(:,i);
%     group{end+1} = 'Non-selective';
% end
tbl = table(group',out(:,1),out(:,2),out(:,3),'VariableNames',{'Group','CGRP-CEA response','Delta Novel response','Delta Selectivity'});
writetable(tbl,'Z:\Chris\matlab\cz\cta-source-data\Fig-4e.csv')
%% summary plot - Novel-preferring units by CGRP response
figure('Position', get(0, 'Screensize'))

subplot(1,3,1)
hold on
axis square
plot([0 0],[-1 1],'k','LineWidth',1)
times = -5:.01:10;
times = times(1:end-1) + mean(diff(times))/2;
idx = find(data.stats.significant & data.stats.preference>0);
[~,idx2] = sort(data.stats.response.cgrp_period(idx),'descend'); idx = idx(idx2(1:floor(length(idx2)*CGRP_cutoff)));
A1 = data.psth.smooth.novel(idx,501:end);
A2 = data.psth.smooth.novel_retrieval(idx,501:end);
B1 = data.psth.smooth.water(idx,501:end);
B2 = data.psth.smooth.water_retrieval(idx,501:end);
fill([times fliplr(times)],[mean(B1)+std(B1)/sqrt(size(B1,1)) fliplr(mean(B1)-std(B1)/sqrt(size(B1,1)))],'k','LineStyle','none','FaceAlpha',0.1)
fill([times fliplr(times)],[mean(B2)+std(B2)/sqrt(size(B2,1)) fliplr(mean(B2)-std(B2)/sqrt(size(B2,1)))],'k','LineStyle','none','FaceAlpha',0.1)
fill([times fliplr(times)],[mean(A1)+std(A1)/sqrt(size(A1,1)) fliplr(mean(A1)-std(A1)/sqrt(size(A1,1)))],'r','LineStyle','none','FaceAlpha',0.1)
fill([times fliplr(times)],[mean(A2)+std(A2)/sqrt(size(A2,1)) fliplr(mean(A2)-std(A2)/sqrt(size(A2,1)))],'r','LineStyle','none','FaceAlpha',0.1)
c = plot(times,mean(B1),'k','LineWidth',1);
d = plot(times,mean(B2),'k--','LineWidth',1);
a = plot(times,mean(A1),'r','LineWidth',1);
b = plot(times,mean(A2),'r--','LineWidth',1);
ylim([-.2 1])
yticks(-.2:.2:1)
xticks(-5:5:10)
xlim([-5 10])
xlabel('Time (s)')
ylabel('Spiking (σ)')
title('Novel-preferring Units (High-CGRP)')
legend([a,b,c,d],{'Novel/Pairing','Novel/Retrieval','Water/Pairing','Water/Retrieval'},'box','off')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')

subplot(1,3,2)
hold on
axis square
plot([0 0],[-1 1],'k','LineWidth',1)
times = -5:.01:10;
times = times(1:end-1) + mean(diff(times))/2;
idx = find(data.stats.significant & data.stats.preference>0);
[~,idx2] = sort(data.stats.response.cgrp_period(idx),'descend'); idx = idx(idx2(floor(length(idx2)*CGRP_cutoff)+1:end));
A1 = data.psth.smooth.novel(idx,501:end);
A2 = data.psth.smooth.novel_retrieval(idx,501:end);
B1 = data.psth.smooth.water(idx,501:end);
B2 = data.psth.smooth.water_retrieval(idx,501:end);
fill([times fliplr(times)],[mean(B1)+std(B1)/sqrt(size(B1,1)) fliplr(mean(B1)-std(B1)/sqrt(size(B1,1)))],'k','LineStyle','none','FaceAlpha',0.1)
fill([times fliplr(times)],[mean(B2)+std(B2)/sqrt(size(B2,1)) fliplr(mean(B2)-std(B2)/sqrt(size(B2,1)))],'k','LineStyle','none','FaceAlpha',0.1)
fill([times fliplr(times)],[mean(A1)+std(A1)/sqrt(size(A1,1)) fliplr(mean(A1)-std(A1)/sqrt(size(A1,1)))],'r','LineStyle','none','FaceAlpha',0.1)
fill([times fliplr(times)],[mean(A2)+std(A2)/sqrt(size(A2,1)) fliplr(mean(A2)-std(A2)/sqrt(size(A2,1)))],'r','LineStyle','none','FaceAlpha',0.1)
c = plot(times,mean(B1),'k','LineWidth',1);
d = plot(times,mean(B2),'k--','LineWidth',1);
a = plot(times,mean(A1),'r','LineWidth',1);
b = plot(times,mean(A2),'r--','LineWidth',1);
ylim([-.2 1])
yticks(-.2:.2:1)
xticks(-5:5:10)
xlim([-5 10])
xlabel('Time (s)')
ylabel('Spiking (σ)')
legend([a,b,c,d],{'Novel/Pairing','Novel/Retrieval','Water/Pairing','Water/Retrieval'},'box','off')
title('Novel-preferring Units (Low-CGRP)')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')

subplot(1,3,3)
hold on
axis square
for i = 0:.1:2.9
    plot([i i],[.980 .995],'Color',[0 .5 1],'LineWidth',2)
end
times = -1:.01:4;
times = times(1:end-1) + mean(diff(times))/2;
idx = find(data.stats.significant & data.stats.preference>0);
[~,idx2] = sort(data.stats.response.cgrp_period(idx),'descend'); idx = idx(idx2(floor(length(idx2)*CGRP_cutoff)+1:end));
A1 = data.psth.raw.cgrp(idx,101:600);
fill([times fliplr(times)],[mean(A1)+std(A1)/sqrt(size(A1,1)) fliplr(mean(A1)-std(A1)/sqrt(size(A1,1)))],'k','LineStyle','none','FaceAlpha',0.15)
b = plot(times,mean(A1),'k','LineWidth',1);
idx = find(data.stats.significant & data.stats.preference>0);
[~,idx2] = sort(data.stats.response.cgrp_period(idx),'descend'); idx = idx(idx2(1:floor(length(idx2)*CGRP_cutoff)));
A1 = data.psth.raw.cgrp(idx,101:600);
fill([times fliplr(times)],[mean(A1)+std(A1)/sqrt(size(A1,1)) fliplr(mean(A1)-std(A1)/sqrt(size(A1,1)))],'r','LineStyle','none','FaceAlpha',0.15)
a = plot(times,mean(A1),'r','LineWidth',1);
ylim([-.2 1])
yticks(-.2:.2:1)
legend([a,b],{'High-CGRP','Low-CGRP'},'box','off')
xticks(-1:1:4)
xlim([-1 4])
xlabel('Time (s)')
ylabel('Spiking (σ)')
title('CGRP Train Responses')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')

saveas(gcf,['plots-png/CGRP-chrmine-retrieval/CGRP-chrmine-retrieval-2-',regions,'-',num2str(min_amp),'uV-',num2str(window(2)),'sec-',num2str(FDR*100),'pctFDR-',unittype{end}],'png')
set(gcf,'renderer','Painters')
saveas(gcf,['plots-eps/CGRP-chrmine-retrieval/CGRP-chrmine-retrieval-2-',regions,'-',num2str(min_amp),'uV-',num2str(window(2)),'sec-',num2str(FDR*100),'pctFDR-',unittype{end}],'epsc')
%% PCA plots - population-wide

cmap = struct;

idx1 = find(data.stats.significant & data.stats.preference>0);
idx2 = find(data.stats.significant & data.stats.preference<0);
idx3 = find(~data.stats.significant);
idx = sort([idx1,idx2]);
%[~,idx_cgrp] = sort(data.stats.response.cgrp(idx),'descend'); idx = idx_cgrp(1:floor(length(idx)*CGRP_cutoff));

PCAdata.novel = data.psth.smooth.novel(idx,501:end) - mean([data.psth.smooth.novel(idx,501:600),data.psth.smooth.water(idx,501:600)],2);
PCAdata.water = data.psth.smooth.water(idx,501:end) - mean([data.psth.smooth.novel(idx,501:600),data.psth.smooth.water(idx,501:600)],2);
PCAdata.cgrp  = data.psth.smooth.cgrp(idx,101:600)  - mean(data.psth.smooth.cgrp(idx,101:200),2);
PCAdata.novel_retrieval = data.psth.smooth.novel_retrieval(idx,501:end) - mean([data.psth.smooth.novel(idx,501:600),data.psth.smooth.water(idx,501:600)],2);
PCAdata.water_retrieval = data.psth.smooth.water_retrieval(idx,501:end) - mean([data.psth.smooth.novel(idx,501:600),data.psth.smooth.water(idx,501:600)],2);

% PCAdata2 = struct;
% for i = 1:size(PCAdata.novel,1)
%     for j = 1:size(PCAdata.novel,2)/10
%         PCAdata2.novel(i,j) = mean(PCAdata.novel(i,(j-1)*10+1:j*10));
%     end
% end
% for i = 1:size(PCAdata.water,1)
%     for j = 1:size(PCAdata.water,2)/10
%         PCAdata2.water(i,j) = mean(PCAdata.water(i,(j-1)*10+1:j*10));
%     end
% end
% for i = 1:size(PCAdata.cgrp,1)
%     for j = 1:size(PCAdata.cgrp,2)/10
%         PCAdata2.cgrp(i,j) = mean(PCAdata.cgrp(i,(j-1)*10+1:j*10));
%     end
% end
% for i = 1:size(PCAdata.novel_retrieval,1)
%     for j = 1:size(PCAdata.novel_retrieval,2)/10
%         PCAdata2.novel_retrieval(i,j) = mean(PCAdata.novel_retrieval(i,(j-1)*10+1:j*10));
%     end
% end
% for i = 1:size(PCAdata.water_retrieval,1)
%     for j = 1:size(PCAdata.water_retrieval,2)/10
%         PCAdata2.water_retrieval(i,j) = mean(PCAdata.water_retrieval(i,(j-1)*10+1:j*10));
%     end
% end
% PCAdata = PCAdata2;

t=max(abs([PCAdata.novel,PCAdata.water]),[],2);
PCAdata.novel = PCAdata.novel./t;
PCAdata.water = PCAdata.water./t;
%t=max(abs([PCAdata.novel_retrieval,PCAdata.water_retrieval]),[],2);
PCAdata.novel_retrieval = PCAdata.novel_retrieval./t;
PCAdata.water_retrieval = PCAdata.water_retrieval./t;
t=max(abs([PCAdata.cgrp]),[],2);
PCAdata.cgrp = PCAdata.cgrp./t;

cmap.novel = cbrewer('seq','Reds',size(PCAdata.novel,2)*1.5,'spline'); cmap.novel = cmap.novel(size(PCAdata.novel,2)/2+1:size(PCAdata.novel,2)*1.5,:);
cmap.water = cbrewer('seq','Blues',size(PCAdata.water,2)*1.5,'spline'); cmap.water = cmap.water(size(PCAdata.water,2)/2+1:size(PCAdata.water,2)*1.5,:);
cmap.cgrp = cbrewer('seq','Greens',size(PCAdata.cgrp,2)*1.5,'spline'); cmap.cgrp = cmap.cgrp(size(PCAdata.cgrp,2)/2+1:size(PCAdata.cgrp,2)*1.5,:); cmap.cgrp(cmap.cgrp<0) = 0;

x = [PCAdata.novel(:,501:1000) PCAdata.water(:,501:1000)]; x = x-mean(x);
% x = [PCAdata.novel(:,51:100) PCAdata.water(:,51:100)]; x = x-mean(x);
[coeff,score,latent,tsquared,explained,mu] = pca(x,'Centered',false);

figure('Position', get(0, 'Screensize'))

water = struct; novel = struct;
water_retrieval = struct; novel_retrieval = struct;

ax1 = subplot(1,3,1);
axis square
hold on
pc1 = []; pc2 = [];
vec = struct;
for i = 1:size(PCAdata.water,2)
    x = PCAdata.water(:,i); x = x-mean(x);
    y = x'*score; pc1(i) = y(1); pc2(i) = y(2);
    scatter(pc1(i),pc2(i),64,cmap.water(i,:),'filled')
    if i>1
        plot([pc1(i-1) pc1(i)],[pc2(i-1) pc2(i)],'Color',cmap.water(i,:),'LineWidth',2)
    end
end
water.pc1 = pc1;
water.pc2 = pc2;
[~,idx]=sort(pc1,'descend');
vec.water(1,1) = mean(pc1(1:100));
vec.water(1,2) = mean(pc2(1:100));
vec.water(2,1) = mean(pc1(idx(1:100)));
vec.water(2,2) = mean(pc2(idx(1:100)));
for i = 1:size(PCAdata.novel,2)
    x = PCAdata.novel(:,i); x = x-mean(x);
    y = x'*score; pc1(i) = y(1); pc2(i) = y(2);
    scatter(pc1(i),pc2(i),64,cmap.novel(i,:),'filled')
    if i>1
        plot([pc1(i-1) pc1(i)],[pc2(i-1) pc2(i)],'Color',cmap.novel(i,:),'LineWidth',2)
    end
end
novel.pc1 = pc1;
novel.pc2 = pc2;
[~,idx]=sort(pc1,'descend');
vec.novel(1,1) = mean(pc1(1:100));
vec.novel(1,2) = mean(pc2(1:100));
vec.novel(2,1) = mean(pc1(idx(1:100)));
vec.novel(2,2) = mean(pc2(idx(1:100)));
pc1 = []; pc2 = [];
for i = 1:size(PCAdata.cgrp,2)
    x = PCAdata.cgrp(:,i); x = x-mean(x);
    y = x'*score; pc1(i) = y(1); pc2(i) = y(2);
    scatter(pc1(i),pc2(i),64,cmap.cgrp(i,:),'filled')
    if i>1
        plot([pc1(i-1) pc1(i)],[pc2(i-1) pc2(i)],'Color',cmap.cgrp(i,:),'LineWidth',2)
    end
end
[~,idx]=sort(pc1,'descend');
vec.cgrp(1,1) = mean(pc1(1:100));
vec.cgrp(1,2) = mean(pc2(1:100));
vec.cgrp(2,1) = mean(pc1(idx(1:100)));
vec.cgrp(2,2) = mean(pc2(idx(1:100)));
x1 = xlim;
y1 = ylim;
xticks([])
yticks([])
xlabel('PC1')
ylabel('PC2')
set(gca,'FontSize',16,'LineWidth',1)
title('Pairing','FontWeight','normal')
hold off

ax2 = subplot(1,3,2);
axis square
hold on
pc1 = []; pc2 = [];
for i = 1:size(PCAdata.water_retrieval,2)
    x = PCAdata.water_retrieval(:,i); x = x-mean(x);
    y = x'*score; pc1(i) = y(1); pc2(i) = y(2);
    scatter(pc1(i),pc2(i),64,cmap.water(i,:),'filled')
    if i>1
        plot([pc1(i-1) pc1(i)],[pc2(i-1) pc2(i)],'Color',cmap.water(i,:),'LineWidth',2)
    end
end
water_retrieval.pc1 = pc1;
water_retrieval.pc2 = pc2;
[~,idx]=sort(pc1,'descend');
vec.water_retrieval(1,1) = mean(pc1(1:100));
vec.water_retrieval(1,2) = mean(pc2(1:100));
vec.water_retrieval(2,1) = mean(pc1(idx(1:100)));
vec.water_retrieval(2,2) = mean(pc2(idx(1:100)));
for i = 1:size(PCAdata.novel_retrieval,2)
    x = PCAdata.novel_retrieval(:,i); x = x-mean(x);
    y = x'*score; pc1(i) = y(1); pc2(i) = y(2);
    scatter(pc1(i),pc2(i),64,cmap.novel(i,:),'filled')
    if i>1
        plot([pc1(i-1) pc1(i)],[pc2(i-1) pc2(i)],'Color',cmap.novel(i,:),'LineWidth',2)
    end
end
novel_retrieval.pc1 = pc1;
novel_retrieval.pc2 = pc2;
[~,idx]=sort(pc1,'descend');
vec.novel_retrieval(1,1) = mean(pc1(1:100));
vec.novel_retrieval(1,2) = mean(pc2(1:100));
vec.novel_retrieval(2,1) = mean(pc1(idx(1:100)));
vec.novel_retrieval(2,2) = mean(pc2(idx(1:100)));
y2 = ylim;
x2 = xlim;
xticks([])
yticks([])
xlabel('PC1')
ylabel('PC2')
set(gca,'FontSize',16,'LineWidth',1)
title('Retrieval','FontWeight','normal')
hold off

ax3 = subplot(1,3,3);
hold on
axis square
plot(vec.novel(:,1), vec.novel(:,2),'Color',cmap.novel(75,:),'LineWidth',4)
plot(vec.novel_retrieval(:,1), vec.novel_retrieval(:,2),'Color',cmap.novel(75,:),'LineWidth',4,'LineStyle','--')
plot(vec.water(:,1), vec.water(:,2),'Color',cmap.water(75,:),'LineWidth',4)
plot(vec.water_retrieval(:,1), vec.water_retrieval(:,2),'Color',cmap.water(75,:),'LineWidth',4,'LineStyle','--')
plot(vec.cgrp(:,1), vec.cgrp(:,2),'Color',cmap.cgrp(25,:),'LineWidth',4)
scatter(vec.novel(2,1), vec.novel(2,2),128,cmap.novel(75,:),'filled')
scatter(vec.novel_retrieval(2,1), vec.novel_retrieval(2,2),128,cmap.novel(75,:),'filled')
scatter(vec.water(2,1), vec.water(2,2),128,cmap.water(75,:),'filled')
scatter(vec.water_retrieval(2,1), vec.water_retrieval(2,2),128,cmap.water(75,:),'filled')
scatter(vec.cgrp(2,1), vec.cgrp(2,2),128,cmap.cgrp(25,:),'filled')
xticks([])
yticks([])
xlabel('PC1')
ylabel('PC2')
legend({'Novel/Pairing','Novel/Retrieval','Water/Pairing','Water/Retrieval','CGRP'},'location','northwest','box','off')
set(gca,'FontSize',16,'LineWidth',1)
title('Comparison','FontWeight','normal')
hold off

x = [min([x1,x2]) max([x1 x2])];
y = [min([y1,y2]) max([y1 y2])];
set(ax1,'xlim',x,'ylim',y)
set(ax2,'xlim',x,'ylim',y)
set(ax3,'xlim',x,'ylim',y)

saveas(gcf,['plots-png/CGRP-chrmine-retrieval/CGRP-chrmine-retrieval-3-',regions,'-',num2str(min_amp),'uV-',num2str(window(2)),'sec-',num2str(FDR*100),'pctFDR-',unittype{end}],'png')
set(gcf,'renderer','Painters')
saveas(gcf,['plots-eps/CGRP-chrmine-retrieval/CGRP-chrmine-retrieval-3-',regions,'-',num2str(min_amp),'uV-',num2str(window(2)),'sec-',num2str(FDR*100),'pctFDR-',unittype{end}],'epsc')
%%

cmap = struct;
cmap.before = cbrewer('seq','Greys',size(PCAdata.novel,2)*1.5,'spline'); cmap.before = cmap.before(size(PCAdata.novel,2)/2+1:size(PCAdata.novel,2)*1.5,:);
cmap.novel = cbrewer('seq','Reds',size(PCAdata.novel,2)*1.5,'spline'); cmap.novel = cmap.novel(size(PCAdata.novel,2)/2+1:size(PCAdata.novel,2)*1.5,:);
cmap.water = cbrewer('seq','Blues',size(PCAdata.water,2)*1.5,'spline'); cmap.water = cmap.water(size(PCAdata.water,2)/2+1:size(PCAdata.water,2)*1.5,:);

figure('Position', get(0, 'Screensize'))

% subplot(2,2,2)
% hold on
% times = -5:.01:10;
% times = times(1:end-1) + median(diff(times))/2;
% for i = 1:length(times)
%     scatter(times(i),water.pc1(i),128,cmap.water(i,:),'filled')
% end
% for i = 1:length(times)
%     scatter(times(i),novel.pc1(i),128,cmap.novel(i,:),'filled')
% end
% y = ylim;
% plot([0 0],y,'k','LineWidth',1)
% ylim(y)
% yticks([])
% xticks(-5:5:15)
% xlabel('Time (s)')
% ylabel('PC1')
% xlim([-5 10])
% set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
% hold off

subplot(1,4,1)
axis square
hold on
times = -5:.01:10;
times = times(1:end-1) + median(diff(times))/2;
for i = 1:length(times)
    scatter(times(i),water.pc2(i),64,cmap.water(i,:),'filled')
end
for i = 1:length(times)
    scatter(times(i),novel.pc2(i),64,cmap.novel(i,:),'filled')
end
y = ylim;
plot([0 0],y,'k','LineWidth',1)
ylim(y)
yticks([])
xticks(-5:5:15)
xlabel('Time (s)')
ylabel('PC2')
title('Pairing')
xlim([-5 10])
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

subplot(1,4,2)
axis square
hold on
times = -5:.01:10;
times = times(1:end-1) + median(diff(times))/2;
for i = 1:length(times)
    scatter(times(i),water_retrieval.pc2(i),64,cmap.water(i,:),'filled')
end
for i = 1:length(times)
    scatter(times(i),novel_retrieval.pc2(i),64,cmap.novel(i,:),'filled')
end
plot([0 0],y,'k','LineWidth',1)
ylim(y)
yticks([])
xticks(-5:5:15)
xlabel('Time (s)')
ylabel('PC2')
title('Retrieval')
xlim([-5 10])
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

subplot(1,4,3)
axis square
hold on
times = -5:.01:10;
times = times(1:end-1) + median(diff(times))/2;
for i = 1:length(times)
    scatter(times(i),novel.pc2(i),64,cmap.before(i,:),'filled')
end
for i = 1:length(times)
    scatter(times(i),novel_retrieval.pc2(i),64,cmap.novel(i,:),'filled')
end
y = ylim;
plot([0 0],y,'k','LineWidth',1)
ylim(y)
yticks([])
xticks(-5:5:15)
xlabel('Time (s)')
ylabel('PC2')
title('Pairing vs. Retrieval (Flavor)')
xlim([-5 10])
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

subplot(1,4,4)
axis square
hold on
times = -5:.01:10;
times = times(1:end-1) + median(diff(times))/2;
for i = 1:length(times)
    scatter(times(i),novel.pc2(i)-water.pc2(i),64,cmap.before(i,:),'filled')
end
for i = 1:length(times)
    scatter(times(i),novel_retrieval.pc2(i)-water_retrieval.pc2(i),64,cmap.novel(i,:),'filled')
end
y = ylim;
plot([0 0],y,'k','LineWidth',1)
ylim(y)
yticks([])
xticks(-5:5:15)
xlabel('Time (s)')
ylabel('PC2')
title('Pairing vs. Retrieval (Flavor-Water)')
xlim([-5 10])
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

saveas(gcf,['plots-png/CGRP-chrmine-retrieval/CGRP-chrmine-retrieval-4-',regions,'-',num2str(min_amp),'uV-',num2str(window(2)),'sec-',num2str(FDR*100),'pctFDR-',unittype{end}],'png')
set(gcf,'renderer','Painters')
saveas(gcf,['plots-eps/CGRP-chrmine-retrieval/CGRP-chrmine-retrieval-4-',regions,'-',num2str(min_amp),'uV-',num2str(window(2)),'sec-',num2str(FDR*100),'pctFDR-',unittype{end}],'epsc')