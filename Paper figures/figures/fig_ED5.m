function fig_ED5(data_path)

disp('Generating panels for Extended Data Figure 5...')
%% Fig ED5b

openfig([data_path,'\Neuropixels\CGRP stim acute\cgrp-stim-acute-trajectories.fig']);
set(gcf,'Position', get(0, 'Screensize'))
%% Fig ED5c

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 5c','FontWeight','bold')
load([data_path,'\Neuropixels\CGRP stim acute\cgrp-stim-acute-summary.mat'],'data')

counter = 1;
binedges = -(60*5):0.01:(60*10);
edges = -1:.01:2;
times = [-(60*5):0.01:(60*10)]+.005; times = times(1:end-1);
PSTH = [];
locdata = struct;
for i = 1:length(data)
    for j = 1:length(data(i).spikes)
        t = histcounts(data(i).spikes{j}-data(i).laser(1),binedges);
        idx1 = find(times<=-1,1,'last');
        idx2 = find(times<=data(i).laser(end)-data(i).laser(1)+2,1,'last');
        mu = mean(t(idx1:idx2));
        sigma = std(t(idx1:idx2));
        for k = 1:60
            spks = data(i).spikes{j} - data(i).laser((k-1)*10+1);
            PSTH(counter,k,:) = (histcounts(spks,edges)-mu)./sigma;
        end
        if ismember(data(i).location{j},{'PA','AAA','LA'})
            data(i).location{j} = 'Other';
        end
        locdata.region{counter} = data(i).location{j};
        locdata.x(counter) = data(i).x(j);
        locdata.y(counter) = data(i).y(j);
        locdata.z(counter) = data(i).z(j);
        counter = counter + 1;
    end
end

X = squeeze(mean(PSTH,2)); X = X - mean(X(:,1:100),2);
clusters_n = 4;
rng(12345)
GMModel=fitgmdist(X(:,101:200),clusters_n,'CovarianceType','diagonal','RegularizationValue',1e-5,'SharedCovariance',false,'Replicates',100);
clusterX = cluster(GMModel,X(:,101:200));

clusters = cell(0,0);
for i = 1:clusters_n
    clusters{i} = X(clusterX==i,:);
end
[~,idx] = sort(cellfun(@(x) mean(mean(x(:,101:200),2)),clusters),'ascend');
clusters_idx = idx;
clusters = cell(0,0);
for i = 1:clusters_n
    clusters{i} = X(clusterX==clusters_idx(i),:);
end
heatmap_clusters = [];
for i = 1:clusters_n
    heatmap_clusters(end+1:end+size(clusters{i},1),:) = clusters{i};
    heatmap_clusters(end+1:end+50,:) = NaN;
end

cmap = flipud(cbrewer('div','RdBu',1000,'spline')); cmap(cmap<0) = 0;
cmap_clusters = [225 196 225; ...
    216.75 216.75 216.75; ...
    208 231 213; ...
    54 161 86]/255;
subplot(1,2,1);
hold on
axis square
heatmap(heatmap_clusters(1:end-50,:),[],[],[],'Colormap',cmap,'ColorLevels',1000,'MaxColorValue',.5,'MinColorValue',-.5,'NaNColor',[1 1 1]);
s0 = 0.5;
for i = 1:clusters_n
    s1 = s0+ sum(clusterX==clusters_idx(i));
    plot([4.5 4.5],[s0 s1],'Color',cmap_clusters(i,:),'LineWidth',12)
    s0 = s1+50;
end
yticks([])
xticks([0:100:300]+.5)
xticklabels({'-1','0','1','2'})
ylabel('CGRP stim response types')
xlabel('Time (s)')
colormap(gca,cmap)
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

subplot(1,2,2)
hold on
axis square
for i = 1:clusters_n
    plot(edges(1:end-1)+.005,mean(clusters{i}),'Color',cmap_clusters(i,:),'LineWidth',1)
end
for i = 0:.1:.9
    plot([i i],[1.90 1.95],'Color',[54 161 86]/255,'LineWidth',2)
end
xlabel('Time (s)'); xlim([-1 2]); xticks(-1:1:2);
ylabel('Spiking (σ)'); ylim([-0.5 2]); yticks(-.5:.5:2);
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off
%% Fig ED5d

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 5d','FontWeight','bold')

summary = []; n = [];
regions_short = {'CEAc','CEAl','CEAm','BMAa','COAa','IA','BLAp','PAA','MEA','BLAa','Other','BLAv','BMAp','COAp'};
for i = 1:length(regions_short)
    n(i) = sum(cellfun(@(x) isequal(x,regions_short{i}),locdata.region));
    for j = 1:clusters_n
        summary(i,j) = sum(clusterX(cellfun(@(x) isequal(x,regions_short{i}),locdata.region))==clusters_idx(j))/n(i);
    end
end
hold on
axis square
b = bar(1:length(summary),summary(:,3:4)*100,'stacked','LineStyle','none');
b(1).FaceColor = cmap_clusters(3,:);
b(2).FaceColor = cmap_clusters(4,:);
xlim([0 length(summary)+1]); xticks(1:length(summary)); xticklabels(regions_short)
ylabel('CGRP-activated neurons (%)'); ylim([0 50]); yticks(0:25:50);
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off
%% Fig ED5e

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 5e','FontWeight','bold')

load([data_path,'\Neuropixels\CGRP stim acute\allen_ccfv3_modified_cz_v0.mat'],'atlas','RegionLibrary')
amygdala_regions = find(contains(RegionLibrary.reduced.region(1:200),'amyg'));
cmap_HSV = hsv(length(amygdala_regions));
cmap_HSV = cmap_HSV*.75;
cmap_HSV = circshift(cmap_HSV,-4,1);
brainoutline = atlas>0;
[ML,AP,DV] = bregma2ccf(locdata.x, locdata.y, locdata.z);
plane = 260;
B = bwboundaries(squeeze(brainoutline(:,plane,:)));
atlasmaskplot2 = atlas==1306;
B2a = bwboundaries(squeeze(atlasmaskplot2(:,plane,:)));
atlasmaskplot2 = atlas==1115;
B2b = bwboundaries(squeeze(atlasmaskplot2(:,plane,:)));
atlasmaskplot2 = atlas==1;
B2d = bwboundaries(squeeze(atlasmaskplot2(:,plane,:)));
subplot(1,2,1);
hold on
axis square
for i = 1:length(amygdala_regions)
    idx = RegionLibrary.reduced{amygdala_regions(i),1}+1;
    mask = squeeze(atlas(:,plane,:)) == idx;
    B3 = bwboundaries(mask);
    for j = 1:size(B3)
        fill(B3{j}(:,1),B3{j}(:,2),cmap_HSV(i,:),'linestyle','none','facealpha',0.5)
    end
end
for i = 1:length(B2a)
    fill(B2a{i}(:,1),B2a{i}(:,2),[1 1 1],'LineWidth',1)
end
for i = 1:length(B2b)
    fill(B2b{i}(:,1),B2b{i}(:,2),[.85 .85 .85],'LineWidth',1)
end
for i = 1:length(B2d)
    fill(B2d{i}(:,1),B2d{i}(:,2),[.85 .85 .85],'LineWidth',1)
end
for i = 1:length(B)
    plot(B{i}(:,1),B{i}(:,2),'k','LineWidth',1)
end
for j = 1:clusters_n
    x = ML(clusterX==clusters_idx(j))+1*rand(size(ML(clusterX==clusters_idx(j))))-.5;
    y = DV(clusterX==clusters_idx(j));
    if j == 1
        scatter(x,y,50,'filled','MarkerFaceColor',[0.5 0 0.5],'MarkerFaceAlpha',0.25)
    elseif j == 2
        scatter(x,y,50,'filled','MarkerFaceColor',[.85 .85 .85])
    elseif j == 3
        scatter(x,y,50,'filled','MarkerFaceColor',[54 161 86]/255,'MarkerFaceAlpha',0.25)
    else
        scatter(x,y,50,'filled','MarkerFaceColor',[54 161 86]/255)
    end
end
plot([277 397 397 277 277],[187 187 307 307 187],'k--','LineWidth',1)
xlim([277 397])
ylim([187 307])
xticks([]); yticks([]);
set(gca,'Ydir','reverse','Xdir','reverse')
title('Coronal projection')
axis off
set(gca,'FontSize',12)
hold off

subplot(1,2,2);
plane = 340;
B = bwboundaries(squeeze(brainoutline(plane,:,:)));
atlasmaskplot2 = atlas==1306;
B2a = bwboundaries(squeeze(atlasmaskplot2(plane,:,:)));
atlasmaskplot2 = atlas==1115;
B2b = bwboundaries(squeeze(atlasmaskplot2(plane,:,:)));
atlasmaskplot2 = atlas==1;
B2d = bwboundaries(squeeze(atlasmaskplot2(plane,:,:)));
hold on
axis square
for i = 1:length(amygdala_regions)
    idx = RegionLibrary.reduced{amygdala_regions(i),1}+1;
    mask = squeeze(atlas(plane,:,:)) == idx;
    B3 = bwboundaries(mask);
    for j = 1:size(B3)
        fill(B3{j}(:,1),B3{j}(:,2),cmap_HSV(i,:),'linestyle','none','facealpha',0.5)
    end
end
for i = 1:length(B2a)
    fill(B2a{i}(:,1),B2a{i}(:,2),[1 1 1],'LineWidth',1)
end
for i = 1:length(B2b)
    fill(B2b{i}(:,1),B2b{i}(:,2),[.85 .85 .85],'LineWidth',1)
end
for i = 1:length(B2d)
    fill(B2d{i}(:,1),B2d{i}(:,2),[.85 .85 .85],'LineWidth',1)
end
for i = 1:length(B)
    plot(B{i}(:,1),B{i}(:,2),'k','LineWidth',1)
end
for j = 1:clusters_n
    x = AP(clusterX==clusters_idx(j))+1*rand(size(AP(clusterX==clusters_idx(j))))-.5;
    y = DV(clusterX==clusters_idx(j));
    if j == 1
        scatter(x,y,50,'filled','MarkerFaceColor',[0.5 0 0.5],'MarkerFaceAlpha',0.25)
    elseif j == 2
        scatter(x,y,50,'filled','MarkerFaceColor',[.85 .85 .85])
    elseif j == 3
        scatter(x,y,50,'filled','MarkerFaceColor',[54 161 86]/255,'MarkerFaceAlpha',0.25)
    else
        scatter(x,y,50,'filled','MarkerFaceColor',[54 161 86]/255)
    end
end
plot([212 332 332 212 212],[187 187 307 307 187],'k--','LineWidth',1)
xlim([212 332]); yticks([]);
ylim([187 307]); xticks([]);
set(gca,'Ydir','reverse')
set(gca,'FontSize',12)
axis off
title('Sagittal projection')

drawnow
function [ML,AP,DV] = bregma2ccf(x,y,z)
bregma = [5739 5400 332];
ML = nan(size(x));
AP = nan(size(y));
DV = nan(size(z));
for i = 1:length(x)
    ML(i) = (bregma(1) + x(i))/25 + 1;
    AP(i) = (bregma(2) - y(i))/25 + 1;
    DV(i) = (bregma(3) - z(i))/25 + 1;
end