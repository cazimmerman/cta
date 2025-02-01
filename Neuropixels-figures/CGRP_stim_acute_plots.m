%% setp
clear all; close all;
addpath(genpath('Z:\Chris\matlab\heatmap\'));
addpath(genpath('Z:\Chris\matlab\cz\neuropixels-utils\'));
min_amp = 20;
unittype = {'good','mua'};
regions = 'amygdala'; % Options: 'CEA','notCEA','amygdala','all'
%% load data
fpath = 'Z:\Chris\matlab\cz\neuropixels-cgrp-acute\data\';
fnames = dir([fpath,'*.mat']);

cmap_RdBu = flipud(cbrewer('div','RdBu',1000,'spline')); cmap_RdBu(cmap_RdBu<0) = 0;
cmap_Grey = cbrewer('seq','Greys',1000,'spline');
regions_short = {'CEAc','CEAl','CEAm','BMAa','COAa','IA','BLAp','PAA','MEA','BLAa','Other','BLAv','BMAp','COAp'};

library = struct;
library.CEA = {'CEAc','CEAl','CEAm','CEAv','CEAast'};
library.amygdala = {'CEAc','CEAl','CEAm','CEAv','CEAast','BMAa','BMAp','BLAa','BLAp','BLAv','COAa','COAp','IA','MEA','LA','PAA','AAA','PA'};
library.all = cell(0,0);

disp('Loading data...')
data_out = struct;
traj = struct;
for i = 1:length(fnames)
    clear data
    load([fpath,fnames(i).name]);
    data_out(i).laser = data.final.laser(:,1);
    data_out(i).spikes = [data.final.ephys.spikes];
    data_out(i).region = [data.final.ephys.region];
    data_out(i).label = [data.final.ephys.label];
    data_out(i).x = [data.final.ephys.x];
    data_out(i).y = [data.final.ephys.y];
    data_out(i).z = [data.final.ephys.z];
    data_out(i).firing_rate = [data.final.ephys.stats.firing_rate];
    data_out(i).amplitude = [data.final.ephys.stats.amplitude];
    data_out(i).trajectory_full = data.final.trajectory.full;
    data_out(i).trajectory_recording = data.final.trajectory.recording;
    
    traj(i).full = data.final.trajectory.full;
    traj(i).recording = data.final.trajectory.recording;
    library.all = unique([library.all data_out(i).region]);
end
data = data_out;
save('cgrp-stim-acute.mat','data')
%%
counter = 1;
binedges = -(60*5):0.01:(60*10);
edges = -1:.01:2;
times = [-(60*5):0.01:(60*10)]+.005; times = times(1:end-1);
PSTH = [];
locdata = struct;
for i = 1:length(X)
    for j = 1:length(data_out(i).spikes)
        if data(i).amplitude(j)>=min_amp && ismember(data(i).label{j},unittype) ...
           && ismember(data(i).region{j},library.(regions))% data_in{i}.meta.isi_viol(j)<=0.1
       
            t = histcounts(data(i).spikes{j}-data(i).laser(1),binedges);
            idx1 = find(times<=-1,1,'last');
            idx2 = find(times<=data(i).laser(end)-data(i).laser(1)+2,1,'last');
            mu = mean(t(idx1:idx2));
            sigma = std(t(idx1:idx2));
            for k = 1:60
                spks = data(i).spikes{j} - data(i).laser((k-1)*10+1);
                PSTH(counter,k,:) = (histcounts(spks,edges)-mu)./sigma;
            end
            if ismember(data(i).region{j},{'PA','AAA','LA'})
                data(i).region{j} = 'Other';
            end
            locdata.region{counter} = data(i).region{j};
            locdata.x(counter) = data(i).x(j);
            locdata.y(counter) = data(i).y(j);
            locdata.z(counter) = data(i).z(j);
            counter = counter + 1;
            
        end
    end
end
%% summary stats

x = unique(locdata.region);
y = countcats(categorical(locdata.region));
[y,idx] = sort(y,'descend');
x = x(idx);
disp('Yield summary:')
t = table(x',y','VariableNames',{'Region','Yield'});
disp(t)
%% GMM

X = squeeze(mean(PSTH,2)); X = X - mean(X(:,1:100),2);
T = mean(PSTH(:,:,1:100),3); T = T - mean(T(:,1:20),2);
clusters_n = 4;
rng(12345)
disp(['Fitting GMM with ',num2str(clusters_n),' clusters...'])
GMModel=fitgmdist(X(:,101:200),clusters_n,'CovarianceType','diagonal','RegularizationValue',1e-5,'SharedCovariance',false,'Replicates',100);
save('GMModel.mat','GMModel')
clusterX = cluster(GMModel,X(:,101:200));
a=histcounts(clusterX,0.5:1:clusters_n+.5)./length(clusterX)*100;

clusters = cell(0,0);
for i = 1:clusters_n
    clusters{i} = X(clusterX==i,:);
end
[~,idx] = sort(cellfun(@(x) mean(mean(x(:,101:200),2)),clusters),'ascend');
clusters_idx = idx;%idx([2,1,3,4]);
clusters = cell(0,0);
clusters_tonic = cell(0,0);
for i = 1:clusters_n
    ii = clusters_idx(i);
    clusters{i} = X(clusterX==ii,:);
    clusters_tonic{i} = T(clusterX==ii,:);
end

clusters_heatmap_evoked = [];
clusters_heatmap_tonic = [];
clusters_ylabels = [];
for i = 1:clusters_n
    clusters_heatmap_evoked(end+1:end+size(clusters{i},1),:) = clusters{i};
    clusters_heatmap_evoked(end+1:end+50,:) = NaN;
    clusters_heatmap_tonic(end+1:end+size(clusters{i},1),:) = clusters_tonic{i};
    clusters_heatmap_tonic(end+1:end+50,:) = NaN;
    if i == 1
        clusters_ylabels(i) = size(clusters{i},1)/2;
    else
        clusters_ylabels(i) = sum(cellfun(@(x) size(x,1),clusters(1:i-1))) + 20*(i-1) + size(clusters{i},1)/2;
    end
end

clusters_cmap = [.75 .75 1; ...
    .75 .75 .75;...
    1 .75 .75; ...
    1 .25 .25];
%% try hierarchical clustering
% treedata = X(:,101:200);
% tree = linkage(treedata,'ward','chebychev');
% D = pdist(treedata);
% leafOrder = optimalleaforder(tree,D);
% figure('Position', get(0,'Screensize'))
% subplot(1,2,1)
% hold on
% hLines = dendrogram(tree,length(treedata)+1,'ColorThreshold',4,'Reorder',leafOrder);
% xtickangle(90)
% hold off
% subplot(1,2,2)
% hold on
% heatmap(X(leafOrder,:),[],[],[],'Colormap',cmap_RdBu,'ColorLevels',1000,'MaxColorValue',.5,'MinColorValue',-.5);
% hold off
%% evoked summary plot

close all
figure('Position', get(0, 'Screensize'))
load('Z:\Chris\matlab\cz\neuropixels-utils\Allen_CCFv3_RegionLibrary.mat','RegionLibrary')
AMYG = find(contains(RegionLibrary.reduced.region,'amyg'));
AMYG = find(contains(RegionLibrary.reduced.region(1:200),'amyg'));
cmap_HSV = hsv(length(AMYG)); 
cmap_HSV = cmap_HSV*.75;
cmap_HSV = circshift(cmap_HSV,-4,1);

ax1 = subplot(2,3,1);
hold on
axis square
heatmap(clusters_heatmap_evoked(1:end-50,:),[],[],[],'Colormap',cmap_RdBu,'ColorLevels',1000,'MaxColorValue',.5,'MinColorValue',-.5,'NaNColor',[1 1 1]);
plot([100 100]+.5,[0 length(clusters_heatmap_evoked(1:end-50,:))]+0.5,'k','LineWidth',1)
plot([200 200]+.5,[0 length(clusters_heatmap_evoked(1:end-50,:))]+0.5,'k','LineWidth',1)
s0 = 0.5;
for i = 1:clusters_n
    s1 = s0+ sum(clusterX==clusters_idx(i));
    plot([4.5 4.5],[s0 s1],'Color',clusters_cmap(i,:),'LineWidth',10)
    s0 = s1+50;
end
yticks(clusters_ylabels)
yticklabels(1:clusters_n)
xticks([0:100:300]+.5)
xticklabels({'-1','0','1','2'})
ylabel('Cluster')
xlabel('Time (s)')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
h=gca; h.YAxis.TickLength = [0 0];
title('Evoked activity')
colormap(gca,cmap_RdBu)
c1 = colorbar('Location','westoutside','FontSize',16);
c1.Position = c1.Position+[-.03 0 0 0];
c1.Ticks = ([-.5 0 .5]);
c1.Label.String = 'Spiking (σ)';
hold off

subplot(2,3,2)
hold on
axis square
for i = 1:clusters_n
    plot(edges(1:end-1)+.005,mean(clusters{i}),'Color',clusters_cmap(i,:),'LineWidth',1)
end
plot([0 0],[1.90 1.95],'Color',[0 .5 1],'LineWidth',2)
plot([0.1 0.1],[1.90 1.95],'Color',[0 .5 1],'LineWidth',2)
plot([0.2 0.2],[1.90 1.95],'Color',[0 .5 1],'LineWidth',2)
plot([0.3 0.3],[1.90 1.95],'Color',[0 .5 1],'LineWidth',2)
plot([0.4 0.4],[1.90 1.95],'Color',[0 .5 1],'LineWidth',2)
plot([0.5 0.5],[1.90 1.95],'Color',[0 .5 1],'LineWidth',2)
plot([0.6 0.6],[1.90 1.95],'Color',[0 .5 1],'LineWidth',2)
plot([0.7 0.7],[1.90 1.95],'Color',[0 .5 1],'LineWidth',2)
plot([0.8 0.8],[1.90 1.95],'Color',[0 .5 1],'LineWidth',2)
plot([0.9 0.9],[1.90 1.95],'Color',[0 .5 1],'LineWidth',2)
ylim([-0.5 2])
xlim([-1 2])
yticks(-.5:.5:2)
xticks(-1:1:2)
xlabel('Time (s)')
ylabel('Spiking (σ)')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

subplot(2,3,3)
summary = [];
n = [];
for i = 1:length(regions_short)
    n(i) = sum(cellfun(@(x) isequal(x,regions_short{i}),locdata.region));
    for j = 1:clusters_n
        summary(i,j) = sum(clusterX(cellfun(@(x) isequal(x,regions_short{i}),locdata.region))==clusters_idx(j))/n(i);
    end
end
hold on
axis square
b = bar(1:length(summary),summary(:,3:4)*100,'stacked','LineWidth',1);
%text([1:5,6.5:7.5,9:11,12.5:13.5,15,16.5],sum(summary(:,3:5)*100,2),num2str(n'),'vert','bottom','horiz','center');
b(1).FaceColor = clusters_cmap(3,:);
b(2).FaceColor = clusters_cmap(4,:);
b2 = bar(1:length(summary),-summary(:,1)*100,'LineWidth',1);
b2(1).FaceColor = clusters_cmap(1,:);
xticks(1:length(summary))
xticklabels(regions_short)
xtickangle(90)
yticks(-50:25:50)
yticklabels({'50','25','0','25','50'})
ylim([-50 50])
xlim([0 length(summary)+1])
ylabel('Units (%)')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

summary = [];
n = [];
for i = 1:length(regions_short)
    n(i) = sum(cellfun(@(x) isequal(x,regions_short{i}),locdata.region));
    for j = 1:clusters_n
        summary(i,j) = sum(clusterX(cellfun(@(x) isequal(x,regions_short{i}),locdata.region))==clusters_idx(j))/n(i);
    end
end

t = table(regions_short',summary(:,4),summary(:,3),summary(:,2),summary(:,1),'VariableNames',{'Region','Strong act','Weak act','Unmod','Inh'});
writetable(t,'Z:\Chris\matlab\cz\cta-source-data\Fig-ED5d.csv')

load('Z:\Chris\matlab\cz\neuropixels-utils\Allen_CCFv3_Fos_AST.mat','atlas')
brainoutline = atlas>0;
load('Z:\Chris\matlab\cz\neuropixels-utils\Fos_Density_CGRP.mat','CGRP')
[ML,AP,DV] = bregma2ccf(locdata.x, locdata.y, locdata.z);
plane = 260;
B = bwboundaries(squeeze(brainoutline(:,plane,:)));
atlasmaskplot2 = atlas==1306;
B2a = bwboundaries(squeeze(atlasmaskplot2(:,plane,:)));
atlasmaskplot2 = atlas==1115;
B2b = bwboundaries(squeeze(atlasmaskplot2(:,plane,:)));
atlasmaskplot2 = atlas==1;
B2d = bwboundaries(squeeze(atlasmaskplot2(:,plane,:)));
subplot(2,2,3);
hold on
axis square
data = [CGRP.Novel(:,plane,:),CGRP.Familiar(:,plane,:)];
data = squeeze(mean(data,2));
colormap(gcf,cmap_Grey)
heatmap(rot90(data,-1),[],[],[],'UseFigureColormap',true,'ColorLevels',1000,'MaxColorValue',1,'MinColorValue',0,'NaNColor',[1 1 1]);
for i = 1:length(AMYG)
    idx = RegionLibrary.reduced{AMYG(i),1}+1;
    mask = squeeze(atlas(:,plane,:)) == idx;
    B3 = bwboundaries(mask);
    for j = 1:size(B3)
        plot(B3{j}(:,1),B3{j}(:,2),'Color',cmap_HSV(i,:),'LineWidth',2)
    end
end
for i = 1:length(B2a)
    fill(B2a{i}(:,1),B2a{i}(:,2),[1 1 1],'LineWidth',2)
end
for i = 1:length(B2b)
    fill(B2b{i}(:,1),B2b{i}(:,2),[.85 .85 .85],'LineWidth',2)
end
for i = 1:length(B2d)
    fill(B2d{i}(:,1),B2d{i}(:,2),[.85 .85 .85],'LineWidth',2)
end
for i = 1:length(B)
    plot(B{i}(:,1),B{i}(:,2),'k','LineWidth',2)
end
for j = 1:clusters_n
    x = ML(clusterX==clusters_idx(j))+1*rand(size(ML(clusterX==clusters_idx(j))))-.5;
    y = DV(clusterX==clusters_idx(j));
    scatter(x,y,12,'filled','MarkerFaceColor',clusters_cmap(j,:))
end
plot([340 340],[187 307],'k','LineWidth',1)
plot([277 397 397 277 277],[187 187 307 307 187],'k','LineWidth',1)
xlim([277 397])
ylim([187 307])
set(gca,'Ydir','reverse')
title('Coronal projection')
set(gca,'FontSize',16)
c1 = colorbar('Location','westoutside','FontSize',16);
c1.Position = c1.Position+[-.015 0 0 0];
c1.Ticks = [0 .5 1];
c1.Label.String = 'Fos (% cells per mm^3)';
hold off

subplot(2,2,4);
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
data = [CGRP.Novel(plane,:,:);CGRP.Familiar(plane,:,:)];
data = squeeze(mean(data,1));
heatmap(fliplr(rot90(data,-1)),[],[],[],'UseFigureColormap',true,'ColorLevels',1000,'MaxColorValue',1,'MinColorValue',0,'NaNColor',[1 1 1]);
for i = 1:length(AMYG)
    idx = RegionLibrary.reduced{AMYG(i),1}+1;
    mask = squeeze(atlas(plane,:,:)) == idx;
    B3 = bwboundaries(mask);
    for j = 1:size(B3)
        plot(B3{j}(:,1),B3{j}(:,2),'Color',cmap_HSV(i,:),'LineWidth',2)
    end
end
for i = 1:length(B2a)
    fill(B2a{i}(:,1),B2a{i}(:,2),[1 1 1],'LineWidth',2)
end
for i = 1:length(B2b)
    fill(B2b{i}(:,1),B2b{i}(:,2),[.85 .85 .85],'LineWidth',2)
end
for i = 1:length(B2d)
    fill(B2d{i}(:,1),B2d{i}(:,2),[.85 .85 .85],'LineWidth',2)
end
for i = 1:length(B)
    plot(B{i}(:,1),B{i}(:,2),'k','LineWidth',2)
end
for j = 1:clusters_n
    x = AP(clusterX==clusters_idx(j))+1*rand(size(AP(clusterX==clusters_idx(j))))-.5;
    y = DV(clusterX==clusters_idx(j));
    scatter(x,y,12,'filled','MarkerFaceColor',clusters_cmap(j,:))
end
plot([260 260],[187 307],'k','LineWidth',1)
plot([212 332 332 212 212],[187 187 307 307 187],'k','LineWidth',1)
xlim([212 332])
ylim([187 307])
set(gca,'Ydir','reverse')
set(gca,'FontSize',16)
title('Sagittal projection')
c1 = colorbar('Location','westoutside','FontSize',16);
c1.Position = c1.Position+[-.015 0 0 0];
c1.Ticks = [0 .5 1];
c1.Label.String = 'Fos (% cells per mm^3)';
hold off
colormap(ax1,cmap_RdBu)

saveas(gcf,['plots-png/summary-plot-acute-CGRP-1-',num2str(min_amp),'uV-',unittype{end}],'png')
set(gcf,'renderer','painters')
saveas(gcf,['plots-eps/summary-plot-acute-CGRP-1-',num2str(min_amp),'uV-',unittype{end}],'epsc')
%% baseline summary plot

figure('Position', get(0, 'Screensize'))

p = [];
for i = 1:size(T,1)
    p(i) = ranksum(T(i,1:20),T(i,end-19:end));
    %p(i) = fitlm(1:60,T(i,:)).Coefficients.pValue(2);
    %[a,b] = corrcoef(1:60,T(i,:)); p(i) = b(2,1);
end

ax1 = subplot(2,3,1);
hold on
axis square
bd{1} = T(p<=0.05 & mean(T(:,end-19:end),2)'>0,:);
[~,ii] = sort(mean(bd{1}(:,41:60),2)); bd{1} = bd{1}(ii,:);
bd{2} = T(p>0.05,:);
[~,ii] = sort(mean(bd{2}(:,41:60),2)); bd{2} = bd{2}(ii,:);
bd{3} = T(p<=0.05 & mean(T(:,end-19:end),2)'<0,:);
[~,ii] = sort(mean(bd{3}(:,41:60),2)); bd{3} = bd{3}(ii,:);
heatmap([bd{3}; nan(20,60); bd{2}; nan(20,60); bd{1}],[],[],[],'Colormap',cmap_RdBu,'ColorLevels',1000,'MaxColorValue',.5,'MinColorValue',-.5,'NaNColor',[1 1 1]);
xticks([0:10:60]+0.5)
xticklabels({'0','10','20','30','40','50','60'})
xlabel('Trial')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
h=gca; h.YAxis.TickLength = [0 0];
title('Baseline activity')
colormap(gca,cmap_RdBu)
c1 = colorbar('Location','westoutside','FontSize',16);
c1.Position = c1.Position+[-.015 0 0 0];
c1.Ticks = ([-.5 0 .5]);
c1.Label.String = 'Spiking (σ)';
hold off

subplot(2,3,2)
hold on
axis square
summary = [];
n = [];
for i = 1:clusters_n
    n(i) = sum(clusterX==clusters_idx(i));
    summary(i) = sum(p<=0.05 & [clusterX==clusters_idx(i)]' & mean(T(:,end-19:end),2)'>0)/n(i);
end
b = bar([1:clusters_n],summary*100,'FaceColor',[1 .75 .75],'LineWidth',1);
%text([1:5],summary*100+1,num2str(n'),'vert','bottom','horiz','center');
for i = 1:clusters_n
    n(i) = sum(clusterX==clusters_idx(i));
    summary(i) = sum(p<=0.05 & [clusterX==clusters_idx(i)]' & mean(T(:,end-19:end),2)'<0)/n(i);
end
b = bar([1:clusters_n],-summary*100,'FaceColor',[.75 .75 1],'LineWidth',1);
xticks([1:clusters_n])
xlabel('Cluster')
yticks(-20:10:20)
yticklabels({'20','10','0','10','20'})
ylim([-20 20])
xlim([0 clusters_n+1])
ylabel('Units (%)')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')

subplot(2,3,3)
hold on
axis square
for i = 1:length(regions_short)
    n(i) = sum(cellfun(@(x) isequal(x,regions_short{i}),locdata.region));
    summary(i) = sum(p<=0.05 & cellfun(@(x) isequal(x,regions_short{i}),locdata.region) & mean(T(:,end-19:end),2)'>0)/n(i);
end
bar(1:length(summary),summary*100,'FaceColor',[1 .75 .75],'LineWidth',1);
%text([1:5,6.5:7.5,9:11,12.5:13.5,15,16.5],summary*100+1,num2str(n'),'vert','bottom','horiz','center');
for i = 1:length(regions_short)
    n(i) = sum(cellfun(@(x) isequal(x,regions_short{i}),locdata.region));
    summary(i) = sum(p<=0.05 & cellfun(@(x) isequal(x,regions_short{i}),locdata.region) & mean(T(:,end-19:end),2)'<0)/n(i);
end
bar(1:length(summary),-summary*100,'FaceColor',[.75 .75 1],'LineWidth',1);
xticks(1:length(summary))
xticklabels(regions_short)
xtickangle(90)
yticks(-30:10:30)
yticklabels({'30','20','10','0','10','20','30'})
ylim([-30 30])
xlim([0 length(summary)+1])
ylabel('Units (%)')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

subplot(2,2,3);
load('Z:\Chris\matlab\cz\neuropixels-utils\Allen_CCFv3_Fos_AST.mat','atlas')
brainoutline = atlas>0;
load('Z:\Chris\matlab\cz\neuropixels-utils\Fos_Density_CGRP.mat','CGRP')
[ML,AP,DV] = bregma2ccf(locdata.x, locdata.y, locdata.z);
plane = 260;
B = bwboundaries(squeeze(brainoutline(:,plane,:)));
atlasmaskplot2 = atlas==1306;
B2a = bwboundaries(squeeze(atlasmaskplot2(:,plane,:)));
atlasmaskplot2 = atlas==1115;
B2b = bwboundaries(squeeze(atlasmaskplot2(:,plane,:)));
atlasmaskplot2 = atlas==1;
B2d = bwboundaries(squeeze(atlasmaskplot2(:,plane,:)));
hold on
axis square
data = [CGRP.Novel(:,plane,:),CGRP.Familiar(:,plane,:)];
data = squeeze(mean(data,2));
colormap(gcf,cmap_Grey)
heatmap(rot90(data,-1),[],[],[],'UseFigureColormap',true,'ColorLevels',1000,'MaxColorValue',1,'MinColorValue',0,'NaNColor',[1 1 1]);
for i = 1:length(AMYG)
    idx = RegionLibrary.reduced{AMYG(i),1}+1;
    mask = squeeze(atlas(:,plane,:)) == idx;
    B3 = bwboundaries(mask);
    for j = 1:size(B3)
        plot(B3{j}(:,1),B3{j}(:,2),'Color',cmap_HSV(i,:),'LineWidth',2)
    end
end
for i = 1:length(B2a)
    fill(B2a{i}(:,1),B2a{i}(:,2),[1 1 1],'LineWidth',2)
end
for i = 1:length(B2b)
    fill(B2b{i}(:,1),B2b{i}(:,2),[.85 .85 .85],'LineWidth',2)
end
for i = 1:length(B2d)
    fill(B2d{i}(:,1),B2d{i}(:,2),[.85 .85 .85],'LineWidth',2)
end
for i = 1:length(B)
    plot(B{i}(:,1),B{i}(:,2),'k','LineWidth',2)
end
[~,idx] = sort(p,'descend');
T = T(idx,:);
AP = AP(idx);
ML = ML(idx);
DV = DV(idx);
for j = 1:size(T,1)
    x = ML(j)+1*rand(size(ML(j)))-.5;
    y = DV(j);
    idx = round((mean(T(j,end-19:end))+.25)/.5*1000);
    idx(idx<1) = 1; idx(idx>1000) = 1000;
    scatter(x,y,12,'filled','MarkerFaceColor',cmap_RdBu(idx,:))
end
plot([340 340],[187 307],'k','LineWidth',1)
plot([277 397 397 277 277],[187 187 307 307 187],'k','LineWidth',1)
xlim([277 397])
ylim([187 307])
set(gca,'Ydir','reverse')
title('Coronal projection')
set(gca,'FontSize',16)
c1 = colorbar('Location','westoutside','FontSize',16);
c1.Position = c1.Position+[-.015 0 0 0];
c1.Ticks = [0 .5 1];
c1.Label.String = 'Fos (% cells per mm^3)';
hold off

subplot(2,2,4);
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
data = [CGRP.Novel(plane,:,:);CGRP.Familiar(plane,:,:)];
data = squeeze(mean(data,1));
heatmap(fliplr(rot90(data,-1)),[],[],[],'UseFigureColormap',true,'ColorLevels',1000,'MaxColorValue',1,'MinColorValue',0,'NaNColor',[1 1 1]);
for i = 1:length(AMYG)
    idx = RegionLibrary.reduced{AMYG(i),1}+1;
    mask = squeeze(atlas(plane,:,:)) == idx;
    B3 = bwboundaries(mask);
    for j = 1:size(B3)
        plot(B3{j}(:,1),B3{j}(:,2),'Color',cmap_HSV(i,:),'LineWidth',2)
    end
end
for i = 1:length(B2a)
    fill(B2a{i}(:,1),B2a{i}(:,2),[1 1 1],'LineWidth',2)
end
for i = 1:length(B2b)
    fill(B2b{i}(:,1),B2b{i}(:,2),[.85 .85 .85],'LineWidth',2)
end
for i = 1:length(B2d)
    fill(B2d{i}(:,1),B2d{i}(:,2),[.85 .85 .85],'LineWidth',2)
end
for i = 1:length(B)
    plot(B{i}(:,1),B{i}(:,2),'k','LineWidth',2)
end
for j = 1:size(T,1)
    x = AP(j)+1*rand(size(AP(j)))-.5;
    y = DV(j);
    idx = round((mean(T(j,end-19:end))+.25)/.5*1000);
    idx(idx<1) = 1; idx(idx>1000) = 1000;
    scatter(x,y,12,'filled','MarkerFaceColor',cmap_RdBu(idx,:))
end
plot([260 260],[187 307],'k','LineWidth',1)
plot([212 332 332 212 212],[187 187 307 307 187],'k','LineWidth',1)
xlim([212 332])
ylim([187 307])
set(gca,'Ydir','reverse')
set(gca,'FontSize',16)
c1 = colorbar('Location','westoutside','FontSize',16);
c1.Position = c1.Position+[-.015 0 0 0];
c1.Ticks = [0 .5 1];
c1.Label.String = 'Fos (% cells per mm^3)';
title('Sagittal projection')
colormap(ax1,cmap_RdBu)

saveas(gcf,['plots-png/summary-plot-acute-CGRP-2-',num2str(min_amp),'uV-',unittype{end}],'png')
set(gcf,'renderer','painters')
saveas(gcf,['plots-eps/summary-plot-acute-CGRP-2-',num2str(min_amp),'uV-',unittype{end}],'epsc')
%% firing rate summary

X = struct;
traj = struct;
for i = 1:length(fnames)
    clear data
    load([fpath,fnames(i).name]);
    X(i).region = [data.final.ephys.region];
    X(i).firing_rate = [data.final.ephys.stats.firing_rate];
    X(i).amplitude = [data.final.ephys.stats.amplitude];
end

FR = [];
AMP = [];
REGION = cell(0,0);
for i = 1:length(fnames)
    FR = [FR X(i).firing_rate];
    AMP = [AMP abs(X(i).amplitude)];
    REGION = [REGION X(i).region];
end
xvals = [1:5,6.5:7.5,9:11,12.5:13.5,15,16.5,18,19.5];
figure('Position', get(0, 'Screensize'))
subplot(1,2,1)
hold on
axis square
for i = 1:length(regions_short)
    vals = FR(find(cellfun(@(x) isequal(x,regions_short{i}),REGION)));
    j = prctile(vals,[10 25 50 75 90]);
    plot([xvals(i) xvals(i)],[j(1) j(5)],'LineWidth',2,'Color',cmap_HSV(i,:))
    plot([xvals(i) xvals(i)],[j(2) j(4)],'LineWidth',8,'Color',cmap_HSV(i,:))
    plot([xvals(i)-.3 xvals(i)+.3],[j(3) j(3)],'LineWidth',2,'Color',cmap_HSV(i,:))
end
xticks([1:5,6.5:7.5,9:11,12.5:13.5,15,16.5,18,19.5])
xticklabels(regions_short)
xtickangle(90)
yticks(0:5:20)
ylim([0 20])
xlim([0 20.5])
ylabel('Firing rate (Hz)')
set(gca,'FontSize',24,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off
subplot(1,2,2)
hold on
axis square
for i = 1:length(regions_short)
    vals = AMP(find(cellfun(@(x) isequal(x,regions_short{i}),REGION)));
    j = prctile(vals,[10 25 50 75 90]);
    plot([xvals(i) xvals(i)],[j(1) j(5)],'LineWidth',2,'Color',cmap_HSV(i,:))
    plot([xvals(i) xvals(i)],[j(2) j(4)],'LineWidth',8,'Color',cmap_HSV(i,:))
    plot([xvals(i)-.3 xvals(i)+.3],[j(3) j(3)],'LineWidth',2,'Color',cmap_HSV(i,:))
end
xticks([1:5,6.5:7.5,9:11,12.5:13.5,15,16.5,18,19.5])
xticklabels(regions_short)
xtickangle(90)
yticks(20:20:100)
ylim([20 100])
xlim([0 20.5])
ylabel('Amplitude (µV)')
set(gca,'FontSize',24,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

saveas(gcf,['plots-png/summary-plot-acute-CGRP-3-',num2str(min_amp),'uV-',unittype{end}],'png')
set(gcf,'renderer','painters')
saveas(gcf,['plots-eps/summary-plot-acute-CGRP-3-',num2str(min_amp),'uV-',unittype{end}],'epsc')