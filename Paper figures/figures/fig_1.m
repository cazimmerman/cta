function fig_1(data_path)

disp('Generating panels for Figure 1...')
%% Fig 1b

figure('Position', get(0, 'Screensize'))
sgtitle('Figure 1b','FontWeight','bold')
T = readtable([data_path,'/source data/Fig-1b.csv']);

subplot(1,2,1)
hold on
axis square
idx = find(cellfun(@(x,y) isequal(x,'Novel').*isequal(y,'Saline'),T.Flavour,T.Injection));
data = [T.Day1(idx), T.Day2(idx), T.Day3(idx)]*100;
plot(1:3,data,'color',[.85 .85 .85],'linewidth',1)
a = plot(1:3,mean(data),'color','k','linewidth',1);
errorbar(1:3,mean(data),std(data)./sqrt(size(data,1)),'color','k','linestyle','none','capsize',0,'linewidth',1)
idx = find(cellfun(@(x,y) isequal(x,'Novel').*isequal(y,'LiCl'),T.Flavour,T.Injection));
data = [T.Day1(idx), T.Day2(idx), T.Day3(idx)]*100;
plot(1:3,data,'color',[252 216 213]/255,'linewidth',1)
b = plot(1:3,mean(data),'color',[229 45 38]/255,'linewidth',1);
errorbar(1:3,mean(data),std(data)./sqrt(size(data,1)),'color',[229 45 38]/255,'linestyle','none','capsize',0,'linewidth',1)
xlabel('Retrieval day'); xlim([0.5 3.5]); xticks(1:3);
ylabel('Flavour preference (%)'); ylim([0 100]); yticks(0:50:100);
title('Novel flavour')
legend([a,b],{'Saline','LiCl'})
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

subplot(1,2,2)
hold on
axis square
idx = find(cellfun(@(x,y) isequal(x,'Familiar').*isequal(y,'Saline'),T.Flavour,T.Injection));
data = [T.Day1(idx), T.Day2(idx), T.Day3(idx)]*100;
plot(1:3,data,'color',[.85 .85 .85],'linewidth',1)
a = plot(1:3,mean(data),'color','k','linewidth',1);
errorbar(1:3,mean(data),std(data)./sqrt(size(data,1)),'color','k','linestyle','none','capsize',0,'linewidth',1)
idx = find(cellfun(@(x,y) isequal(x,'Familiar').*isequal(y,'LiCl'),T.Flavour,T.Injection));
data = [T.Day1(idx), T.Day2(idx), T.Day3(idx)]*100;
plot(1:3,data,'color',[252 216 213]/255,'linewidth',1)
b = plot(1:3,mean(data),'color',[229 45 38]/255,'linewidth',1);
errorbar(1:3,mean(data),std(data)./sqrt(size(data,1)),'color',[229 45 38]/255,'linestyle','none','capsize',0,'linewidth',1)
xlabel('Retrieval day'); xlim([0.5 3.5]); xticks(1:3);
ylabel('Flavour preference (%)'); ylim([0 100]); yticks(0:50:100);
title('Familiar flavour')
legend([a,b],{'Saline','LiCl'})
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

T = readtable([data_path,'/GLMMs/CFA-behavior-GLMM-output.csv']);
p = multicmp(T{4,2:4},'up',0.05);
StatsTbl = table({'1b, top'},{'Day 1: LiCl vs. Saline'},{'GLMM marginal effect'},{'3 days'},{[8 8]},T{3,2},p(1), ...
    'VariableNames',{'Figure panel','Group','Statistical test','Multiple comparisons','Sample size','Test statistic','P-value'});
StatsTbl(end+1,:) = table({'1b, top'},{'Day 2: LiCl vs. Saline'},{'GLMM marginal effect'},{'3 days'},{[8 8]},T{3,3},p(2));
StatsTbl(end+1,:) = table({'1b, top'},{'Day 3: LiCl vs. Saline'},{'GLMM marginal effect'},{'3 days'},{[8 8]},T{3,4},p(3));
p = multicmp(T{4,5:7},'up',0.05);
StatsTbl(end+1,:) = table({'1b, bottom'},{'Day 1: LiCl vs. Saline'},{'GLMM marginal effect'},{'3 days'},{[8 8]},T{3,5},p(1));
StatsTbl(end+1,:) = table({'1b, bottom'},{'Day 2: LiCl vs. Saline'},{'GLMM marginal effect'},{'3 days'},{[8 8]},T{3,6},p(2));
StatsTbl(end+1,:) = table({'1b, bottom'},{'Day 3: LiCl vs. Saline'},{'GLMM marginal effect'},{'3 days'},{[8 8]},T{3,7},p(3));
%% Fig 1e

figure('Position', get(0, 'Screensize'))
sgtitle('Figure 1e','FontWeight','bold')
load([data_path,'/FOS imaging/FOS-GLMM-statistics.mat'],'data')

hold on
axis square
x = struct;
idx = find(data.GLMMoutput.Eq2.modelstats.significant);
x.Consume = data.GLMMoutput.Eq2.flavor.Zstat(idx,1);
x.Malaise = data.GLMMoutput.Eq2.flavor.Zstat(idx,2);
x.Retrieval = data.GLMMoutput.Eq2.flavor.Zstat(idx,3);
violinplot(x);
xlim([0.5 3.5]); xticks(1:3);
ylabel('Novel – Familiar ΔFOS (Z)'); ylim([-5 5]); yticks(-5:2.5:5);
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

stat = []; p = [];
[~,p(1),stat(1)] = kstest2(x.Consume,x.Malaise);
[~,p(2),stat(2)] = kstest2(x.Consume,x.Retrieval);
[~,p(3),stat(3)] = kstest2(x.Malaise,x.Retrieval);
p = multicmp(p,'up',0.05);
StatsTbl(end+1,:) = table({'1e'},{'Consume vs. Malaise'},{'Kolmogorov-Smirnov'},{'3 pairs of time points'},{[length(x.Consume)]},stat(1),p(1));
StatsTbl(end+1,:) = table({'1e'},{'Consume vs. Retrieval'},{'Kolmogorov-Smirnov'},{'3 pairs of time points'},{[length(x.Consume)]},stat(2),p(2));
StatsTbl(end+1,:) = table({'1e'},{'Malaise vs. Retrieval'},{'Kolmogorov-Smirnov'},{'3 pairs of time points'},{[length(x.Consume)]},stat(3),p(3));
%% Fig 1f

figure('Position', get(0, 'Screensize'))
sgtitle('Figure 1f','FontWeight','bold')
load([data_path,'/FOS imaging/FOS-GLMM-statistics.mat'],'data')
T = readtable([data_path,'/source data/Fig-1e.csv']);

ax1 = subplot(1,2,1);
cmap = flipud(cbrewer('div','RdBu',1000,'spline')); cmap(cmap<0) = 0;
significant = find(data.GLMMoutput.Eq2.modelstats.significant);
treedata = data.GLMMoutput.Eq2.flavor.Zstat(significant,:);
treenames = T.Abbreviation;
tree = linkage(treedata,'ward','chebychev');
D = pdist(treedata);
leafOrder = optimalleaforder(tree,D);
[hLines,~,outperm] = dendrogram(tree,length(significant)+1,'Labels',[],'ColorThreshold',4.7,'Reorder',leafOrder,'Orientation','left');
for i = 1:length(hLines)
    zz(i,:) = hLines(i).Color;
    if isequal(round(hLines(i).Color,1),[1.0000 0.6000 0])
        hLines(i).Color = [229 45 38]/255;
    elseif ~isequal(hLines(i).Color,[0 0 0])
        hLines(i).Color = [.85 .85 .85];
    end
end
set(hLines,'LineWidth',1)
xticks([]);
ylim([0.5 length(treedata)+0.5]); yticks(1:length(leafOrder)); yticklabels(treenames(leafOrder));;
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0, 0],'TickDir','out')
ax1.YAxis.FontSize = 6;
hold off

ax2 = subplot(1,2,2);
hold on
heatmap(treedata(leafOrder,1:3),[],[],[],'Colormap',cmap,'ColorLevels',1000,'MaxColorValue',5,'MinColorValue',-5,'NaNColor',[1 1 1]);
N_clusters = length(unique(zz,'rows'))-1;
idx = cluster(tree,'maxclust',N_clusters);
idx = idx(outperm);
list = unique(idx,'stable');
counter = [];
for i = 1:length(list)-1
    counter(i) = sum(idx==list(i));
    plot([0 3]+0.5,[sum(counter) sum(counter)]+0.5,'k','LineWidth',1)
end
yticks(1:length(leafOrder)); yticklabels(treenames(leafOrder));
xticks(1:3); xticklabels({'Consume','Malaise','Retrieval'});
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0, 0],'TickDir','out')
ax2.YAxis.FontSize = 6;
hold off
%% Fig 1h

if exist([data_path,'/FOS imaging/kernel-density-estimates'],'dir') ~= 7
    disp('Unzipping file: .../FOS imaging/kernel-density-estimates.zip')
    unzip([data_path,'/FOS imaging/kernel-density-estimates.zip'],[data_path,'/FOS imaging/'])
end

fname = '/FOS imaging/kernel-density-estimates/kde-consumption-novel.npy';
KDE.Consume.Novel = readNPY([data_path,fname]);
fname = '/FOS imaging/kernel-density-estimates/kde-consumption-familiar.npy';
KDE.Consume.Familiar = readNPY([data_path,fname]);
fname = '/FOS imaging/kernel-density-estimates/kde-malaise-novel.npy';
KDE.Malaise.Novel = readNPY([data_path,fname]);
fname = '/FOS imaging/kernel-density-estimates/kde-malaise-familiar.npy';
KDE.Malaise.Familiar = readNPY([data_path,fname]);
fname = '/FOS imaging/kernel-density-estimates/kde-retrieval-novel.npy';
KDE.Retrieval.Novel = readNPY([data_path,fname]);
fname = '/FOS imaging/kernel-density-estimates/kde-retrieval-familiar.npy';
KDE.Retrieval.Familiar = readNPY([data_path,fname]);

load([data_path,'/FOS imaging/modified-atlas/allen_ccfv3_modified_cz.mat'],'atlas')

atlasmask  = atlas>=1028 | atlas<=1;
KDE.Consume.Novel(atlasmask) = NaN;
KDE.Consume.Familiar(atlasmask) = NaN;
KDE.Malaise.Novel(atlasmask) = NaN;
KDE.Malaise.Familiar(atlasmask) = NaN;
KDE.Retrieval.Novel(atlasmask) = NaN;
KDE.Retrieval.Familiar(atlasmask) = NaN;

cmap = flipud(cbrewer('div','RdBu',1000,'spline')); cmap(cmap<0) = 0;
lims2 = [-0.5 0.5];


figure('Position', get(0, 'Screensize'))
sgtitle('Figure 1h','FontWeight','bold')
load([data_path,'/FOS imaging/FOS-GLMM-statistics.mat'],'data'); data_in = data;
T = readtable([data_path,'/source data/Fig-1e.csv']);

idx = find(data.GLMMoutput.Eq2.modelstats.significant);
idx_amygdala = idx(find(T.Cluster==1));
planes_nov = [210,265,325];

for pl = 1:2
    plane = planes_nov(pl+1);
    subplot(2,3,1+(pl-1)*3)
    hold on
    axis equal
    
    mask = ismember(atlas,1115:1305); atlas(mask) = 1115;
    mask = ismember(atlas,1306:1317); atlas(mask) = 1306;
    mask = ismember(atlas,1028:1114); atlas(mask) = 1;
    brainoutline = atlas>0;
    B = bwboundaries(fliplr(squeeze(brainoutline(:,plane,:))));
    atlasmaskplot2 = atlas==1306;
    B2a = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    atlasmaskplot2 = atlas==1115;
    B2b = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    atlasmaskplot2 = atlas==1;
    B2d = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    
    data = [KDE.Consume.Novel(:,plane,:),KDE.Consume.Familiar(:,plane,:)];
    data = squeeze(data(:,1,:))-squeeze(data(:,2,:));
    colormap(gcf,cmap)
    heatmap(rot90(data),[],[],[],'UseFigureColormap',true,'ColorLevels',1000,'MaxColorValue',max(lims2),'MinColorValue',min(lims2),'NaNColor',[1 1 1]);
    for i = 1:size(data_in.regions.index,1)
        if ~ismember(i,idx_amygdala)
            idx = data_in.regions.index(i)+1;
            mask = squeeze(atlas(:,plane,:)) == idx;
            B3 = bwboundaries(fliplr(mask));
            for j = 1:length(B3)
                plot(B3{j}(:,1),B3{j}(:,2),'Color',[.85 .85 .85],'LineWidth',1)
            end
        end
    end
    for i = 1:size(data_in.regions.index,1)
        if ismember(i,idx_amygdala)
            idx = data_in.regions.index(i)+1;
            mask = squeeze(atlas(:,plane,:)) == idx;
            B3 = bwboundaries(fliplr(mask));
            for j = 1:length(B3)
                plot(B3{j}(:,1),B3{j}(:,2),'Color',[228 45 38]/255,'LineWidth',1)
            end
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
    title('Consume')
    set(gca,'FontSize',12)
    xlim([.5 228.5])
    axis off
    hold off
    
    subplot(2,3,2+(pl-1)*3)
    hold on
    axis equal
    
    mask = ismember(atlas,1115:1305); atlas(mask) = 1115;
    mask = ismember(atlas,1306:1317); atlas(mask) = 1306;
    mask = ismember(atlas,1028:1114); atlas(mask) = 1;
    brainoutline = atlas>0;
    B = bwboundaries(fliplr(squeeze(brainoutline(:,plane,:))));
    atlasmaskplot2 = atlas==1306;
    B2a = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    atlasmaskplot2 = atlas==1115;
    B2b = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    atlasmaskplot2 = atlas==1;
    B2d = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    data = [KDE.Malaise.Novel(:,plane,:),KDE.Malaise.Familiar(:,plane,:)];
    data = squeeze(data(:,1,:))-squeeze(data(:,2,:));
    colormap(gcf,cmap)
    heatmap(rot90(data),[],[],[],'UseFigureColormap',true,'ColorLevels',1000,'MaxColorValue',max(lims2),'MinColorValue',min(lims2),'NaNColor',[1 1 1]);
    for i = 1:size(data_in.regions.index,1)
        if ~ismember(i,idx_amygdala)
            idx = data_in.regions.index(i)+1;
            mask = squeeze(atlas(:,plane,:)) == idx;
            B3 = bwboundaries(fliplr(mask));
            for j = 1:length(B3)
                plot(B3{j}(:,1),B3{j}(:,2),'Color',[.85 .85 .85],'LineWidth',1)
            end
        end
    end
    for i = 1:size(data_in.regions.index,1)
        if ismember(i,idx_amygdala)
            idx = data_in.regions.index(i)+1;
            mask = squeeze(atlas(:,plane,:)) == idx;
            B3 = bwboundaries(fliplr(mask));
            for j = 1:length(B3)
                plot(B3{j}(:,1),B3{j}(:,2),'Color',[228 45 38]/255,'LineWidth',1)
            end
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
    title(['Malaise'])
    set(gca,'FontSize',12)
    xlim([.5 228.5])
    axis off
    hold off
    
    subplot(2,3,3+(pl-1)*3)
    hold on
    axis equal
    
    mask = ismember(atlas,1115:1305); atlas(mask) = 1115;
    mask = ismember(atlas,1306:1317); atlas(mask) = 1306;
    mask = ismember(atlas,1028:1114); atlas(mask) = 1;
    brainoutline = atlas>0;
    B = bwboundaries(fliplr(squeeze(brainoutline(:,plane,:))));
    atlasmaskplot2 = atlas==1306;
    B2a = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    atlasmaskplot2 = atlas==1115;
    B2b = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    atlasmaskplot2 = atlas==1;
    B2d = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
    data = [KDE.Retrieval.Novel(:,plane,:),KDE.Retrieval.Familiar(:,plane,:)];
    data = squeeze(data(:,1,:))-squeeze(data(:,2,:));
    colormap(gcf,cmap)
    heatmap(rot90(data),[],[],[],'UseFigureColormap',true,'ColorLevels',1000,'MaxColorValue',max(lims2),'MinColorValue',min(lims2),'NaNColor',[1 1 1]);
    for i = 1:size(data_in.regions.index,1)
        if ~ismember(i,idx_amygdala)
            idx = data_in.regions.index(i)+1;
            mask = squeeze(atlas(:,plane,:)) == idx;
            B3 = bwboundaries(fliplr(mask));
            for j = 1:length(B3)
                plot(B3{j}(:,1),B3{j}(:,2),'Color',[.85 .85 .85],'LineWidth',1)
            end
        end
    end
    for i = 1:size(data_in.regions.index,1)
        if ismember(i,idx_amygdala)
            idx = data_in.regions.index(i)+1;
            mask = squeeze(atlas(:,plane,:)) == idx;
            B3 = bwboundaries(fliplr(mask));
            for j = 1:length(B3)
                plot(B3{j}(:,1),B3{j}(:,2),'Color',[228 45 38]/255,'LineWidth',1)
            end
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
    title(['Retrieval'])
    set(gca,'FontSize',12)
    xlim([.5 228.5])
    axis off
    hold off
    
end
%% Fig 1i

figure('Position', get(0, 'Screensize'))
sgtitle('Figure 1i','FontWeight','bold')
T1 = readtable([data_path,'/FOS imaging/region_info.csv']);
T2 = readtable([data_path,'/FOS imaging/sample_info.csv']);

axis square
hold on
counts_norm = [T1{find(cellfun(@(x) isequal(x,'CentralAmygdalarNucleus'),T1.region)),4:end}./T1{1,4:end}./1.987]*100;
phaselist = {'Consumption','Malaise','Retrieval'};
for j = 1:length(phaselist)
    
    PHASEidx = find(cellfun(@(x) isequal(phaselist{j},x),T2.Timepoint)&cellfun(@(x) isequal('Novel',x),T2.Novel));
    x = counts_norm(PHASEidx);
    xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
    scatter(ones(1,length(x))*j+rand(1,length(x))*.1-.05+.125,x,100,'filled','MarkerFaceColor',[252 216 213]/255,'MarkerEdgeColor','w')
    a = plot([-0.125 0.125]+j+.125,[xx(3) xx(3)],'color',[229 45 38]/255,'LineWidth',1);
    plot([j j]+.125,[xx(2) xx(4)],'color',[229 45 38]/255,'LineWidth',1)
    plot([j j]+.125,[xx(1) xx(5)],'color',[229 45 38]/255,'LineWidth',1)
    PHASEidx = find(cellfun(@(x) isequal(phaselist{j},x),T2.Timepoint)&cellfun(@(x) isequal('Familiar',x),T2.Novel));
    x = counts_norm(PHASEidx);
    xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
    scatter(ones(1,length(x))*j+rand(1,length(x))*.1-.05-.125,x,100,'filled','MarkerFaceColor',[216 231 243]/255,'MarkerEdgeColor','w')
    b = plot([-0.125 0.125]+j-.125,[xx(3) xx(3)],'color',[55 136 193]/255,'LineWidth',1);
    plot([j j]-.125,[xx(2) xx(4)],'color',[55 136 193]/255,'LineWidth',1)
    plot([j j]-.125,[xx(1) xx(5)],'color',[55 136 193]/255,'LineWidth',1)
end
xlim([0.5 3.5]); xticks(1:3); xticklabels({'Consume','Malaise','Retrieval'});
ylabel('CEA FOS (% per mm^3)'); ylim([0 1]); yticks(0:.25:1);
legend([b,a],{'Familiar','Novel'})
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

load([data_path,'/FOS imaging/FOS-GLMM-statistics.mat'],'data')
StatsTbl(end+1,:) = table({'1i'},{'Consume: Novel vs. Familiar'},{'GLMM marginal effect'},{'3 time points'},{[12 12]},data.GLMMoutput.CEA.Eq2.flavor.Zstat(1),data.GLMMoutput.CEA.Eq2.flavor.pvalues_corrected(1));
StatsTbl(end+1,:) = table({'1i'},{'Malaise: Novel vs. Familiar'},{'GLMM marginal effect'},{'3 time points'},{[12 12]},data.GLMMoutput.CEA.Eq2.flavor.Zstat(2),data.GLMMoutput.CEA.Eq2.flavor.pvalues_corrected(2));
StatsTbl(end+1,:) = table({'1i'},{'Retrieval: Novel vs. Familiar'},{'GLMM marginal effect'},{'3 time points'},{[12 12]},data.GLMMoutput.CEA.Eq2.flavor.Zstat(3),data.GLMMoutput.CEA.Eq2.flavor.pvalues_corrected(3));

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