function fig_ED1(data_path)

disp('Generating panels for Extended Data Figure 1...')
%% Fig ED1a

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 1a','FontWeight','bold')
load([data_path,'/Fos imaging/Fos-GLMM-statistics.mat'],'data')

subplot(1,2,1)
counts_norm = [data.GLMMinput.counts./data.GLMMinput.offset./data.regions.size']*100;
idx_nov = find([data.GLMMoutput.Eq2.flavor.pvalues_corrected(:,1)<=0.05].*[data.GLMMoutput.Eq2.flavor.estimates(:,1)>0]);
T = readtable([data_path,'/source data/Fig-1e.csv']);
hold on
axis square
for j = 1:length(idx_nov)
    idx = find(cellfun(@(x) isequal('Consumption',x),data.GLMMinput.phase)&cellfun(@(x) isequal('Novel',x),data.GLMMinput.flavor));
    x = counts_norm(idx,idx_nov(j));
    xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
    scatter(ones(1,length(x))*j+rand(1,length(x))*.1-.05+.125,x,100,'filled','MarkerFaceColor',[252 216 213]/255,'MarkerEdgeColor','w')
    a = plot([-0.125 0.125]+j+.125,[xx(3) xx(3)],'color',[229 45 38]/255,'LineWidth',1);
    plot([j j]+.125,[xx(2) xx(4)],'color',[229 45 38]/255,'LineWidth',1)
    plot([j j]+.125,[xx(1) xx(5)],'color',[229 45 38]/255,'LineWidth',1)
    
    idx = find(cellfun(@(x) isequal('Consumption',x),data.GLMMinput.phase)&cellfun(@(x) isequal('Familiar',x),data.GLMMinput.flavor));
    x = counts_norm(idx,idx_nov(j));
    xx = prctile(x,[10 25 50 75 90]);
    xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
    scatter(ones(1,length(x))*j+rand(1,length(x))*.1-.05-.125,x,100,'filled','MarkerFaceColor',[216 231 243]/255,'MarkerEdgeColor','w')
    b = plot([-0.125 0.125]+j-.125,[xx(3) xx(3)],'color',[55 136 193]/255,'LineWidth',1);
    plot([j j]-.125,[xx(2) xx(4)],'color',[55 136 193]/255,'LineWidth',1)
    plot([j j]-.125,[xx(1) xx(5)],'color',[55 136 193]/255,'LineWidth',1)
    TickLabels(j) = T.Abbreviation(find(cellfun(@(x) isequal(x,data.regions.name(idx_nov(j),:)),T.Name)),:);
end
xticks(1:length(idx_nov))
xticklabels(TickLabels)
ylim([0 1])
yticks(0:.5:1)
xlim([0.5 length(idx_nov)+.5])
ylabel('Fos (% per mm^3)')
legend([b,a],{'Familiar','Novel'})
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

StatsTbl = table({'ED 1a'},{[TickLabels{1},': Novel vs. Familiar']},{'GLMM marginal effect'},{'3 timepoints'},{[12 12]},data.GLMMoutput.Eq2.flavor.Zstat(idx_nov(1),1),data.GLMMoutput.Eq2.flavor.pvalues_corrected(idx_nov(1),1), ...
    'VariableNames',{'Figure panel','Group','Statistical test','Multiple comparisons','Sample size','Test statistic','P-value'});
for i = 2:length(idx_nov)
    StatsTbl(end+1,:) = table({'ED 1a'},{[TickLabels{i},': Novel vs. Familiar']},{'GLMM marginal effect'},{'3 timepoints'},{[12 12]},data.GLMMoutput.Eq2.flavor.Zstat(idx_nov(i),1),data.GLMMoutput.Eq2.flavor.pvalues_corrected(idx_nov(i),1));
end

subplot(1,2,2)
load([data_path,'/Fos imaging/modified-atlas/allen_ccfv3_modified_cz.mat'],'atlas','RegionLibrary')
hold on
for i = 1:length(idx_nov)
    t = atlas==(RegionLibrary.reduced.index(idx_nov(i))+1);
    p = patch(isosurface(t));
    p.FaceColor = [229 45 38]/255;
    p.FaceAlpha = 0.4;
    p.LineStyle = 'none';
end
p = patch(isosurface(atlas>0));
p.FaceAlpha = 0.05;
p.FaceColor = [44 31 22]/255;
p.LineStyle = 'none';
h = get(gca,'DataAspectRatio');
set(gca,'DataAspectRatio',[1 1 h(3)],'CameraViewAngleMode','Manual');
set(gca,'Ydir','reverse','Zdir','reverse')
axis off
xlim('auto')
ylim('auto')
zlim('auto')
view(-20,20)
hold off
%% Fig ED1b

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 1b','FontWeight','bold')
load([data_path,'/Fos imaging/Fos-GLMM-statistics.mat'],'data')

subplot(1,2,1)
counts_norm = [data.GLMMinput.counts./data.GLMMinput.offset./data.regions.size']*100;
idx_fam = find([data.GLMMoutput.Eq2.flavor.pvalues_corrected(:,1)<=0.05].*[data.GLMMoutput.Eq2.flavor.estimates(:,1)<0]);
T = readtable([data_path,'/source data/Fig-1e.csv']);
hold on
axis square
for j = 1:length(idx_fam)
    idx = find(cellfun(@(x) isequal('Consumption',x),data.GLMMinput.phase)&cellfun(@(x) isequal('Novel',x),data.GLMMinput.flavor));
    x = counts_norm(idx,idx_fam(j));
    xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
    scatter(ones(1,length(x))*j+rand(1,length(x))*.1-.05+.125,x,100,'filled','MarkerFaceColor',[252 216 213]/255,'MarkerEdgeColor','w')
    a = plot([-0.125 0.125]+j+.125,[xx(3) xx(3)],'color',[229 45 38]/255,'LineWidth',1);
    plot([j j]+.125,[xx(2) xx(4)],'color',[229 45 38]/255,'LineWidth',1)
    plot([j j]+.125,[xx(1) xx(5)],'color',[229 45 38]/255,'LineWidth',1)
    
    idx = find(cellfun(@(x) isequal('Consumption',x),data.GLMMinput.phase)&cellfun(@(x) isequal('Familiar',x),data.GLMMinput.flavor));
    x = counts_norm(idx,idx_fam(j));
    xx = prctile(x,[10 25 50 75 90]);
    xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
    scatter(ones(1,length(x))*j+rand(1,length(x))*.1-.05-.125,x,100,'filled','MarkerFaceColor',[216 231 243]/255,'MarkerEdgeColor','w')
    b = plot([-0.125 0.125]+j-.125,[xx(3) xx(3)],'color',[55 136 193]/255,'LineWidth',1);
    plot([j j]-.125,[xx(2) xx(4)],'color',[55 136 193]/255,'LineWidth',1)
    plot([j j]-.125,[xx(1) xx(5)],'color',[55 136 193]/255,'LineWidth',1)
    TickLabels(j) = T.Abbreviation(find(cellfun(@(x) isequal(x,data.regions.name(idx_fam(j),:)),T.Name)),:);
end
xticks(1:length(idx_fam))
xticklabels(TickLabels)
ylim([0 2])
yticks(0:1:2)
xlim([0.5 length(idx_fam)+.5])
ylabel('Fos (% per mm^3)')
legend([b,a],{'Familiar','Novel'})
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

for i = 1:length(idx_fam)
    StatsTbl(end+1,:) = table({'ED 1b'},{[TickLabels{i},': Novel vs. Familiar']},{'GLMM marginal effect'},{'3 timepoints'},{[12 12]},data.GLMMoutput.Eq2.flavor.Zstat(idx_fam(i),1),data.GLMMoutput.Eq2.flavor.pvalues_corrected(idx_fam(i),1));
end

subplot(1,2,2)
load([data_path,'/Fos imaging/modified-atlas/allen_ccfv3_modified_cz.mat'],'atlas','RegionLibrary')
hold on
for i = 1:length(idx_fam)
    t = atlas==(RegionLibrary.reduced.index(idx_fam(i))+1);
    p = patch(isosurface(t));
    p.FaceColor = [55 136 193]/255;
    p.FaceAlpha = 0.4;
    p.LineStyle = 'none';
end
p = patch(isosurface(atlas>0));
p.FaceAlpha = 0.05;
p.FaceColor = [44 31 22]/255;
p.LineStyle = 'none';
h = get(gca,'DataAspectRatio');
set(gca,'DataAspectRatio',[1 1 h(3)],'CameraViewAngleMode','Manual');
set(gca,'Ydir','reverse','Zdir','reverse')
axis off
xlim('auto')
ylim('auto')
zlim('auto')
view(-20,20)
hold off
%% Fig ED1c

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 1c','FontWeight','bold')

fname = '/Fos imaging/kernel-density-estimates/kde-consumption-novel.npy';
KDE.Consume.Novel = readNPY([data_path,fname]);
fname = '/Fos imaging/kernel-density-estimates/kde-consumption-familiar.npy';
KDE.Consume.Familiar = readNPY([data_path,fname]);
fname = '/Fos imaging/kernel-density-estimates/kde-malaise-novel.npy';
KDE.Malaise.Novel = readNPY([data_path,fname]);
fname = '/Fos imaging/kernel-density-estimates/kde-malaise-familiar.npy';
KDE.Malaise.Familiar = readNPY([data_path,fname]);
fname = '/Fos imaging/kernel-density-estimates/kde-retrieval-novel.npy';
KDE.Retrieval.Novel = readNPY([data_path,fname]);
fname = '/Fos imaging/kernel-density-estimates/kde-retrieval-familiar.npy';
KDE.Retrieval.Familiar = readNPY([data_path,fname]);
load([data_path,'/Fos imaging/modified-atlas/allen_ccfv3_modified_cz.mat'],'atlas','RegionLibrary')
atlasmask  = atlas>=1028 | atlas<=1;
KDE.Consume.Novel(atlasmask) = NaN;
KDE.Consume.Familiar(atlasmask) = NaN;
KDE.Malaise.Novel(atlasmask) = NaN;
KDE.Malaise.Familiar(atlasmask) = NaN;
KDE.Retrieval.Novel(atlasmask) = NaN;
KDE.Retrieval.Familiar(atlasmask) = NaN;
cmap = flipud(cbrewer('div','RdBu',1000,'spline')); cmap(cmap<0) = 0;
lims1 = [0 1];
lims2 = [-0.5 0.5];
cmap1 = cbrewer('seq','Greys',1000,'spline');
cmap2 = flipud(cbrewer('div','RdBu',1000,'spline')); cmap2(cmap2<0) = 0;
planes = 99:20:479;

for pl = 1:length(planes)
    plane = planes(pl);
    subplot(2,10,pl)
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
    data = squeeze(mean(data,2));
    colormap(gcf,cmap1)
    heatmap(rot90(data),[],[],[],'UseFigureColormap',true,'ColorLevels',1000,'MaxColorValue',max(lims1),'MinColorValue',min(lims1),'NaNColor',[1 1 1]);
    for i = 1:size(RegionLibrary.reduced,1)
        idx = RegionLibrary.reduced{i,1}+1;
        mask = squeeze(atlas(:,plane,:)) == idx;
        B3 = bwboundaries(fliplr(mask));
        for j = 1:length(B3)
            plot(B3{j}(:,1),B3{j}(:,2),'Color',[.85 .85 .85],'LineWidth',1)
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
    set(gca,'FontSize',12)
    xlim([.5 228.5])
    axis off
    hold off
end
%% Fig ED1d

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 1d','FontWeight','bold')

for pl = 1:length(planes)
    plane = planes(pl);
    subplot(2,10,pl)
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
    colormap(gcf,cmap2)
    heatmap(rot90(data),[],[],[],'UseFigureColormap',true,'ColorLevels',1000,'MaxColorValue',max(lims2),'MinColorValue',min(lims2),'NaNColor',[1 1 1]);
    for i = 1:size(RegionLibrary.reduced,1)
        idx = RegionLibrary.reduced{i,1}+1;
        mask = squeeze(atlas(:,plane,:)) == idx;
        B3 = bwboundaries(fliplr(mask));
        for j = 1:length(B3)
            plot(B3{j}(:,1),B3{j}(:,2),'Color',[.85 .85 .85],'LineWidth',1)
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
    set(gca,'FontSize',12)
    xlim([.5 228.5])
    axis off
    hold off
end
%% Fig ED1e

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 1e','FontWeight','bold')

for pl = 1:length(planes)
    plane = planes(pl);
    subplot(2,10,pl)
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
    data = squeeze(mean(data,2));
    colormap(gcf,cmap1)
    heatmap(rot90(data),[],[],[],'UseFigureColormap',true,'ColorLevels',1000,'MaxColorValue',max(lims1),'MinColorValue',min(lims1),'NaNColor',[1 1 1]);
    for i = 1:size(RegionLibrary.reduced,1)
        idx = RegionLibrary.reduced{i,1}+1;
        mask = squeeze(atlas(:,plane,:)) == idx;
        B3 = bwboundaries(fliplr(mask));
        for j = 1:length(B3)
            plot(B3{j}(:,1),B3{j}(:,2),'Color',[.85 .85 .85],'LineWidth',1)
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
    set(gca,'FontSize',12)
    xlim([.5 228.5])
    axis off
    hold off
end
%% Fig ED1f

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 1f','FontWeight','bold')

for pl = 1:length(planes)
    plane = planes(pl);
    subplot(2,10,pl)
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
    colormap(gcf,cmap2)
    heatmap(rot90(data),[],[],[],'UseFigureColormap',true,'ColorLevels',1000,'MaxColorValue',max(lims2),'MinColorValue',min(lims2),'NaNColor',[1 1 1]);
    for i = 1:size(RegionLibrary.reduced,1)
        idx = RegionLibrary.reduced{i,1}+1;
        mask = squeeze(atlas(:,plane,:)) == idx;
        B3 = bwboundaries(fliplr(mask));
        for j = 1:length(B3)
            plot(B3{j}(:,1),B3{j}(:,2),'Color',[.85 .85 .85],'LineWidth',1)
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
    set(gca,'FontSize',12)
    xlim([.5 228.5])
    axis off
    hold off
end
%% Fig ED1g

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 1g','FontWeight','bold')

for pl = 1:length(planes)
    plane = planes(pl);
    subplot(2,10,pl)
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
    data = squeeze(mean(data,2));
    colormap(gcf,cmap1)
    heatmap(rot90(data),[],[],[],'UseFigureColormap',true,'ColorLevels',1000,'MaxColorValue',max(lims1),'MinColorValue',min(lims1),'NaNColor',[1 1 1]);
    for i = 1:size(RegionLibrary.reduced,1)
        idx = RegionLibrary.reduced{i,1}+1;
        mask = squeeze(atlas(:,plane,:)) == idx;
        B3 = bwboundaries(fliplr(mask));
        for j = 1:length(B3)
            plot(B3{j}(:,1),B3{j}(:,2),'Color',[.85 .85 .85],'LineWidth',1)
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
    set(gca,'FontSize',12)
    xlim([.5 228.5])
    axis off
    hold off
end
%% Fig ED1h

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 1h','FontWeight','bold')

for pl = 1:length(planes)
    plane = planes(pl);
    subplot(2,10,pl)
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
    colormap(gcf,cmap2)
    heatmap(rot90(data),[],[],[],'UseFigureColormap',true,'ColorLevels',1000,'MaxColorValue',max(lims2),'MinColorValue',min(lims2),'NaNColor',[1 1 1]);
    for i = 1:size(RegionLibrary.reduced,1)
        idx = RegionLibrary.reduced{i,1}+1;
        mask = squeeze(atlas(:,plane,:)) == idx;
        B3 = bwboundaries(fliplr(mask));
        for j = 1:length(B3)
            plot(B3{j}(:,1),B3{j}(:,2),'Color',[.85 .85 .85],'LineWidth',1)
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
    set(gca,'FontSize',12)
    xlim([.5 228.5])
    axis off
    hold off
end

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