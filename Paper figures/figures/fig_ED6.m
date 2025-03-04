function fig_ED6(data_path)

disp('Generating panels for Extended Data Figure 6...')
%% Fig ED6b

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 6b','FontWeight','bold')
load([data_path,'/Fos imaging/Fos-GLMM-statistics.mat'],'data');

hold on
axis square
counts_norm = [data.GLMMinput.counts./data.GLMMinput.offset./data.regions.size']*100;
idx = find(data.regions.name=='Parabrachial nucleus');
phases = {'Consumption','Malaise','CGRP stim','Retrieval'};
for j = 1:length(phases)
    PHASEidx = find(cellfun(@(x) isequal(phases{j},x),data.GLMMinput.phase));
    x = counts_norm(PHASEidx,idx);
    xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
    scatter(ones(1,length(x))*j+rand(1,length(x))*.1-.05,x,64,'filled','MarkerFaceColor',[.85 .85 .85],'MarkerEdgeColor','w')
    plot([-0.125 0.125]+j,[xx(3) xx(3)],'k','LineWidth',1)
    plot([j j],[xx(2) xx(4)],'k','LineWidth',1)
    plot([j j],[xx(1) xx(5)],'k','LineWidth',1)
end
xlim([.5 4.5]); xticks(1:4); xticklabels({'Consume','Malaise','CGRP stim','Retrieval'})
ylabel('PB Fos (% per mm^3)'); ylim([0 2.5]); yticks(0:.5:2.5);
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off
%% Fig ED6c

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 6c','FontWeight','bold')
load([data_path,'/Fos imaging/Fos-GLMM-statistics.mat'],'data')

T = readtable([data_path,'/source data/Fig-1e.csv']);
idx = find(data.GLMMoutput.Eq2.modelstats.significant);
idx_amygdala = idx(find(T.Cluster==1))';
idx_other = setdiff(find(data.GLMMoutput.Eq2.modelstats.significant),[find(data.regions.name=='Parabrachial nucleus');idx_amygdala]);

subplot(1,2,1)
hold on
axis square
A = data.GLMMoutput.Eq3.coefficients.Zstat(idx_amygdala,2)';
B = data.GLMMoutput.Eq3.coefficients.Zstat(idx_amygdala,3)';
a = scatter(A,B,64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([-4 3])
ylim([-4 3])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[y_fit] = polyval(p2,x,S);
fitresult = fit(A',B','poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[252 216 213]/255,'LineStyle','none'),
plot(x,y_fit,'color',[229 45 38]/255,'LineWidth',1)
scatter(A,B,100,[229 45 38]/255,'filled','MarkerEdgeColor','w')
xlim([-4 3])
ylim([-4 3])
xticks(xlim)
yticks(ylim)
xlabel(['Malaise',char(10),'Average Fos (Z)'])
ylabel(['CGRP stim',char(10),'Average Fos (Z)'])
title('Amygdala network')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

[r,p] = corr(A',B');
StatsTbl = table({'ED 6c, top'},{'Malaise vs. CGRP stim'},{'Pearson correlation'},{'N/A'},{[length(A)]},r,p, ...
    'VariableNames',{'Figure panel','Group','Statistical test','Multiple comparisons','Sample size','Test statistic','P-value'});

subplot(1,2,2)
hold on
axis square
A = data.GLMMoutput.Eq3.coefficients.Zstat(idx_other,2)';
B = data.GLMMoutput.Eq3.coefficients.Zstat(idx_other,3)';
a = scatter(A,B,64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([-13 2])
ylim([-13 2])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[y_fit] = polyval(p2,x,S);
fitresult = fit(A',B','poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[.85 .85 .85],'LineStyle','none'),
plot(x,y_fit,'k','LineWidth',1)
scatter(A,B,100,'k','filled','MarkerEdgeColor','w')
xlim([-13 2])
ylim([-13 2])
xticks(xlim)
yticks(ylim)
xlabel(['Malaise',char(10),'Average Fos (Z)'])
ylabel(['CGRP stim',char(10),'Average Fos (Z)'])
title('Other regions')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

[r,p] = corr(A',B');
StatsTbl(end+1,:) = table({'ED 6c, bottom'},{'Malaise vs. CGRP stim'},{'Pearson correlation'},{'N/A'},{[length(A)]},r,p);
%% Fig ED6d

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 6d','FontWeight','bold')
load([data_path,'/Fos imaging/Fos-GLMM-statistics.mat'],'data')

T = readtable([data_path,'/source data/Fig-1e.csv']);
idx = find(data.GLMMoutput.Eq2.modelstats.significant);
idx_amygdala = idx(find(T.Cluster==1))';
idx_other = setdiff(find(data.GLMMoutput.Eq2.modelstats.significant),[find(data.regions.name=='Parabrachial nucleus');idx_amygdala]);

subplot(1,2,1)
hold on
axis square
A = data.GLMMoutput.Eq4.flavor.Zstat(idx_amygdala,1)';
B = data.GLMMoutput.Eq4.flavor.Zstat(idx_amygdala,2)';
a = scatter(A,B,64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([-1 3])
ylim([0 5])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[y_fit] = polyval(p2,x,S);
fitresult = fit(A',B','poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[252 216 213]/255,'LineStyle','none'),
plot(x,y_fit,'color',[229 45 38]/255,'LineWidth',1)
scatter(A,B,100,[229 45 38]/255,'filled','MarkerEdgeColor','w')
xlim([-1 3])
ylim([0 5])
xticks(xlim)
yticks(ylim)
xlabel(['Malaise',char(10),'Novel – Familiar ΔFos (Z)'])
ylabel(['CGRP stim',char(10),'Novel – Familiar ΔFos (Z)'])
title('Amygdala network')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

[r,p] = corr(A',B');
StatsTbl(end+1,:) = table({'ED 6D, top'},{'Malaise vs. CGRP stim'},{'Pearson correlation'},{'N/A'},{[length(A)]},r,p);

subplot(1,2,2)
hold on
axis square
A = data.GLMMoutput.Eq4.flavor.Zstat(idx_other,1)';
B = data.GLMMoutput.Eq4.flavor.Zstat(idx_other,2)';
a = scatter(A,B,64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([-3 3])
ylim([-4 5])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[y_fit] = polyval(p2,x,S);
fitresult = fit(A',B','poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[.85 .85 .85],'LineStyle','none'),
plot(x,y_fit,'k','LineWidth',1)
scatter(A,B,100,'k','filled','MarkerEdgeColor','w')
xlim([-3 3])
ylim([-4 5])
xticks(xlim)
yticks(ylim)
xlabel(['Malaise',char(10),'Novel – Familiar ΔFos (Z)'])
ylabel(['CGRP stim',char(10),'Novel – Familiar ΔFos (Z)'])
title('Other regions')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

[r,p] = corr(A',B');
StatsTbl(end+1,:) = table({'ED 6D, bottom'},{'Malaise vs. CGRP stim'},{'Pearson correlation'},{'N/A'},{[length(A)]},r,p);
%% Fig ED6e

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 6e','FontWeight','bold')

fname = '/Fos imaging/kernel-density-estimates/kde-cgrp-stim-novel.npy';
KDE.CGRP.Novel = readNPY([data_path,fname]);
fname = '/Fos imaging/kernel-density-estimates/kde-cgrp-stim-familiar.npy';
KDE.CGRP.Familiar = readNPY([data_path,fname]);
load([data_path,'/Fos imaging/modified-atlas/allen_ccfv3_modified_cz.mat'],'atlas','RegionLibrary')
atlasmask  = atlas>=1028 | atlas<=1;
KDE.CGRP.Novel(atlasmask) = NaN;
KDE.CGRP.Familiar(atlasmask) = NaN;
KDE.Malaise.Novel(atlasmask) = NaN;
KDE.Malaise.Familiar(atlasmask) = NaN;
KDE.Retrieval.Novel(atlasmask) = NaN;
KDE.Retrieval.Familiar(atlasmask) = NaN;
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
    data = [KDE.CGRP.Novel(:,plane,:),KDE.CGRP.Familiar(:,plane,:)];
    data = squeeze(mean(data,2));
    colormap(gcf,cmap1)
    heatmap(rot90(data),[],[],[],'UseFigureColormap',true,'ColorLevels',1000,'MaxColorValue',200,'MinColorValue',0,'NaNColor',[1 1 1]);
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
%% Fig ED6f

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 6f','FontWeight','bold')

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
    data = [KDE.CGRP.Novel(:,plane,:),KDE.CGRP.Familiar(:,plane,:)];
    data = squeeze(data(:,1,:))-squeeze(data(:,2,:));
    colormap(gcf,cmap2)
    heatmap(rot90(data),[],[],[],'UseFigureColormap',true,'ColorLevels',1000,'MaxColorValue',100,'MinColorValue',-100,'NaNColor',[1 1 1]);
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