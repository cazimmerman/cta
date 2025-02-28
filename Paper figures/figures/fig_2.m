function fig_2(data_path)

disp('Generating panels for Figure 2...')
%% Fig 2b

figure('Position', get(0, 'Screensize'))
sgtitle('Figure 2b','FontWeight','bold')
load([data_path,'\Photometry\cgrp-photometry-licl.mat'],'data')

hold on
axis square
fill([data.time fliplr(data.time)],[mean(data.GCaMP)+std(data.GCaMP)/sqrt(size(data.GCaMP,1)) fliplr(mean(data.GCaMP)-std(data.GCaMP)/sqrt(size(data.GCaMP,1)))],'k','LineStyle','none','FaceAlpha',0.1)
plot(data.time,mean(data.GCaMP),'k','LineWidth',1)
xlabel('Time (min)'); xlim([-10 30]); xticks(-10:10:30);
ylabel('CGRP activity (σ)'); ylim([-1 5]); yticks(-1:1:5);
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off
%% Fig 2d

figure('Position', get(0, 'Screensize'))
sgtitle('Figure 2d','FontWeight','bold')
T = readtable([data_path,'\source data\Fig-2d.csv']);

hold on
axis square
idx = find(cellfun(@(x) isequal(x,'Familiar'),T.Group));
x = T.Preference(idx)*100;
xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
scatter(ones(1,length(x))+rand(1,length(x))*.1-0.05,x,100,'filled','MarkerFaceColor',[216 231 243]/255,'MarkerEdgeColor','w')
plot([-0.125 0.125]+1,[xx(3) xx(3)],'color',[55 136 193]/255,'LineWidth',1)
plot([1 1],[xx(2) xx(4)],'color',[55 136 193]/255,'LineWidth',1)
idx = find(cellfun(@(x) isequal(x,'Novel'),T.Group));
x = T.Preference(idx)*100;
xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
scatter(ones(1,length(x))*2+rand(1,length(x))*.1-0.05,x,100,'filled','MarkerFaceColor',[252 216 213]/255,'MarkerEdgeColor','w')
plot([-0.125 0.125]+2,[xx(3) xx(3)],'color',[229 45 38]/255,'LineWidth',1)
plot([1 1]*2,[xx(2) xx(4)],'color',[229 45 38]/255,'LineWidth',1)
xlim([0.5 2.5]); xticks(1:2); xticklabels({'Familiar','Novel'});
ylabel('Flavor preference (%)'); ylim([0 100]); yticks(0:50:100);
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

[p,~,stats] = ranksum(T.Preference(find(cellfun(@(x) isequal(x,'Novel'),T.Group))),T.Preference(find(cellfun(@(x) isequal(x,'Familiar'),T.Group))),'method','exact');
StatsTbl = table({'2d'},{'Novel vs. Familiar'},{'Wilcoxon rank-sum'},{'N/A'},{[length(find(cellfun(@(x) isequal(x,'Novel'),T.Group))) length(find(cellfun(@(x) isequal(x,'Familiar'),T.Group)))]},stats.ranksum,p(1), ...
    'VariableNames',{'Figure panel','Group','Statistical test','Multiple comparisons','Sample size','Test statistic','P-value'});
%% Fig 2e

figure('Position', get(0, 'Screensize'))
sgtitle('Figure 2e','FontWeight','bold')
T = readtable([data_path,'\source data\Fig-2e.csv']);

hold on
axis square
idx = find(cellfun(@(x) isequal(x,'YFP'),T.Group));
x = T.Preference(idx)*100;
xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
scatter(ones(1,length(x))+rand(1,length(x))*.1-0.05,x,100,'filled','MarkerFaceColor',[.85 .85 .85],'MarkerEdgeColor','w')
plot([-0.125 0.125]+1,[xx(3) xx(3)],'color','k','LineWidth',1)
plot([1 1],[xx(2) xx(4)],'color','k','LineWidth',1)
plot([1 1],[xx(1) xx(5)],'color','k','LineWidth',1)
idx = find(cellfun(@(x) isequal(x,'ChRmine'),T.Group));
x = T.Preference(idx)*100;
xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
scatter(ones(1,length(x))*2+rand(1,length(x))*.1-0.05,x,100,'filled','MarkerFaceColor',[252 216 213]/255,'MarkerEdgeColor','w')
plot([-0.125 0.125]+2,[xx(3) xx(3)],'color',[229 45 38]/255,'LineWidth',1)
plot([1 1]*2,[xx(2) xx(4)],'color',[229 45 38]/255,'LineWidth',1)
plot([1 1]*2,[xx(1) xx(5)],'color',[229 45 38]/255,'LineWidth',1)
xlim([0.5 2.5]); xticks(1:2); xticklabels({'YFP','ChRmine'});
ylabel('Flavor preference (%)'); ylim([0 100]); yticks(0:50:100);
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

[p,~,stats] = ranksum(T.Preference(find(cellfun(@(x) isequal(x,'ChRmine'),T.Group))),T.Preference(find(cellfun(@(x) isequal(x,'YFP'),T.Group))),'method','exact');
StatsTbl(end+1,:) = table({'2e'},{'ChRmine vs. YFP'},{'Wilcoxon rank-sum'},{'N/A'},{[length(find(cellfun(@(x) isequal(x,'ChRmine'),T.Group))) length(find(cellfun(@(x) isequal(x,'YFP'),T.Group)))]},stats.ranksum,p(1));
%% Fig 2f

figure('Position', get(0, 'Screensize'))
sgtitle('Figure 2f','FontWeight','bold')
T = readtable([data_path,'\source data\Fig-2f.csv']);

hold on
axis square
idx = find(cellfun(@(x) isequal(x,'YFP'),T.Group));
x = T.Preference(idx)*100;
xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
scatter(ones(1,length(x))+rand(1,length(x))*.1-0.05,x,100,'filled','MarkerFaceColor',[.85 .85 .85],'MarkerEdgeColor','w')
plot([-0.125 0.125]+1,[xx(3) xx(3)],'color','k','LineWidth',1)
plot([1 1],[xx(2) xx(4)],'color','k','LineWidth',1)
plot([1 1],[xx(1) xx(5)],'color','k','LineWidth',1)
idx = find(cellfun(@(x) isequal(x,'eOPN3'),T.Group));
x = T.Preference(idx)*100;
xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
scatter(ones(1,length(x))*2+rand(1,length(x))*.1-0.05,x,100,'filled','MarkerFaceColor',[252 216 213]/255,'MarkerEdgeColor','w')
plot([-0.125 0.125]+2,[xx(3) xx(3)],'color',[229 45 38]/255,'LineWidth',1)
plot([1 1]*2,[xx(2) xx(4)],'color',[229 45 38]/255,'LineWidth',1)
plot([1 1]*2,[xx(1) xx(5)],'color',[229 45 38]/255,'LineWidth',1)
xlim([0.5 2.5]); xticks(1:2); xticklabels({'YFP','eOPN3'});
ylabel('Flavor preference (%)'); ylim([0 50]); yticks(0:25:50);
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

[p,~,stats] = ranksum(T.Preference(find(cellfun(@(x) isequal(x,'eOPN3'),T.Group))),T.Preference(find(cellfun(@(x) isequal(x,'YFP'),T.Group))),'method','exact');
StatsTbl(end+1,:) = table({'2f'},{'eOPN3 vs. YFP'},{'Wilcoxon rank-sum'},{'N/A'},{[length(find(cellfun(@(x) isequal(x,'eOPN3'),T.Group))) length(find(cellfun(@(x) isequal(x,'YFP'),T.Group)))]},stats.ranksum,p(1));
%% Fig 2h

figure('Position', get(0, 'Screensize'))
sgtitle('Figure 2h','FontWeight','bold')
T1 = readtable([data_path,'\Fos imaging\region_info.csv']);
T2 = readtable([data_path,'\Fos imaging\sample_info.csv']);

axis square
hold on
counts_CEA = T1{find(cellfun(@(x) isequal(x,'CentralAmygdalarNucleus'),T1.region)),4:end};
counts_PB = T1{find(cellfun(@(x) isequal(x,'ParabrachialNucleus'),T1.region)),4:end};
counts_total = T1{1,4:end};

idx_cgrp = find(cellfun(@(x) isequal('CGRP stim',x),T2.Timepoint));
counts_norm = nan(size(counts_CEA));
counts_norm(idx_cgrp) = counts_CEA(idx_cgrp)./counts_PB(idx_cgrp)./mean(counts_total(idx_cgrp)./counts_PB(idx_cgrp));
counts_norm = [counts_norm./1.987']*100;
PHASEidx = find(cellfun(@(x) isequal('CGRP stim',x),T2.Timepoint)&cellfun(@(x) isequal('Familiar',x),T2.Novel));
x = counts_norm(PHASEidx);
xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
scatter(ones(1,length(x))+rand(1,length(x))*.1-0.05,x,100,'filled','MarkerFaceColor',[216 231 243]/255,'MarkerEdgeColor','w')
plot([-0.125 0.125]+1,[xx(3) xx(3)],'color',[55 136 193]/255,'LineWidth',1)
plot([1 1],[xx(2) xx(4)],'color',[55 136 193]/255,'LineWidth',1)
plot([1 1],[xx(1) xx(5)],'color',[55 136 193]/255,'LineWidth',1)
PHASEidx = find(cellfun(@(x) isequal('CGRP stim',x),T2.Timepoint)&cellfun(@(x) isequal('Novel',x),T2.Novel));
x = counts_norm(PHASEidx);
xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
scatter(ones(1,length(x))*2+rand(1,length(x))*.1-0.05,x,100,'filled','MarkerFaceColor',[252 216 213]/255,'MarkerEdgeColor','w')
plot([-0.125 0.125]+2,[xx(3) xx(3)],'color',[229 45 38]/255,'LineWidth',1)
plot([1 1]*2,[xx(2) xx(4)],'color',[229 45 38]/255,'LineWidth',1)
plot([1 1]*2,[xx(1) xx(5)],'color',[229 45 38]/255,'LineWidth',1)
xlim([0.5 2.5]); xticks(1:2); xticklabels({'Familiar','Novel'});
ylabel('CEA Fos (% per mm^3)'); ylim([0 2]); yticks(0:1:2);
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

load([data_path,'\Fos imaging\Fos-GLMM-statistics.mat'],'data')
N1 = length(find(cellfun(@(x) isequal('CGRP stim',x),T2.Timepoint)&cellfun(@(x) isequal('Novel',x),T2.Novel)));
N2 = length(find(cellfun(@(x) isequal('CGRP stim',x),T2.Timepoint)&cellfun(@(x) isequal('Familiar',x),T2.Novel)));
StatsTbl(end+1,:) = table({'2h'},{'Novel vs. Familiar'},{'GLMM marginal effect'},{'N/A'},{[N1 N2]},data.GLMMoutput.CEA.Eq4.flavor.Zstat(2),data.GLMMoutput.CEA.Eq4.flavor.pvalues_raw(2));
%% Fig 2i

figure('Position', get(0, 'Screensize'))
sgtitle('Figure 2i','FontWeight','bold')
load([data_path,'\Fos imaging\Fos-GLMM-statistics.mat'],'data')

T = readtable([data_path,'\source data\Fig-1e.csv']);
idx = find(data.GLMMoutput.Eq2.modelstats.significant);
amygdala_regions = idx(find(T.Cluster==1))';

idx_licl = find(cellfun(@(x) isequal('Malaise',x),data.GLMMinput.phase));
idx_cgrp = find(cellfun(@(x) isequal('CGRP stim',x),data.GLMMinput.phase));
idx_nov = find(cellfun(@(x) isequal('Novel',x),data.GLMMinput.flavor));
idx_fam = find(cellfun(@(x) isequal('Familiar',x),data.GLMMinput.flavor));
idx_amygdala = amygdala_regions;
idx_other = setdiff(find(data.GLMMoutput.Eq2.modelstats.significant),[find(data.regions.name=='Parabrachial nucleus');amygdala_regions]);

counts_norm = nan(size(data.GLMMinput.counts));
total = data.GLMMinput.offset./data.GLMMinput.pbn;
idx = sort([idx_nov,idx_fam]);
counts_norm(intersect(idx_licl,idx),:) = data.GLMMinput.counts(intersect(idx_licl,idx),:)./data.GLMMinput.pbn(intersect(idx_licl,idx))./mean(total(intersect(idx_licl,idx)));
counts_norm(intersect(idx_cgrp,idx),:) = data.GLMMinput.counts(intersect(idx_cgrp,idx),:)./data.GLMMinput.pbn(intersect(idx_cgrp,idx))./mean(total(intersect(idx_cgrp,idx)));
counts_norm = [counts_norm./data.regions.size']*100;

subplot(1,2,1)
hold on
axis square
A = mean(counts_norm(idx_licl,idx_amygdala));
B = mean(counts_norm(idx_cgrp,idx_amygdala));
a = scatter(A,B,64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([0 1.2])
ylim([0 ceil(max(B)*10)/10])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[y_fit] = polyval(p2,x,S);
fitresult = fit(A',B','poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[252 216 213]/255,'LineStyle','none'),
plot(x,y_fit,'color',[229 45 38]/255,'LineWidth',1)
scatter(A,B,100,[229 45 38]/255,'filled','MarkerEdgeColor','w')
xlim([0 1.2])
ylim([0 1.8])
xticks([0 1.2])
yticks([0 1.8])
xlabel(['Malaise',char(10),'Average Fos (% per mm^3)'])
ylabel(['CGRP stim',char(10),'Average Fos (% per mm^3)'])
title('Amygdala network')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

[r,p] = corr(A',B');
StatsTbl(end+1,:) = table({'2i, top'},{'Malaise vs. CGRP stim'},{'Pearson correlation'},{'N/A'},{[length(A)]},r,p);

subplot(1,2,2)
hold on
axis square
A = mean(counts_norm(idx_licl,idx_other));
B = mean(counts_norm(idx_cgrp,idx_other));
a = scatter(A,B,64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([0 ceil(max(A)*10)/10])
ylim([0 ceil(max(B)*10)/10])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[y_fit] = polyval(p2,x,S);
fitresult = fit(A',B','poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[.85 .85 .85],'LineStyle','none'),
plot(x,y_fit,'k','LineWidth',1)
scatter(A,B,100,'k','filled','MarkerEdgeColor','w')
xlim([0 1.6])
ylim([0 1.2])
xticks([0 1.6])
yticks([0 1.2])
xlabel(['Malaise',char(10),'Average Fos (% per mm^3)'])
ylabel(['CGRP stim',char(10),'Average Fos (% per mm^3)'])
title('Other regions')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

[r,p] = corr(A',B');
StatsTbl(end+1,:) = table({'2i, bottom'},{'Malaise vs. CGRP stim'},{'Pearson correlation'},{'N/A'},{[length(A)]},r,p);
%% Fig 2j

figure('Position', get(0, 'Screensize'))
sgtitle('Figure 2j','FontWeight','bold')
load([data_path,'\Fos imaging\Fos-GLMM-statistics.mat'],'data')

T = readtable([data_path,'\source data\Fig-1e.csv']);
idx = find(data.GLMMoutput.Eq2.modelstats.significant);
amygdala_regions = idx(find(T.Cluster==1))';

idx_licl = find(cellfun(@(x) isequal('Malaise',x),data.GLMMinput.phase));
idx_cgrp = find(cellfun(@(x) isequal('CGRP stim',x),data.GLMMinput.phase));
idx_nov = find(cellfun(@(x) isequal('Novel',x),data.GLMMinput.flavor));
idx_fam = find(cellfun(@(x) isequal('Familiar',x),data.GLMMinput.flavor));
idx_amygdala = amygdala_regions;
idx_other = setdiff(find(data.GLMMoutput.Eq2.modelstats.significant),[find(data.regions.name=='Parabrachial nucleus');amygdala_regions]);

counts_norm = nan(size(data.GLMMinput.counts));
total = data.GLMMinput.offset./data.GLMMinput.pbn;
idx = sort([idx_nov,idx_fam]);
counts_norm(intersect(idx_licl,idx),:) = data.GLMMinput.counts(intersect(idx_licl,idx),:)./data.GLMMinput.pbn(intersect(idx_licl,idx))./mean(total(intersect(idx_licl,idx)));
counts_norm(intersect(idx_cgrp,idx),:) = data.GLMMinput.counts(intersect(idx_cgrp,idx),:)./data.GLMMinput.pbn(intersect(idx_cgrp,idx))./mean(total(intersect(idx_cgrp,idx)));
counts_norm = [counts_norm./data.regions.size']*100;

subplot(1,2,1)
hold on
axis square
A = mean(counts_norm(intersect(idx_licl,idx_nov),idx_amygdala))-mean(counts_norm(intersect(idx_licl,idx_fam),idx_amygdala));
B = mean(counts_norm(intersect(idx_cgrp,idx_nov),idx_amygdala))-mean(counts_norm(intersect(idx_cgrp,idx_fam),idx_amygdala));
a = scatter(A,B,64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([floor(min(A)*10)/10 ceil(max(A)*10)/10])
ylim([-0.1 ceil(max(B)*10)/10])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[y_fit] = polyval(p2,x,S);
fitresult = fit(A',B','poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[252 216 213]/255,'LineStyle','none'),
plot(x,y_fit,'color',[229 45 38]/255,'LineWidth',1)
scatter(A,B,100,[229 45 38]/255,'filled','MarkerEdgeColor','w')
xlim([-.1 .3])
ylim([-.1 .9])
xticks([-.1 .3])
yticks([-.1 .9])
xlabel(['Malaise',char(10),'Novel – Familiar ΔFos (% per mm^3)'])
ylabel(['CGRP stim',char(10),'Novel – Familiar ΔFos (% per mm^3)'])
title('Amygdala network')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

[r,p] = corr(A',B');
StatsTbl(end+1,:) = table({'2j, top'},{'Malaise vs. CGRP stim'},{'Pearson correlation'},{'N/A'},{[length(A)]},r,p);

subplot(1,2,2)
hold on
axis square
A = nanmean(counts_norm(intersect(idx_licl,idx_nov),idx_other))-nanmean(counts_norm(intersect(idx_licl,idx_fam),idx_other));
B = nanmean(counts_norm(intersect(idx_cgrp,idx_nov),idx_other))-nanmean(counts_norm(intersect(idx_cgrp,idx_fam),idx_other));
a = scatter(A,B,64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([floor(min(A)*10)/10 ceil(max(A)*10)/10])
ylim([floor(min(B)*10)/10 0.5])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[y_fit] = polyval(p2,x,S);
fitresult = fit(A',B','poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[.85 .85 .85],'LineStyle','none'),
plot(x,y_fit,'k','LineWidth',1)
scatter(A,B,100,'k','filled','MarkerEdgeColor','w')
xlim([-.4 .6])
ylim([-.1 .5])
xticks([-.4 .6])
yticks([-.1 .5])
xlabel(['Malaise',char(10),'Novel – Familiar ΔFos (% per mm^3)'])
ylabel(['CGRP stim',char(10),'Novel – Familiar ΔFos (% per mm^3)'])
title('Other regions')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

[r,p] = corr(A',B');
StatsTbl(end+1,:) = table({'2j, bottom'},{'Malaise vs. CGRP stim'},{'Pearson correlation'},{'N/A'},{[length(A)]},r,p);
%% Fig 2k

fname = '\Fos imaging\kernel-density-estimates\kde-cgrp-stim-novel.npy';
KDE.CGRP.Novel = readNPY([data_path,fname]);
fname = '\Fos imaging\kernel-density-estimates\kde-cgrp-stim-familiar.npy';
KDE.CGRP.Familiar = readNPY([data_path,fname]);

load([data_path,'\Fos imaging\modified-atlas\allen_ccfv3_modified_cz.mat'],'atlas','RegionLibrary')

atlasmask  = atlas>=1028 | atlas<=1;
KDE.CGRP.Novel(atlasmask) = NaN;
KDE.CGRP.Familiar(atlasmask) = NaN;

cmap = flipud(cbrewer('div','RdBu',1000,'spline')); cmap(cmap<0) = 0;

figure('Position', get(0, 'Screensize'))
sgtitle('Figure 2k','FontWeight','bold')
load([data_path,'\Fos imaging\Fos-GLMM-statistics.mat'],'data'); data_in = data;
T = readtable([data_path,'\source data\Fig-1e.csv']);

idx = find(data.GLMMoutput.Eq2.modelstats.significant);
idx_amygdala = idx(find(T.Cluster==1));
planes_nov = [210,265,325];
for pl = 1:2
    plane = planes_nov(pl+1);
    subplot(2,1,pl)
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
    colormap(gcf,cmap)
    heatmap(rot90(data),[],[],[],'UseFigureColormap',true,'ColorLevels',1000,'MaxColorValue',100,'MinColorValue',-100,'NaNColor',[1 1 1]);
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
    set(gca,'FontSize',12)
    xlim([.5 228.5])
    axis off
    hold off
end
%% Fig 2m

figure('Position', get(0, 'Screensize'))
sgtitle('Figure 2m','FontWeight','bold')
T = readtable([data_path,'\source data\Fig-2m.csv']);

hold on
axis square
for j = 1:3
    idx = cellfun(@(x) isequal(x,'Novel'),T.Flavor);
    x = T{idx,j+1}*100;
    xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
    scatter(ones(1,length(x))*j+rand(1,length(x))*.1-.05+.125,x,100,'filled','MarkerFaceColor',[252 216 213]/255,'MarkerEdgeColor','w')
    a = plot([-0.125 0.125]+j+.125,[xx(3) xx(3)],'color',[229 45 38]/255,'LineWidth',1);
    plot([j j]+.125,[xx(2) xx(4)],'color',[229 45 38]/255,'LineWidth',1)
    plot([j j]+.125,[xx(1) xx(5)],'color',[229 45 38]/255,'LineWidth',1)
    idx = cellfun(@(x) isequal(x,'Familiar'),T.Flavor);
    x = T{idx,j+1}*100;
    xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
    scatter(ones(1,length(x))*j+rand(1,length(x))*.1-.05-.125,x,100,'filled','MarkerFaceColor',[216 231 243]/255,'MarkerEdgeColor','w')
    b = plot([-0.125 0.125]+j-.125,[xx(3) xx(3)],'color',[55 136 193]/255,'LineWidth',1);
    plot([j j]-.125,[xx(2) xx(4)],'color',[55 136 193]/255,'LineWidth',1)
    plot([j j]-.125,[xx(1) xx(5)],'color',[55 136 193]/255,'LineWidth',1)
end
xlim([0.5 3.5]); xticks(1:3); xticklabels({'{\it Sst}','{\it Prkcd}','{\it Calcrl}'});
ylabel('{\it Fos}^+ CEA cells (%)'); ylim([0 100]); yticks(0:50:100);
legend([b,a],{'Familiar','Novel'})
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

idxNov = find(cellfun(@(x) isequal(x,'Novel'),T.Flavor));
idxFam = find(cellfun(@(x) isequal(x,'Familiar'),T.Flavor));
[p,~,stats] = ranksum(T.Sst_(idxNov),T.Sst_(idxFam),'method','exact');
StatsTbl(end+1,:) = table({'2m'},{'Sst: Novel vs. Familiar'},{'Wilcoxon rank-sum'},{'N/A'},{[length(idxNov) length(idxFam)]},stats.ranksum,p);
[p,~,stats] = ranksum(T.Prkcd_(idxNov),T.Prkcd_(idxFam),'method','exact');
StatsTbl(end+1,:) = table({'2m'},{'Prkcd: Novel vs. Familiar'},{'Wilcoxon rank-sum'},{'N/A'},{[length(idxNov) length(idxFam)]},stats.ranksum,p);
[p,~,stats] = ranksum(T.Calcrl_(idxNov),T.Calcrl_(idxFam),'method','exact');
StatsTbl(end+1,:) = table({'2m'},{'Calcrl: Novel vs. Familiar'},{'Wilcoxon rank-sum'},{'N/A'},{[length(idxNov) length(idxFam)]},stats.ranksum,p);
%% Fig 2n

figure('Position', get(0, 'Screensize'))
sgtitle('Figure 2n','FontWeight','bold')
T = readtable([data_path,'\source data\Fig-2n.csv']);

hold on
axis square

for j = 1:8
    idx = cellfun(@(x) isequal(x,'Novel'),T.Flavor);
    x = T{idx,j+1}*100;
    xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
    scatter(ones(1,length(x))*j+rand(1,length(x))*.1-.05+.125,x,100,'filled','MarkerFaceColor',[252 216 213]/255,'MarkerEdgeColor','w')
    a = plot([-0.125 0.125]+j+.125,[xx(3) xx(3)],'color',[229 45 38]/255,'LineWidth',1);
    plot([j j]+.125,[xx(2) xx(4)],'color',[229 45 38]/255,'LineWidth',1)
    plot([j j]+.125,[xx(1) xx(5)],'color',[229 45 38]/255,'LineWidth',1)
    idx = cellfun(@(x) isequal(x,'Familiar'),T.Flavor);
    x = T{idx,j+1}*100;
    xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
    scatter(ones(1,length(x))*j+rand(1,length(x))*.1-.05-.125,x,100,'filled','MarkerFaceColor',[216 231 243]/255,'MarkerEdgeColor','w')
    b = plot([-0.125 0.125]+j-.125,[xx(3) xx(3)],'color',[55 136 193]/255,'LineWidth',1);
    plot([j j]-.125,[xx(2) xx(4)],'color',[55 136 193]/255,'LineWidth',1)
    plot([j j]-.125,[xx(1) xx(5)],'color',[55 136 193]/255,'LineWidth',1)
end
xlim([0.5 8.5]); xticks(1:8); xticklabels({'none','{\it Sst}','{\it Prkcd}','{\it Calcrl}','{\it Sst/Prkcd}','{\it Sst/Calcrl}','{\it Prkcd/Calcrl}','all'});
ylabel('{\it Fos}^+ CEA cells (%)'); ylim([0 50]); yticks(0:25:50);
legend([b,a],{'Familiar','Novel'})
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

idxNov = find(cellfun(@(x) isequal(x,'Novel'),T.Flavor));
idxFam = find(cellfun(@(x) isequal(x,'Familiar'),T.Flavor));
[p,~,stats] = ranksum(T.None(idxNov),T.None(idxFam),'method','exact');
StatsTbl(end+1,:) = table({'2n'},{'None: Novel vs. Familiar'},{'Wilcoxon rank-sum'},{'N/A'},{[length(idxNov) length(idxFam)]},stats.ranksum,p);
[p,~,stats] = ranksum(T.Sst_only(idxNov),T.Sst_only(idxFam),'method','exact');
StatsTbl(end+1,:) = table({'2n'},{'Sst: Novel vs. Familiar'},{'Wilcoxon rank-sum'},{'N/A'},{[length(idxNov) length(idxFam)]},stats.ranksum,p);
[p,~,stats] = ranksum(T.Prkcd_only(idxNov),T.Prkcd_only(idxFam),'method','exact');
StatsTbl(end+1,:) = table({'2n'},{'Prkcd: Novel vs. Familiar'},{'Wilcoxon rank-sum'},{'N/A'},{[length(idxNov) length(idxFam)]},stats.ranksum,p);
[p,~,stats] = ranksum(T.Calcrl_only(idxNov),T.Calcrl_only(idxFam),'method','exact');
StatsTbl(end+1,:) = table({'2n'},{'Calcrl: Novel vs. Familiar'},{'Wilcoxon rank-sum'},{'N/A'},{[length(idxNov) length(idxFam)]},stats.ranksum,p);
[p,~,stats] = ranksum(T.Sst_Prkcd_only(idxNov),T.Sst_Prkcd_only(idxFam),'method','exact');
StatsTbl(end+1,:) = table({'2n'},{'Sst/Prkcd: Novel vs. Familiar'},{'Wilcoxon rank-sum'},{'N/A'},{[length(idxNov) length(idxFam)]},stats.ranksum,p);
[p,~,stats] = ranksum(T.Sst_Calcrl_only(idxNov),T.Sst_Calcrl_only(idxFam),'method','exact');
StatsTbl(end+1,:) = table({'2n'},{'Sst/Calcrl: Novel vs. Familiar'},{'Wilcoxon rank-sum'},{'N/A'},{[length(idxNov) length(idxFam)]},stats.ranksum,p);
[p,~,stats] = ranksum(T.Prkcd_Calcrl_only(idxNov),T.Prkcd_Calcrl_only(idxFam),'method','exact');
StatsTbl(end+1,:) = table({'2n'},{'Prkcd/Calcrl: Novel vs. Familiar'},{'Wilcoxon rank-sum'},{'N/A'},{[length(idxNov) length(idxFam)]},stats.ranksum,p);
[p,~,stats] = ranksum(T.all(idxNov),T.all(idxFam),'method','exact');
StatsTbl(end+1,:) = table({'2n'},{'All: Novel vs. Familiar'},{'Wilcoxon rank-sum'},{'N/A'},{[length(idxNov) length(idxFam)]},stats.ranksum,p);

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