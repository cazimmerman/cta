function fig_ED2(data_path)

disp('Generating panels for Extended Data Figure 2...')
%% Fig ED2b

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 2b','FontWeight','bold')
T = readtable([data_path,'/source data/Fig-ED2b.csv']);

hold on
axis square
idx = find(cellfun(@(x) isequal(x,'YFP'),T.Group));
x = T.Preference(idx)*100;
xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
scatter(ones(1,length(x))+rand(1,length(x))*.1-0.05,x,100,'filled','MarkerFaceColor',[.85 .85 .85],'MarkerEdgeColor','w')
plot([-0.125 0.125]+1,[xx(3) xx(3)],'color','k','LineWidth',1)
plot([1 1],[xx(2) xx(4)],'color','k','LineWidth',1)
plot([1 1],[xx(1) xx(5)],'color','k','LineWidth',1)
idx = find(cellfun(@(x) isequal(x,'hM3D'),T.Group));
x = T.Preference(idx)*100;
xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
scatter(ones(1,length(x))*2+rand(1,length(x))*.1-0.05,x,100,'filled','MarkerFaceColor',[252 216 213]/255,'MarkerEdgeColor','w')
plot([-0.125 0.125]+2,[xx(3) xx(3)],'color',[229 45 38]/255,'LineWidth',1)
plot([1 1]*2,[xx(2) xx(4)],'color',[229 45 38]/255,'LineWidth',1)
plot([1 1]*2,[xx(1) xx(5)],'color',[229 45 38]/255,'LineWidth',1)
xlim([0.5 2.5]); xticks(1:2); xticklabels({'YFP','hM3D'});
ylabel('Flavor preference (%)'); ylim([0 100]); yticks(0:50:100);
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

[p,~,stats] = ranksum(T.Preference(find(cellfun(@(x) isequal(x,'hM3D'),T.Group))),T.Preference(find(cellfun(@(x) isequal(x,'YFP'),T.Group))),'method','exact');
StatsTbl = table({'ED 2b'},{'hM3D vs. YFP'},{'Wilcoxon rank-sum'},{'N/A'},{[length(find(cellfun(@(x) isequal(x,'hM3D'),T.Group))) length(find(cellfun(@(x) isequal(x,'YFP'),T.Group)))]},stats.ranksum,p(1), ...
    'VariableNames',{'Figure panel','Group','Statistical test','Multiple comparisons','Sample size','Test statistic','P-value'});
%% Fig ED2d

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 2d','FontWeight','bold')
T1 = readtable([data_path,'/Fos imaging/region_info.csv']);
T2 = readtable([data_path,'/Fos imaging/sample_info.csv']);

axis square
hold on
counts_norm = [T1{find(cellfun(@(x) isequal(x,'LateralSeptalComplex'),T1.region)),4:end}./T1{1,4:end}./3.553]*100;
idx = find(cellfun(@(x) isequal('Malaise + LS-YFP',x),T2.Timepoint));
x = counts_norm(idx);
xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
scatter(ones(1,length(x))+rand(1,length(x))*.1-0.05,x,100,'filled','MarkerFaceColor',[.85 .85 .85],'MarkerEdgeColor','w')
plot([-0.125 0.125]+1,[xx(3) xx(3)],'color','k','LineWidth',1)
plot([1 1],[xx(2) xx(4)],'color','k','LineWidth',1)
plot([1 1],[xx(1) xx(5)],'color','k','LineWidth',1)
idx = find(cellfun(@(x) isequal('Malaise + LS-hM3D',x),T2.Timepoint));
x = counts_norm(idx);
xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
scatter(ones(1,length(x))*2+rand(1,length(x))*.1-0.05,x,100,'filled','MarkerFaceColor',[252 216 213]/255,'MarkerEdgeColor','w')
plot([-0.125 0.125]+2,[xx(3) xx(3)],'color',[229 45 38]/255,'LineWidth',1)
plot([1 1]*2,[xx(2) xx(4)],'color',[229 45 38]/255,'LineWidth',1)
plot([1 1]*2,[xx(1) xx(5)],'color',[229 45 38]/255,'LineWidth',1)
xlim([0.5 2.5]); xticks(1:2); xticklabels({'YFP','hM3D'});
ylabel('LS Fos (% per mm^3)'); ylim([0 5]); yticks(0:2.5:5);
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

load([data_path,'/Fos imaging/Fos-GLMM-statistics.mat'],'data')
StatsTbl(end+1,:) = table({'ED 2d'},{'hM3D vs. YFP'},{'GLMM coefficient estimate'},{'N/A'},{[12 12]},data.GLMMoutput.LS.Eq5.coefficients.Zstat(2),data.GLMMoutput.LS.Eq5.coefficients.pvalues_raw(2));
%% Fig ED2e

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 2e','FontWeight','bold')
T1 = readtable([data_path,'/Fos imaging/region_info.csv']);
T2 = readtable([data_path,'/Fos imaging/sample_info.csv']);

axis square
hold on
counts_norm = [T1{find(cellfun(@(x) isequal(x,'CentralAmygdalarNucleus'),T1.region)),4:end}./T1{1,4:end}./1.987]*100;
idx = find(cellfun(@(x) isequal('Malaise + LS-YFP',x),T2.Timepoint));
x = counts_norm(idx);
xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
scatter(ones(1,length(x))+rand(1,length(x))*.1-0.05,x,100,'filled','MarkerFaceColor',[.85 .85 .85],'MarkerEdgeColor','w')
plot([-0.125 0.125]+1,[xx(3) xx(3)],'color','k','LineWidth',1)
plot([1 1],[xx(2) xx(4)],'color','k','LineWidth',1)
plot([1 1],[xx(1) xx(5)],'color','k','LineWidth',1)
idx = find(cellfun(@(x) isequal('Malaise + LS-hM3D',x),T2.Timepoint));
x = counts_norm(idx);
xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
scatter(ones(1,length(x))*2+rand(1,length(x))*.1-0.05,x,100,'filled','MarkerFaceColor',[252 216 213]/255,'MarkerEdgeColor','w')
plot([-0.125 0.125]+2,[xx(3) xx(3)],'color',[229 45 38]/255,'LineWidth',1)
plot([1 1]*2,[xx(2) xx(4)],'color',[229 45 38]/255,'LineWidth',1)
plot([1 1]*2,[xx(1) xx(5)],'color',[229 45 38]/255,'LineWidth',1)
xlim([0.5 2.5]); xticks(1:2); xticklabels({'YFP','hM3D'});
ylabel('CEA Fos (% per mm^3)'); ylim([0 .5]); yticks(0:.25:.5);
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

load([data_path,'/Fos imaging/Fos-GLMM-statistics.mat'],'data')
StatsTbl(end+1,:) = table({'ED 2e'},{'hM3D vs. YFP'},{'GLMM coefficient estimate'},{'N/A'},{[12 12]},data.GLMMoutput.CEA.Eq5.coefficients.Zstat(2),data.GLMMoutput.CEA.Eq5.coefficients.pvalues_raw(2));
%% Fig ED2f

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 2f','FontWeight','bold')
load([data_path,'/Fos imaging/Fos-GLMM-statistics.mat'],'data')

hold on
axis square
counts_norm = [data.GLMMinput.counts./data.GLMMinput.offset./data.regions.size']*100;

T = readtable([data_path,'/source data/Fig-1e.csv']);
regions = struct;
regions.significant = find(data.GLMMoutput.Eq2.modelstats.significant);
regions.septum = find(contains(data.regions.name,'ept'));
regions.amygdala = regions.significant(find(T.Cluster==1))';
YFP = find(cellfun(@(x) isequal('Malaise + LS-YFP',x),data.GLMMinput.phase));
hM3D = find(cellfun(@(x) isequal('Malaise + LS-hM3D',x),data.GLMMinput.phase));

idx0 = setdiff(regions.significant,[regions.amygdala;regions.septum]);
A = mean(counts_norm(hM3D,idx0));
B = mean(counts_norm(YFP,idx0));
a = scatter(A,B,1,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A,B,1);
x = 0:.01:4;
delete(a)
[y_fit] = polyval(p2,x,S);
fitresult = fit(A',B','poly1');
p2 = predint(fitresult,x,0.95,'functional');
p2(p2<0) = 0; p2(p2>4) = 4;
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[.85 .85 .85],'LineStyle','none'),
plot(x,y_fit,'Color','k','LineWidth',1)

A = mean(counts_norm(hM3D,regions.amygdala));
B = mean(counts_norm(YFP,regions.amygdala));
a = scatter(A,B,1,'r','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A,B,1);
x = 0:.01:4;
delete(a)
[y_fit] = polyval(p2,x,S);
fitresult = fit(A',B','poly1');
p2 = predint(fitresult,x,0.95,'functional');
p2(p2<0) = 0; p2(p2>4) = 4;
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[252 216 213]/255,'LineStyle','none'),
plot(x,y_fit,'Color',[229 45 38]/255,'LineWidth',1)

A = mean(counts_norm(hM3D,regions.septum));
B = mean(counts_norm(YFP,regions.septum));
a = scatter(A,B,1,'r','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A,B,1);
x = 0:.01:4;
delete(a)
[y_fit] = polyval(p2,x,S);
fitresult = fit(A',B','poly1');
p2 = predint(fitresult,x,0.95,'functional');
p2(p2<0) = 0; p2(p2>4) = 4;
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[216 231 243]/255,'LineStyle','none'),
plot(x,y_fit,'Color',[55 136 193]/255,'LineWidth',1)

a = scatter(mean(counts_norm(hM3D,idx0)),mean(counts_norm(YFP,idx0)),100,'k','filled','MarkerEdgeColor','w');
b = scatter(mean(counts_norm(hM3D,regions.amygdala)),mean(counts_norm(YFP,regions.amygdala)),100,[229 45 38]/255,'filled','MarkerEdgeColor','w');
c = scatter(mean(counts_norm(hM3D,regions.septum)),mean(counts_norm(YFP,regions.septum)),100,[55 136 193]/255,'filled','MarkerEdgeColor','w');

xlim([0 4]); xticks(0:1:4);
ylim([0 4]); yticks(0:1:4);
ylabel('YFP Fos^ (% per mm^3)')
xlabel('hM3D Fos^ (% per mm^3)')
legend([b,c,a],{'Amygdala network','Septal complex','Other regions'})
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

G = cell(0,0);
G(1:length(regions.significant)) = {'Other regions'};
[~,idx1]=ismember(regions.amygdala,regions.significant); G(idx1) = {'Amygdala network'};
[~,idx2]=ismember(regions.septum,regions.significant); G(idx2) = {'Septal complex'};
figure
[~,~,~,stats] = aoctool(mean(counts_norm(hM3D,regions.significant)),mean(counts_norm(YFP,regions.significant)),G,.05,'hM3D','YFP','display','off');
% calculate T values using the function called by multcomp
gmeans = stats.slopes; gcov = stats.slopecov;
t = isnan(gmeans);
if any(t)
    gcov(t,:) = 0;
    gcov(:,t) = 0;
end
ng = length(gmeans);
M = nchoosek(1:ng, 2);
g1 = M(:,1);
g2 = M(:,2);
mn = gmeans(g1) - gmeans(g2);
i12 = sub2ind(size(gcov), g1, g2);
gvar = diag(gcov);
se = sqrt(gvar(g1) + gvar(g2) - 2 * gcov(i12));
t = mn./se; t = [-t(1) -t(2) t(3)]; % correct directions of comparisons
[c,~,~,~] = multcompare(stats,0.05,'on','bonferroni','slope','display','off');
close(gcf)
StatsTbl(end+1,:) = table({'ED 2f'},{'Amygdala network vs. Septal complex'},{'One-way ANCOVA (slope)'},{'3 pairs of region groups'},{[length(idx1) length(idx2)]},t(3),c(3,end));
StatsTbl(end+1,:) = table({'ED 2f'},{'Amygdala network vs. Other regions'},{'One-way ANCOVA (slope)'},{'3 pairs of region groups'},{[length(idx1) length(setdiff(regions.significant,[regions.amygdala;regions.septum]))]},t(1),c(1,end));
StatsTbl(end+1,:) = table({'ED 2f'},{'Septal complex vs. Other regions'},{'One-way ANCOVA (slope)'},{'3 pairs of region groups'},{[length(idx2) length(setdiff(regions.significant,[regions.amygdala;regions.septum]))]},t(2),c(2,end));
%% Fig ED2g

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 2g','FontWeight','bold')

fname = '/Fos imaging/kernel-density-estimates/kde-malaise-LS-YFP.npy';
KDE.YFP = readNPY([data_path,fname]);
fname = '/Fos imaging/kernel-density-estimates/kde-malaise-LS-hM3D.npy';
KDE.hM3D = readNPY([data_path,fname]);

load([data_path,'/Fos imaging/modified-atlas/allen_ccfv3_modified_cz.mat'],'atlas')

atlasmask  = atlas>=1028 | atlas<=1;
KDE.YFP(atlasmask) = NaN;
KDE.hM3D(atlasmask) = NaN;

cmap = flipud(cbrewer('div','RdGy',1000,'spline')); cmap(cmap<0) = 0; cmap(cmap>1) = 1;

load([data_path,'/Fos imaging/Fos-GLMM-statistics.mat'],'data'); data_in = data;
T = readtable([data_path,'/source data/Fig-1e.csv']);
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
    
    data = squeeze(KDE.hM3D(:,plane,:))-squeeze(KDE.YFP(:,plane,:));
    colormap(gcf,cmap)
    heatmap(rot90(data),[],[],[],'Colormap',cmap,'ColorLevels',1000,'MaxColorValue',.5,'MinColorValue',-.5,'NaNColor',[1 1 1]);
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