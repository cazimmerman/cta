function fig_ED3(data_path)

disp('Generating panels for Extended Data Figure 3...')
%% Fig ED3a

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 3a','FontWeight','bold')
load([data_path,'/FOS imaging/FOS-GLMM-statistics.mat'],'data')

hold on
axis square
x = struct;
idx = intersect(find(data.GLMMoutput.Eq2.modelstats.significant),find(data.regions.parent=='Cerebral Cortex'));
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
StatsTbl = table({'ED 3a'},{'Consume vs. Malaise'},{'Kolmogorov-Smirnov'},{'3 pairs of time points'},{[length(x.Consume)]},stat(1),p(1), ...
    'VariableNames',{'Figure panel','Group','Statistical test','Multiple comparisons','Sample size','Test statistic','P-value'});
StatsTbl(end+1,:) = table({'ED 3a'},{'Consume vs. Retrieval'},{'Kolmogorov-Smirnov'},{'3 pairs of time points'},{[length(x.Consume)]},stat(2),p(2));
StatsTbl(end+1,:) = table({'ED 3a'},{'Malaise vs. Retrieval'},{'Kolmogorov-Smirnov'},{'3 pairs of time points'},{[length(x.Consume)]},stat(3),p(3));
%% Fig ED3b

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 3b','FontWeight','bold')
load([data_path,'/FOS imaging/FOS-GLMM-statistics.mat'],'data')

hold on
axis square
x = struct;
idx = intersect(find(data.GLMMoutput.Eq2.modelstats.significant),find([data.regions.parent=='Cerebral Nuclei'|data.regions.parent=='Thalamus'|data.regions.parent=='Hypothalamus']));
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
StatsTbl(end+1,:) = table({'ED 3b'},{'Consume vs. Malaise'},{'Kolmogorov-Smirnov'},{'3 pairs of time points'},{[length(x.Consume)]},stat(1),p(1));
StatsTbl(end+1,:) = table({'ED 3b'},{'Consume vs. Retrieval'},{'Kolmogorov-Smirnov'},{'3 pairs of time points'},{[length(x.Consume)]},stat(2),p(2));
StatsTbl(end+1,:) = table({'ED 3b'},{'Malaise vs. Retrieval'},{'Kolmogorov-Smirnov'},{'3 pairs of time points'},{[length(x.Consume)]},stat(3),p(3));
%% Fig ED3c

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 3c','FontWeight','bold')
load([data_path,'/FOS imaging/FOS-GLMM-statistics.mat'],'data')

hold on
axis square
x = struct;
idx = intersect(find(data.GLMMoutput.Eq2.modelstats.significant),find([data.regions.parent=='Pons'|data.regions.parent=='Medulla'|data.regions.parent=='Midbrain']));
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
StatsTbl(end+1,:) = table({'ED 3c'},{'Consume vs. Malaise'},{'Kolmogorov-Smirnov'},{'3 pairs of time points'},{[length(x.Consume)]},stat(1),p(1));
StatsTbl(end+1,:) = table({'ED 3c'},{'Consume vs. Retrieval'},{'Kolmogorov-Smirnov'},{'3 pairs of time points'},{[length(x.Consume)]},stat(2),p(2));
StatsTbl(end+1,:) = table({'ED 3c'},{'Malaise vs. Retrieval'},{'Kolmogorov-Smirnov'},{'3 pairs of time points'},{[length(x.Consume)]},stat(3),p(3));
%% Fig ED3e-n

load([data_path,'/FOS imaging/modified-atlas/allen_ccfv3_modified_cz.mat'],'atlas')
load([data_path,'/FOS imaging/FOS-GLMM-statistics.mat'],'data')
T = readtable([data_path,'/source data/Fig-1e.csv']);
significant = find(data.GLMMoutput.Eq2.modelstats.significant);
panels = {'e','f','g','h','i','j','k','l','m','n'};
for i = 1:max(T.Cluster)
    figure('Position', get(0, 'Screensize'))
    sgtitle(['Extended Data Figure 3',panels{i}],'FontWeight','bold')
    subplot(1,2,1)
    hold on
    regions = significant(T.Cluster==i);
    for j = 1:length(regions)
        t = atlas==data.regions.index(regions(j))+1;
        p = patch(isosurface(t));
        if i == 1
            p.FaceColor = [229 45 38]/255;
        else
            p.FaceColor = [44 31 22]/255;
        end
        p.FaceAlpha = 0.4;
        p.LineStyle = 'none';
    end
    p = patch(isosurface(atlas>0));
    p.FaceAlpha = 0.05;
    p.FaceColor = [44 31 22]/255;
    p.LineStyle = 'none';
    h = get(gca,'DataAspectRatio');
    title(['Cluster ',num2str(i),char(10),num2str(length(regions)),' brain regions'])
    set(gca,'DataAspectRatio',[1 1 h(3)],'CameraViewAngleMode','Manual');
    set(gca,'Ydir','reverse','Zdir','reverse')
    axis off; xlim('auto'); ylim('auto'); zlim('auto'); view(-20,20)
    set(gca,'FontSize',12)
    hold off
    
    subplot(1,2,2)
    axis square
    hold on
    if i == 1
        c1 = [229 45 38]/255;
        c0 = [252 216 213]/255;
    else
        c1 = [0 0 0];
        c0 = [.85 .85 .85];
    end
    x = data.GLMMoutput.Eq2.flavor.Zstat(regions,1);
    xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
    scatter(ones(1,length(x))+rand(1,length(x))*.1-.05,x,100,'filled','MarkerFaceColor',c0,'MarkerEdgeColor','w')
    plot([-0.125 0.125]+1,[xx(3) xx(3)],'color',c1,'LineWidth',1)
    plot([1 1],[xx(2) xx(4)],'color',c1,'LineWidth',1)
    x = data.GLMMoutput.Eq2.flavor.Zstat(regions,2);
    xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
    scatter(ones(1,length(x))*2+rand(1,length(x))*.1-.05,x,100,'filled','MarkerFaceColor',c0,'MarkerEdgeColor','w')
    plot([-0.125 0.125]+2,[xx(3) xx(3)],'color',c1,'LineWidth',1)
    plot([2 2],[xx(2) xx(4)],'color',c1,'LineWidth',1)
    x = data.GLMMoutput.Eq2.flavor.Zstat(regions,3);
    xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
    scatter(ones(1,length(x))*3+rand(1,length(x))*.1-.05,x,100,'filled','MarkerFaceColor',c0,'MarkerEdgeColor','w')
    plot([-0.125 0.125]+3,[xx(3) xx(3)],'color',c1,'LineWidth',1)
    plot([3 3],[xx(2) xx(4)],'color',c1,'LineWidth',1)
    xlim([0.5 3.5]); xticks(1:3); xticklabels({'Consume','Malaise','Retrieval'});
    ylabel('Novel – Familiar ΔFOS (Z)'); ylim([-5 5]); yticks(-5:2.5:5);
    set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
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