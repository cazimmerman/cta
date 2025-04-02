function fig_ED4(data_path)

disp('Generating panels for Extended Data Figure 4...')
%% Fig ED4a

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 4a','FontWeight','bold')
load([data_path,'/FOS imaging/FOS-GLMM-statistics.mat'],'data');

counts_norm = [data.GLMMinput.counts./data.GLMMinput.totalcounts./data.regions.size']*100;
significant = find(data.GLMMoutput.Eq2.modelstats.significant);
treedata = data.GLMMoutput.Eq2.flavor.Zstat(significant,:);
tree = linkage(treedata,'ward','chebychev');
D = pdist(treedata);
leafOrder = optimalleaforder(tree,D);
[~,~,outperm] = dendrogram(tree,length(significant)+1,'Labels',[],'ColorThreshold',4.7,'Reorder',leafOrder,'Orientation','left');
delete(gca)
cmap = flipud(cbrewer('div','RdBu',1000,'spline')); cmap(cmap<0) = 0;

A = counts_norm(find(cellfun(@(x) isequal('Consumption',x),data.GLMMinput.timepoint)),significant(leafOrder));
subplot(1,3,1);
axis square
hold on
heatmap(corr(A),[],[],[],'Colormap',cmap,'ColorLevels',1000,'MaxColorValue',1,'MinColorValue',-1,'NaNColor',[1 1 1]);
idx = cluster(tree,'maxclust',10);
idx = idx(outperm);
list = unique(idx,'stable');
counter = [];
for i = 1:length(list)-1
    counter(i) = sum(idx==list(i));
    plot([0 size(A,2)]+0.5,[sum(counter) sum(counter)]+0.5,'k','LineWidth',1)
    plot([sum(counter) sum(counter)]+0.5,[0 size(A,2)]+0.5,'k','LineWidth',1)
end
title('Consumption')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0, 0],'TickDir','out')
set(gca,'Xdir','reverse')
axis off
hold off

A = counts_norm(find(cellfun(@(x) isequal('Malaise',x),data.GLMMinput.timepoint)),significant(leafOrder));
subplot(1,3,2);
hold on
axis square
heatmap(corr(A),[],[],[],'Colormap',cmap,'ColorLevels',1000,'MaxColorValue',1,'MinColorValue',-1,'NaNColor',[1 1 1]);
idx = cluster(tree,'maxclust',10);
idx = idx(outperm);
list = unique(idx,'stable');
counter = [];
for i = 1:length(list)-1
    counter(i) = sum(idx==list(i));
    plot([0 size(A,2)]+0.5,[sum(counter) sum(counter)]+0.5,'k','LineWidth',1)
    plot([sum(counter) sum(counter)]+0.5,[0 size(A,2)]+0.5,'k','LineWidth',1)
end
title('Malaise')
set(gca,'Xdir','reverse')
axis off
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0, 0],'TickDir','out')
hold off

A = counts_norm(find(cellfun(@(x) isequal('Retrieval',x),data.GLMMinput.timepoint)),significant(leafOrder));
subplot(1,3,3);
axis square
hold on
heatmap(corr(A),[],[],[],'Colormap',cmap,'ColorLevels',1000,'MaxColorValue',1,'MinColorValue',-1,'NaNColor',[1 1 1]);
idx = cluster(tree,'maxclust',10);
idx = idx(outperm);
list = unique(idx,'stable');
counter = [];
for i = 1:length(list)-1
    counter(i) = sum(idx==list(i));
    plot([0 size(A,2)]+0.5,[sum(counter) sum(counter)]+0.5,'k','LineWidth',1)
    plot([sum(counter) sum(counter)]+0.5,[0 size(A,2)]+0.5,'k','LineWidth',1)
end
title('Retrieval')
set(gca,'Xdir','reverse')
axis off
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0, 0],'TickDir','out')
hold off
%% Fig ED4b

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 4b','FontWeight','bold')
load([data_path,'/FOS imaging/FOS-GLMM-statistics.mat'],'data');

counts_norm = [data.GLMMinput.counts./data.GLMMinput.totalcounts./data.regions.size']*100;
significant = find(data.GLMMoutput.Eq2.modelstats.significant);
treedata = data.GLMMoutput.Eq2.flavor.Zstat(significant,:);
tree = linkage(treedata,'ward','chebychev');
D = pdist(treedata);
leafOrder = optimalleaforder(tree,D);
[~,~,outperm] = dendrogram(tree,length(significant)+1,'Labels',[],'ColorThreshold',4.7,'Reorder',leafOrder,'Orientation','left');
idx = cluster(tree,'maxclust',10);
idx = idx(outperm);
A = counts_norm(:,significant(leafOrder));
CorrStruct = nan(12,3);
for i = 1
    in = find(idx==10);
    for j = 1:length(in)
        for k = 1
            t1 = []; t2 = []; t3 = [];
            test = setdiff(find(idx==10),in(j));
            for l = 1:length(test)
                t1(1) = corr(A(find(cellfun(@(x) isequal('Consumption',x),data.GLMMinput.timepoint)),in(j)),A(find(cellfun(@(x) isequal('Consumption',x),data.GLMMinput.timepoint)),test(l)));
                t2(1) = corr(A(find(cellfun(@(x) isequal('Malaise',x),data.GLMMinput.timepoint)),in(j)),A(find(cellfun(@(x) isequal('Malaise',x),data.GLMMinput.timepoint)),test(l)));
                t3(1) = corr(A(find(cellfun(@(x) isequal('Retrieval',x),data.GLMMinput.timepoint)),in(j)),A(find(cellfun(@(x) isequal('Retrieval',x),data.GLMMinput.timepoint)),test(l)));
            end
            CorrStruct(1,j) = mean(t1);
            CorrStruct(2,j) = mean(t2);
            CorrStruct(3,j) = mean(t3);
        end
    end
end
delete(gca)
hold on
axis square
x = CorrStruct(1,:);
xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
scatter(ones(1,length(x))+rand(1,length(x))*.1-.05,x,100,'filled','MarkerFaceColor',[252 216 213]/255,'MarkerEdgeColor','w')
plot([-0.125 0.125]+1,[xx(3) xx(3)],'Color',[229 45 38]/255,'LineWidth',1)
plot([1 1],[xx(2) xx(4)],'Color',[229 45 38]/255,'LineWidth',1)
p = []; stat = [];
[p(1),~,s] = signrank(x,0,'method','exact'); stat(1) = s.signedrank;
x = CorrStruct(2,:);
xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
scatter(ones(1,length(x))*2+rand(1,length(x))*.1-.05,x,100,'filled','MarkerFaceColor',[252 216 213]/255,'MarkerEdgeColor','w')
plot([-0.125 0.125]+2,[xx(3) xx(3)],'Color',[229 45 38]/255,'LineWidth',1)
plot([2 2],[xx(2) xx(4)],'Color',[229 45 38]/255,'LineWidth',1)
[p(2),~,s] = signrank(x,0,'method','exact'); stat(2) = s.signedrank;
x = CorrStruct(3,:);
xx = [NaN mean(x)-std(x)./sqrt(length(x)) mean(x) mean(x)+std(x)./sqrt(length(x)) NaN];
scatter(ones(1,length(x))*3+rand(1,length(x))*.1-.05,x,100,'filled','MarkerFaceColor',[252 216 213]/255,'MarkerEdgeColor','w')
plot([-0.125 0.125]+3,[xx(3) xx(3)],'Color',[229 45 38]/255,'LineWidth',1)
plot([3 3],[xx(2) xx(4)],'Color',[229 45 38]/255,'LineWidth',1)
[p(3),~,s] = signrank(x,0,'method','exact'); stat(3) = s.signedrank;
xlim([0.5 3.5]); xticks([1:3]); xticklabels({'Consume','Malaise','Retrieval'})
ylim([0 1]); yticks(0:.25:1); ylabel(['Animal-by-animal FOS correlation among',char(10),'individual amygdala network regions'])
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

p = multicmp(p,'up',0.05);
StatsTbl = table({'ED 4b'},{'Consume'},{'Wilcoxon signed-rank'},{'3 time points'},{[length(x)]},stat(1),p(1), ...
    'VariableNames',{'Figure panel','Group','Statistical test','Multiple comparisons','Sample size','Test statistic','P-value'});
StatsTbl(end+1,:) = table({'ED 4b'},{'Malaise'},{'Wilcoxon signed-rank'},{'3 time points'},{[length(x)]},stat(2),p(2));
StatsTbl(end+1,:) = table({'ED 4b'},{'Retrieval'},{'Wilcoxon signed-rank'},{'3 time points'},{[length(x)]},stat(3),p(3));
%% Fig ED4c

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 4c','FontWeight','bold')
load([data_path,'/FOS imaging/FOS-GLMM-statistics.mat'],'data');

idx = cluster(tree,'maxclust',10);
idx = idx(outperm);
list = unique(idx,'stable');
significant = find(data.GLMMoutput.Eq2.modelstats.significant);
treedata = data.GLMMoutput.Eq2.flavor.Zstat(significant,:);

AME = [];
for i = 1:length(list)
    AME(i,:) = mean(treedata(leafOrder(find(idx==list(i))),:));
end
B = corr(counts_norm(find(cellfun(@(x) isequal('Consumption',x),data.GLMMinput.timepoint)),significant(leafOrder))); t = [];
list = unique(idx,'stable');
for i = 1:length(list)
    for j = 1:length(list)
        t(i,j) = mean(B(find(idx==list(i)),find(idx==list(j))),'all');
    end
end
AmygCor(:,1) = t(:,find(list==idx(find(data.regions.name(significant(leafOrder))=="Central amygdalar nucleus, medial part"))));
B = corr(counts_norm(find(cellfun(@(x) isequal('Malaise',x),data.GLMMinput.timepoint)),significant(leafOrder))); t = [];
for i = 1:length(list)
    for j = 1:length(list)
        t(i,j) = mean(B(find(idx==list(i)),find(idx==list(j))),'all');
    end
end
AmygCor(:,2) = t(:,find(list==idx(find(data.regions.name(significant(leafOrder))=="Central amygdalar nucleus, medial part"))));
B = corr(counts_norm(find(cellfun(@(x) isequal('Retrieval',x),data.GLMMinput.timepoint)),significant(leafOrder))); t = [];
for i = 1:length(list)
    for j = 1:length(list)
        t(i,j) = mean(B(find(idx==list(i)),find(idx==list(j))),'all');
    end
end
AmygCor(:,3) = t(:,find(list==idx(find(data.regions.name(significant(leafOrder))=="Central amygdalar nucleus, medial part"))));

hold on
axis square
idx = setdiff(1:10,find(list==idx(find(data.regions.name(significant(leafOrder))=="Central amygdalar nucleus, medial part"))));
A = AME(idx,:);
B = AmygCor(idx,:);
a = scatter(A(:),B(:),100,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A(:)',B(:),1);
xlim([-4 4])
ylim([-.4 .4])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[y_fit] = polyval(p2,x,S);
fitresult = fit(A(:),B(:),'poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[.85 .85 .85],'LineStyle','none'),
plot(x,y_fit,'k','LineWidth',1)
for i = length(A):-1:1
    scatter(A(i,1),B(i,1),100,cmap(round((A(i,1)+5)*100),:),'o','filled','MarkerEdgeColor','k')
    scatter(A(i,2),B(i,2),100,cmap(round((A(i,2)+5)*100),:),'^','filled','MarkerEdgeColor','k')
    scatter(A(i,3),B(i,3),100,cmap(round((A(i,3)+5)*100),:),'square','filled','MarkerEdgeColor','k')
end
xlabel(['Novel – Familiar ΔFOS (Z),',char(10),'other clusters']); xlim([-4 4]); xticks(-4:2:4)
ylabel(['Animal-by-animal FOS correlation,',char(10),'amygdala network vs. other clusters']); ylim([-.4 .4]); yticks(-.4:.2:.4)
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

[r,p] = corr(A(:),B(:));
StatsTbl(end+1,:) = table({'ED 4c'},{'Amygdala cluster vs. Other clusters'},{'Pearson correlation'},{'N/A'},{[length(A(:))]},r,p);
%% Fig ED4d

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 4d','FontWeight','bold')
load([data_path,'/FOS imaging/FOS-GLMM-statistics.mat'],'data');
T1 = readtable([data_path,'/FOS imaging/region_info.csv']);

counts_norm_CEA = [T1{find(cellfun(@(x) isequal(x,'CentralAmygdalarNucleus'),T1.region)),4:end}./T1{1,4:end}./1.987]'*100;
counts_norm = [data.GLMMinput.counts./data.GLMMinput.totalcounts./data.regions.size']*100;

regionB = find(data.regions.name=='Agranular insular area, posterior part');
subplot(2,3,1)
hold on
axis square
A = counts_norm_CEA(find(cellfun(@(x) isequal('Consumption',x),data.GLMMinput.timepoint)));
B = counts_norm(find(cellfun(@(x) isequal('Consumption',x),data.GLMMinput.timepoint)),regionB);
a = scatter(A,B,64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([0 1])
ylim([0 1])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[y_fit] = polyval(p2,x,S);
fitresult = fit(A,B,'poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[.85 .85 .85],'LineStyle','none')
plot(x,y_fit,'k','LineWidth',1)
scatter(A,B,100,'k','filled','MarkerEdgeColor','w')
xlabel('CEA FOS (% per mm^3)'); xlim([0 1]); xticks(xlim);
ylabel('AIp FOS (% per mm^3)'); ylim([0 1]); yticks(ylim);
title('Consumption')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

[r,p] = corr(A,B);
StatsTbl(end+1,:) = table({'ED 4d, left'},{'CEA vs. AIp'},{'Pearson correlation'},{'N/A'},{[length(A(:))]},r,p);

regionB = find(data.regions.name=='Bed nuclei of the stria terminalis');
subplot(2,3,4)
hold on
axis square
A = counts_norm_CEA(find(cellfun(@(x) isequal('Consumption',x),data.GLMMinput.timepoint)));
B = counts_norm(find(cellfun(@(x) isequal('Consumption',x),data.GLMMinput.timepoint)),regionB);
a = scatter(A,B,64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([0 1])
ylim([0 1])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[y_fit] = polyval(p2,x,S);
fitresult = fit(A,B,'poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[.85 .85 .85],'LineStyle','none')
plot(x,y_fit,'k','LineWidth',1)
scatter(A,B,100,'k','filled','MarkerEdgeColor','w')
xlabel('CEA FOS (% per mm^3)'); xlim([0 1]); xticks(xlim);
ylabel('BST FOS (% per mm^3)'); ylim([0 1]); yticks(ylim);
title('Consumption')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

[r,p] = corr(A,B);
StatsTbl(end+1,:) = table({'ED 4d, left'},{'CEA vs. BST'},{'Pearson correlation'},{'N/A'},{[length(A(:))]},r,p);

regionB = find(data.regions.name=='Agranular insular area, posterior part');
subplot(2,3,2)
hold on
axis square
A = counts_norm_CEA(find(cellfun(@(x) isequal('Malaise',x),data.GLMMinput.timepoint)));
B = counts_norm(find(cellfun(@(x) isequal('Malaise',x),data.GLMMinput.timepoint)),regionB);
a = scatter(A,B,64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([0 1])
ylim([0 1])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[y_fit] = polyval(p2,x,S);
fitresult = fit(A,B,'poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[.85 .85 .85],'LineStyle','none')
plot(x,y_fit,'k','LineWidth',1)
scatter(A,B,100,'k','filled','MarkerEdgeColor','w')
xlabel('CEA FOS (% per mm^3)'); xlim([0 1]); xticks(xlim);
ylabel('AIp FOS (% per mm^3)'); ylim([0 1]); yticks(ylim);
title('Malaise')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

[r,p] = corr(A,B);
StatsTbl(end+1,:) = table({'ED 4d, middle'},{'CEA vs. AIp'},{'Pearson correlation'},{'N/A'},{[length(A(:))]},r,p);

regionB = find(data.regions.name=='Bed nuclei of the stria terminalis');
subplot(2,3,5)
hold on
axis square
A = counts_norm_CEA(find(cellfun(@(x) isequal('Malaise',x),data.GLMMinput.timepoint)));
B = counts_norm(find(cellfun(@(x) isequal('Malaise',x),data.GLMMinput.timepoint)),regionB);
a = scatter(A,B,64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([0 1])
ylim([0 1])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[y_fit] = polyval(p2,x,S);
fitresult = fit(A,B,'poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[.85 .85 .85],'LineStyle','none')
plot(x,y_fit,'k','LineWidth',1)
scatter(A,B,100,'k','filled','MarkerEdgeColor','w')
xlabel('CEA FOS (% per mm^3)'); xlim([0 1]); xticks(xlim);
ylabel('BST FOS (% per mm^3)'); ylim([0 1]); yticks(ylim);
title('Malaise')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

[r,p] = corr(A,B);
StatsTbl(end+1,:) = table({'ED 4d, middle'},{'CEA vs. BST'},{'Pearson correlation'},{'N/A'},{[length(A(:))]},r,p);

regionB = find(data.regions.name=='Agranular insular area, posterior part');
subplot(2,3,3)
hold on
axis square
A = counts_norm_CEA(find(cellfun(@(x) isequal('Retrieval',x),data.GLMMinput.timepoint)));
B = counts_norm(find(cellfun(@(x) isequal('Retrieval',x),data.GLMMinput.timepoint)),regionB);
a = scatter(A,B,64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([0 1])
ylim([0 1])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[y_fit] = polyval(p2,x,S);
fitresult = fit(A,B,'poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[.85 .85 .85],'LineStyle','none')
plot(x,y_fit,'k','LineWidth',1)
scatter(A,B,100,'k','filled','MarkerEdgeColor','w')
xlabel('CEA FOS (% per mm^3)'); xlim([0 1]); xticks(xlim);
ylabel('AIp FOS (% per mm^3)'); ylim([0 1]); yticks(ylim);
title('Retrieval')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

[r,p] = corr(A,B);
StatsTbl(end+1,:) = table({'ED 4d, right'},{'CEA vs. AIp'},{'Pearson correlation'},{'N/A'},{[length(A(:))]},r,p);

regionB = find(data.regions.name=='Bed nuclei of the stria terminalis');
subplot(2,3,6)
hold on
axis square
A = counts_norm_CEA(find(cellfun(@(x) isequal('Retrieval',x),data.GLMMinput.timepoint)));
B = counts_norm(find(cellfun(@(x) isequal('Retrieval',x),data.GLMMinput.timepoint)),regionB);
a = scatter(A,B,64,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([0 1])
ylim([0 1])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[y_fit] = polyval(p2,x,S);
fitresult = fit(A,B,'poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[.85 .85 .85],'LineStyle','none')
plot(x,y_fit,'k','LineWidth',1)
scatter(A,B,100,'k','filled','MarkerEdgeColor','w')
xlabel('CEA FOS (% per mm^3)'); xlim([0 1]); xticks(xlim);
ylabel('BST FOS (% per mm^3)'); ylim([0 1]); yticks(ylim);
title('Retrieval')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

[r,p] = corr(A,B);
StatsTbl(end+1,:) = table({'ED 4d, right'},{'CEA vs. BST'},{'Pearson correlation'},{'N/A'},{[length(A(:))]},r,p);

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