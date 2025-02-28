function fig_4(data_path)

disp('Generating panels for Figure 4...')
%% Fig 4b

figure('Position', get(0, 'Screensize'))
sgtitle('Figure 4b','FontWeight','bold')
load([data_path,'\Neuropixels\CGRP stim CFA\cgrp-stim-retrieval-summary.mat'],'data')

idx = find(data.stats.significant & data.stats.preference>0);
[~,idx2] = sort(data.stats.response.cgrp_period(idx),'descend'); idx = idx(idx2);
A1 = data.psth.smooth.novel(idx,:);
A2 = data.psth.smooth.water(idx,:);
A4 = data.psth.smooth.novel_retrieval(idx,:);
A5 = data.psth.smooth.water_retrieval(idx,:);
A6 = data.psth.cgrp_period(idx,:);

idx = find(data.stats.significant & data.stats.preference<0);
[~,idx2] = sort(data.stats.response.cgrp_period(idx),'descend'); idx = idx(idx2);
B1 = data.psth.smooth.novel(idx,:);
B2 = data.psth.smooth.water(idx,:);
B4 = data.psth.smooth.novel_retrieval(idx,:);
B5 = data.psth.smooth.water_retrieval(idx,:);
B6 = data.psth.cgrp_period(idx,:);

idx = find(~data.stats.significant);
[~,idx2] = sort(data.stats.response.cgrp_period(idx),'descend'); idx = idx(idx2);
C1 = data.psth.smooth.novel(idx,:);
C2 = data.psth.smooth.water(idx,:);
C4 = data.psth.smooth.novel_retrieval(idx,:);
C5 = data.psth.smooth.water_retrieval(idx,:);
C6 = data.psth.cgrp_period(idx,:);

pl1 = [A1(:,501:2000);nan(15,1500);B1(:,501:2000);nan(15,1500);C1(:,501:2000)];
pl2 = [A2(:,501:2000);nan(15,1500);B2(:,501:2000);nan(15,1500);C2(:,501:2000)];
pl = [pl1,nan(size(pl1,1),75),pl2];
cmap = flipud(cbrewer('div','RdBu',1000,'spline')); cmap(cmap<0) = 0;

subplot(1,3,1)
hold on
heatmap(flipud(pl),[],[],[],'Colormap',cmap,'ColorLevels',1000,'MaxColorValue',.5,'MinColorValue',-.5,'NaNColor',[1 1 1]);
nomod = size(C1,1); nopref = size(B1,1); flavorpref = size(A1,1);
yticks([nomod./2 nomod+nopref./2+15 nomod+nopref+flavorpref./2+30]+0.5)
ytickangle(90)
yticklabels({'Non-selective','Water','Flavor'})
xticks([750 750+1500+75])
xticklabels({'Flavor','Water'})
title('Conditioning day')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0, 0],'TickDir','out')
hold off

pl6 = [A6;nan(15,75);B6;nan(15,75);C6];
subplot(1,3,2)
hold on
heatmap(flipud(pl6),[],[],[],'Colormap',cmap,'ColorLevels',1000,'MaxColorValue',.5,'MinColorValue',-.5,'NaNColor',[1 1 1]);
xticks([0:15:75]+.5)
title('Conditioning day')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0, 0],'TickDir','out')
hold off

pl1 = [A4(:,501:2000);nan(15,1500);B4(:,501:2000);nan(15,1500);C4(:,501:2000)];
pl2 = [A5(:,501:2000);nan(15,1500);B5(:,501:2000);nan(15,1500);C5(:,501:2000)];
pl = [pl1,nan(size(pl1,1),75),pl2];
subplot(1,3,3)
hold on
heatmap(flipud(pl),[],[],[],'Colormap',cmap,'ColorLevels',1000,'MaxColorValue',.5,'MinColorValue',-.5,'NaNColor',[1 1 1]);
xticks([750 750+1500+75])
xticklabels({'Flavor','Water'})
title('Retrieval day')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0, 0],'TickDir','out')
hold off
%% Fig 4c

figure('Position', get(0, 'Screensize'))
sgtitle('Figure 4c','FontWeight','bold')
load([data_path,'\Neuropixels\CGRP stim CFA\cgrp-stim-retrieval-summary.mat'],'data')

subplot(2,2,1)
hold on
axis square
times = -5:.01:10;
times = times(1:end-1) + mean(diff(times))/2;
idx = find(data.stats.significant & data.stats.preference>0);
A1 = data.psth.smooth.novel(idx,501:end);
A2 = data.psth.smooth.novel_retrieval(idx,501:end);
fill([times fliplr(times)],[mean(A1)+std(A1)/sqrt(size(A1,1)) fliplr(mean(A1)-std(A1)/sqrt(size(A1,1)))],[.85 .85 .85],'LineStyle','none')
a = plot(times,mean(A1),'k','LineWidth',1);
fill([times fliplr(times)],[mean(A2)+std(A2)/sqrt(size(A2,1)) fliplr(mean(A2)-std(A2)/sqrt(size(A2,1)))],[252 216 213]/255,'LineStyle','none')
b = plot(times,mean(A2),'color',[229 45 38]/255,'LineWidth',1);
ylabel('Spiking (σ)'); ylim([-.2 .8]); yticks(-.2:.2:.8);
xlabel('Time (s)'); xlim([-5 10]); xticks(-5:5:10)
title('All Flavor-preferring')
legend([a,b],{'Conditioning day','Retrieval day'})
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

subplot(2,2,2)
hold on
axis square
times = -5:.01:10;
times = times(1:end-1) + mean(diff(times))/2;
idx = find(data.stats.significant & data.stats.preference>0);
[~,idx2] = sort(data.stats.response.cgrp_period(idx),'descend'); idx = idx(idx2(1:floor(length(idx2)*0.1)));
A1 = data.psth.smooth.novel(idx,501:end);
A2 = data.psth.smooth.novel_retrieval(idx,501:end);
fill([times fliplr(times)],[mean(A1)+std(A1)/sqrt(size(A1,1)) fliplr(mean(A1)-std(A1)/sqrt(size(A1,1)))],[.85 .85 .85],'LineStyle','none')
a = plot(times,mean(A1),'k','LineWidth',1);
fill([times fliplr(times)],[mean(A2)+std(A2)/sqrt(size(A2,1)) fliplr(mean(A2)-std(A2)/sqrt(size(A2,1)))],[252 216 213]/255,'LineStyle','none')
b = plot(times,mean(A2),'color',[229 45 38]/255,'LineWidth',1);
ylabel('Spiking (σ)'); ylim([-.2 .8]); yticks(-.2:.2:.8);
xlabel('Time (s)'); xlim([-5 10]); xticks(-5:5:10)
title('High CGRP response')
legend([a,b],{'Conditioning day','Retrieval day'})
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

subplot(2,2,3)
hold on
axis square
times = -5:.01:10;
times = times(1:end-1) + mean(diff(times))/2;
idx = find(data.stats.significant & data.stats.preference>0);
[~,idx2] = sort(data.stats.response.cgrp_period(idx),'descend'); idx = idx(idx2(floor(length(idx2)*0.1)+1:end));
A1 = data.psth.smooth.novel(idx,501:end);
A2 = data.psth.smooth.novel_retrieval(idx,501:end);
fill([times fliplr(times)],[mean(A1)+std(A1)/sqrt(size(A1,1)) fliplr(mean(A1)-std(A1)/sqrt(size(A1,1)))],[.85 .85 .85],'LineStyle','none')
a = plot(times,mean(A1),'k','LineWidth',1);
fill([times fliplr(times)],[mean(A2)+std(A2)/sqrt(size(A2,1)) fliplr(mean(A2)-std(A2)/sqrt(size(A2,1)))],[252 216 213]/255,'LineStyle','none')
b = plot(times,mean(A2),'color',[229 45 38]/255,'LineWidth',1);
ylabel('Spiking (σ)'); ylim([-.2 .8]); yticks(-.2:.2:.8);
xlabel('Time (s)'); xlim([-5 10]); xticks(-5:5:10)
title('Low CGRP response')
legend([a,b],{'Conditioning day','Retrieval day'})
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

subplot(2,2,4)
hold on
axis square
for i = 0:.1:2.9
    plot([i i],[.980 .995],'Color',[54 161 86]/255,'LineWidth',2)
end
times = -1:.01:4;
times = times(1:end-1) + mean(diff(times))/2;
idx = find(data.stats.significant & data.stats.preference>0);
[~,idx2] = sort(data.stats.response.cgrp_period(idx),'descend'); idx = idx(idx2(floor(length(idx2)*0.1)+1:end));
A1 = data.psth.raw.cgrp(idx,101:600);
fill([times fliplr(times)],[mean(A1)+std(A1)/sqrt(size(A1,1)) fliplr(mean(A1)-std(A1)/sqrt(size(A1,1)))],[.85 .85 .85],'LineStyle','none')
b = plot(times,mean(A1),'k','LineWidth',1);
idx = find(data.stats.significant & data.stats.preference>0);
[~,idx2] = sort(data.stats.response.cgrp_period(idx),'descend'); idx = idx(idx2(1:floor(length(idx2)*0.1)));
A1 = data.psth.raw.cgrp(idx,101:600);
fill([times fliplr(times)],[mean(A1)+std(A1)/sqrt(size(A1,1)) fliplr(mean(A1)-std(A1)/sqrt(size(A1,1)))],[252 216 213]/255,'LineStyle','none')
a = plot(times,mean(A1),'color',[229 45 38]/255,'LineWidth',1);
ylim([-.2 1])
yticks(-.2:.2:1)
legend([a,b],{'High CGRP response','Low CGRP response'})
xticks(-1:1:4)
xlim([-1 4])
xlabel('Time (s)')
ylabel('Spiking (σ)')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off
%% Fig 4d

figure('Position', get(0, 'Screensize'))
sgtitle('Figure 4d','FontWeight','bold')
load([data_path,'\Neuropixels\CGRP stim CFA\cgrp-stim-retrieval-summary.mat'],'data')

idx = find(data.stats.significant & data.stats.preference>0);
A = data.stats.response.cgrp_period(idx);
ranksum(data.stats.response.novel(idx),data.stats_retrieval.response.novel(idx));
B = data.stats_retrieval.response.novel_subtracted(idx)-data.stats.response.novel_subtracted(idx);
subplot(2,3,1)
hold on
axis square
a = scatter(A,B,100,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([-.5 2])
ylim([-.5 1.5])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[y_fit] = polyval(p2,x,S);
fitresult = fit(A',B','poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[252 216 213]/255,'LineStyle','none'),
plot(x,y_fit,'color',[229 45 38]/255,'LineWidth',1)
scatter(A,B,100,[229 45 38]/255,'filled','MarkerEdgeColor','w')
xlim([x(1) x(end)])
yticks(ylim); xticks(xlim);
xlabel('CGRP response (σ)')
ylabel('ΔNovel flavor (σ)')
title('Flavor-preferring')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

[r,p] = corr(A',B');
StatsTbl = table({'4d, left'},{'CGRP vs. ΔNovel flavor'},{'Pearson correlation'},{'N/A'},{[length(A)]},r,p, ...
    'VariableNames',{'Figure panel','Group','Statistical test','Multiple comparisons','Sample size','Test statistic','P-value'});

B = data.stats_retrieval.preference_subtracted(idx)-data.stats.preference_subtracted(idx);
subplot(2,3,4)
hold on
axis square
a = scatter(A,B,100,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([-.5 2])
ylim([-.5 1.5])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[y_fit] = polyval(p2,x,S);
fitresult = fit(A',B','poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[252 216 213]/255,'LineStyle','none'),
plot(x,y_fit,'color',[229 45 38]/255,'LineWidth',1)
scatter(A,B,100,[229 45 38]/255,'filled','MarkerEdgeColor','w')
xlim([x(1) x(end)])
yticks(ylim); xticks(xlim);
xlabel('CGRP response (σ)')
ylabel('ΔSelectivity (σ)')
title('Flavor-preferring')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

[r,p] = corr(A',B');
StatsTbl(end+1,:) = table({'4d, left'},{'CGRP vs. ΔSelectivity'},{'Pearson correlation'},{'N/A'},{[length(A)]},r,p);

idx = find(data.stats.significant & data.stats.preference<0);
A = data.stats.response.cgrp_period(idx);
B = data.stats_retrieval.response.novel_subtracted(idx)-data.stats.response.novel_subtracted(idx);
subplot(2,3,2)
hold on
axis square
a = scatter(A,B,100,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([-.5 2])
ylim([-.5 1.5])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[y_fit] = polyval(p2,x,S);
fitresult = fit(A',B','poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[216 231 243]/255,'LineStyle','none'),
plot(x,y_fit,'color',[55 136 193]/255,'LineWidth',1)
scatter(A,B,100,[55 136 193]/255,'filled','MarkerEdgeColor','w')
xlim([x(1) x(end)])
yticks(ylim); xticks(xlim);
xlabel('CGRP response (σ)')
ylabel('ΔNovel flavor (σ)')
title('Water-preferring')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

[r,p] = corr(A',B');
StatsTbl(end+1,:) = table({'4d, middle'},{'CGRP vs. ΔNovel flavor'},{'Pearson correlation'},{'N/A'},{[length(A)]},r,p);

B = data.stats_retrieval.preference_subtracted(idx)-data.stats.preference_subtracted(idx);
subplot(2,3,5)
hold on
axis square
a = scatter(A,B,100,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([-.5 2])
ylim([-.5 1.5])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[y_fit] = polyval(p2,x,S);
fitresult = fit(A',B','poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[216 231 243]/255,'LineStyle','none'),
plot(x,y_fit,'color',[55 136 193]/255,'LineWidth',1)
scatter(A,B,100,[55 136 193]/255,'filled','MarkerEdgeColor','w')
xlim([x(1) x(end)])
yticks(ylim); xticks(xlim);
xlabel('CGRP response (σ)')
ylabel('ΔSelectivity (σ)')
title('Water-preferring')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

[r,p] = corr(A',B');
StatsTbl(end+1,:) = table({'4d, middle'},{'CGRP vs. ΔSelectivity'},{'Pearson correlation'},{'N/A'},{[length(A)]},r,p);

idx = find(~data.stats.significant);
A = data.stats.response.cgrp_period(idx);
B = data.stats_retrieval.response.novel_subtracted(idx)-data.stats.response.novel_subtracted(idx);
subplot(2,3,3)
hold on
axis square
a = scatter(A,B,100,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([-.5 2])
ylim([-.5 1.5])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[y_fit] = polyval(p2,x,S);
fitresult = fit(A',B','poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[.85 .85 .85],'LineStyle','none'),
plot(x,y_fit,'k','LineWidth',1)
scatter(A,B,100,'k','filled','MarkerEdgeColor','w')
xlim([x(1) x(end)])
yticks(ylim); xticks(xlim);
xlabel('CGRP response (σ)')
ylabel('ΔNovel flavor (σ)')
title('Non-selective')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

[r,p] = corr(A',B');
StatsTbl(end+1,:) = table({'4d, right'},{'CGRP vs. ΔNovel flavor'},{'Pearson correlation'},{'N/A'},{[length(A)]},r,p);

B = data.stats_retrieval.preference_subtracted(idx)-data.stats.preference_subtracted(idx);
subplot(2,3,6)
hold on
axis square
a = scatter(A,B,100,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([-.5 2])
ylim([-.5 1.5])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[y_fit] = polyval(p2,x,S);
fitresult = fit(A',B','poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[.85 .85 .85],'LineStyle','none'),
plot(x,y_fit,'k','LineWidth',1)
scatter(A,B,100,'k','filled','MarkerEdgeColor','w')
xlim([x(1) x(end)])
yticks(ylim); xticks(xlim);
xlabel('CGRP response (σ)')
ylabel('ΔSelectivity (σ)')
title('Non-selective')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

[r,p] = corr(A',B');
StatsTbl(end+1,:) = table({'4d, right'},{'CGRP vs. ΔSelectivity'},{'Pearson correlation'},{'N/A'},{[length(A)]},r,p);
%% Fig 4e

figure('Position', get(0, 'Screensize'))
sgtitle('Figure 4e','FontWeight','bold')
load([data_path,'\Neuropixels\CGRP-CEA stim CFA\cgrp-cea-stim-retrieval-summary.mat'],'data')

idx = find(data.stats.significant & data.stats.preference>0);
A = data.stats.response.cgrp_period(idx);
ranksum(data.stats.response.novel(idx),data.stats_retrieval.response.novel(idx));
B = data.stats_retrieval.response.novel_subtracted(idx)-data.stats.response.novel_subtracted(idx);
subplot(1,2,1)
hold on
axis square
a = scatter(A,B,100,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([-.5 1.5])
ylim([-.5 1])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[y_fit] = polyval(p2,x,S);
fitresult = fit(A',B','poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[252 216 213]/255,'LineStyle','none'),
plot(x,y_fit,'color',[229 45 38]/255,'LineWidth',1)
scatter(A,B,100,[229 45 38]/255,'filled','MarkerEdgeColor','w')
xlim([x(1) x(end)])
yticks(ylim); xticks(xlim);
xlabel('CGRP^{CEA} response (σ)')
ylabel('ΔNovel flavor (σ)')
title('Flavor-preferring')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

[r,p] = corr(A',B');
StatsTbl(end+1,:) = table({'4e'},{'CGRP vs. ΔNovel flavor'},{'Pearson correlation'},{'N/A'},{[length(A)]},r,p);

B = data.stats_retrieval.preference_subtracted(idx)-data.stats.preference_subtracted(idx);
subplot(1,2,2)
hold on
axis square
a = scatter(A,B,100,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A',B,1);
xlim([-.5 1.5])
ylim([-.5 1])
x = xlim;
x = x(1):.01:x(2);
delete(a)
[y_fit] = polyval(p2,x,S);
fitresult = fit(A',B','poly1');
p2 = predint(fitresult,x,0.95,'functional');
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[252 216 213]/255,'LineStyle','none'),
plot(x,y_fit,'color',[229 45 38]/255,'LineWidth',1)
scatter(A,B,100,[229 45 38]/255,'filled','MarkerEdgeColor','w')
xlim([x(1) x(end)])
yticks(ylim); xticks(xlim);
xlabel('CGRP^{CEA} response (σ)')
ylabel('ΔSelectivity (σ)')
title('Flavor-preferring')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

[r,p] = corr(A',B');
StatsTbl(end+1,:) = table({'4e'},{'CGRP vs. ΔSelectivity'},{'Pearson correlation'},{'N/A'},{[length(A)]},r,p);
%% Fig 4f

figure('Position', get(0, 'Screensize'))
sgtitle('Figure 4f','FontWeight','bold')
load([data_path,'\Neuropixels\Familiarization\familiarization-summary.mat'],'data')

subplot(1,2,1)
hold on
axis square
times = -5:.01:10;
times = times(1:end-1) + mean(diff(times))/2;
idx = find(data.stats.significant & data.stats.preference>0);
A1 = data.psth.smooth.novel(idx,501:end);
A2 = data.psth.smooth.novel_retrieval(idx,501:end);
fill([times fliplr(times)],[mean(A1)+std(A1)/sqrt(size(A1,1)) fliplr(mean(A1)-std(A1)/sqrt(size(A1,1)))],[.85 .85 .85],'LineStyle','none')
a = plot(times,mean(A1),'k','LineWidth',1);
fill([times fliplr(times)],[mean(A2)+std(A2)/sqrt(size(A2,1)) fliplr(mean(A2)-std(A2)/sqrt(size(A2,1)))],[216 231 243]/255,'LineStyle','none')
b = plot(times,mean(A2),'color',[55 136 193]/255,'LineWidth',1);
ylim([-.1 .3])
yticks(-.1:.1:.3)
xticks(-5:5:10)
xlim([-5 10])
xlabel('Time (s)')
ylabel('Spiking (σ)')
title('All Flavor-preferring')
legend([a,b],{'Novel day','Familiar day'})
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

subplot(1,2,2)
axis square
hold on
idx = find(data.stats.significant & data.stats.preference>0);
A = data.stats.response.novel(idx);
B = data.stats_retrieval.response.novel(idx);
simpleboxplot(1,A,'k')
simpleboxplot(2,B,[55 136 193]/255)
ylabel('Flavor response (σ)'); ylim([-.1 .2]); yticks([-.1 .2])
xticks([1,2])
xticklabels({'Novel day','Familiar day'})
xlim([0.25 2.75])
title('All Flavor-preferring')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

[p,~,stat] = signrank(A,B,'method','exact');
StatsTbl(end+1,:) = table({'4f'},{'Novel vs. Familiar'},{'Wilcoxon signed-rank'},{'N/A'},{[length(A)]},stat.signedrank,p);

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