function fig_ED8(data_path)

disp('Generating panels for Extended Data Figure 8...')
%% Fig ED8a

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 8a','FontWeight','bold')
load([data_path,'/Neuropixels/CGRP stim CFA/cgrp-stim-conditioning-behavior.mat'],'data')

axis square
hold on
fill([data.time fliplr(data.time)],[mean(data.water)+std(data.water)/sqrt(size(data.water,1)) fliplr(mean(data.water)-std(data.water)/sqrt(size(data.water,1)))]*0.02,[216 231 243]/255,'LineStyle','none');
a = plot(data.time,mean(data.water)*0.02,'Color',[55 136 193]/255,'LineWidth',1);
fill([data.time fliplr(data.time)],[mean(data.flavor)+std(data.flavor)/sqrt(size(data.flavor,1)) fliplr(mean(data.flavor)-std(data.flavor)/sqrt(size(data.flavor,1)))]*0.02,[252 216 213]/255,'LineStyle','none');
b = plot(data.time,mean(data.flavor)*0.02,'Color',[229 45 38]/255,'LineWidth',1);
xlabel('Time (min)'); xlim([0 25]); xticks(0:5:25);
ylabel('Cumulative intake (ml)'); ylim([0 0.6]); yticks(0:.2:.6);
legend([a,b],{'Water','Novel flavor'})
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off
%% Fig ED8b

figure('Position', get(0, 'Screensize'))
sgtitle('Extended Data Figure 8b','FontWeight','bold')
load([data_path,'/Neuropixels/CGRP stim CFA/cgrp-stim-conditioning-summary.mat'],'data')

subplot(1,3,1)
axis square
hold on
times = -5:.01:10; times = times(1:end-1) + mean(diff(times))/2;
idx = find(~data.stats.significant);
fill([times fliplr(times)],[mean(data.psth.smooth.novel(idx,501:end))+std(data.psth.smooth.novel(idx,501:end))/sqrt(size(data.psth.smooth.novel(idx,501:end),1)) fliplr(mean(data.psth.smooth.novel(idx,501:end))-std(data.psth.smooth.novel(idx,501:end))/sqrt(size(data.psth.smooth.novel(idx,501:end),1)))],[.85 .85 .85],'LineStyle','none')
a = plot(times,mean(data.psth.smooth.novel(idx,501:end)),'color','k','LineWidth',1);
idx = find(data.stats.significant & data.stats.preference<0);
fill([times fliplr(times)],[mean(data.psth.smooth.novel(idx,501:end))+std(data.psth.smooth.novel(idx,501:end))/sqrt(size(data.psth.smooth.novel(idx,501:end),1)) fliplr(mean(data.psth.smooth.novel(idx,501:end))-std(data.psth.smooth.novel(idx,501:end))/sqrt(size(data.psth.smooth.novel(idx,501:end),1)))],[216 231 243]/255,'LineStyle','none')
b = plot(times,mean(data.psth.smooth.novel(idx,501:end)),'color',[55 136 193]/255,'LineWidth',1);
idx = find(data.stats.significant & data.stats.preference>0);
fill([times fliplr(times)],[mean(data.psth.smooth.novel(idx,501:end))+std(data.psth.smooth.novel(idx,501:end))/sqrt(size(data.psth.smooth.novel(idx,501:end),1)) fliplr(mean(data.psth.smooth.novel(idx,501:end))-std(data.psth.smooth.novel(idx,501:end))/sqrt(size(data.psth.smooth.novel(idx,501:end),1)))],[252 216 213]/255,'LineStyle','none')
c = plot(times,mean(data.psth.smooth.novel(idx,501:end)),'color',[229 45 38]/255,'LineWidth',1);
xlabel('Time (s)'); xlim([-5 10]); xticks(-5:5:10);
ylabel('Spiking (σ)'); ylim([-.1 .2]); yticks(-.1:.1:.2);
title('Novel flavor response')
legend([c,b,a],{'Flavor-pref','Water-pref','Non-selective'})
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

subplot(1,3,2)
axis square
hold on
times = -5:.01:10; times = times(1:end-1) + mean(diff(times))/2;
idx = find(~data.stats.significant);
fill([times fliplr(times)],[mean(data.psth.smooth.water(idx,501:end))+std(data.psth.smooth.water(idx,501:end))/sqrt(size(data.psth.smooth.water(idx,501:end),1)) fliplr(mean(data.psth.smooth.water(idx,501:end))-std(data.psth.smooth.water(idx,501:end))/sqrt(size(data.psth.smooth.water(idx,501:end),1)))],[.85 .85 .85],'LineStyle','none')
a = plot(times,mean(data.psth.smooth.water(idx,501:end)),'color','k','LineWidth',1);
idx = find(data.stats.significant & data.stats.preference<0);
fill([times fliplr(times)],[mean(data.psth.smooth.water(idx,501:end))+std(data.psth.smooth.water(idx,501:end))/sqrt(size(data.psth.smooth.water(idx,501:end),1)) fliplr(mean(data.psth.smooth.water(idx,501:end))-std(data.psth.smooth.water(idx,501:end))/sqrt(size(data.psth.smooth.water(idx,501:end),1)))],[216 231 243]/255,'LineStyle','none')
b = plot(times,mean(data.psth.smooth.water(idx,501:end)),'color',[55 136 193]/255,'LineWidth',1);
idx = find(data.stats.significant & data.stats.preference>0);
fill([times fliplr(times)],[mean(data.psth.smooth.water(idx,501:end))+std(data.psth.smooth.water(idx,501:end))/sqrt(size(data.psth.smooth.water(idx,501:end),1)) fliplr(mean(data.psth.smooth.water(idx,501:end))-std(data.psth.smooth.water(idx,501:end))/sqrt(size(data.psth.smooth.water(idx,501:end),1)))],[252 216 213]/255,'LineStyle','none')
c = plot(times,mean(data.psth.smooth.water(idx,501:end)),'color',[229 45 38]/255,'LineWidth',1);
xlabel('Time (s)'); xlim([-5 10]); xticks(-5:5:10);
ylabel('Spiking (σ)'); ylim([-.1 .2]); yticks(-.1:.1:.2);
title('Water response')
legend([c,b,a],{'Flavor-pref','Water-pref','Non-selective'})
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

subplot(1,3,3)
hold on
axis square
a=pie([sum(data.stats.significant & data.stats.preference>0) sum(data.stats.significant & data.stats.preference<0) sum(~data.stats.significant)],{'Flavor-pref','Water-pref','Non-selective'});
a(1).FaceColor = [229 45 38]/255; a(1).EdgeColor = [1 1 1];
a(3).FaceColor = [55 136 193]/255; a(3).EdgeColor = [1 1 1];
a(5).FaceColor = [.85 .85 .85]; a(5).EdgeColor = [1 1 1];
xticks([]); yticks([])
title('Flavor selectivity')
axis off
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

drawnow