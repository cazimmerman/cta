function plotphotometryPKA(data,port)

%% housekeeping
warning('off')
cmap = flipud(cbrewer('div','RdBu',1000,'spline')); cmap(cmap<0) = 0;
warning('on')

idx = min([data.corrected.TTL2(1,1),data.corrected.TTL5(1,1)]);
idx = [find(data.corrected.time>(idx-60*5),1) find(data.corrected.time==max(data.corrected.time),1)];
x = (data.corrected.time(idx(1):idx(2))-data.corrected.time(idx(1)))./60-5;
y = data.corrected.AKAR465(idx(1):idx(2))*1000;

data.corrected.PSTH_rear465 = data.corrected.PSTH_rear465./data.corrected.PSTH_rear405;
data.corrected.PSTH_front465 = data.corrected.PSTH_front465./data.corrected.PSTH_front405;
if port
    data.corrected.PSTH_water465 = data.corrected.PSTH_rear465 - mean(data.corrected.PSTH_rear465(:,1:1000),2);
    data.corrected.PSTH_flavor465 = data.corrected.PSTH_front465 - mean(data.corrected.PSTH_front465(:,1:1000),2);
    data.corrected.PSTH_water405 = data.corrected.PSTH_rear405 - mean(data.corrected.PSTH_rear405(:,1:1000),2);
    data.corrected.PSTH_flavor405 = data.corrected.PSTH_front405 - mean(data.corrected.PSTH_front405(:,1:1000),2);
    for i = 1:length(data.corrected.TTL2)
        if i < length(data.corrected.TTL2) || data.corrected.TTL2(i,1) < data.corrected.TTL5(end,1)
            if i < length(data.corrected.TTL2)
                data.corrected.NextReward_water(i) = min([data.corrected.TTL2(i+1,1) data.corrected.TTL5(find(data.corrected.TTL5(:,1)>data.corrected.TTL2(i,1),1,'first'))]) - data.corrected.TTL2(i,1);
            else
                data.corrected.NextReward_water(i) = data.corrected.TTL5(find(data.corrected.TTL5(:,1)>data.corrected.TTL2(i,1),1,'first')) - data.corrected.TTL2(i,1);
            end
        end
    end
    for i = 1:length(data.corrected.TTL5)
        if i < length(data.corrected.TTL5) || data.corrected.TTL5(i,1) < data.corrected.TTL2(end,1)
            if i < length(data.corrected.TTL5)
                data.corrected.NextReward_flavor(i) = min([data.corrected.TTL5(i+1,1) data.corrected.TTL2(find(data.corrected.TTL2(:,1)>data.corrected.TTL5(i,1),1,'first'))]) - data.corrected.TTL5(i,1);
            else
                data.corrected.NextReward_flavor(i) = data.corrected.TTL2(find(data.corrected.TTL2(:,1)>data.corrected.TTL5(i,1),1,'first')) - data.corrected.TTL5(i,1);
            end
        end
    end
else
    data.corrected.PSTH_flavor465 = data.corrected.PSTH_rear465 - mean(data.corrected.PSTH_rear465(:,1:1000),2);
    data.corrected.PSTH_water465 = data.corrected.PSTH_front465 - mean(data.corrected.PSTH_front465(:,1:1000),2);
    data.corrected.PSTH_flavor405 = data.corrected.PSTH_rear405 - mean(data.corrected.PSTH_rear405(:,1:1000),2);
    data.corrected.PSTH_water405 = data.corrected.PSTH_front405 - mean(data.corrected.PSTH_front405(:,1:1000),2);
    for i = 1:length(data.corrected.TTL2)
        if i < length(data.corrected.TTL2) || data.corrected.TTL2(i,1) < data.corrected.TTL5(end,1)
            if i < length(data.corrected.TTL2)
                data.corrected.NextReward_flavor(i) = min([data.corrected.TTL2(i+1,1) data.corrected.TTL5(find(data.corrected.TTL5(:,1)>data.corrected.TTL2(i,1),1,'first'))]) - data.corrected.TTL2(i,1);
            else
                data.corrected.NextReward_flavor(i) = data.corrected.TTL5(find(data.corrected.TTL5(:,1)>data.corrected.TTL2(i,1),1,'first')) - data.corrected.TTL2(i,1);
            end
        end
    end
    for i = 1:length(data.corrected.TTL5)
        if i < length(data.corrected.TTL5) || data.corrected.TTL5(i,1) < data.corrected.TTL2(end,1)
            if i < length(data.corrected.TTL5)
                data.corrected.NextReward_water(i) = min([data.corrected.TTL5(i+1,1) data.corrected.TTL2(find(data.corrected.TTL2(:,1)>data.corrected.TTL5(i,1),1,'first'))]) - data.corrected.TTL5(i,1);
            else
                data.corrected.NextReward_water(i) = data.corrected.TTL2(find(data.corrected.TTL2(:,1)>data.corrected.TTL5(i,1),1,'first')) - data.corrected.TTL5(i,1);
            end
        end
    end
end
sigma465 = [data.corrected.PSTH_flavor465(:,1:1000);data.corrected.PSTH_water465(:,1:1000)]; sigma465 = std(sigma465(:));
data.corrected.PSTH_flavor465 = data.corrected.PSTH_flavor465./sigma465;
data.corrected.PSTH_water465 = data.corrected.PSTH_water465./sigma465;
sigma405 = [data.corrected.PSTH_flavor405(:,1:1000);data.corrected.PSTH_water405(:,1:1000)]; sigma405 = std(sigma405(:));
data.corrected.PSTH_flavor405 = data.corrected.PSTH_flavor405./sigma405;
data.corrected.PSTH_water405 = data.corrected.PSTH_water405./sigma405;

data.corrected.PSTH_flavor = data.corrected.PSTH_flavor465;%-data.corrected.PSTH_flavor405;
data.corrected.PSTH_water = data.corrected.PSTH_water465;%-data.corrected.PSTH_water405;
%% generate plots

figure
subplot(2,2,1)
hold on
axis square
%lims = [floor(min(y)) ceil(max(y))];
lims = [floor(min(y)) ceil(max(y))];
if port
    for i = 1:length(data.corrected.TTL2)
        t = (data.corrected.TTL2(i,1)-data.corrected.time(idx(1)))/60-5;
        X = fill([t t+15/60 t+15/60 t],[lims(1) lims(1) lims(2) lims(2)],[216 231 243]/255,'LineStyle','none');
    end
    for i = 1:length(data.corrected.TTL5)
        t = (data.corrected.TTL5(i,1)-data.corrected.time(idx(1)))/60-5;
        Y = fill([t t+15/60 t+15/60 t],[lims(1) lims(1) lims(2) lims(2)],[252 216 213]/255,'LineStyle','none');
    end
else
    for i = 1:length(data.corrected.TTL2)
        t = (data.corrected.TTL2(i,1)-data.corrected.time(idx(1)))/60-5;
        Y = fill([t t+15/60 t+15/60 t],[lims(1) lims(1) lims(2) lims(2)],[252 216 213]/255,'LineStyle','none');
    end
    for i = 1:length(data.corrected.TTL5)
        t = (data.corrected.TTL5(i,1)-data.corrected.time(idx(1)))/60-5;
        X = fill([t t+15/60 t+15/60 t],[lims(1) lims(1) lims(2) lims(2)],[216 231 243]/255,'LineStyle','none');
    end
end
plot(x,y,'k','LineWidth',1)
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 1]);
xlabel('Time (min) from first reward')
ylabel('PKA activity (mV)')
ylim(lims)
xlim([-5 max(x)])
xticks(-5:5:60)
title('Full session')
legend([Y,X],{['Flavor port (',num2str(size(data.corrected.PSTH_water,1)),')'],['Water port (',num2str(size(data.corrected.PSTH_flavor,1)),')']},'Location','NorthEast','box','on')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
if isequal(data.tdt.fname,'Z:\Chris\data\photometry\CZ_photometry\DataTanks\CZ_DataTanks\20240615_pka702')
    xlim([-5 44.5]); ylim([378 421]);
end
hold off

lims = [-5 10];
subplot(2,2,2)
hold on
axis square
plot([0 0],[-5 15],'k--','LineWidth',1)
times = linspace(-10,60,size(data.corrected.PSTH_flavor,2));
fill([times fliplr(times)],[mean(data.corrected.PSTH_water)+std(data.corrected.PSTH_water)./sqrt(size(data.corrected.PSTH_water,1)) fliplr(mean(data.corrected.PSTH_water)-std(data.corrected.PSTH_water)./sqrt(size(data.corrected.PSTH_water,1)))],[216 231 243]/255,'LineStyle','none')
X=plot(times,mean(data.corrected.PSTH_water),'Color',[55 136 193]/255,'LineWidth',1);
fill([times fliplr(times)],[mean(data.corrected.PSTH_flavor)+std(data.corrected.PSTH_flavor)./sqrt(size(data.corrected.PSTH_flavor,1)) fliplr(mean(data.corrected.PSTH_flavor)-std(data.corrected.PSTH_flavor)./sqrt(size(data.corrected.PSTH_flavor,1)))],[252 216 213]/255,'LineStyle','none')
Y=plot(times,mean(data.corrected.PSTH_flavor),'Color',[229 45 38]/255,'LineWidth',1);

pl = data.corrected.NextReward_flavor;
pl2 = prctile(pl,[25 50 75]);
scatter(pl,[lims(2)-.1*diff(lims)]*ones(size(pl))+.4*rand(size(pl))-.2,16,'MarkerFaceColor',[252 216 213]/255,'MarkerEdgeColor','w','LineWidth',1)
plot([pl2(1) pl2(3)],[lims(2)-.1*diff(lims) lims(2)-.1*diff(lims)],'Color',[229 45 38]/255,'LineWidth',1)
plot([pl2(2) pl2(2)],[lims(2)-.085*diff(lims) lims(2)-.115*diff(lims)],'Color',[229 45 38]/255,'LineWidth',1)

pl = data.corrected.NextReward_water;
pl2 = prctile(pl,[25 50 75]);
scatter(pl,[lims(2)-.15*diff(lims)]*ones(size(pl))+.4*rand(size(pl))-.2,16,'MarkerFaceColor',[216 231 243]/255,'MarkerEdgeColor','w','LineWidth',1)
plot([pl2(1) pl2(3)],[lims(2)-.15*diff(lims) lims(2)-.15*diff(lims)],'Color',[55 136 193]/255,'LineWidth',1)
plot([pl2(2) pl2(2)],[lims(2)-.135*diff(lims) lims(2)-.165*diff(lims)],'Color',[55 136 193]/255,'LineWidth',1)
text(25,lims(2)-.035*diff(lims),'Next reward','FontSize',16,'HorizontalAlignment','Center')

xlabel('Time (s) from reward')
ylabel('PKA activity (σ)')
xlim([-10 60])
ylim(lims)
legend([Y,X],{['Flavor port (',num2str(size(data.corrected.PSTH_flavor,1)),')'],['Water port (',num2str(size(data.corrected.PSTH_water,1)),')']},'Location','SouthWest','box','on')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
title('Average reward PETHs')
hold off

subplot(2,2,3)
hold on
axis square
heatmap(data.corrected.PSTH_water,[],[],[],'Colormap',cmap,'ColorLevels',1000,'MaxColorValue',max(abs(lims)),'MinColorValue',-max(abs(lims)),'NaNColor',[1 1 1]);
plot([find(times>=0,1),find(times>=0,1)]+0.5,[0 size(data.corrected.PSTH_water,1)]+0.5,'k--','LineWidth',1)
fill([0 size(data.corrected.PSTH_water,2) size(data.corrected.PSTH_water,2) 0]+0.5,[0 0 size(data.corrected.PSTH_water,1) size(data.corrected.PSTH_water,1)]+0.5,'k','FaceAlpha',0,'LineWidth',1)
xlabel('Time (s) from reward ')
ylabel('Trial number')
title('Individual reward PETHs: Water port')
yticks([0:5:30]+0.5)
yticklabels({'0','5','10','15','20','25','30'})
xticks([find(times>=-10,1) find(times>=0,1) find(times>=10,1) find(times>=20,1) find(times>=30,1) find(times>=40,1) find(times>=50,1) find(times>=60,1)]+0.5)
xticklabels({'-10','0','10','20','30','40','50','60'})
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
h = colorbar('FontSize',16,'LineWidth',1);
ylabel(h,'PKA activity (σ)');
h.Limits = lims;
h.Ticks = [lims(1) 0 lims(2)];
h.TickLength = 0;
h.Box = 'on';
hold off

subplot(2,2,4)
hold on
axis square
heatmap(data.corrected.PSTH_flavor,[],[],[],'Colormap',cmap,'ColorLevels',1000,'MaxColorValue',max(abs(lims)),'MinColorValue',-max(abs(lims)),'NaNColor',[1 1 1]);
plot([find(times>=0,1),find(times>=0,1)]+0.5,[0 size(data.corrected.PSTH_flavor,1)]+0.5,'k--','LineWidth',1)
fill([0 size(data.corrected.PSTH_flavor,2) size(data.corrected.PSTH_flavor,2) 0]+0.5,[0 0 size(data.corrected.PSTH_flavor,1) size(data.corrected.PSTH_flavor,1)]+0.5,'k','FaceAlpha',0,'LineWidth',1)
yticks([0:5:30]+0.5)
yticklabels({'0','5','10','15','20','25','30'})
xticks([find(times>=-10,1) find(times>=0,1) find(times>=10,1) find(times>=20,1) find(times>=30,1) find(times>=40,1) find(times>=50,1) find(times>=60,1)]+0.5)
xticklabels({'-10','0','10','20','30','40','50','60'})
xlabel('Time (s) from reward')
ylabel('Trial number')
title('Individual reward PETHs: Flavor port')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
h = colorbar('FontSize',16,'LineWidth',1);
ylabel(h,'PKA activity (σ)');
h.Limits = lims;
h.Ticks = [lims(1) 0 lims(2)];
h.TickLength = 0;
h.Box = 'on';
hold off
%% save plots
sgtitle(data.tdt.fname(end-14:end),'interpreter','none','FontWeight','bold','FontSize',24)
saveas(gcf,['figures/',data.tdt.fname(end-14:end)],'png')
% set(gcf,'renderer','Painters')
% saveas(gcf,['figures/',data.tdt.fname(end-14:end)],'epsc')
end