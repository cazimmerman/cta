function plotphotometryPKA_final(mice,sessions,ports)

warning('off')
cmap = flipud(cbrewer('div','RdBu',1000,'spline')); cmap(cmap<0) = 0;
warning('on')

expts = {' (Water)',' (Novel)',' (Familiar)',' (Familiar)'};
% for i = 1:length(mice)
%     Rear465 = [];
%     Front465 = [];
%     Rear405 = [];
%     Front405 = [];
%     for j = 1:length(sessions)
%         clear data
%         load(['data/',sessions{j},'_',mice{i}],'data')
%         Rear465  = [Rear465;  data.corrected.PSTH_rear465 - mean(data.corrected.PSTH_rear465(:,1:1000),2)];
%         Front465 = [Front465; data.corrected.PSTH_front465 - mean(data.corrected.PSTH_front465(:,1:1000),2)];
%         Rear405  = [Rear405;  data.corrected.PSTH_rear405 - mean(data.corrected.PSTH_rear405(:,1:1000),2)];
%         Front405 = [Front405; data.corrected.PSTH_front405 - mean(data.corrected.PSTH_front405(:,1:1000),2)];
%     end
%     X = [Rear465(:,1:1000); Front465(:,1:1000)]; sigma465(i) = std(X(:));
%     X = [Rear405(:,1:1000); Front405(:,1:1000)]; sigma405(i) = std(X(:));
% end

figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 1]);
lims = [-2 5];

load('summary.mat')
[~,mouse_order]=sort(Flavor(2,:),'ascend');

N = length(mice);
Flavor = [];
Water = [];
Flavor_NovelDay = nan(length(mice),30);
Water_NovelDay = nan(length(mice),30);
Next_FlavorFlavor = [];
Next_WaterWater = [];
for j = 1:length(sessions)
    DATA = struct;
    DATA.Full_AKAR465 = nan(N,1E6);
    DATA.Full_Time = [];
    DATA.NextReward_water = [];
    DATA.NextReward_flavor = [];
    for k = 1:N
        clear data
        load(['data/',sessions{j},'_',mice{k}],'data')
        idx = min([data.corrected.TTL2(1,1),data.corrected.TTL5(1,1)]);
        idx = [find(data.corrected.time>(idx-60*5),1) find(data.corrected.time>(idx),1) find(data.corrected.time==max(data.corrected.time),1)];
        if idx(3)-idx(1)==269999
            idx(3) = idx(3)+1;
        end
        
        data.corrected.PSTH_rear465 = data.corrected.PSTH_rear465./data.corrected.PSTH_rear405;
        data.corrected.PSTH_front465 = data.corrected.PSTH_front465./data.corrected.PSTH_front405;
        
        data.corrected.PSTH_rear465 = data.corrected.PSTH_rear465 - mean(data.corrected.PSTH_rear465(:,501:900),2);
        data.corrected.PSTH_front465 = data.corrected.PSTH_front465 - mean(data.corrected.PSTH_front465(:,501:900),2);
        sigma465 = [data.corrected.PSTH_rear465(:,501:900);data.corrected.PSTH_front465(:,501:900)]; sigma465 = std(sigma465(:));
        data.corrected.PSTH_rear465 = data.corrected.PSTH_rear465./sigma465;
        data.corrected.PSTH_front465 = data.corrected.PSTH_front465./sigma465;
        data.corrected.PSTH_rear405 = data.corrected.PSTH_rear405 - mean(data.corrected.PSTH_rear405(:,501:900),2);
        data.corrected.PSTH_front405 = data.corrected.PSTH_front405 - mean(data.corrected.PSTH_front405(:,501:900),2);
        sigma405 = [data.corrected.PSTH_rear405(:,501:900);data.corrected.PSTH_front405(:,501:900)]; sigma405 = std(sigma405(:));
        data.corrected.PSTH_rear405 = data.corrected.PSTH_rear405./sigma405;
        data.corrected.PSTH_front405 = data.corrected.PSTH_front405./sigma405;
        
        %         data.corrected.PSTH_rear465 = data.corrected.PSTH_rear465 - mean(data.corrected.PSTH_rear465(:,1:1000),2);
        %         data.corrected.PSTH_front465 = data.corrected.PSTH_front465 - mean(data.corrected.PSTH_front465(:,1:1000),2);
        %         data.corrected.PSTH_rear465 = data.corrected.PSTH_rear465./sigma465(i);
        %         data.corrected.PSTH_front465 = data.corrected.PSTH_front465./sigma465(i);
        %         data.corrected.PSTH_rear405 = data.corrected.PSTH_rear405 - mean(data.corrected.PSTH_rear405(:,1:1000),2);
        %         data.corrected.PSTH_front405 = data.corrected.PSTH_front405 - mean(data.corrected.PSTH_front405(:,1:1000),2);
        %         data.corrected.PSTH_rear405 = data.corrected.PSTH_rear405./sigma405(i);
        %         data.corrected.PSTH_front405 = data.corrected.PSTH_front405./sigma405(i);
        
        if j == 2
            x = []; y = [];
            for i = 1:29
                if data.corrected.TTL2(i+1,1) < data.corrected.TTL5(find(data.corrected.TTL5(:,1)>data.corrected.TTL2(i,1),1,'first'),1) & (data.corrected.TTL2(i+1,1)-data.corrected.TTL2(i,1))<60 & ...
                   ~any(data.corrected.TTL5(data.corrected.TTL5(:,1)>data.corrected.TTL2(i,1),1)-data.corrected.TTL2(i,1) < 60)
                    x = [x i];
                end
                if data.corrected.TTL5(i+1,1) < data.corrected.TTL2(find(data.corrected.TTL2(:,1)>data.corrected.TTL5(i,1),1,'first'),1) & (data.corrected.TTL5(i+1,1)-data.corrected.TTL5(i,1))<60 & ...
                    ~any(data.corrected.TTL2(data.corrected.TTL2(:,1)>data.corrected.TTL5(i,1),1)-data.corrected.TTL5(i,1) < 60)
                    y = [y i];
                end
            end
        end
        
        if ports(k)
            DATA.PSTH_Flavor(k,:) = mean(data.corrected.PSTH_front465);%-mean(data.corrected.PSTH_front405);
            DATA.PSTH_Water(k,:) = mean(data.corrected.PSTH_rear465);%-mean(data.corrected.PSTH_rear405);
            if j == 2
                PSTH_FlavorFlavor(k,:) = mean(data.corrected.PSTH_front465(y,:));
                PSTH_WaterWater(k,:) = mean(data.corrected.PSTH_rear465(x,:));
                Flavor_NovelDay(k,:) = mean(data.corrected.PSTH_front465(:,1501:2500),2);
                Water_NovelDay(k,:) = mean(data.corrected.PSTH_rear465(:,1501:2500),2);
            end
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
            if j == 2
                Next_FlavorFlavor = [Next_FlavorFlavor data.corrected.NextReward_water(x)];
                Next_WaterWater = [Next_WaterWater data.corrected.NextReward_flavor(y)];
            end
        else
            DATA.PSTH_Flavor(k,:) = mean(data.corrected.PSTH_rear465);%-mean(data.corrected.PSTH_rear405);
            DATA.PSTH_Water(k,:) = mean(data.corrected.PSTH_front465);%-mean(data.corrected.PSTH_front405);
            if j == 2
                PSTH_FlavorFlavor(k,:) = mean(data.corrected.PSTH_rear465(x,:));
                PSTH_WaterWater(k,:) = mean(data.corrected.PSTH_front465(y,:));
                Flavor_NovelDay(k,:) = mean(data.corrected.PSTH_rear465(:,1501:2500),2);
                Water_NovelDay(k,:) = mean(data.corrected.PSTH_front465(:,1501:2500),2);
            end
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
            if j == 2
                Next_FlavorFlavor = [Next_FlavorFlavor data.corrected.NextReward_water(y)];
                Next_WaterWater = [Next_WaterWater data.corrected.NextReward_flavor(x)];
            end
        end
        
        DATA.NextReward_water = [DATA.NextReward_water data.corrected.NextReward_water];
        DATA.NextReward_flavor = [DATA.NextReward_flavor data.corrected.NextReward_flavor];
    end
    
    Flavor(j,:) = mean(DATA.PSTH_Flavor(:,1501:2500),2);
    Water(j,:) = mean(DATA.PSTH_Water(:,1501:2500),2);
    
    subplot(3,length(sessions),j)
    hold on
    axis square
    plot([0 0],[-5 15],'k--','LineWidth',1)
    DATA.Full_Times = linspace(-10,60,size(DATA.PSTH_Water,2));
    fill([DATA.Full_Times fliplr(DATA.Full_Times)],[mean(DATA.PSTH_Water)+std(DATA.PSTH_Water)./sqrt(N) fliplr(mean(DATA.PSTH_Water)-std(DATA.PSTH_Water)./sqrt(N))],[216 231 243]/255,'LineStyle','none')
    X=plot(DATA.Full_Times,mean(DATA.PSTH_Water),'Color',[55 136 193]/255,'LineWidth',1);
    fill([DATA.Full_Times fliplr(DATA.Full_Times)],[mean(DATA.PSTH_Flavor)+std(DATA.PSTH_Flavor)./sqrt(N) fliplr(mean(DATA.PSTH_Flavor)-std(DATA.PSTH_Flavor)./sqrt(N))],[252 216 213]/255,'LineStyle','none')
    Y=plot(DATA.Full_Times,mean(DATA.PSTH_Flavor),'Color',[229 45 38]/255,'LineWidth',1);
    pl = DATA.NextReward_flavor;
    pl2 = prctile(pl,[25 50 75]);
    scatter(pl,[lims(2)-.1*diff(lims)]*ones(size(pl))+.4*rand(size(pl))-.2,16,'MarkerFaceColor',[252 216 213]/255,'MarkerEdgeColor','w','LineWidth',1)
    plot([pl2(1) pl2(3)],[lims(2)-.1*diff(lims) lims(2)-.1*diff(lims)],'Color',[229 45 38]/255,'LineWidth',1)
    plot([pl2(2) pl2(2)],[lims(2)-.085*diff(lims) lims(2)-.115*diff(lims)],'Color',[229 45 38]/255,'LineWidth',1)
    pl = DATA.NextReward_water;
    pl2 = prctile(pl,[25 50 75]);
    scatter(pl,[lims(2)-.20*diff(lims)]*ones(size(pl))+.4*rand(size(pl))-.2,16,'MarkerFaceColor',[216 231 243]/255,'MarkerEdgeColor','w','LineWidth',1)
    plot([pl2(1) pl2(3)],[lims(2)-.20*diff(lims) lims(2)-.20*diff(lims)],'Color',[55 136 193]/255,'LineWidth',1)
    plot([pl2(2) pl2(2)],[lims(2)-.185*diff(lims) lims(2)-.215*diff(lims)],'Color',[55 136 193]/255,'LineWidth',1)
    %text(25,lims(2)-.035*diff(lims),'Next reward','FontSize',16,'HorizontalAlignment','Center')
    xlabel('Time (s) from reward')
    ylabel('PKA activity (σ)')
    xlim([-10 30])
    xticks(-10:10:30)
    yticks(lims(1):lims(2))
    ylim(lims)
    set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
    title([sessions{j}(1:4),'-',sessions{j}(5:6),'-',sessions{j}(7:8),expts{j}])
    hold off
    
    subplot(3,length(sessions),j+length(sessions)*2)
    hold on
    axis square
    heatmap(DATA.PSTH_Water(mouse_order,:),[],[],[],'Colormap',cmap,'ColorLevels',1000,'MaxColorValue',max(abs(lims)),'MinColorValue',-max(abs(lims)),'NaNColor',[1 1 1]);
    plot([find(DATA.Full_Times>=0,1),find(DATA.Full_Times>=0,1)]+0.5,[0 N]+0.5,'k--','LineWidth',1)
    %fill([0 size(DATA.PSTH_Water,2) size(DATA.PSTH_Water,2) 0]+0.5,[0 0 N N]+0.5,'k','FaceAlpha',0,'LineWidth',1)
    yticks([])
    %yticklabels(cellfun(@(x) x(end-5:end),mice,'UniformOutput',false))
    xticks([find(DATA.Full_Times>=-10,1) find(DATA.Full_Times>=0,1) find(DATA.Full_Times>=10,1) find(DATA.Full_Times>=20,1)  find(DATA.Full_Times>=30,1)]+0.5)
    xticklabels({'-10','0','10','20','30'})
    xlim([.5 4.0015e+03])
    xlabel('Time (s) from reward')
    ylabel('Individual mice')
    title('Water port')
    set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
    h = colorbar('FontSize',16,'LineWidth',1);
    ylabel(h,'PKA activity (σ)');
    h.Limits = lims;
    h.Ticks = [lims(1) 0 lims(2)];
    h.TickLength = 0;
    h.Box = 'on';
    hold off
    
    subplot(3,length(sessions),j+length(sessions))
    hold on
    axis square
    heatmap(DATA.PSTH_Flavor(mouse_order,:),[],[],[],'Colormap',cmap,'ColorLevels',1000,'MaxColorValue',max(abs(lims)),'MinColorValue',-max(abs(lims)),'NaNColor',[1 1 1]);
    plot([find(DATA.Full_Times>=0,1),find(DATA.Full_Times>=0,1)]+0.5,[0 N]+0.5,'k--','LineWidth',1)
    %fill([0 size(DATA.PSTH_Flavor,2) size(DATA.PSTH_Flavor,2) 0]+0.5,[0 0 N N]+0.5,'k','FaceAlpha',0,'LineWidth',1)
    xlabel('Time (s) from reward')
    ylabel('Individual mice')
    title('Flavor port')
    yticks([])
    %yticklabels(cellfun(@(x) x(end-5:end),mice,'UniformOutput',false))
    xticks([find(DATA.Full_Times>=-10,1) find(DATA.Full_Times>=0,1) find(DATA.Full_Times>=10,1) find(DATA.Full_Times>=20,1)  find(DATA.Full_Times>=30,1)]+0.5)
    xticklabels({'-10','0','10','20','30'})
    xlim([.5 4.0015e+03])
    set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
    h = colorbar('FontSize',16,'LineWidth',1);
    ylabel(h,'PKA activity (σ)');
    h.Limits = lims;
    h.Ticks = [lims(1) 0 lims(2)];
    h.TickLength = 0;
    h.Box = 'on';
    hold off
end
saveas(gcf,['figures/summary1'],'png')
set(gcf,'renderer','Painters')
saveas(gcf,['figures/summary1'],'epsc')

figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 1]);

save('summary.mat','Water','Flavor')
subplot(1,3,1)
hold on
axis square
plot(0:(length(sessions)-1),Water','Color',[216 231 243]/255,'LineWidth',1)
plot(0:(length(sessions)-1),Flavor','Color',[252 216 213]/255,'LineWidth',1)
A=plot(0:(length(sessions)-1),mean(Water,2),'Color',[55 136 193]/255,'LineWidth',1);
errorbar(0:(length(sessions))-1,mean(Water,2),std(Water,0,2)./sqrt(length(Water)),'Color',[55 136 193]/255,'CapSize',0,'LineWidth',1,'LineStyle','none')
B=plot(0:length(sessions)-1,mean(Flavor,2),'Color',[229 45 38]/255,'LineWidth',1);
errorbar(0:length(sessions)-1,mean(Flavor,2),std(Flavor,0,2)./sqrt(length(Flavor)),'Color',[229 45 38]/255,'CapSize',0,'LineWidth',1,'LineStyle','none')
xticks(0:(length(sessions))-1)
ylim([-2 8])
yticks(-2:2:8)
xlim([-.25 length(sessions)-.75])
xticklabels({'Water','Novel','Familiar','Familiar','Familiar'})
legend([B,A],{['Flavor port'],['Water port']},'Location','NorthEast','box','on')
xlabel('Day')
ylabel('PKA activity (σ)')

GLMM = readmatrix('PKAphotometry_GLMM.csv');
p = multicmp(GLMM(6,2:5),'up',0.05);
text(0.05,0.95,['Water: ',num2str(p(1),2),char(10),...
                'Novel: ',num2str(p(2),2),char(10),...
                'Familiar 1: ',num2str(p(3),2),char(10),...
                'Familiar 2: ',num2str(p(4),2)],'Units','Normalized','FontSize',12,'VerticalAlignment','top')

title('Response across days')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')

subplot(1,3,2)
hold on
axis square
plot(1:30,Water_NovelDay','Color',[216 231 243]/255,'LineWidth',1)
plot(1:30,Flavor_NovelDay','Color',[252 216 213]/255,'LineWidth',1)
A=plot(1:30,mean(Water_NovelDay',2),'Color',[55 136 193]/255,'LineWidth',1);
errorbar(1:30,mean(Water_NovelDay',2),std(Water_NovelDay',0,2)./sqrt(length(Water_NovelDay')),'Color',[55 136 193]/255,'CapSize',0,'LineWidth',1,'LineStyle','none')
B=plot(1:30,mean(Flavor_NovelDay',2),'Color',[229 45 38]/255,'LineWidth',1);
errorbar(1:30,mean(Flavor_NovelDay',2),std(Flavor_NovelDay',0,2)./sqrt(length(Flavor_NovelDay')),'Color',[229 45 38]/255,'CapSize',0,'LineWidth',1,'LineStyle','none')
xticks(0:5:30)
ylim([-30 30])
yticks(-30:10:30)
xlim([0 31])
%xticklabels({'Water','Novel','Familiar d1','Familiar d2','Familiar d3'})
legend([B,A],{['Novel flavor'],['Water']},'Location','Southwest','box','on')
xlabel('Trial')
ylabel('PKA activity (σ)')
title('Response across trials (Novel day)')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')


lims = [-3 6];
subplot(1,3,3)
hold on
axis square
plot([0 0],[-5 15],'k--','LineWidth',1)
DATA.Full_Times = linspace(-10,60,size(DATA.PSTH_Water,2));
fill([DATA.Full_Times fliplr(DATA.Full_Times)],[mean(PSTH_WaterWater)+std(PSTH_WaterWater)./sqrt(N) fliplr(mean(PSTH_WaterWater)-std(PSTH_WaterWater)./sqrt(N))],[216 231 243]/255,'LineStyle','none')
X=plot(DATA.Full_Times,mean(PSTH_WaterWater),'Color',[55 136 193]/255,'LineWidth',1);
fill([DATA.Full_Times fliplr(DATA.Full_Times)],[mean(PSTH_FlavorFlavor)+std(PSTH_FlavorFlavor)./sqrt(N) fliplr(mean(PSTH_FlavorFlavor)-std(PSTH_FlavorFlavor)./sqrt(N))],[252 216 213]/255,'LineStyle','none')
Y=plot(DATA.Full_Times,mean(PSTH_FlavorFlavor),'Color',[229 45 38]/255,'LineWidth',1);
pl = Next_FlavorFlavor;
pl2 = prctile(pl,[25 50 75]);
scatter(pl,[lims(2)-.05*diff(lims)]*ones(size(pl))+.4*rand(size(pl))-.2,16,'MarkerFaceColor',[252 216 213]/255,'MarkerEdgeColor','w','LineWidth',1)
plot([pl2(1) pl2(3)],[lims(2)-.05*diff(lims) lims(2)-.05*diff(lims)],'Color',[229 45 38]/255,'LineWidth',1)
plot([pl2(2) pl2(2)],[lims(2)-.035*diff(lims) lims(2)-.065*diff(lims)],'Color',[229 45 38]/255,'LineWidth',1)
pl = Next_WaterWater;
pl2 = prctile(pl,[25 50 75]);
scatter(pl,[lims(2)-.15*diff(lims)]*ones(size(pl))+.4*rand(size(pl))-.2,16,'MarkerFaceColor',[216 231 243]/255,'MarkerEdgeColor','w','LineWidth',1)
plot([pl2(1) pl2(3)],[lims(2)-.15*diff(lims) lims(2)-.15*diff(lims)],'Color',[55 136 193]/255,'LineWidth',1)
plot([pl2(2) pl2(2)],[lims(2)-.135*diff(lims) lims(2)-.165*diff(lims)],'Color',[55 136 193]/255,'LineWidth',1)
xlabel('Time (s) from reward')
ylabel('PKA activity (σ)')
xlim([-10 60])
yticks(lims(1):lims(2))
ylim(lims)
legend([Y,X],{['Novel flavor ×2'],['Water ×2']},'Location','Southwest','box','on')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
title('Repeated reward types (Novel day)')
hold off

saveas(gcf,['figures/summary2'],'png')
set(gcf,'renderer','Painters')
saveas(gcf,['figures/summary2'],'epsc')

figure
set(gcf,'Units','Normalized','OuterPosition',[0 0 1 1]);

subplot(2,2,1)
load('Z:\Chris\matlab\cz\photometry\data\20240614_pka701.mat')
port = 0;
idx = min([data.corrected.TTL2(1,1),data.corrected.TTL5(1,1)]);
idx = [find(data.corrected.time>(idx-60*5),1) find(data.corrected.time==max(data.corrected.time),1)];
x = (data.corrected.time(idx(1):idx(2))-data.corrected.time(idx(1)))./60-5;
y = data.corrected.AKAR465(idx(1):idx(2))./data.corrected.AKAR405(idx(1):idx(2));
hold on
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
ylabel('PKA activity (AU)')
ylim([2.35 2.75])
xlim([-5 20])
xticks(-5:5:20)
title('PKA701')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

subplot(2,2,2)
load('Z:\Chris\matlab\cz\photometry\data\20240614_pka713.mat')
port = 1;
idx = min([data.corrected.TTL2(1,1),data.corrected.TTL5(1,1)]);
idx = [find(data.corrected.time>(idx-60*5),1) find(data.corrected.time==max(data.corrected.time),1)];
x = (data.corrected.time(idx(1):idx(2))-data.corrected.time(idx(1)))./60-5;
y = data.corrected.AKAR465(idx(1):idx(2))./data.corrected.AKAR405(idx(1):idx(2));
hold on
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
ylabel('PKA activity (AU)')
ylim([2.15 2.8])
xlim([-5 20])
xticks(-5:5:20)
title('PKA713')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

subplot(2,2,3)
load('Z:\Chris\matlab\cz\photometry\data\20240614_pka705.mat')
port = 1;
idx = min([data.corrected.TTL2(1,1),data.corrected.TTL5(1,1)]);
idx = [find(data.corrected.time>(idx-60*5),1) find(data.corrected.time==max(data.corrected.time),1)];
x = (data.corrected.time(idx(1):idx(2))-data.corrected.time(idx(1)))./60-5;
y = data.corrected.AKAR465(idx(1):idx(2))./data.corrected.AKAR405(idx(1):idx(2));
hold on
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
ylabel('PKA activity (AU)')
ylim([2.45 3])
xlim([-5 20])
xticks(-5:5:20)
title('PKA705')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

subplot(2,2,4)
load('Z:\Chris\matlab\cz\photometry\data\20240614_pka703.mat')
port = 1;
idx = min([data.corrected.TTL2(1,1),data.corrected.TTL5(1,1)]);
idx = [find(data.corrected.time>(idx-60*5),1) find(data.corrected.time==max(data.corrected.time),1)];
x = (data.corrected.time(idx(1):idx(2))-data.corrected.time(idx(1)))./60-5;
y = data.corrected.AKAR465(idx(1):idx(2))./data.corrected.AKAR405(idx(1):idx(2));
hold on
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
ylabel('PKA activity (AU)')
ylim([2.3 3.15])
xlim([-5 20])
xticks(-5:5:20)
title('PKA703')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

saveas(gcf,['figures/examples'],'png')
set(gcf,'renderer','Painters')
saveas(gcf,['figures/examples'],'epsc')
end