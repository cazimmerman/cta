function fig_5(data_path)

disp('Generating panels for Figure 5...')
%% Fig 5d

figure('Position', get(0, 'Screensize'))
sgtitle('Figure 5d','FontWeight','bold')

lims = [-2 5];
for j = 1:4
    subplot(1,4,j)
    load([data_path,'\Photometry\pka-photometry-day',num2str(j-1),'.mat'],'data')
    N = length(data.raw.port);
    hold on
    axis square
    fill([data.psth.time fliplr(data.psth.time)],[mean(data.psth.PortB)+std(data.psth.PortB)./sqrt(N) fliplr(mean(data.psth.PortB)-std(data.psth.PortB)./sqrt(N))],[216 231 243]/255,'LineStyle','none')
    X=plot(data.psth.time,mean(data.psth.PortB),'Color',[55 136 193]/255,'LineWidth',1);
    fill([data.psth.time fliplr(data.psth.time)],[mean(data.psth.PortA)+std(data.psth.PortA)./sqrt(N) fliplr(mean(data.psth.PortA)-std(data.psth.PortA)./sqrt(N))],[252 216 213]/255,'LineStyle','none')
    Y=plot(data.psth.time,mean(data.psth.PortA),'Color',[229 45 38]/255,'LineWidth',1);
    pl = data.next_reward.PortA;
    pl2 = prctile(pl,[25 50 75]);
    scatter(pl,[lims(2)-.1*diff(lims)]*ones(size(pl))+.4*rand(size(pl))-.2,25,'MarkerFaceColor',[252 216 213]/255,'MarkerEdgeColor','w','LineWidth',1)
    plot([pl2(1) pl2(3)],[lims(2)-.1*diff(lims) lims(2)-.1*diff(lims)],'Color',[229 45 38]/255,'LineWidth',1)
    plot([pl2(2) pl2(2)],[lims(2)-.085*diff(lims) lims(2)-.115*diff(lims)],'Color',[229 45 38]/255,'LineWidth',1)
    pl = data.next_reward.PortB;
    pl2 = prctile(pl,[25 50 75]);
    scatter(pl,[lims(2)-.20*diff(lims)]*ones(size(pl))+.4*rand(size(pl))-.2,25,'MarkerFaceColor',[216 231 243]/255,'MarkerEdgeColor','w','LineWidth',1)
    plot([pl2(1) pl2(3)],[lims(2)-.20*diff(lims) lims(2)-.20*diff(lims)],'Color',[55 136 193]/255,'LineWidth',1)
    plot([pl2(2) pl2(2)],[lims(2)-.185*diff(lims) lims(2)-.215*diff(lims)],'Color',[55 136 193]/255,'LineWidth',1)
    xlabel('Time (s)')
    ylabel('PKA activity (σ)')
    xlim([-10 30])
    xticks(-10:10:30)
    yticks(lims(1):lims(2))
    ylim(lims)
    title(['Day ',num2str(j-1)])
    legend([Y,X],{'Port A','Port B'})
    set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
    hold off
end
%% Fig 5e

figure('Position', get(0, 'Screensize'))
sgtitle('Figure 5e','FontWeight','bold')

cmap = flipud(cbrewer('div','RdBu',1000,'spline')); cmap(cmap<0) = 0;
load([data_path,'\Photometry\pka-photometry-day1.mat'],'data')
[~,idx]=sort(mean(data.psth.PortA(:,1501:2500),2));
lims = [-2 5];
for j = 1:4
    load([data_path,'\Photometry\pka-photometry-day',num2str(j-1),'.mat'],'data')
    
    subplot(2,4,j)
    hold on
    axis square
    heatmap(data.psth.PortA(idx,:),[],[],[],'Colormap',cmap,'ColorLevels',1000,'MaxColorValue',max(abs(lims)),'MinColorValue',-max(abs(lims)),'NaNColor',[1 1 1]);
    yticks([])
    xticks([find(data.psth.time>=-10,1) find(data.psth.time>=0,1) find(data.psth.time>=10,1) find(data.psth.time>=20,1)  find(data.psth.time>=30,1)]+0.5)
    xticklabels({'-10','0','10','20','30'})
    xlim([.5 4.0015e+03])
    xlabel('Time (s)')
    ylabel('Port A')
    title(['Day ',num2str(j-1)])
    set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
    hold off
    
    subplot(2,4,j+4)
    hold on
    axis square
    heatmap(data.psth.PortB(idx,:),[],[],[],'Colormap',cmap,'ColorLevels',1000,'MaxColorValue',max(abs(lims)),'MinColorValue',-max(abs(lims)),'NaNColor',[1 1 1]);
    yticks([])
    xticks([find(data.psth.time>=-10,1) find(data.psth.time>=0,1) find(data.psth.time>=10,1) find(data.psth.time>=20,1)  find(data.psth.time>=30,1)]+0.5)
    xticklabels({'-10','0','10','20','30'})
    xlim([.5 4.0015e+03])
    xlabel('Time (s)')
    ylabel('Port B')
    title(['Day ',num2str(j-1)])
    set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
    hold off
end
%% Fig 5f

figure('Position', get(0, 'Screensize'))
sgtitle('Figure 5f','FontWeight','bold')

cmap = flipud(cbrewer('div','RdBu',1000,'spline')); cmap(cmap<0) = 0;
PortA = []; PortB = [];
for j = 1:4
    load([data_path,'\Photometry\pka-photometry-day',num2str(j-1),'.mat'],'data')
    PortA(j,:) = mean(data.psth.PortA(:,1501:2500),2);
    PortB(j,:) = mean(data.psth.PortB(:,1501:2500),2);
end

hold on
axis square
plot(0:3,PortB','Color',[216 231 243]/255,'LineWidth',1)
plot(0:3,PortA','Color',[252 216 213]/255,'LineWidth',1)
A=plot(0:3,mean(PortB,2),'Color',[55 136 193]/255,'LineWidth',1);
errorbar(0:3,mean(PortB,2),std(PortB,0,2)./sqrt(length(PortB)),'Color',[55 136 193]/255,'CapSize',0,'LineWidth',1,'LineStyle','none')
B=plot(0:3,mean(PortA,2),'Color',[229 45 38]/255,'LineWidth',1);
errorbar(0:3,mean(PortA,2),std(PortA,0,2)./sqrt(length(PortA)),'Color',[229 45 38]/255,'CapSize',0,'LineWidth',1,'LineStyle','none')
xticks(0:3)
ylim([-2 8])
yticks(-2:2:8)
xlim([-.5 3.5])
legend([B,A],{'Port A','Port B'})
xlabel('Day')
ylabel('PKA activity (σ)')
set(gca,'FontSize',12,'LineWidth',1,'TickLength',[0.015, 0],'TickDir','out')
hold off

T = readtable([data_path,'\GLMMs\PKA-photometry-GLMM-output.csv']);
p = multicmp(T{4,2:5},'up',0.05);
StatsTbl = table({'5f'},{'Day 0: Port A vs. Port B'},{'GLMM marginal effect'},{'4 days'},{[size(PortA,2)]},-T{3,2},p(1), ...
    'VariableNames',{'Figure panel','Group','Statistical test','Multiple comparisons','Sample size','Test statistic','P-value'});
StatsTbl(end+1,:) = table({'5f'},{'Day 1: Port A vs. Port B'},{'GLMM marginal effect'},{'4 days'},{[size(PortA,2)]},-T{3,3},p(2));
StatsTbl(end+1,:) = table({'5f'},{'Day 2: Port A vs. Port B'},{'GLMM marginal effect'},{'4 days'},{[size(PortA,2)]},-T{3,4},p(3));
StatsTbl(end+1,:) = table({'5f'},{'Day 3: Port A vs. Port B'},{'GLMM marginal effect'},{'4 days'},{[size(PortA,2)]},-T{3,5},p(4));
%% Fig 5g

figure('Position', get(0, 'Screensize'))
sgtitle('Figure 5g','FontWeight','bold')
idx_amygdala = [80 81 82];

load([data_path,'\Fos imaging\modified-atlas\allen_ccfv3_modified_cz.mat'],'atlas','RegionLibrary')
load([data_path,'\Photometry\pka-photometry-day1.mat'],'data')
PortA = mean(data.psth.PortA(:,1501:2500),2);
PortB = mean(data.psth.PortB(:,1501:2500),2);
data.raw.fiber_location(data.raw.fiber_location==0)=nan;
data.raw.fiber_location = round(data.raw.fiber_location);

cmap = flipud(cbrewer('div','RdBu',1000,'spline')); cmap(cmap<0) = 0;
cmap_amygdala = [192 0 72;...
    168 0 0;...
    192 72 0]/255;
brainoutline = atlas>0;

subplot(1,2,1);
plane = round(nanmean(reshape(squeeze(data.raw.fiber_location(:,1)),[numel(squeeze(data.raw.fiber_location(:,1))),1])));
B = bwboundaries(squeeze(brainoutline(:,plane,:)));
atlasmaskplot2 = atlas==1306;
B2a = bwboundaries(squeeze(atlasmaskplot2(:,plane,:)));
atlasmaskplot2 = atlas==1115;
B2b = bwboundaries(squeeze(atlasmaskplot2(:,plane,:)));
atlasmaskplot2 = atlas==1;
B2d = bwboundaries(squeeze(atlasmaskplot2(:,plane,:)));
hold on
axis square
axis equal
for i = 1:3
    idx = RegionLibrary.reduced{idx_amygdala(i),1};
    mask = squeeze(atlas(:,plane,:)) == idx;
    B3 = bwboundaries(mask);
    for j = 1:size(B3)
        fill(B3{j}(:,1),B3{j}(:,2),cmap_amygdala(i,:),'LineStyle','none','FaceAlpha',0.5)
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
plot([309 364 364 309 309],[200 200 255 255 200],'k--','LineWidth',1)
[~,idx]=sort(PortA);
for j = 1:length(idx)
    i = idx(j);
    C = round((PortA(i)+5)*100);
    if C >1000
        C = 1000;
    elseif C<1
        C = 1;
    end
    rectangle('Position',[data.raw.fiber_location(i,3)-2 data.raw.fiber_location(i,2) 4 4],'Curvature',[1 1],'linewidth',1,'FaceColor',cmap(C,:))
end
xlim([309 364])
ylim([200 255])
set(gca,'Ydir','reverse')
set(gca,'Xdir','reverse')
title('Novel flavor')
axis off
set(gca,'FontSize',12)
hold off

subplot(1,2,2);
plane = round(nanmean(reshape(squeeze(data.raw.fiber_location(:,1)),[numel(squeeze(data.raw.fiber_location(:,1))),1])));
B = bwboundaries(squeeze(brainoutline(:,plane,:)));
atlasmaskplot2 = atlas==1306;
B2a = bwboundaries(squeeze(atlasmaskplot2(:,plane,:)));
atlasmaskplot2 = atlas==1115;
B2b = bwboundaries(squeeze(atlasmaskplot2(:,plane,:)));
atlasmaskplot2 = atlas==1;
B2d = bwboundaries(squeeze(atlasmaskplot2(:,plane,:)));
hold on
axis square
axis equal
for i = 1:3
    idx = RegionLibrary.reduced{idx_amygdala(i),1};
    mask = squeeze(atlas(:,plane,:)) == idx;
    B3 = bwboundaries(mask);
    for j = 1:size(B3)
        fill(B3{j}(:,1),B3{j}(:,2),cmap_amygdala(i,:),'LineStyle','none','FaceAlpha',0.5)
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
plot([309 364 364 309 309],[200 200 255 255 200],'k--','LineWidth',1)
[~,idx]=sort(PortB);
for j = 1:length(idx)
    i = idx(j);
    C = round((PortB(i)+5)*100);
    if C >1000
        C = 1000;
    elseif C<1
        C = 1;
    end
    rectangle('Position',[data.raw.fiber_location(i,3)-2 data.raw.fiber_location(i,2) 4 4],'Curvature',[1 1],'linewidth',1,'FaceColor',cmap(C,:))
end
xlim([309 364])
ylim([200 255])
set(gca,'Ydir','reverse')
set(gca,'Xdir','reverse')
title('Water')
axis off
set(gca,'FontSize',12)
hold off

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