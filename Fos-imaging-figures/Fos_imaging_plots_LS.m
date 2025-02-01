%% Setup
addpath(genpath('Z:\Chris\matlab\violin-plot\'));
addpath(genpath('Z:\Chris\matlab\heatmap\'));
addpath(genpath('Z:\Chris\matlab\numpy-matlab\'));
clear all; close all; clc
path.main = 'Z:\Chris\data\clearmap2\';
path.file = '\rawdata\resolution_3.6x\region_summary_statistics_classified_350size_0intensity.csv';
library = 'no_cerebellum_NEW';
load(['regions_allen_',library,'.mat'],'RegionLibrary')
load('significant_regions.mat','regions')
%% Load
fname  = cell(0,0);
flavor = cell(0,0);
phase  = cell(0,0);
batch  = cell(0,0);
sex    = cell(0,0);
code   = cell(0,0);

%%%%% BATCH 0 %%%%%
fname{end+1} = 'zimmerman_26\zimmerman_26-579\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-YFP'; batch{end+1} = 'Batch08'; sex{end+1} = 'Male';   code{end+1} = '08-099';
fname{end+1} = 'zimmerman_26\zimmerman_26-598\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-YFP'; batch{end+1} = 'Batch08'; sex{end+1} = 'Male';   code{end+1} = '08-100';
fname{end+1} = 'zimmerman_26\zimmerman_26-599\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-YFP'; batch{end+1} = 'Batch08'; sex{end+1} = 'Male';   code{end+1} = '08-101';
fname{end+1} = 'zimmerman_26\zimmerman_26-224\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-YFP'; batch{end+1} = 'Batch08'; sex{end+1} = 'Female'; code{end+1} = '08-102';
fname{end+1} = 'zimmerman_26\zimmerman_26-578\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-YFP'; batch{end+1} = 'Batch08'; sex{end+1} = 'Female'; code{end+1} = '08-103';
fname{end+1} = 'zimmerman_26\zimmerman_26-580\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-YFP'; batch{end+1} = 'Batch08'; sex{end+1} = 'Female'; code{end+1} = '08-104';

fname{end+1} = 'zimmerman_27\zimmerman_27-220\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-hM3D'; batch{end+1} = 'Batch08'; sex{end+1} = 'Male';   code{end+1} = '08-105';
fname{end+1} = 'zimmerman_27\zimmerman_27-221\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-hM3D'; batch{end+1} = 'Batch08'; sex{end+1} = 'Male';   code{end+1} = '08-106';
fname{end+1} = 'zimmerman_27\zimmerman_27-581\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-hM3D'; batch{end+1} = 'Batch08'; sex{end+1} = 'Male';   code{end+1} = '08-107';
fname{end+1} = 'zimmerman_27\zimmerman_27-595\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-hM3D'; batch{end+1} = 'Batch08'; sex{end+1} = 'Female'; code{end+1} = '08-108';
fname{end+1} = 'zimmerman_27\zimmerman_27-596\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-hM3D'; batch{end+1} = 'Batch08'; sex{end+1} = 'Female'; code{end+1} = '08-109';
fname{end+1} = 'zimmerman_27\zimmerman_27-597\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-hM3D'; batch{end+1} = 'Batch08'; sex{end+1} = 'Female'; code{end+1} = '08-110';

%%%%% BATCH 1 %%%%%
fname{end+1} = 'zimmerman_28\zimmerman_28-276\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-YFP'; batch{end+1} = 'Batch09'; sex{end+1} = 'Male';   code{end+1} = '09-111';
fname{end+1} = 'zimmerman_28\zimmerman_28-277\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-YFP'; batch{end+1} = 'Batch09'; sex{end+1} = 'Male';   code{end+1} = '09-112';
fname{end+1} = 'zimmerman_28\zimmerman_28-297\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-YFP'; batch{end+1} = 'Batch09'; sex{end+1} = 'Male';   code{end+1} = '09-113';
fname{end+1} = 'zimmerman_28\zimmerman_28-298\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-YFP'; batch{end+1} = 'Batch09'; sex{end+1} = 'Female'; code{end+1} = '09-114';
fname{end+1} = 'zimmerman_28\zimmerman_28-299\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-YFP'; batch{end+1} = 'Batch09'; sex{end+1} = 'Female'; code{end+1} = '09-115';
fname{end+1} = 'zimmerman_28\zimmerman_28-300\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-YFP'; batch{end+1} = 'Batch09'; sex{end+1} = 'Female'; code{end+1} = '09-116';

fname{end+1} = 'zimmerman_29\zimmerman_29-295\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-hM3D'; batch{end+1} = 'Batch09'; sex{end+1} = 'Male';   code{end+1} = '09-117';
fname{end+1} = 'zimmerman_29\zimmerman_29-296\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-hM3D'; batch{end+1} = 'Batch09'; sex{end+1} = 'Male';   code{end+1} = '09-118';
fname{end+1} = 'zimmerman_29\zimmerman_29-235\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-hM3D'; batch{end+1} = 'Batch09'; sex{end+1} = 'Female'; code{end+1} = '09-119';
fname{end+1} = 'zimmerman_29\zimmerman_29-236\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-hM3D'; batch{end+1} = 'Batch09'; sex{end+1} = 'Female'; code{end+1} = '09-120';
fname{end+1} = 'zimmerman_29\zimmerman_29-278\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-hM3D'; batch{end+1} = 'Batch09'; sex{end+1} = 'Female'; code{end+1} = '09-121';
fname{end+1} = 'zimmerman_29\zimmerman_29-294\imaging_request_1'; flavor{end+1} = 'Novel'; phase{end+1} = 'Malaise + LS-hM3D'; batch{end+1} = 'Batch09'; sex{end+1} = 'Female'; code{end+1} = '09-122';

warning('off')
for i = 1:length(fname)
    try
        raw_data{i} = readtable([path.main,fname{i},path.file]);
    catch
        disp(['Missing file: ',path.main,fname{i},path.file]);
    end
end
warning('on')
counts = [];
for j = 1:length(fname)
    counts(j,:) = raw_data{j}{1,:};
    %counts(j,:) = raw_data{j}{2,:};
    %disp([size(raw_data{i},2) sum(strcmp('CentralAmygdalarNucleus_VentralPart',raw_data{i}.Properties.VariableNames))])
end

[~,idx] = sort(fname);
A = cellfun(@(x) x(1:12),fname,'UniformOutput',false)';
B = cellfun(@(x) x(27:29),fname,'UniformOutput',false)';
C = cellfun(@(x) x(end),fname,'UniformOutput',false)';
D1 = counts(:,find(cellfun(@(x) isequal("CentralAmygdalarNucleus",x),RegionLibrary.full{:,2})));
D2 = counts(:,find(cellfun(@(x) isequal("LateralSeptalComplex",x),RegionLibrary.full{:,2})));
D3 = counts(:,find(cellfun(@(x) isequal("MedialSeptalComplex",x),RegionLibrary.full{:,2})));
D4 = counts(:,find(cellfun(@(x) isequal("IntercalatedAmygdalarNucleus",x),RegionLibrary.full{:,2})));
D5 = counts(:,find(cellfun(@(x) isequal("MedialAmygdalarNucleus",x),RegionLibrary.full{:,2})));
E = counts(:,1) - counts(:,find(cellfun(@(x) isequal("Cerebellum",x),RegionLibrary.full{:,2})));
tbl = table([A(idx)],[B(idx)],[C(idx)],[flavor(idx)]',[phase(idx)]',[sex(idx)]',[batch(idx)]',D1,D2,D3,E,D4,D5);
tbl.Properties.VariableNames = {'Cohort','Sample','Imaging session','Flavor','Phase','Sex','Batch code','CEA Fos','LS Fos','MS Fos','Total Fos','IA Fos','MEA Fos'};
writetable(tbl,'sample_info_LS_hM3D.csv')

fpath = 'Z:\Chris\data\clearmap2\kde-visualization\cta-final\';

fname = 'kde_20um_cta_licl_LS_hM3D_weighted_decr75pc.npy';
KDE.hM3D = readNPY([fpath,fname]);
fname = 'kde_20um_cta_licl_LS_YFP_weighted_decr75pc.npy';
KDE.YFP = readNPY([fpath,fname]);

fpath = 'Z:\Chris\data\clearmap2\utilities\allen-atlas-cz\';
fname = 'annotation_2017_25um_sagittal_16bit_hierarch_labels_fillmissing_cz_v2.tif';
filename = [fpath,fname];
tstack  = Tiff(filename);
[I,J] = size(tstack.read());
K = length(imfinfo(filename));
data = zeros(K,I,J);
data(1,:,:)  = tstack.read();
for n = 2:K
    tstack.nextDirectory()
    data(n,:,:) = tstack.read();
end
atlas = data;
mask = ismember(atlas,1115:1305); atlas(mask) = 1115;
mask = ismember(atlas,1306:1317); atlas(mask) = 1306;
mask = ismember(atlas,1028:1114); atlas(mask) = 1;

atlasmask  = atlas>=1028 | atlas<=1;
KDE.Stim_Resilient(atlasmask) = NaN;
KDE.Stim_Susceptible(atlasmask) = NaN;
KDE.Stim(atlasmask) = NaN;
KDE.NoStim(atlasmask) = NaN;

method = 'diff';
lims1 = [0 1];
if isequal(method,'log')
    lims2 = [-2 2];
    lambda = 1E-2;
elseif isequal(method,'diff')
    lims2 = [-.5 .5];
else
    disp('Invalid comparison method.')
    return
end
%% Plots

YFP = find(cellfun(@(x) isequal(x,'YFP'),tbl.Phase));
hM3D = find(cellfun(@(x) isequal(x,'hM3D'),tbl.Phase));
cmap2 = flipud(cbrewer('div','RdGy',1000,'spline')); cmap2(cmap2<0) = 0; cmap2(cmap2>1) = 1;

%close all
figure('Position', get(0, 'Screensize'))

subplot(2,4,1)
plane = 180;
brainoutline = atlas>0;
B = bwboundaries(fliplr(squeeze(brainoutline(:,plane,:))));
atlasmaskplot2 = atlas==1306;
B2a = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
atlasmaskplot2 = atlas==1115;
B2b = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
atlasmaskplot2 = atlas==1;
B2d = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
hold on
axis equal
data = squeeze(KDE.hM3D(:,plane,:))-squeeze(KDE.YFP(:,plane,:));
heatmap(rot90(data),[],[],[],'Colormap',cmap2,'ColorLevels',1000,'MaxColorValue',2,'MinColorValue',-2,'NaNColor',[1 1 1]);
for i = 1:size(RegionLibrary.reduced,1)
    idx = RegionLibrary.reduced{i,1}+1;
    mask = squeeze(atlas(:,plane,:)) == idx;
    B3 = bwboundaries(fliplr(mask));
    if ismember(i,regions.amygdala)
        clr = [229 45 38]/255;
    elseif  i ==76
        clr = [55 136 193]/255;
    else
        clr = [.85 .85 .85];
    end
    for j = 1:length(B3)
        plot(B3{j}(:,1),B3{j}(:,2),'Color',clr,'LineWidth',1)
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
set(gca,'FontSize',16)
axis off
c2 = colorbar('Location','westoutside','FontSize',16);
c2.Position = c2.Position+[-.01 0 0 0];
c2.Ticks = [-2 0 2];
c2.Label.String = '\DeltaFos (% per mm^3)';
xlim([.5 228.5])
hold off

subplot(2,4,2)
plane = 265;
brainoutline = atlas>0;
B = bwboundaries(fliplr(squeeze(brainoutline(:,plane,:))));
atlasmaskplot2 = atlas==1306;
B2a = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
atlasmaskplot2 = atlas==1115;
B2b = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
atlasmaskplot2 = atlas==1;
B2d = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
hold on
axis equal
data = squeeze(KDE.hM3D(:,plane,:))-squeeze(KDE.YFP(:,plane,:));
heatmap(rot90(data),[],[],[],'Colormap',cmap2,'ColorLevels',1000,'MaxColorValue',.5,'MinColorValue',-.5,'NaNColor',[1 1 1]);
for i = 1:size(RegionLibrary.reduced,1)
    idx = RegionLibrary.reduced{i,1}+1;
    mask = squeeze(atlas(:,plane,:)) == idx;
    B3 = bwboundaries(fliplr(mask));
    if ismember(i,regions.amygdala)
        clr = [229 45 38]/255;
    elseif  i ==76
        clr = [55 136 193]/255;
    else
        clr = [.85 .85 .85];
    end
    for j = 1:length(B3)
        plot(B3{j}(:,1),B3{j}(:,2),'Color',clr,'LineWidth',1)
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
set(gca,'FontSize',16)
c2 = colorbar('Location','westoutside','FontSize',16);
c2.Position = c2.Position+[-.01 0 0 0];
c2.Ticks = [-.5 0 .5];
c2.Label.String = '\DeltaFos (% per mm^3)';
xlim([.5 228.5])
axis off
hold off

subplot(2,4,3)
plane = 325;
brainoutline = atlas>0;
B = bwboundaries(fliplr(squeeze(brainoutline(:,plane,:))));
atlasmaskplot2 = atlas==1306;
B2a = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
atlasmaskplot2 = atlas==1115;
B2b = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
atlasmaskplot2 = atlas==1;
B2d = bwboundaries(fliplr(squeeze(atlasmaskplot2(:,plane,:))));
hold on
axis equal
data = squeeze(KDE.hM3D(:,plane,:))-squeeze(KDE.YFP(:,plane,:));
heatmap(rot90(data),[],[],[],'Colormap',cmap2,'ColorLevels',1000,'MaxColorValue',.5,'MinColorValue',-.5,'NaNColor',[1 1 1]);
for i = 1:size(RegionLibrary.reduced,1)
    idx = RegionLibrary.reduced{i,1}+1;
    mask = squeeze(atlas(:,plane,:)) == idx;
    B3 = bwboundaries(fliplr(mask));
    if ismember(i,regions.amygdala)
        clr = [229 45 38]/255;
    elseif  i ==76
        clr = [55 136 193]/255;
    else
        clr = [.85 .85 .85];
    end
    for j = 1:length(B3)
        plot(B3{j}(:,1),B3{j}(:,2),'Color',clr,'LineWidth',1)
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
set(gca,'FontSize',16)
c2 = colorbar('Location','westoutside','FontSize',16);
c2.Position = c2.Position+[-.01 0 0 0];
c2.Ticks = [-.5 0 .5];
c2.Label.String = '\DeltaFos (% per mm^3)';
xlim([.5 228.5])
axis off
hold off

subplot(2,3,4)
hold on
axis square
data = counts(:,RegionLibrary.reduced{:,1}+1)./E*100./RegionLibrary.reduced{:,5}';
plot([0 4],[0 4],'k--','LineWidth',1)
regions.septum = find(contains(RegionLibrary.reduced.region,'ept'));
idx0 = setdiff(regions.significant,[regions.amygdala;regions.septum]);

A = mean(data(hM3D,idx0));
B = mean(data(YFP,idx0));
a = scatter(A,B,1,'k','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A,B,1);
x = 0:.01:4;
delete(a)
[r,p]=corr(A',B');
[y_fit] = polyval(p2,x,S);
fitresult = fit(A',B','poly1');
p2 = predint(fitresult,x,0.95,'functional');
p2(p2<0) = 0; p2(p2>4) = 4;
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[.85 .85 .85],'LineStyle','none'),
plot(x,y_fit,'Color','k','LineWidth',1)

A = mean(data(hM3D,regions.amygdala));
B = mean(data(YFP,regions.amygdala));
a = scatter(A,B,1,'r','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A,B,1);
x = 0:.01:4;
delete(a)
[r,p]=corr(A',B');
[y_fit] = polyval(p2,x,S);
fitresult = fit(A',B','poly1');
p2 = predint(fitresult,x,0.95,'functional');
p2(p2<0) = 0; p2(p2>4) = 4;
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[252 216 213]/255,'LineStyle','none'),
plot(x,y_fit,'Color',[229 45 38]/255,'LineWidth',1)

A = mean(data(hM3D,regions.septum));
B = mean(data(YFP,regions.septum));
a = scatter(A,B,1,'r','filled','MarkerEdgeColor','w');
[p2,S] = polyfit(A,B,1);
x = 0:.01:4;
delete(a)
[r,p]=corr(A',B');
[y_fit] = polyval(p2,x,S);
fitresult = fit(A',B','poly1');
p2 = predint(fitresult,x,0.95,'functional');
p2(p2<0) = 0; p2(p2>4) = 4;
fill([x fliplr(x)],[p2(:,1)' flipud(p2(:,2))'],[216 231 243]/255,'LineStyle','none'),
plot(x,y_fit,'Color',[55 136 193]/255,'LineWidth',1)

a = scatter(mean(data(hM3D,idx0)),mean(data(YFP,idx0)),128,'k','filled','MarkerEdgeColor','w');
b = scatter(mean(data(hM3D,regions.amygdala)),mean(data(YFP,regions.amygdala)),128,[229 45 38]/255,'filled','MarkerEdgeColor','w');
c = scatter(mean(data(hM3D,regions.septum)),mean(data(YFP,regions.septum)),128,[55 136 193]/255,'filled','MarkerEdgeColor','w');
xlim([0 4]); xticks(0:1:4);
ylim([0 4]); yticks(0:1:4);
ylabel('YFP Fos^ (% per mm^3)')
xlabel('hM3D Fos^ (% per mm^3)')
legend([b,c,a],{'Amygdala net.','Septum','Other regions'},'location','northwest')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

a = [mean(data(hM3D,regions.amygdala)), mean(data(hM3D,regions.septum)), mean(data(hM3D,idx0))]';
b = [mean(data(YFP,regions.amygdala)), mean(data(YFP,regions.septum)), mean(data(YFP,idx0))]';
t = table(a,b,'VariableNames',{'hM3D','YFP'});
writetable(t,'Z:\Chris\matlab\cz\cta-source-data\Fig-ED2f.csv')

subplot(2,3,5)
hold on
axis square
data = tbl.('CEA Fos')./tbl.('Total Fos')*100/1.987;
bar(1,mean(data(YFP)),'k','linestyle','none')
bar(2,mean(data(hM3D)),'FaceColor',[229 45 38]/255,'linestyle','none')
errorbar(1,mean(data(YFP)),std(data(YFP))/sqrt(length(YFP)),'k','linewidth',2,'CapSize',50)
errorbar(2,mean(data(hM3D)),std(data(hM3D))/sqrt(length(hM3D)),'Color',[229 45 38]/255,'linewidth',2,'CapSize',50)
scatter(ones(size(YFP))+0.2*rand(size(YFP))-.1,data(YFP),128,[.85 .85 .85],'filled','MarkerEdgeColor','w','linewidth',1)
scatter(2*ones(size(hM3D))+0.2*rand(size(hM3D))-.1,data(hM3D),128,[252 216 213]/255,'filled','MarkerEdgeColor','w','linewidth',1)
xlim([.25 2.75])
xticks([1,2])
xticklabels({'YFP','hM3D'})
ylabel('CEA Fos (% per mm^3)')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

t = table(data,'VariableNames',{'CEA Fos'});
t = [tbl.Phase,t];
writetable(t,'Z:\Chris\matlab\cz\cta-source-data\Fig-ED2e.csv')

subplot(2,3,6)
hold on
axis square
data = tbl.('LS Fos')./tbl.('Total Fos')*100/3.553;
bar(1,mean(data(YFP)),'k','linestyle','none')
bar(2,mean(data(hM3D)),'FaceColor',[229 45 38]/255,'linestyle','none')
errorbar(1,mean(data(YFP)),std(data(YFP))/sqrt(length(YFP)),'k','linewidth',2,'CapSize',50)
errorbar(2,mean(data(hM3D)),std(data(hM3D))/sqrt(length(hM3D)),'Color',[229 45 38]/255,'linewidth',2,'CapSize',50)
scatter(ones(size(YFP))+0.2*rand(size(YFP))-.1,data(YFP),128,[.85 .85 .85],'filled','MarkerEdgeColor','w','linewidth',1)
scatter(2*ones(size(hM3D))+0.2*rand(size(hM3D))-.1,data(hM3D),128,[252 216 213]/255,'filled','MarkerEdgeColor','w','linewidth',1)
xlim([.25 2.75])
xticks([1,2])
xticklabels({'YFP','hM3D'})
ylabel('LS Fos (% per mm^3)')
set(gca,'FontSize',16,'LineWidth',1,'TickLength',[0.025, 0],'TickDir','out')
hold off

t = table(data,'VariableNames',{'LS Fos'});
t = [tbl.Phase,t];
writetable(t,'Z:\Chris\matlab\cz\cta-source-data\Fig-ED2d.csv')

saveas(gcf,'LS-Fos','png')
set(gcf,'renderer','Painters')
saveas(gcf,'LS-Fos','epsc')
%% Stats
YFP = find(cellfun(@(x) isequal(x,'YFP'),tbl.Phase));
hM3D = find(cellfun(@(x) isequal(x,'hM3D'),tbl.Phase));
data = counts(:,RegionLibrary.reduced{:,1}+1)./E*100./RegionLibrary.reduced{:,5}';
yfp = mean(data(YFP,regions.significant));
hm3d = mean(data(hM3D,regions.significant));
[~,idx1]=ismember(regions.amygdala,regions.significant);
[~,idx2]=ismember(regions.septum,regions.significant);
C = zeros(length(regions.significant),1); C(idx1) = 1; C(idx2) = 2;
X = cell(0,0);
for i = 1:length(C)
    if C(i)==0
        X{i} = 'Other regions';
    elseif C(i) ==1
        X{i} = 'Amygdala net.';
    else
        X{i} = 'Septum';
    end
end
[h,atab,ctab,stats] = aoctool(hm3d,yfp,X,.05,'YFP','hM3D');
[c,m,h,gnames] = multcompare(stats,0.05,"on","bonferroni","slope");
[mn,se] = tValue(stats.slopes,stats.slopecov); % copied from multcompare
t = mn./se;
disp([c t])