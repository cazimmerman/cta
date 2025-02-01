%%
clear all; close all; clc

tools.libary     = 'Z:\Chris\data\clearmap2\utilities\allen-atlas-cz\allen_CCFv3_ontology_cz\allen_CCFv3_library_cz.mat';
tools.ontology   = 'Z:\Chris\data\clearmap2\utilities\allen-atlas-cz\allen_CCFv3_ontology_cz\';

fpath.in         = 'Z:\Chris\data\clearmap2\utilities\allen-atlas\';
fpath.out	     = 'Z:\Chris\data\clearmap2\utilities\allen-atlas-cz\';

fname.scan.in    = 'average_template_25_sagittal_forDVscans.tif';
fname.scan.out   = 'average_template_25_sagittal_forDVscans_cz.tif';
fname.atlas.in   = 'annotation_2017_25um_sagittal_16bit_hierarch_labels_fillmissing.tif';
fname.atlas.out1 = 'annotation_2017_25um_sagittal_16bit_hierarch_labels_fillmissing_cz.tif';
fname.atlas.out2 = 'annotation_2017_25um_sagittal_16bit_hierarch_labels_fillmissing_cz_v2.tif';
fname.atlas.out2b = 'annotation_2017_25um_sagittal_16bit_hierarch_labels_fillmissing_cz_v2_Neuroglancer.tif';
fname.atlas.out3 = 'annotation_2017_25um_sagittal_16bit_hierarch_labels_fillmissing_cz_v2_kde.tif';
%% symmetrical scan

filename = [fpath.in,fname.scan.in];
tstack  = Tiff(filename);
[I,J] = size(tstack.read());
K = length(imfinfo(filename));
data = zeros(K,I,J);
data(1,:,:)  = tstack.read();
for n = 2:K
    tstack.nextDirectory()
    data(n,:,:) = tstack.read();
end
x1 = data(1:456/2,:,:);
x1 = [x1; x1(end:-1:1,:,:)];

filename_out = [fpath.out,fname.scan.out];
copyfile(filename,filename_out)
tstack_out = Tiff(filename_out,'r+');
K = length(imfinfo(filename_out));
tstack_out.write(uint16(squeeze(x1(1,:,:))))
for n = 2:K
    tstack_out.nextDirectory()
    tstack_out.write(uint16(squeeze(x1(n,:,:))))
end
%% symmetrical data

filename = [fpath.in,fname.atlas.in];
tstack  = Tiff(filename);
[I,J] = size(tstack.read());
K = length(imfinfo(filename));
data = zeros(K,I,J);
data(1,:,:)  = tstack.read();
for n = 2:K
    tstack.nextDirectory()
    data(n,:,:) = tstack.read();
end
x1 = data(1:456/2,:,:);
x1 = [x1; x1(end:-1:1,:,:)];

filename_out = [fpath.out,fname.atlas.out1];
copyfile(filename,filename_out)
tstack_out = Tiff(filename_out,'r+');
K = length(imfinfo(filename_out));
tstack_out.write(uint16(squeeze(x1(1,:,:))))
for n = 2:K
    tstack_out.nextDirectory()
    tstack_out.write(uint16(squeeze(x1(n,:,:))))
end
%% symmetrical data + region reassignments

filename = [fpath.out,fname.atlas.out1];
tstack  = Tiff(filename);
[I,J] = size(tstack.read());
K = length(imfinfo(filename));
data = zeros(K,I,J);
data(1,:,:)  = tstack.read();
for n = 2:K
    tstack.nextDirectory()
    data(n,:,:) = tstack.read();
end

load(tools.libary)
for i = 1:size(RegionLibrary.reduced,1)
    val = RegionLibrary.reduced{i,1};
    fname1 = [tools.ontology,'region_',num2str(val),'.p'];
    fid = py.open(fname1,'rb');
    fdata = int64(py.pickle.load(fid));
    if ~isempty(fdata)
        for j = 1:length(fdata)
            data(data==fdata(j)+1) = val+1;
        end
    end
end

target = [];
reassign = cell(0,0);
% Lateral geniculate complex
target(end+1) = 671;
reassign{end+1} = [714,716:718];
% Paraventricular hypothalamic nucelus
target(end+1) = 731;
reassign{end+1} = 793:797;
% Periventricular hypothalamic nucleus
target(end+1) = 742;
reassign{end+1} = [743,760,761];
% Midbrain reticular nucleus
target(end+1) = 838;
reassign{end+1} = 837;
% Pontine reticular nucleus
target(end+1) = 944;
reassign{end+1} = 926:927;
% Spinal nucleus of the trigeminal
target(end+1) = 967;
reassign{end+1} = 968:974;
% "Isocortex-un1/un2/un3/un4/un5/un6" ---> Isocortex
target(end+1) = 5;
reassign{end+1} = [381:386,565];
% "SSp-un/un2/un3" ---> Primary somatosensory cortex
target(end+1) = 37;
reassign{end+1} = 93:101;
% "HPF-un" ---> Hippocampal formation
target(end+1) = 462;
reassign{end+1} = 563;
% "Isocortex-un7/un8" ---> Thalamus
target(end+1) = 650;
reassign{end+1} = [724,725];
% "IB-un" ---> Hypothalamus
target(end+1) = 726;
reassign{end+1} = 817;
% "MB-un" ---> Midbrain (looks like pons at posterior end)
target(end+1) = 818;
reassign{end+1} = 894;
% "root-un1" ---> Hippocampal formation
target(end+1) = 462;
reassign{end+1} = 1340;
% "root-un2" ---> Medulla (looks like cerebellum at anterior end)
target(end+1) = 948;
reassign{end+1} = 1341;
% "root-un3" ---> Pons (looks like midbrain at anterior end)
target(end+1) = 896;
reassign{end+1} = 1342;
% "root-un4" ---> Hypothalamus
target(end+1) = 726;
reassign{end+1} = 1343;
% "root-un5" ---> Hypothalamus (looks like midbrain at posterior end)
target(end+1) = 726;
reassign{end+1} = 1344;
% "root-un6" ---> Midbrain
target(end+1) = 818;
reassign{end+1} = 1345;
% "root-un7" ---> Midbrain
target(end+1) = 818;
reassign{end+1} = 1346;

for i = 1:length(target)
    mask = ismember(data,reassign{i}+1);
    data(mask) = target(i)+1;
end

% remove anterior olfactory areas
mask = (data>=388 & data<=462) | data>=1115; mask(:,112:end,:) = 0; data(mask) = 1;
data = Add_CEAast(data);

filename_out = [fpath.out,fname.atlas.out2];
copyfile(filename,filename_out)
tstack_out = Tiff(filename_out,'r+');
K = length(imfinfo(filename_out));
tstack_out.write(uint16(squeeze(data(1,:,:))))
for n = 2:K
    tstack_out.nextDirectory()
    tstack_out.write(uint16(squeeze(data(n,:,:))))
end

%% symmetrical data + region reassignments + brain outline (for Neuroglancer viz)

filename = [fpath.out,fname.atlas.out2];
tstack  = Tiff(filename);
[I,J] = size(tstack.read());
K = length(imfinfo(filename));
data = zeros(K,I,J);
data(1,:,:)  = tstack.read();
for n = 2:K
    tstack.nextDirectory()
    data(n,:,:) = tstack.read();
end

mask = data==0;
ind = find(mask);
[idx1,idx2,idx3] = ind2sub(size(mask), ind);
mask = zeros(size(data));
mask2 = data>1;
for i = 1:length(idx1)
    try
    if sum(mask2([idx1(i)-1:idx1(i)+1],[idx2(i)-1:idx2(i)+1],[idx3(i)-1:idx3(i)+1]),'all')>0
        mask(idx1(i),idx2(i),idx3(i)) = 1;
    end
    catch
    end
end
mask = logical(mask);
data(mask) = 1;

filename_out = [fpath.out,fname.atlas.out2b];
copyfile(filename,filename_out)
tstack_out = Tiff(filename_out,'r+');
K = length(imfinfo(filename_out));
tstack_out.write(uint16(squeeze(data(1,:,:))))
for n = 2:K
    tstack_out.nextDirectory()
    tstack_out.write(uint16(squeeze(data(n,:,:))))
end

%% symmetrical data + region reassignments + root/cerebellum/fibers/ventricles removed (for KDE)

filename = [fpath.out,fname.atlas.out2];
tstack  = Tiff(filename);
[I,J] = size(tstack.read());
K = length(imfinfo(filename));
data = zeros(K,I,J);
data(1,:,:)  = tstack.read();
for n = 2:K
    tstack.nextDirectory()
    data(n,:,:) = tstack.read();
end

mask = data<=1 | data>=1028; data(mask) = 0;

filename_out = [fpath.out,fname.atlas.out3];
copyfile(filename,filename_out)
tstack_out = Tiff(filename_out,'r+');
K = length(imfinfo(filename_out));
tstack_out.write(uint16(squeeze(data(1,:,:))))
for n = 2:K
    tstack_out.nextDirectory()
    tstack_out.write(uint16(squeeze(data(n,:,:))))
end
%% function
function atlas = Add_CEAast(atlas)

atlas(atlas==581) = 1E5;
[r,c,v] = ind2sub(size(atlas),find(atlas == 1E5));
idx = find(c<241);
idx = unique([idx; find((v<150|v>269))]);
for i = 1:length(idx)
    atlas(r(idx(i)),c(idx(i)),v(idx(i))) = 581;
end
idx = [];
[r,c,v] = ind2sub(size(atlas),find(atlas == 1E5));
val = 237;
AMYG = 606:608;
for i = 1:size(atlas,2)
    [r1,c1] = find(ismember(squeeze(atlas(:,i,:)),AMYG));
    if ~isempty(r1)
        idx = [idx ;find(c==i & v>median(c1))];
        val = median(c1);
    else
        idx = [idx ;find(c==i & v>val)];
    end
end
for i = 1:length(idx)
    atlas(r(idx(i)),c(idx(i)),v(idx(i))) = 3E5;
end
% x=bwconncomp(atlas==1E5,6);
% [~,idx]=sort(cellfun(@length,x.PixelIdxList),'descend');
% atlas(x.PixelIdxList{idx(1)}) = 2E5;
% atlas(x.PixelIdxList{idx(2)}) = 2E5;
% x=bwconncomp(atlas==3E5,6);
% [~,idx]=sort(cellfun(@length,x.PixelIdxList),'descend');
% atlas(x.PixelIdxList{idx(1)}) = 4E5;
% atlas(x.PixelIdxList{idx(2)}) = 4E5;
% atlas(atlas==1E5) = 581;
% atlas(atlas==3E5) = 581;
atlas(atlas==1E5) = 605;
atlas(atlas==3E5) = 605;
end