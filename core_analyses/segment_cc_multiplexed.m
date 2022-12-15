function [nsynapses,imdata,syndata] = segment_cc_multiplexed(fov,params)
% Automated synapse segmentation from multi-ExR data
% This script has not been updated/maintained since June 2022
% Please contact Margaret at mschro@mit.edu if you want to collaborate on
% this, but it was deprioritized before manuscript submission

%Extract parameters
parentfolder = params.parentfolder;
xystep = params.xystep;
zstep = params.zstep;
method = params.method;
postinds = params.postinds;%channel indices of post-synaptic stains
preinds = params.preinds;%channel indices of pre-synaptic stains
morphclose = params.morphclose;
minchannels = params.minchannels;
savechunks = params.savechunks;
chunkfoldername = params.chunkfoldername;
refch = params.refch;
matchinginds = params.matchinginds;
ctrlshift_xy = floor(params.controlshift/xystep);
ctrlshift_z = 0;%floor(params.controlshift/zstep);
buffersize = params.buffersize;

fnames = dir([parentfolder '*' fov '*.tif']); %binary images

if params.dothresh %if images need to be binarized
    if params.dofilt %if we're passing them through median filter
        [Iall,Icombo] = load_and_combine_and_mask_multiplexed_channels(parentfolder,fov,fnames,refch,params.thresh_pct,params.medfilt);
    end
else
    [Iall,Icombo] = load_and_combine_multiplexed_channels([parentfolder 'masks/'],fov,fnames,refch);
end

if method == 'minimum'
    Icombo(Icombo<minchannels) = 0;
    Icombo(Icombo>=minchannels) = 1;
elseif method == 'prepost' 
    Icombo = id_synapses_pre_and_post_coloc(preinds,postinds,Iall,xystep);
else
    disp(['Please specify method for synapse identification, either prepost or minimum'])
end

closingrad = params.closingrad*(1/xystep); %radius for morphological closing, in nm
if morphclose
        %I_c=imclose(I_f0,strel('sphere',closingrad)); 
    I_c=imclose(Icombo,strel('disk',floor(closingrad))); %morphological closing
else
    I_c=Icombo; %no morphological closing
end
    
%apply lower filter on both after closing
%set lower limit on ROI size to minimize noise: .1 um in each dimension
lowerlim = ceil((params.lowerlim^3)*(1/xystep)*(1/xystep)*(1/zstep));
I_f = bwareaopen(I_c, lowerlim); %filtration, removes binary objects LESS than this size

CC_combo = bwconncomp(I_f,26);%find all connected components in the combined image
nobjects = CC_combo.NumObjects;

%size upper filtration
upperlim = (params.upperlim^3)*(1/xystep)*(1/xystep)*(1/zstep); %upper limit of 1um^3
    
centroidstab = regionprops3(CC_combo,'Centroid');
centroids = centroidstab.Centroid;
imsize = size(I_c);

%create segmentation image
I_s = create_synapse_seg_mask(CC_combo,buffersize,upperlim,xystep,zstep,imsize);

%CCfinal = bwconncomp(I_s,26); %connected components on filtered merged image
Ifinal = I_s.*I_f;
%CCfinal = CC_combo(setdiff(1:nobjects,toremove));
CCfinal = bwconncomp(Ifinal,26);
nsynapses=CCfinal.NumObjects;
CCseg = bwconncomp(I_s,26);

%save down the segmentation as tiff stack in same folder
imtosave1 = uint8(255*I_s);
imtosave2 = uint8(255*Ifinal);
tempstr = split(fnames(1).name,'_');
tempstr = split(tempstr{1},'-');
name=tempstr{1};%[tempstr{2} '-' tempstr{3}];
imwrite(imtosave1(:,:,1),[parentfolder 'masks/Seg-' name '.tif'],'tiff')
imwrite(imtosave2(:,:,1),[parentfolder 'masks/Bin_filt-' name '.tif'],'tiff')
for jj = 2:size(imtosave1,3)
    imwrite(imtosave1(:,:,jj),[parentfolder 'masks/Seg-' name '.tif'],'WriteMode','append')
    imwrite(imtosave2(:,:,jj),[parentfolder 'masks/Bin_filt-' name '.tif'],'WriteMode','append')
end

%% now extract measures of interest from segmented synapses and save the segmented synapses, if desired

%get correlation and overall overlap between channels for whole image
for chidx = 1:size(Iall,2)
    Itemp = Iall{chidx};
    Isyn{chidx} = Itemp.*I_s;
    Isyn_trans{chidx} = imtranslate(Isyn{chidx},[ctrlshift_xy 0 0]);
end

%calculate pairwise correlation coefficients and overlaps between all channels
imdata = calculate_pairwise_correlation_and_overlap(Iall,Isyn,Isyn_trans,ctrlshift_xy);

if savechunks
    save_chunks(parentfolder,fov,chunkfoldername,CCseg,refch,name);
end

centroids = regionprops3(CCfinal,'Centroid');
centroids = centroids.Centroid;
volumes = regionprops3(CCfinal,'Volume');
volumes = volumes.Volume;
extents = regionprops3(CCfinal,'Extent');
extents = extents.Extent;
eqdis = regionprops3(CCfinal,'EquivDiameter');
eqdis = eqdis.EquivDiameter;
orients = regionprops3(CCfinal,'Orientation');
orients = orients.Orientation;
SAs = regionprops3(CCfinal,'SurfaceArea');
SAs = SAs.SurfaceArea;

for obj = 1:nsynapses %extract relevant data for each synapse, all channels
    syndata(obj).vol = volumes(obj);
    syndata(obj).extent = extents(obj);
    syndata(obj).equiv_diameter = eqdis(obj);
    syndata(obj).centroid = centroids(obj);
    syndata(obj).orient = orients(obj);
    syndata(obj).SA = SAs(obj);
end

end


