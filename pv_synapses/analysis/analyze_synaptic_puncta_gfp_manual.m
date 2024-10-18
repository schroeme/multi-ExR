function [data] = analyze_synaptic_puncta_gfp_manual(params)

%Extract parameters
maskfolder = params.outfolder;
targetfolder = params.targetoutfolder;
parentfolder = params.parentfolder;

%Get filenames for existing masks
gfpp_segs = dir([maskfolder filesep '*bbox_seg_GFPpos.tif']); %bounding boxes
all_masks = dir([maskfolder filesep '*pp_seg_all.tif']); %not bounding boxes
fnames = dir([parentfolder filesep '*.tif']); %preprocessed images
target_masks = dir([targetfolder filesep '*.tif']);
nfovs = length(gfpp_segs);

nsyn_estimate = ceil(nfovs*params.nsyn_estimate*1.1); 

npunctas_all = zeros(nsyn_estimate,1);
vols_all = zeros(nsyn_estimate,1);
meanints_all = zeros(nsyn_estimate,1);
names = cell(nsyn_estimate,1);
meanints_refnorm_all = zeros(nsyn_estimate,1);
counter = 0;
selected_tf = zeros(nsyn_estimate,1);

%loop through each field of view in the folder
for fidx = 1:nfovs
    
    %load up masks
    gfpp_seg = mat2gray(loadtiff([maskfolder gfpp_segs(fidx).name]));
    all_seg = mat2gray(loadtiff([maskfolder all_masks(fidx).name]));
    target_mask = double(loadtiff([targetfolder target_masks(fidx).name]));

    %load up background-subtracted image (NOT masked, we need this for mean
    %intensity)
    img = loadtiff([parentfolder fnames(fidx).name]);
    img = double(img);
    ch2img = img(:,:,2:3:end-1);

    CC_all = bwconncomp(all_seg,26);

    %get bounding boxes for all synapses
    bboxes_all = regionprops3(CC_all,'BoundingBox');
    nobjects_all = CC_all.NumObjects;

    buffer_half = params.nbuffer/2;
    zbuffer_half = 2;

    %Next, quantify signal within bounding boxes
    for ss = 1:nobjects_all
        counter = counter+1;
        bbox = bboxes_all(ss,1).BoundingBox;
        ulf_x = bbox(1);
        ulf_y = bbox(2);
        ulf_z = bbox(3);
        width_x = bbox(4);
        width_y = bbox(5);
        width_z = bbox(6);

        %Note: somehow these bounding boxes are larger than 
        % Calculate the start and end indices for each dimension
        xStart = max(1, floor(ulf_x) - buffer_half);
        xEnd = min(size(target_mask,2), floor(ulf_x + width_x) + buffer_half);
        yStart = max(1, floor(ulf_y) - buffer_half);
        yEnd = min(size(target_mask,1), floor(ulf_y + width_y) + buffer_half);
        zStart = max(1, floor(ulf_z) - zbuffer_half);
        zEnd = min(size(target_mask,3), floor(ulf_z + width_z) + zbuffer_half);

        %Grab the bounding box area in target channel
        ch2chunk = ch2img(yStart:yEnd, xStart:xEnd, zStart:zEnd);

        %Grab the bounding box area in the reference channel, binarized
        refchunk = all_seg(yStart:yEnd, xStart:xEnd, zStart:zEnd);

        %use previously generated mask for this
        chunk_filt = target_mask(yStart:yEnd, xStart:xEnd, zStart:zEnd);

        %grab the corresponding area in the "selected" GFP+ mask
        selected_filt = gfpp_seg(yStart:yEnd, xStart:xEnd, zStart:zEnd);

        if mean(selected_filt(:)) ==1
            selected_tf(counter,1) = 1;
        else
            selected_tf(counter,1)=0;
        end

        ch2masked = ch2chunk .* chunk_filt;
        CC2 = bwconncomp(ch2masked,26);%find all connected components in the combined image
        n2 = CC2.NumObjects;

        if n2 > 0
            mean_int = sum(ch2masked(:))/nnz(ch2masked(:));
            mean_int_refnorm = sum(ch2masked(:))/nnz(refchunk(:));
            voltab=regionprops3(CC2,'Volume');
            volstemp=sum(voltab.Volume,1);
            vols = volstemp * params.vol_converter;
            npuncta = n2;
        else
            mean_int = 0;
            vols = 0;
            npuncta = n2;
            mean_int_refnorm = 0;
        end

        name = [fnames(fidx).name(1:end-4) '_syn' num2str(ss)];
        names{counter} = name;
        meanints_all(counter,1) = mean_int;
        npunctas_all(counter,1) = npuncta;
        vols_all(counter,1) = vols;
        meanints_refnorm_all(counter,1) = mean_int_refnorm;
        if mod(counter,10000) == 0
            disp(['Analyzed ' num2str(counter) ' synapses!']);
        end

    end

end
    
data.names = names;
data.total_vol = vols_all;
data.mean_int = meanints_all;
data.npuncta = npunctas_all;
data.selected_tf = selected_tf;
data.mean_int_refnorm = meanints_refnorm_all;
