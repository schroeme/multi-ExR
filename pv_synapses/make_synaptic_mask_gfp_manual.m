function [fnames,nsyn_GFP,nsyn_all,fov_vol_all] = make_synaptic_mask_gfp_manual(params)

%Extract parameters
parentfolder = params.parentfolder;
fnames = dir([parentfolder filesep '*.tif']);
mask_fnames = dir([params.maskfolder filesep '*.tif']);
nfovs = length(fnames);

nsyn_GFP = zeros(nfovs,1);
nsyn_all = zeros(nfovs,1);
fov_vol_all = zeros(nfovs,1);

%loop through each field of view in the folder
for fidx = 1:nfovs
    img = loadtiff([parentfolder fnames(fidx).name]);
    img = double(img);

    ch1img = img(:,:,1:3:end-2);
    ch2img = img(:,:,2:3:end-1);
    ch3img = img(:,:,3:3:end);

    %% Create synapse and target segmentations
    if params.gfilt
        ch1img_blur = imgaussfilt3(ch1img,params.sigma);
    else
        ch1img_blur=ch1img;
    end

    %Use channel 1 to segment synapses
    imbin = binarize_intensity_threshold(ch1img_blur,params);
    imbin_target = binarize_intensity_threshold(ch2img,params);

    if strcmp(params.filt,'med')
        imbin_filt = medfilt3(imbin,params.filt_size);
        imbin_target_filt = medfilt3(imbin_target,params.filt_size);
    else
        imbin_filt = imbin;
        imbin_target_filt = imbin_target;
    end

    if params.sizefilt
        mask = bwareaopen(imbin_filt, params.lowerlim); %filtration, removes binary objects LESS than this size
    else
        mask = imbin_filt;
    end

    %for now, DO NOT SIZE FILTER THE TARGET CHANNEL. Otherwise we bias the
    %larger puncta
    if params.target_sizefilt
        target_mask = bwareaopen(imbin_target_filt, params.target_lowerlim);
    else
        target_mask = imbin_target_filt;
    end

    if params.doplot
        figure;
        subplot(1,2,1)
        imagesc(max(ch1img,[],3))
        title('Original MIP')
        subplot(1,2,2)
        imagesc(max(mask,[],3))
        title('Segmentation MIP')
    end

    %% load mask for GFP
    gfp_mask_2d = loadtiff([params.maskfolder mask_fnames(fidx).name]);
    gfp_mask = repmat(gfp_mask_2d, 1, 1, size(ch3img,3)); %with help from chatgpt 4.0 on 8/3/24

    %% Create bounding boxes
    mask_pos = gfp_mask & mask;%GFP+ synapses
    
    %get rid of any small puncta remaining after masking for GFP
    mask_pos = bwareaopen(mask_pos, params.lowerlim);
    
    CC_pos = bwconncomp(mask_pos,26);%find all connected components in the combined image
    nobjects_pos = CC_pos.NumObjects;
    bboxes_pos = regionprops3(CC_pos,'BoundingBox');
    nsyn_GFP(fidx,1) = nobjects_pos;
    fov_vol_all(fidx,1) = size(ch1img,1)*params.xystep*size(ch1img,2)*params.xystep*size(ch1img,3); %volume of the field of view in cubic microns

   %From MATLAB documentation: Smallest cuboid containing the region,
    % returned as a 1-by-6 vector of the form [ulf_x ulf_y ulf_z width_x width_y width_z].
    % ulf_x, ulf_y, and ulf_z specify the upper-left front corner of the cuboid.
    % width_x, width_y, and width_z specify the width of the cuboid along each dimension.

    % Below with help from ChatGPT 4.0 by MES on 7/7/24
    segimg_pos = false(size(mask_pos));

    buffer_half = params.nbuffer/2;
    zbuffer_half = 2;%floor(buffer_half/2);
    for k = 1:nobjects_pos
        bbox = bboxes_pos(k,1).BoundingBox;
        ulf_x = bbox(1);
        ulf_y = bbox(2);
        ulf_z = bbox(3);
        width_x = bbox(4);
        width_y = bbox(5);
        width_z = bbox(6);

        % Calculate the start and end indices for each dimension
        xStart = max(1, floor(ulf_x) - buffer_half);
        xEnd = min(size(mask,2), floor(ulf_x + width_x) + buffer_half);
        yStart = max(1, floor(ulf_y) - buffer_half);
        yEnd = min(size(mask,1), floor(ulf_y + width_y) + buffer_half);
        zStart = max(1, floor(ulf_z) - zbuffer_half);
        zEnd = min(size(mask,3), floor(ulf_z + width_z) + zbuffer_half);
        
        % Set the bounding box region to 1 (white)
        segimg_pos(yStart:yEnd, xStart:xEnd, zStart:zEnd) = true;
    end

    %For now, save a mask of ALL synapses (not just GFP-negative)
    CC_all = bwconncomp(mask,26);
    segimg = false(size(mask));
    bboxes = regionprops3(CC_all,'BoundingBox');
    nobjects = CC_all.NumObjects;
    nsyn_all(fidx,1) = nobjects;

    for k = 1:nobjects
        bbox = bboxes(k,1).BoundingBox;
        ulf_x = bbox(1);
        ulf_y = bbox(2);
        ulf_z = bbox(3);
        width_x = bbox(4);
        width_y = bbox(5);
        width_z = bbox(6);

        % Calculate the start and end indices for each dimension
        xStart = max(1, floor(ulf_x) - buffer_half);
        xEnd = min(size(mask,2), floor(ulf_x + width_x) + buffer_half);
        yStart = max(1, floor(ulf_y) - buffer_half);
        yEnd = min(size(mask,1), floor(ulf_y + width_y) + buffer_half);
        zStart = max(1, floor(ulf_z) - zbuffer_half);
        zEnd = min(size(mask,3), floor(ulf_z + width_z) + zbuffer_half);
        
        % Set the bounding box region to 1 (white)
        segimg(yStart:yEnd, xStart:xEnd, zStart:zEnd) = true;
    end

    % mask_neg = ~segimg_pos & mask;%GFP- synapses
    % mask_neg = bwareaopen(mask_neg, params.lowerlim); %filtration, removes binary objects LESS than this size
    % 
    % CC_neg = bwconncomp(mask_neg,26);%find all connected components in the combined image
    % nobjects_neg = CC_neg.NumObjects;
    % bboxes_neg = regionprops3(CC_neg,'BoundingBox');
    % 
    % segimg_neg = false(size(mask_pos));
    % for k = 1:nobjects_neg
    %     bbox = bboxes_neg(k,1).BoundingBox;
    %     ulf_x = bbox(1);
    %     ulf_y = bbox(2);
    %     ulf_z = bbox(3);
    %     width_x = bbox(4);
    %     width_y = bbox(5);
    %     width_z = bbox(6);
    % 
    %     % Calculate the start and end indices for each dimension
    %     xStart = max(1, floor(ulf_x) - buffer_half);
    %     xEnd = min(size(mask,2), floor(ulf_x + width_x) + buffer_half);
    %     yStart = max(1, floor(ulf_y) - buffer_half);
    %     yEnd = min(size(mask,1), floor(ulf_y + width_y) + buffer_half);
    %     zStart = max(1, floor(ulf_z) - zbuffer_half);
    %     zEnd = min(size(mask,3), floor(ulf_z + width_z) + zbuffer_half);
    % 
    %     % Set the bounding box region to 1 (white)
    %     segimg_neg(yStart:yEnd, xStart:xEnd, zStart:zEnd) = true;
    % end

    if params.savemasks
        options.overwrite=true;
        saveastiff(uint8(mask_pos),[params.outfolder fnames(fidx).name(1:end-4) '_seg_GFPpos.tif'],options);
        saveastiff(uint8(segimg_pos),[params.outfolder fnames(fidx).name(1:end-4) '_bbox_seg_GFPpos.tif'],options);

        saveastiff(uint8(mask),[params.outfolder fnames(fidx).name(1:end-4) '_seg_all.tif'],options);
        saveastiff(uint8(segimg),[params.outfolder fnames(fidx).name(1:end-4) '_bbox_seg_all.tif'],options);

        % saveastiff(uint8(mask_neg),[params.outfolder fnames(fidx).name(1:end-4) '_seg_GFPneg.tif'],options);
        % saveastiff(uint8(segimg_neg),[params.outfolder fnames(fidx).name(1:end-4) '_bbox_seg_GFPneg.tif'],options);
    end

    if params.savetargetmasks
        options.overwrite=true;
        saveastiff(uint8(target_mask),[params.targetoutfolder fnames(fidx).name(1:end-4) '_seg.tif'],options);
    end

end
