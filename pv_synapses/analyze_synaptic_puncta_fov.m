function [data] = analyze_synaptic_puncta_fov(params)
% Analyzes validation data, with same proteins stained in consecutive rounds
% Outputs 
% num_puncta: the number of puncta in each synapse
% punctavol: the total puncta volume for each synapse
% AR: the aspect ratio of the most ellipsoid volume in each synapse
% SNR: signal to noise ratio in synaptic ROI

%Extract parameters
xystep = params.xystep;
zstep = params.zstep; % um/voxel in z
parentfolder = params.parentfolder;
fnames = dir([parentfolder filesep '*.tif']);
nfovs = length(fnames);

nsyn_estimate = ceil(nfovs*params.nsyn_estimate*1.1); 

npunctas_all = zeros(nsyn_estimate,3);
vols_all = zeros(nsyn_estimate,3);
meanints_all = zeros(nsyn_estimate,3);
names = cell(nsyn_estimate,1);
counter = 0;

%loop through each field of view in the folder
for fidx = 1:nfovs
    img = loadtiff([parentfolder fnames(fidx).name]);
    img = double(img);
    nslices=size(img,3);
    nzslices=nslices/params.nchannels;

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
    target_mask = imbin_target_filt;
    
    if params.doplot
        figure;
        subplot(1,2,1)
        imagesc(max(ch1img,[],3))
        title('Original MIP')
        subplot(1,2,2)
        imagesc(max(mask,[],3))
        title('Segmentation MIP')
    end
    CC_signal = bwconncomp(mask,26);%find all connected components in the combined image
    nobjects = CC_signal.NumObjects;

    %% use channel 3 to segment GFP
    if params.gfilt
        ch3img_blur = imgaussfilt3(ch3img,params.gfp_sigma);
    else
        ch3img_blur=ch1img;
    end

    old_multiplier=params.thresh_multiplier;
    params.thresh_multiplier=params.gfp_thresh;
    imbingfp = binarize_intensity_threshold(ch3img_blur,params);
    params.thresh_multiplier=old_multiplier;

    if strcmp(params.filt,'med')
        imbin_filt_gfp = medfilt3(imbingfp,params.filt_size);
    else
        imbin_filt_gfp = imbingfp;
    end

    if params.sizefilt
        mask_gfp = bwareaopen(imbin_filt_gfp, params.lowerlim); %filtration, removes binary objects LESS than this size
    else
        mask_gfp = imbin_filt_gfp;
    end

    if params.savegfp
        options.overwrite=true;
        saveastiff(double(mask_gfp),[params.gfpfolder fnames(fidx).name(1:end-4) '_gfp_seg.tif'],options);
    end

    bboxes = regionprops3(CC_signal,'BoundingBox');
    %From MATLAB documentation: Smallest cuboid containing the region,
    % returned as a 1-by-6 vector of the form [ulf_x ulf_y ulf_z width_x width_y width_z].
    % ulf_x, ulf_y, and ulf_z specify the upper-left front corner of the cuboid.
    % width_x, width_y, and width_z specify the width of the cuboid along each dimension.

    % Below with help from ChatGPT 4.0 by MES on 7/7/24
    segimg = false(size(mask));

    buffer_half = params.nbuffer/2;
    zbuffer_half = 2;%floor(buffer_half/2);
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

    if params.savemasks
        options.overwrite=true;
        saveastiff(double(mask),[params.outfolder fnames(fidx).name(1:end-4) '_seg.tif'],options);
        saveastiff(double(segimg),[params.outfolder fnames(fidx).name(1:end-4) '_bbox_seg.tif'],options);
    end

    if params.savetargetmasks
        options.overwrite=true;
        saveastiff(double(target_mask),[params.targetoutfolder fnames(fidx).name(1:end-4) '_seg.tif'],options);
    end

    %Next, quantify signal within bounding boxes
    for ss = 1:nobjects
        counter = counter+1;
        bbox = bboxes(ss,1).BoundingBox;
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

        %Grab the bounding box area in all channels
        ch1chunk = ch1img(yStart:yEnd, xStart:xEnd, zStart:zEnd);
        ch2chunk = ch2img(yStart:yEnd, xStart:xEnd, zStart:zEnd);
        ch3chunk = ch3img(yStart:yEnd, xStart:xEnd, zStart:zEnd);

        %use previously generated mask fort his
        ch1chunk_filt = mask(yStart:yEnd, xStart:xEnd, zStart:zEnd);
        ch2chunk_filt = target_mask(yStart:yEnd, xStart:xEnd, zStart:zEnd);
        ch3chunk_filt = mask_gfp(yStart:yEnd, xStart:xEnd, zStart:zEnd);

        ch1masked = ch1chunk .* ch1chunk_filt;
        ch2masked = ch2chunk .* ch2chunk_filt;
        ch3masked = ch3chunk .* ch3chunk_filt;

        CC1 = bwconncomp(ch1masked,26);%find all connected components in the combined image
        CC2 = bwconncomp(ch2masked,26);%find all connected components in the combined image
        CC3 = bwconncomp(ch3masked,26);%find all connected components in the combined image

        n1 = CC1.NumObjects;
        n2 = CC2.NumObjects;
        n3 = CC3.NumObjects;

        if n1 > 0
            mean_int(1,1) = sum(ch1masked(:))/nnz(ch1masked(:));
            voltab=regionprops3(CC1,'Volume');
            volstemp=sum(voltab.Volume,1);
            vols(1,1) = volstemp * params.vol_converter;
            npuncta(1,1) = n1;
        else
            mean_int(1,1) = 0;
            vols(1,1) = 0;
            npuncta(1,1) = n1;
        end

        if n2 > 0
            mean_int(1,2) = sum(ch2masked(:))/nnz(ch2masked(:));
            voltab=regionprops3(CC2,'Volume');
            volstemp=sum(voltab.Volume,1);
            vols(1,2) = volstemp * params.vol_converter;
            npuncta(1,2) = n2;
        else
            mean_int(1,2) = 0;
            vols(1,2) = 0;
            npuncta(1,2) = n2;
        end

        if n3 > 0
            mean_int(1,3) = sum(ch3masked(:))/nnz(ch3masked(:));
            voltab=regionprops3(CC3,'Volume');
            volstemp=sum(voltab.Volume,1);
            vols(1,3) = volstemp * params.vol_converter;
            npuncta(1,3) = n3;
        else
            mean_int(1,3) = 0;
            vols(1,3) = 0;
            npuncta(1,3) = n3;
        end

        name = [fnames(fidx).name(1:end-4) '_syn' num2str(ss)];
        names{counter} = name;
        meanints_all(counter,:) = mean_int;
        npunctas_all(counter,:) = npuncta;
        vols_all(counter,:) = vols;
        if mod(counter,10000) == 0
            disp(['Analyzed ' num2str(counter) ' synapses!']);
        end

    end

end
    
data.names = names;
data.total_vol = vols_all;
data.mean_int = meanints_all;
data.npuncta = npunctas_all;
