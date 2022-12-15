function thresholds = determine_intensity_threshold_5xFAD(params,roi,synchannels)
%calculates correlation and volume overlap between different AB channels

folder = params.parentfolder;

%Extract parameters
xystep = params.xystep;
zstep = params.zstep;
nchannels = length(synchannels);

for im1 = 1:nchannels %nested loops because we are looking at colocalization
    chsplits = split(synchannels{im1},"-");
    roundno = chsplits{1};
    chno = chsplits{2};
    im1_files = dir([folder filesep '*' roi '*' roundno '*' chno '*.tif']);

    %mask out the morphology channel to exclude nonspecific signal
    
    if length(im1_files)>0
        img1 = load3DTif_uint16([folder filesep im1_files(1).name]);    

        morph_files = dir([folder filesep '*' roi '*' params.morph_round '*' params.morph_ch '*.tif']);
        morph_img = load3DTif_uint16([folder filesep morph_files(1).name]); 
        thresh_morph = prctile(morph_img(:),99);
        morph_bin = imbinarize(morph_img,thresh_morph);

        %filter the images
        if params.dofilt
            morph_bin = medfilt3(morph_bin,params.filt_size);
        end
        
        se=strel('disk',params.morph_rad);
        morph_bin=imdilate(morph_bin,se);
        
        img1 = double(img1);
        img1 = img1 .* imcomplement(max(morph_bin,[],3));
        img1(img1==0) = mean(img1(:)); %set masked region to mean intensity

        %threshold and binarize the images
        if strcmp(params.thresh_method,'otsu')
            thresh1 = graythresh(img1(:));

            img1_bin = imbinarize(img1,thresh1);
        elseif strcmp(params.thresh_method,'pct')
            thresh1 = prctile(img1(:),params.thresh_pctl);

            img1_bin = imbinarize(img1,thresh1);
        elseif strcmp(params.thresh_method,'zscore')
            meanint = mean(img1(:));
            stdev = std(img1(:));

            thresh1 = meanint + stdev * params.thresh_multiplier;
        end

        thresholds(1,im1) = thresh1;
    else
        thresholds(1,im1)=NaN;
    end
end



