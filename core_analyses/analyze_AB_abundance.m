function data = analyze_AB_abundance(params,roi,AB_chs)
%calculates correlation and volume overlap between different AB channels

folder = params.parentfolder;

%Extract parameters
xystep = params.xystep;
zstep = params.zstep;
nchannels = length(AB_chs);

for im1 = 1:nchannels %nested loops because we are looking at colocalization
    chsplits = split(AB_chs{im1},"-");
    roundno = chsplits{1};
    chno = chsplits{2};
    im1_files = dir([folder filesep '*' roi '*' roundno '*' chno '*.tif']);
    img1 = load3DTif_uint16([folder filesep im1_files(1).name]);    
    
    %threshold and binarize the images
    if strcmp(params.thresh_method,'otsu')
        thresh1 = graythresh(img1(:));
        
        img1_bin = imbinarize(img1,thresh1);
    elseif strcmp(params.thresh_method,'zscore')
        mean1 = mean(img1(:));

        stdev1 = std(img1(:));

        thresh1 = mean1 + params.thresh_multiplier*stdev1;

        img1_bin = imbinarize(img1,thresh1);
    elseif strcmp(params.thresh_method,'pct')
        thresh1 = prctile(img1(:),params.thresh_pctl);
        
        img1_bin = imbinarize(img1,thresh1);
    elseif strcmp(params.thresh_method,'absolute')
        img1_bin = imbinarize(img1,params.thresholds(im1));
    end
    
    %filter the images
    if params.dofilt
        img1_bin = medfilt3(img1_bin,params.filt_size);
    end
    
    if params.subtract_morph
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

        img1_bin = double(img1_bin) - double(max(morph_bin,[],3));
        img1_bin(img1_bin<0) = 0;

    else
        img1_bin = img1_bin;
    end

    %apply a size filter to get rid of bright, high-frequency noise
    img1_bin = bwareaopen(img1_bin,params.lowerlim);
    
    vol(im1,1) = nnz(img1_bin(:)) * params.vol_converter;

    
     %plot overlap of binary images
    if params.doplot
        figure;
        imshow(max(img1_bin,[],3));
        title([roi '- Mask']);
    end
    
    %mask the images to get rid of background noise.
    img1_masked = double(img1 .* img1_bin);
    
    %calculate signal within masked region
    mean_signal(im1,1) = mean(img1_masked(img1_masked>0));
    int_signal(im1,1) = sum(sum(img1_masked(img1_masked>0)));

end

data.vol = vol;
data.mean_signal = mean_signal;
data.int_signal = int_signal;

end


