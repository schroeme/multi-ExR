function data = analyze_camkii(params,roi,camch,AB_chs)
%calculates correlation and volume overlap between different AB channels

folder = params.parentfolder;

%Extract parameters
xystep = params.xystep;
zstep = params.zstep;
%thresholds = params.thresholds;

chsplits = split(camch,"-");
roundno = chsplits{1};
chno = chsplits{2};
im1_files = dir([folder filesep '*' roi '*' roundno '*' chno '*.tif']);

if length(im1_files)>0
    img1 = load3DTif_uint16([folder filesep im1_files(1).name]);    

    %threshold and binarize the images
    if strcmp(params.thresh_method,'absolute')
        img1_bin = imbinarize(img1,thresholds(im1));
    elseif strcmp(params.thresh_method,'pct')
        thresh1 = prctile(img1(:),params.thresh_pctl);
        img1_bin = imbinarize(img1,thresh1);
    end

    %filter the images
    if params.dofilt
        img1_bin = medfilt3(img1_bin,params.filt_size);
    end

    %mask out ABeta channel
    AB_mask = zeros(size(img1_bin));
    for aa = 1:length(AB_chs)
        chsplits = split(AB_chs{aa},"-");
        roundno = chsplits{1};
        chno = chsplits{2};

        AB_files = dir([folder filesep '*' roi '*' roundno '*' chno '*.tif']);
        AB_img = load3DTif_uint16([folder filesep AB_files(1).name]); 

        thresh_AB = prctile(AB_img(:),99.9);
        AB_bin = imbinarize(AB_img,thresh_AB);

        AB_filt=medfilt3(AB_bin,params.filt_size);

        se=strel('disk',params.morph_rad);
        AB_bin=imdilate(AB_filt,se);
        AB_mask = AB_mask + AB_bin;
    end
    
    AB_mask = AB_mask>1;
    img1_bin = double(img1_bin) - double(max(AB_mask,[],3));
    img1_bin(img1_bin<0) = 0;
        
    %subtract morphology channel if desired
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

    %get rid of very small puncta
    img1_bin = bwareaopen(img1_bin,params.lowerlim);
    
    %get rid of very large volumes
    CC = bwconncomp(img1_bin,26);
    volumes = regionprops3(CC,'Volume');
    volumes = volumes.Volume;
    nobjects=CC.NumObjects;
    
    blank = zeros(size(img1_bin));
    for obidx = 1:nobjects
        if volumes(obidx) <= params.upperlim
            blank(CC.PixelIdxList{obidx})=1;
        end
    end
    
    img1_bin = double(blank);
    %calculate overlap and correlation
    nnz1 = nnz(img1_bin);

    %plot binary images
    if params.doplot
        figure;
        imagesc(max(img1_bin,[],3));
    end

    %mask the images to get rid of background noise. we don't want to
    %correlate that
    img1_masked = double(img1 .* img1_bin);
    vol=nnz1 * params.vol_converter;

    %calculate signal within masked region
    mean_signal = mean(img1_masked(img1_masked>0));
    int_signal = sum(sum(img1_masked(img1_masked>0)));
    
    CC = bwconncomp(img1_bin,26);
    volumes = regionprops3(CC,'Volume');
    volumes = volumes.Volume;
    mean_vol = mean(volumes);
    mean_vol = mean_vol * params.vol_converter;
    nobjects=CC.NumObjects;
    
else
    mean_signal = NaN;
    int_signal=NaN;
    vol=NaN;
    nobjects=NaN;
    mean_vol=NaN;
end

data.vol = vol;
data.mean_signal = mean_signal;
data.int_signal = int_signal;
data.nobjects=nobjects;
data.mean_vol=mean_vol;
end

