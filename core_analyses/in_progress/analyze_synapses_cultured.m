function data = analyze_synapses_cultured(params,roi,synchannels)
%calculates CaMKii volume, signal, and # of objects outside of morphology
%and ABeta channels

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
    
    if length(im1_files)>0
        img1 = load3DTif_uint16([folder filesep im1_files(1).name]);    

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
            
            img1 = double(img1);
            img1 = img1 .* imcomplement(max(morph_bin,[],3));
            img1(img1==0) = mean(img1(:)); %set masked region to mean intensity
            
        end

        %threshold and binarize the images
        if strcmp(params.thresh_method,'absolute')
            thresholds = params.thresholds;
            img1_bin = imbinarize(img1,thresholds(im1));
        elseif strcmp(params.thresh_method,'pct')
            thresh1 = prctile(img1(:),params.thresh_pctl);
            img1_bin = imbinarize(img1,thresh1);
        elseif strcmp(params.thresh_method,'zscore')
            meanint = mean(img1(:));
            stdev = std(img1(:));
            thresh1 = meanint + params.thresh_multiplier * stdev;
            img1_bin = imbinarize(img1,thresh1);
        end

        %filter the images
        if params.dofilt
            img1_bin = medfilt3(img1_bin,params.filt_size);
        end

        %get rid of very small synapses less than 100 voxels
        img1_bin = bwareaopen(img1_bin,params.lowerlim);
        
        %get rid of very large volumes
        CC = bwconncomp(img1_bin,26);
        volumes = regionprops3(CC,'Volume');
        volumes = volumes.Volume;
        nobjects_temp=CC.NumObjects;

        blank = zeros(size(img1_bin));
        for obidx = 1:nobjects_temp
            if volumes(obidx) <= params.upperlim
                blank(CC.PixelIdxList{obidx})=1;
            end
        end
        
        img1_bin=blank;
        %calculate overlap and correlation
        nnz1 = nnz(img1_bin);

        %plot binary images
        if params.doplot
            figure;
            imagesc(max(img1_bin,[],3));
            title(synchannels{im1})
        end

        %mask the images to get rid of background noise. we don't want to
        %correlate that
        img1_masked = double(img1 .* img1_bin);
        
        vol(im1,1)=nnz1 * params.vol_converter;

        %calculate signal within masked region
        mean_signal(im1,1) = mean(img1_masked(img1_masked>0));
        int_signal(im1,1) = sum(sum(img1_masked(img1_masked>0)));
        CC = bwconncomp(img1_bin,26);
        nobjects(im1,1)=CC.NumObjects;

        %extract volume, aspect ratio from synapses
        voltab = regionprops3(CC,'Volume');
        vols = voltab.Volume;
        vols = vols * params.vol_converter;

        paxtab = regionprops3(CC,'PrincipalAxisLength');
        pax = paxtab.PrincipalAxisLength;
        if size(pax,2) > 1
            ARs = pax(:,1)./pax(:,2);
        else
            ARs = NaN;
        end

        vol_temp(im1,1) = nanmean(vols);
        vol_std_temp(im1,1) = std(vols);
        AR_temp(im1,1)= nanmean(ARs);
    else
        mean_signal(im1,1) = NaN;
        int_signal(im1,1)=NaN;
        vol(im1,1)=NaN;
        nobjects(im1,1)=NaN;
        vol_temp(im1,1) = NaN;
        AR_temp(im1,1)=NaN;
        vol_std_temp(im1,1)=NaN;
    end
end

data.vol = vol;
data.mean_signal = mean_signal;
data.int_signal = int_signal;
data.nobjects=nobjects;
data.AR = AR_temp;
data.mean_vol = vol_temp;
data.std_vol = vol_std_temp;

end


