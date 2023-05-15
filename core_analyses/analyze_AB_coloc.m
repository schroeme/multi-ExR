function data = analyze_AB_coloc(params,roi,AB_chs)
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
    
    for im2 = im1:nchannels
        chsplits = split(AB_chs{im2},"-");
        roundno = chsplits{1};
        chno = chsplits{2};
        im2_files = dir([folder filesep '*' roi '*' roundno '*' chno '*.tif']);
        img2 = load3DTif_uint16([folder filesep im2_files(1).name]);
        
        %threshold and binarize the images
        if strcmp(params.thresh_method,'otsu')
            thresh1 = graythresh(img1(:));
            thresh2 = graythresh(img2(:));
            
            img1_bin = imbinarize(img1,thresh1);
            img2_bin = imbinarize(img2,thresh2);
        elseif strcmp(params.thresh_method,'zscore')
            mean1 = mean(img1(:));
            mean2 = mean(img2(:));

            stdev1 = std(img1(:));
            stdev2 = std(img2(:));

            thresh1 = mean1 + params.thresh_multiplier*stdev1;
            thresh2 = mean2 + params.thresh_multiplier*stdev2;

            img1_bin = imbinarize(img1,thresh1);
            img2_bin = imbinarize(img2,thresh2);

        elseif strcmp(params.thresh_method,'pct')
            thresh1 = prctile(img1(:),params.thresh_pctl);
            thresh2 = prctile(img2(:),params.thresh_pctl);
            
            img1_bin = imbinarize(img1,thresh1);
            img2_bin = imbinarize(img2,thresh2);
        elseif strcmp(params.thresh_method,'absolute')
            img1_bin = imbinarize(img1,params.thresholds(im1));
            img2_bin = imbinarize(img2,params.thresholds(im2));
        end
        
        %filter the images
        if params.dofilt
            img1_bin = medfilt3(img1_bin,params.filt_size);
            img2_bin = medfilt3(img2_bin,params.filt_size);
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

            img2_bin = double(img2_bin) - double(max(morph_bin,[],3));
            img2_bin(img2_bin<0) = 0;
        else
            img1_bin = img1_bin;
            img2_bin = img2_bin;
        end

        %apply a size filter to get rid of bright, high-frequency noise
        img1_bin = bwareaopen(img1_bin,params.lowerlim);
        img2_bin = bwareaopen(img2_bin,params.lowerlim);
        
        vol(im1,1) = nnz(img1_bin(:)) * params.vol_converter;

        %calculate overlap and correlation
        nnz1 = nnz(img1_bin);
        nnz2 = nnz(img2_bin);
        
        overlap = nnz(img1_bin & img2_bin);
        overlap_pct = overlap*2 / (nnz1+nnz2);
        
         %plot overlap of binary images
        if params.doplot
            figure;
            imshowpair(max(img1_bin,[],3),max(img2_bin,[],3));
            title([roi '- Fraction mutually overlapping = ' num2str(overlap_pct)]);
        end

%         cc1 = bwconncomp(img1_bin);
%         cc2 = bwconncomp(img2_bin);

%         if params.doplot
%             L1 = labelmatrix(cc1);
%             L2 = labelmatrix(cc2);
%             mip1 = max(L1,[],3);
%             mip2 = max(L2,[],3);
% 
%             RGB_label1 = label2rgb(mip1,@turbo,"c","shuffle");
%             figure;
%             imshow(RGB_label1)
% 
%             RGB_label2 = label2rgb(mip2,@turbo,"c","shuffle");
%             figure;
%             imshow(RGB_label2)
%         end
        
        frac_overlap(im1,im2) = overlap_pct;
        
        %mask the images to get rid of background noise. we don't want to
        %correlate that
        img1_masked = double(img1 .* img1_bin);
        img2_masked = double(img2 .* img2_bin);
        
        %calculate pixel-wise correlation of masked images
        pxwise_corr(im1,im2) = corr(img1_masked(:),img2_masked(:));
        
        %calculate signal within masked region
        mean_signal(im1,1) = mean(img1_masked(img1_masked>0));
        int_signal(im1,1) = sum(sum(img1_masked(img1_masked>0)));

    end
end

data.frac_overlap = frac_overlap;
data.pxwise_corr = pxwise_corr;
data.vol = vol;
data.mean_signal = mean_signal;
data.int_signal = int_signal;

end


