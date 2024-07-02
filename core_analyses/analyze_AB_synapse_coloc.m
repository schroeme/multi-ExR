function [frac_overlap,nnz1_all,nnz2_all] = analyze_AB_synapse_coloc(params,roi,chs)
%created by M.E. Schroeder on 7/27/23
%calculates fraction of mutually overlapping volume for specified synaptic
%and beta-amyloid channels WITHIN cropped beta-amyloid nanoclusters
%calculates the fraction of channel 1 volume mutually overlapped with
%channel 2, so pay attention to channel order
%Last updated by mes on 7/27/23

%Extract parameters
nchannels = length(chs);
parentfolder = params.parentfolder;

%do this for a single pair of channels, because we only have 2 pairs of
%interest

chsplits = split(chs{1},"-");
roundno = chsplits{1};
chno = chsplits{2};
im1_files = dir([parentfolder filesep '*' roi '*round0' roundno '*ch0' chno '*.tif']);

ncrops = length(im1_files);

chsplits = split(chs{2},"-");
roundno = chsplits{1};
chno = chsplits{2};
im2_files = dir([parentfolder filesep '*' roi '*round0' roundno '*ch0' chno '*.tif']);

for ss = 1:ncrops %loop through all cropped ROIs
    img1 = loadtiff([parentfolder filesep im1_files(ss).name]);  
    img1 = double(img1);
    img2 = loadtiff([parentfolder filesep im2_files(ss).name]);
    img2 = double(img2);
    
    %threshold and binarize the images
    if strcmp(params.thresh_method,'absolute')
        params.threshold=params.thresholds(im1);
    end
    img1_bin = binarize_intensity_threshold(img1,params);
    if strcmp(params.thresh_method,'absolute')
        params.threshold=params.thresholds(im2);
    end        
    img2_bin = binarize_intensity_threshold(img2,params);

    %filter the images if needed
    if params.dofilt
        if strcmp(params.filt,'med')
            img1_bin = medfilt3(img1_bin,params.filt_size);
            img2_bin = medfilt3(img2_bin,params.filt_size);
        end
    end

    %apply a size filter to get rid of bright, high-frequency noise
    if params.sizefilt
        img1_bin = bwareaopen(img1_bin,params.lowerlim);
        img2_bin = bwareaopen(img2_bin,params.lowerlim);
    end

    %calculate overlap and correlation
    nnz1 = nnz(img1_bin);
    nnz2 = nnz(img2_bin);
    
    overlap = nnz(img1_bin & img2_bin);
    overlap_pct = overlap/nnz1;
    
     %plot overlap of binary images
    if params.doplot
        figure;
        imshowpair(max(img1_bin,[],3),max(img2_bin,[],3));
        title([roi '- Fraction mutually overlapping = ' num2str(overlap_pct)]);
    end
    
    frac_overlap(ss,1) = overlap_pct;
    nnz1_all(ss,1) = nnz1*params.vol_converter;
    nnz2_all(ss,1) = nnz2*params.vol_converter;
end

end



