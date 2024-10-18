function thresholds = determine_intensity_threshold(params,roi,synchannels)
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
    
    if length(im1_files)>0
        disp(im1_files(1).name);
        img1 = load3DTif_uint16([folder filesep im1_files(1).name]);    
        img1 = double(img1);

        %take only the middle 80% of the image (which should definitely
        %have signal
        xlim1 = ceil(size(img1,1)/2) - ceil((size(img1,1)*.8)/2);
        xlim2 = xlim1 + ceil((size(img1,1)*.8));

        ylim1 = xlim1;
        ylim2 = xlim2;

        zlim1 = ceil(size(img1,3)/2) - ceil((size(img1,3)*.8)/2);
        zlim2 = zlim1 + ceil((size(img1,3)*.8));

        img1 = img1(xlim1:xlim2,ylim1:ylim2,zlim1:zlim2);

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



