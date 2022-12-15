function bin = binarize_intensity_threshold(img,params)
%Binarizes an image using an intensity-based threshold

if strcmp(params.thresh_method,'pct')
    threshold = prctile(img(:),params.thresh_pctl);
    bin = imbinarize(img,threshold);
elseif strcmp(params.thresh_method,'zscore')
    meanint = mean(img(:));
    stdint = std(img(:));
    threshold = meanint + params.thresh_multiplier*stdint;
    bin = imbinarize(img,threshold);
elseif strcmp(params.thresh_method,'absolute')
    bin = imbinarize(img,params.threshold);
end
end

