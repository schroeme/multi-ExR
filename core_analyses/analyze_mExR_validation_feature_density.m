function [volume,volume_frac,min_vol,max_vol] = analyze_mExR_validation_feature_density(fov,params)
% Calculates the volume (and as a fraction of totla volume) in each channel

%Extract parameters
xystep = params.xystep;
zstep = params.zstep; % um/voxel in z
nchannels = params.nchannels;
parentfolder = params.parentfolder;
nrounds = params.nrounds;

for rridx = 1:nrounds %loop through each round
    for chidx = 1:nchannels %loop through each channel
            
        fnames = dir([parentfolder '*' fov '*' 'round*' num2str(rridx) '*ch0' num2str(chidx) '*.tif']); 
        img = loadtiff([parentfolder fnames(1).name]);
        
        if strcmp(params.normalization,'minmax')
            img = mat2gray(img);
        end
        
        imbin = binarize_intensity_threshold(img,params);

        closingrad = params.closingrad*(1/xystep); %radius for morphological closing, in nm
        if params.morph_close
                %I_c=imclose(I_f0,strel('sphere',closingrad)); 
            imclosed=imclose(imbin,strel('disk',floor(closingrad))); %morphological closing
        else
            imclosed=imbin; %no morphological closing
        end
        
        %apply lower filter on both after closing
        %set lower limit on ROI size to minimize noise: .1 um in each dimension
        lowerlim = ceil((params.lowerlim^3)*(1/xystep)*(1/xystep)*(1/zstep));
        mask = bwareaopen(imclosed, lowerlim); %filtration, removes binary objects LESS than this size

        CC = bwconncomp(mask,26);%find all connected components in the filtered
        nobjects = CC.NumObjects;

        cc_volumes = regionprops3(CC,'Volume');
        cc_volumes=cc_volumes.Volume;
        
        nobjects(chidx,rridx) = nobjects;
        volume_temp = nnz(mask(:));
        volume(chidx,rridx) = volume_temp * params.vol_converter;
        volume_frac(chidx,rridx) = volume_temp/numel(mask);
        min_vol(chidx,rridx) = min(cc_volumes)/volume_temp;
        max_vol(chidx,rridx) = max(cc_volumes)/volume_temp;

    end
end

end


