minfunction [SNR_overall,SNR_synapses,nsynapses] = analyze_mExR_validation_strip(fov,params)
% Analyzes validation data, with same proteins stained in consecutive rounds
% Outputs 
% SNR: the overall image SNR in each round for each channel
% SNR_synapses: the SNR for each synapse over each round
% nsynapses: the number of synapses in each round


%Extract parameters
xystep = params.xystep;
zstep = params.zstep; % um/voxel in z
nchannels = params.nchannels;
parentfolder = params.parentfolder;
normalization = params.normalization;
subtract_morph = params.subtract_morph;%subtract the morphology channel? no, because we are using it for reg quality analysis
morph_channel = params.morph_channel; %morphology channel used for registration, for subtracting out
rounds = params.rounds;
nrounds = length(rounds);
morphclose = params.morph_close;

%load and binarize the morphology channel
fnames = dir([parentfolder '*' fov '*' 'round001*ch0' morph_channel '*.tif']); %binary images
morph_img = loadtiff([parentfolder fnames(1).name]);
if strcmp(params.normalization,'minmax')
    morph_img = mat2gray(morph_img);
end
morph_bin = binarize_intensity_threshold(morph_img,params);
if strcmp(params.filt,'med')
    morphbin_filt = medfilt3(morph_bin,params.medfilt);
end

closingrad = params.closingrad*(1/xystep); %radius for morphological closing, in nm
if morphclose
        %I_c=imclose(I_f0,strel('sphere',closingrad)); 
    morphclosed=imclose(morphbin_filt,strel('disk',floor(closingrad))); %morphological closing
else
    morphclosed=morphbin_filt; %no morphological closing
end

%apply lower filter on both after closing
%set lower limit on ROI size to minimize noise: .1 um in each dimension
lowerlim = ceil((params.lowerlim^3)*(1/xystep)*(1/xystep)*(1/zstep));
morph_f = bwareaopen(morphclosed, lowerlim); %filtration, removes binary objects LESS than this size

%loop through each channel
for rridx = 1:nrounds
    for chidx = 1:nchannels
        if ~strcmp(num2str(chidx),morph_channel)
            
            round = rounds{rridx};
            fnames = dir([parentfolder '*' fov '*' 'round*' round '_ch0' num2str(chidx) '*.tif']); %binary images
            img = loadtiff([parentfolder fnames(1).name]);
            
            if strcmp(params.normalization,'minmax')
                img = mat2gray(img);
            end
            
            imbin = binarize_intensity_threshold(img,params);
            
            if strcmp(params.filt,'med')
                imbin_filt = medfilt3(imbin,params.medfilt);
            end
            
            if subtract_morph
                imbin_filt = double(imbin_filt) - double(morph_f);
                imbin_filt(imbin_filt<0) = 0;
            end

            closingrad = params.closingrad*(1/xystep); %radius for morphological closing, in nm
            if morphclose
                    %I_c=imclose(I_f0,strel('sphere',closingrad)); 
                imclosed=imclose(imbin_filt,strel('disk',floor(closingrad))); %morphological closing
            else
                imclosed=imbin_filt; %no morphological closing
            end

            %apply lower filter on both after closing
            %set lower limit on ROI size to minimize noise: .1 um in each dimension
            lowerlim = ceil((params.lowerlim^3)*(1/xystep)*(1/xystep)*(1/zstep));
            mask = bwareaopen(imclosed, lowerlim); %filtration, removes binary objects LESS than this size

            CC_signal = bwconncomp(mask,26);%find all connected components in the combined image
            nobjects = CC_signal.NumObjects;
            
            %dilate the image and create a "noise" mask
            se = strel('disk',10);
            dilated = imdilate(mask,se);
            %noise_mask = double(dilated) - double(mask);
            CC_dilated = bwconncomp(dilated,26);
            nobjects_dilated = CC_dilated.NumObjects;
          
%             if nobjects>nobjects_noise
%                 disp(['unequal # of objects between masks'])
%                 disp(fnames(1).name)
%                 nobjects = nobjects_noise;
%             end
            
            if nobjects_dilated > 0
                for oidx = 1:nobjects_dilated %loop through all objects and calculate their SNR
                    pixels = CC_dilated.PixelIdxList{oidx};
                    chunk = img(pixels);
                    bin = binarize_intensity_threshold(chunk,params);
                    inverted = 1-bin;
                    masked_chunk = bin.*chunk;
                    masked_chunk_inv = inverted.*chunk;
                    signal = mean(masked_chunk(masked_chunk>0));
                    noise = mean(masked_chunk_inv(masked_chunk_inv>0));
                    SNR_synapses_temp(oidx,1)=signal/noise;
                end
            else
                SNR_synapses_temp = NaN;
            end

            masked_image = img .* double(mask);
            inverted_mask = 1-double(mask);
            bg_image = img .* inverted_mask;

            SNR = mean(masked_image(masked_image>0))/mean(bg_image(bg_image>0));

            %save down the segmentation as tiff stack in same folder
            imtosave1 = uint8(mask);
            name=fnames(1).name(1:end-4);
            imwrite(imtosave1(:,:,1),[parentfolder 'masks/Seg-' name '.tif'],'tiff')
            for jj = 2:size(imtosave1,3)
                imwrite(imtosave1(:,:,jj),[parentfolder 'masks/Seg-' name '.tif'],'WriteMode','append')
            end

            SNR_overall(rridx,chidx) = SNR;
            nsynapses(rridx,chidx) = nobjects;
            SNR_synapses{rridx,chidx} = SNR_synapses_temp;
            clear SNR_synapses_temp
            
        end
    end
end

end


