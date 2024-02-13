function [SNR,num_puncta,punctavol,punctaint,nsynapses] = analyze_mExR_validation_cropped(fov,params)
% Analyzes validation data, with same proteins stained in consecutive rounds
% Outputs 
% SNR: the overall image SNR in each round for each channel
% SNR_synapses: the SNR for each synapse over each round
% nsynapses: the number of synapses in each round
% punctaint: the mean intensity within each puncta (a.u)


%Extract parameters
xystep = params.xystep;
zstep = params.zstep; % um/voxel in z
channels = params.channels;
nchannels = length(channels);
parentfolder = params.parentfolder;
normalization = params.normalization;
rounds = params.rounds;
nrounds = length(rounds);
morph_channel=params.morph_channel;
morphclose =params.morph_close;

%loop through each channel, each round, and each synapse
for rridx = 1:nrounds
    for chidx = 1:nchannels
        if ~strcmp(num2str(chidx),morph_channel)
            
            chname = channels{chidx};
            round = rounds{rridx};
            fnames = dir([parentfolder '*' fov '*' 'round*' round '_' chname '*.tif']);
            
            nsynapses = length(fnames);
            for ssidx = 1:nsynapses
                img = loadtiff([parentfolder fnames(ssidx).name]);

                if strcmp(params.normalization,'minmax')
                    img = mat2gray(img);
                elseif strcmp(params.normalization,'none')
                    img = double(img);
                end

                imbin = binarize_intensity_threshold(img,params);

                if strcmp(params.filt,'med')
                    imbin_filt = medfilt3(imbin,params.medfilt);
                else
                    imbin_filt=imbin;
                end

                closingrad = params.syn_closingrad*(1/xystep); %radius for morphological closing, in nm
                if morphclose
                        %I_c=imclose(I_f0,strel('sphere',closingrad)); 
                    imclosed=imclose(imbin_filt,strel('disk',floor(closingrad))); %morphological closing
                else
                    imclosed=imbin_filt; %no morphological closing
                end

                %apply lower filter on both after closing
                %set lower limit on ROI size to minimize noise
                lowerlim = ceil((params.lowerlim^3)*(1/xystep)*(1/xystep)*(1/zstep));
                mask = bwareaopen(imclosed, lowerlim); %filtration, removes binary objects LESS than this size

                CC_signal = bwconncomp(mask,26);%find all connected components in the combined image
                nobjects = CC_signal.NumObjects;

                masked_image = img .* double(mask);
                inverted_mask = 1-double(mask);
                bg_image = img .* inverted_mask;

                SNRval = mean(masked_image(masked_image>0))/mean(bg_image(bg_image>0));
                int = mean(masked_image(masked_image>0)); %mean within nonzero pixels of the mask
                
                voltab = regionprops3(CC_signal,'Volume');
                vols = voltab.Volume;

                %create some figures periodically
                if ((params.doplot) && (mod(ssidx,45) == 0))
                    figure();
                    subplot(2,1,1)
                    imagesc(max(img,[],3))
                    title('MIP')
                    subplot(2,1,2)
                    imagesc(max(mask,[],3))
                    title('Mask - MIP')
                end
                %save down the segmentation as tiff stack in same folder
                if params.savechunks
                    imtosave1 = uint8(mask);
                    name=fnames(1).name(1:end-4);
                    imwrite(imtosave1(:,:,1),[parentfolder 'masks/Seg-' name '.tif'],'tiff')
                    for jj = 2:size(imtosave1,3)
                        imwrite(imtosave1(:,:,jj),[parentfolder 'masks/Seg-' name '.tif'],'WriteMode','append')
                    end
                end
                
                SNR_temp(ssidx,1) = SNRval;
                npuncta_temp(ssidx,1) = nobjects;
                vol_temp(ssidx,1) = mean(vols) * params.vol_converter;
                int_temp(ssidx,1) = int;
            
            end
        
            SNR(rridx,chidx) = nanmean(SNR_temp);
            num_puncta(rridx,chidx) = nanmean(npuncta_temp);
            punctavol(rridx,chidx) = nanmean(vol_temp);
            punctaint(rridx,chidx) = nanmean(int_temp);
            clear SNR_temp npuncta_temp vol_temp int_temp
        end
    end
end

end


