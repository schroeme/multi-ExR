function [data] = analyze_synaptic_puncta_cropped(params,fov,chs)
% Analyzes validation data, with same proteins stained in consecutive rounds
% Outputs 
% num_puncta: the number of puncta in each synapse
% punctavol: the total puncta volume for each synapse
% AR: the aspect ratio of the most ellipsoid volume in each synapse
% SNR: signal to noise ratio in synaptic ROI

%Extract parameters
xystep = params.xystep;
zstep = params.zstep; % um/voxel in z
nchannels = length(chs);
parentfolder = params.parentfolder;

%loop through each synaptic chanel and each synapse
for chidx = 1:nchannels
    chsplits = split(chs{chidx},"-");
    roundno = chsplits{1};
    chno = chsplits{2};
    fnames = dir([parentfolder filesep '*' fov '*round0' roundno '*ch0' chno '*.tif']);

    SNR_temp=[];
    npuncta_temp=[];
    vol_temp=[];
    AR_temp=[];
    
    nsynapses = length(fnames);
    for ssidx = 1:nsynapses
        img = loadtiff([parentfolder fnames(ssidx).name]);
        img = double(img);
        %img = mat2gray(img);
       
        if strcmp(params.thresh_method,'absolute')
            params.threshold=params.thresholds(chidx);
        end
        imbin = binarize_intensity_threshold(img,params);

        if strcmp(params.filt,'med')
            imbin_filt = medfilt3(imbin,params.filt_size);
        else
            imbin_filt = imbin;
        end

%         closingrad = params.syn_rad; %radius for morphological closing, in nm
%         if params.morph_close_syn
%                 %I_c=imclose(I_f0,strel('sphere',closingrad)); 
%             imclosed=imclose(imbin_filt,strel('disk',floor(closingrad))); %morphological closing
%         else
%             imclosed=imbin_filt; %no morphological closing
%         end
% 
        %apply lower filter on both after closing
        if params.sizefilt
            mask = bwareaopen(imbin_filt, params.lowerlim); %filtration, removes binary objects LESS than this size
        else
            mask = imbin_filt;
        end
        
        if params.doplot
            figure;
            subplot(1,2,1)
            imagesc(max(img,[],3))
            title(fnames(ssidx).name)
            subplot(1,2,2)
            imagesc(max(mask,[],3))
            title(fnames(ssidx).name)
        end
        CC_signal = bwconncomp(mask,26);%find all connected components in the combined image
        nobjects = CC_signal.NumObjects;
        
        if nobjects > 0 

            masked_image = img .* double(mask);
            mean_int(ssidx,1) = sum(masked_image(:))/nnz(masked_image(:));
            % inverted_mask = 1-double(mask);
            % bg_image = img .* inverted_mask;

%             SNRval = mean(masked_image(masked_image>0))/mean(bg_image(bg_image>0));

            voltab = regionprops3(CC_signal,'Volume');
            vols = voltab.Volume;
            vols = vols * params.vol_converter;

            paxtab = regionprops3(CC_signal,'PrincipalAxisLength');
            pax = paxtab.PrincipalAxisLength;
            ARs = pax(:,1)./pax(:,2);

%             SNR_temp(ssidx,1) = SNRval;
            npuncta_temp(ssidx,1) = nobjects;
            vol_temp(ssidx,1) = sum(vols);
            AR_temp(ssidx,1)= max(ARs);
        else
            
%             SNR_temp(ssidx,1) = 0;
            npuncta_temp(ssidx,1) = nobjects;
            vol_temp(ssidx,1) = 0;
            AR_temp(ssidx,1)= 0;
            mean_int(ssidx,1) =0;
        end

        roi_filename{ssidx,1} = fnames(ssidx).name;

    end
    
    data.num_puncta{chidx} = npuncta_temp;
    data.punctavol{chidx} = vol_temp;
    data.roi_filename{chidx} = roi_filename;
    data.AR{chidx} = AR_temp;
    data.meanint{chidx} = mean_int;

end
    
end
