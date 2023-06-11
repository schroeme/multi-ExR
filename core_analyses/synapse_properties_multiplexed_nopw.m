function [imdata] = synapse_properties_multiplexed_nopw(params,fovid,nsynapses)
% This function extracts multiplexed synapse features (e.g. volume, number
% of puncta for each channel WITHOUT pairwise data (much faster) from cropped synaptic ROIs with multiple channels
% 

%Extract parameters
parent = params.parentfolder;
xystep = params.xystep;
zstep = params.zstep;
disp(fovid);
nchannels = length(params.channels_rounds);

npuncta = zeros(nsynapses,nchannels);
frac_vol_occ = zeros(nsynapses,nchannels);

mean_puncta_vol = zeros(nsynapses,nchannels);
med_puncta_vol = zeros(nsynapses,nchannels);
std_puncta_vol = zeros(nsynapses,nchannels);

mean_puncta_SA = zeros(nsynapses,nchannels);
med_puncta_SA = zeros(nsynapses,nchannels);
std_puncta_SA = zeros(nsynapses,nchannels);
synapse_size = zeros(nsynapses);

for ch_idx = 1:nchannels %loop through each channel
    channel_str = params.channels_rounds{ch_idx};
    splits = strsplit(channel_str,'-');
    roundstr = splits{1};
    chstr = splits{2};
            
    fnames_channel = dir([parent fovid '*round*' roundstr '_ch0' chstr '*.tif']);
%     fnames_round1_ref = dir([parent fovid '*round001_ch04*.tif']);
%     fnames_roundn_ref = dir([parent fovid  '*round*' roundstr '_ch04*.tif']);

%     fnames_thresh = dir([params.threshfolder fovid '*round*' roundstr '_ch0' chstr '*.tif']);%filename for loading in image volume (whole field of view)
%     I_whole = loadtiff([params.threshfolder fnames_thresh(1).name]);
%     I_whole = mat2gray(I_whole);
%     %filter the whole image to eliminate high spatial frequency noise
%     %before calculating the threshold
%     I_whole_filt = medfilt3(I_whole,[9 9 3]);

%     if strcmp(params.thresh_method,'pct')
%         threshold = prctile(I_whole_filt(:),params.thresh_val);
%         I_whole = imbinarize(I_whole,threshold);
%     elseif strcmp(params.thresh_method,'otsu')
%         threshold = graythresh(I_whole_filt(:));
%         threshold = threshold*params.thresh_multiplier;
%         I_whole = imbinarize(I_whole,threshold);
%     end

    for ss_idx = 1:nsynapses %loop through each synapse
        if params.do_binarize

%             fname_round1_ref = fnames_round1_ref(ss_idx).name; %round 1 reference channel, for alignment
%             fname_roundn_ref = fnames_roundn_ref(ss_idx).name;
            fname = fnames_channel(ss_idx).name;

            if length(fname) == 0
                disp(fovid)
                disp(channel_str)
            end
            I = loadtiff([parent fname]);
            I = mat2gray(I);
            Iorig = I;

            if strcmp(params.thresh_method,'pct')
                threshold = prctile(I(:),params.thresh_val);
            elseif strcmp(params.thresh_method,'otsu')
                threshold = graythresh(I(:));
                threshold = threshold*params.thresh_multiplier;           
            elseif strcmp(params.thresh_method,'zscore')
                meanint = mean(I(:));
                stdint = std(I(:));
                threshold = meanint + params.thresh_multiplier * stdint;
            end
            I = imbinarize(I,threshold);
            I = double(I);

            if params.do_med_filt
                I_filt = medfilt3(I,params.filt_size);
            else
                I_filt = I;
            end
        else
            fname = dir([parent filesep 'chunks/Bin-*' fovid '*ch' num2str(ch_idx) '_chunk' num2str(ss_idx) '_*.tif']); %binary from test stain
            I = loadtiff([parent 'chunks/' fname.name]);
            I = mat2gray(I);
        end
        
        I_filt = bwareaopen(I_filt,params.minsize);
        CC_check = bwconncomp(I_filt,26);
        nobjects_check = CC_check.NumObjects;

        %how much of the image is nonzero? should be low for true signal
        frac_nonzero = nnz(I_filt)/numel(I_filt);

        if params.doplot && ss_idx == 10 %only print one synapse 
            figure();
            subplot(1,3,1)
            imagesc(max(Iorig,[],3))
            title([channel_str ' - Original'])
            subplot(1,3,2)
            imagesc(max(I,[],3))
            title('Thresholded')
            subplot(1,3,3)
            imagesc(max(I_filt,[],3))
            title('Filtered')
        end

        if (nobjects_check > 0) && (frac_nonzero < 0.65)
            CC = bwconncomp(I,26);%find all connected components in the binary, filtered image
            npuncta(ss_idx,ch_idx) = CC.NumObjects;
            
            frac_vol_occ(ss_idx,ch_idx) = nnz(I_filt)/numel(I_filt); %fraction of volume occupied
            if ch_idx == 1
                synapse_size(ss_idx) = numel(I_filt);
            end
            
            puncta_vols = regionprops3(CC,'Volume');
            puncta_vols = puncta_vols.Volume;
            puncta_SAs = regionprops3(CC,'SurfaceArea');
            puncta_SAs = puncta_SAs.SurfaceArea;
            
            mean_puncta_vol(ss_idx,ch_idx) = mean(puncta_vols);
            med_puncta_vol(ss_idx,ch_idx) = median(puncta_vols);
            std_puncta_vol(ss_idx,ch_idx) = std(puncta_vols);%./length(puncta_vols);
            
            mean_puncta_SA(ss_idx,ch_idx) = mean(puncta_SAs);
            med_puncta_SA(ss_idx,ch_idx) = median(puncta_SAs);
            std_puncta_SA(ss_idx,ch_idx) = std(puncta_SAs);%./length(puncta_SAs);

            %add in aspect ratio
            paxtab = regionprops3(CC,'PrincipalAxisLength');
            pax = paxtab.PrincipalAxisLength;
            if size(pax,2) > 1
                ARs = pax(:,1)./pax(:,2);
            else
                ARs = NaN;
            end
            max_puncta_AR(ss_idx,ch_idx) = max(ARs);
            
        else
            npuncta(ss_idx,ch_idx)=NaN;
            frac_vol_occ(ss_idx,ch_idx) =NaN;
            
            mean_puncta_vol(ss_idx,ch_idx) = NaN;
            med_puncta_vol(ss_idx,ch_idx) = NaN;
            std_puncta_vol(ss_idx,ch_idx) = NaN;
            
            mean_puncta_SA(ss_idx,ch_idx) = NaN;
            med_puncta_SA(ss_idx,ch_idx) = NaN;
            std_puncta_SA(ss_idx,ch_idx) = NaN;
            max_puncta_AR(ss_idx,ch_idx) = NaN;
                      
        end
    end
end

imdata.npuncta = npuncta;
imdata.size = synapse_size;
imdata.frac_vol_occ = frac_vol_occ;
imdata.mean_puncta_vol = mean_puncta_vol;
imdata.med_puncta_vol = med_puncta_vol;
imdata.std_puncta_vol = std_puncta_vol;

imdata.mean_puncta_SA = mean_puncta_SA;
imdata.med_puncta_SA = med_puncta_SA;
imdata.std_puncta_SA = std_puncta_SA;

imdata.max_puncta_AR=max_puncta_AR;

end

