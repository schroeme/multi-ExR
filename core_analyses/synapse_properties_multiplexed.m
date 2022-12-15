function [imdata,syndata,pwdata] = synapse_properties_multiplexed(params,fovid,nsynapses)
% This function extracts multiplexed synapse features (e.g. volume, number
% of puncta for each channel, pairwise overlaps, pixel-wise correlations,
% etc.) from cropped synaptic ROIs with multiple channels
% 

%Extract parameters
parent = params.parentfolder;
xystep = params.xystep;
zstep = params.zstep;
%ch_names = params.ch_names;
ctrlshift_xy = floor(params.controlshift/xystep);
ctrlshift_z = 0;%floor(params.controlshift/zstep);

fnames_all = dir([parent '*' fovid '*.tif']);
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
    fnames_round1_ref = dir([parent fovid '*round001_ch04*.tif']);
    fnames_roundn_ref = dir([parent fovid  '*round*' roundstr '_ch04*.tif']);

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

            fname_round1_ref = fnames_round1_ref(ss_idx).name; %round 1 reference channel, for alignment
            fname_roundn_ref = fnames_roundn_ref(ss_idx).name;
            fname = fnames_channel(ss_idx).name;

            if length(fname) == 0
                disp(fovid)
                disp(channel_str)
            end
            I = loadtiff([parent fname]);
            I = mat2gray(I);

            Iref_round1 = loadtiff([parent fname_round1_ref]);
            Iref_round1 = mat2gray(Iref_round1);

            Iref_roundn = loadtiff([parent fname_roundn_ref]);
            Iref_roundn = mat2gray(Iref_roundn);

            if (params.doreg) && (~strcmp(roundstr,'1')) %perform an additional registration to overlap to round 1, based on reference channel
                [xyzshift,I_shifted] = get_xyzshift(Iref_round1,I, I, params.pixel_size, params.rmax, [], params.distance);

%                 [optimizer,metric] = imregconfig('multimodal');
% 
%                 tform = imregtform(Iref_roundn,Iref_round1,'affine',optimizer,metric,'PyramidLevels',2);
%                 ref_reg = imwarp(Iref_roundn,tform,'OutputView',imref3d(size(Iref_roundn)));
%                 I_shifted = imwarp(I,tform,'OutputView',imref3d(size(I)));
                I = I_shifted;
                if params.doplot
                    figure();
                    subplot(1,3,1)
                    imagesc(max(Iref_round1,[],3))
                    title('Reference, round 1')
                    subplot(1,3,2)
                    imagesc(max(Iref_roundn,[],3))
                    title('Reference, round n')
                    subplot(1,3,3)
                    %ref_shifted = imtranslate(Iref_roundn,-1*xyzshift);
                    imagesc(max(ref_reg,[],3))
                    title('Shifted channel, round n')
                end
                
            end

            %I = imbinarize(I,threshold); %apply whole field of view threshold calculated above

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

        CC_check = bwconncomp(I_filt,26);
        nobjects_check = CC_check.NumObjects;

        %how much of the image is nonzero? should be low for true signal
        frac_nonzero = nnz(I_filt)/numel(I_filt);

        if (nobjects_check > 0) && (frac_nonzero < 0.65)
            CC = bwconncomp(I,26);%find all connected components in the binary, filtered image
            npuncta(ss_idx,ch_idx) = CC.NumObjects;
            
            frac_vol_occ(ss_idx,ch_idx) = nnz(I)/numel(I); %fraction of volume occupied
            if ch_idx == 1
                synapse_size(ss_idx) = numel(I);
            end
            
            puncta_vols = regionprops3(CC,'Volume');
            puncta_vols = puncta_vols.Volume;
            puncta_SAs = regionprops3(CC,'SurfaceArea');
            puncta_SAs = puncta_SAs.SurfaceArea;
            puncta_idx = regionprops3(CC,'VoxelIdxList');
            puncta_idx = puncta_idx.VoxelIdxList;
            puncta_centroids = regionprops3(CC,'Centroid');
            puncta_centroids = puncta_centroids.Centroid;
            
            mean_puncta_vol(ss_idx,ch_idx) = mean(puncta_vols);
            med_puncta_vol(ss_idx,ch_idx) = median(puncta_vols);
            std_puncta_vol(ss_idx,ch_idx) = std(puncta_vols);%./length(puncta_vols);
            
            mean_puncta_SA(ss_idx,ch_idx) = mean(puncta_SAs);
            med_puncta_SA(ss_idx,ch_idx) = median(puncta_SAs);
            std_puncta_SA(ss_idx,ch_idx) = std(puncta_SAs);%./length(puncta_SAs);
            
            puncta_centroids_all{ss_idx,ch_idx} = puncta_centroids;
            puncta_idx_all{ss_idx,ch_idx} = puncta_idx; 
            
            Itrans = imtranslate(I,[ctrlshift_xy 0 0]);
            Isyn{ch_idx}=I;
            Isyn_trans{ch_idx}=Itrans;
        else
            npuncta(ss_idx,ch_idx) =0;
            
            frac_vol_occ(ss_idx,ch_idx) = 0; %fraction of volume occupied
            
            puncta_vols = 0;
            puncta_SAs = 0;
            puncta_idx = NaN;
            puncta_centroids = NaN;
            
            mean_puncta_vol(ss_idx,ch_idx) = mean(puncta_vols);
            med_puncta_vol(ss_idx,ch_idx) = median(puncta_vols);
            std_puncta_vol(ss_idx,ch_idx) = std(puncta_vols);%./length(puncta_vols);
            
            mean_puncta_SA(ss_idx,ch_idx) = mean(puncta_SAs);
            med_puncta_SA(ss_idx,ch_idx) = median(puncta_SAs);
            std_puncta_SA(ss_idx,ch_idx) = std(puncta_SAs);%./length(puncta_SAs);
            
            puncta_centroids_all{ss_idx,ch_idx} = puncta_centroids;
            puncta_idx_all{ss_idx,ch_idx} = puncta_idx; 
            
            Itrans = imtranslate(I,[ctrlshift_xy 0 0]);
            Isyn{ch_idx}=I;
            Isyn_trans{ch_idx}=I;            
        end
    end
end

for ss_idx = 1:nsynapses
    
    % METRICS AT THE SYNAPSE LEVEL
    % calculate pairwise correlation coefficients and overlaps between all channels
    correlations = zeros(nchannels,nchannels);
    correlations_ctrl = zeros(nchannels,nchannels);

    overlaps = zeros(nchannels,nchannels);
    overlaps_ctrl = zeros(nchannels,nchannels);
    
    mean_frac_overlap = zeros(nchannels,nchannels);
    med_frac_overlap = zeros(nchannels,nchannels);
    std_frac_overlap = zeros(nchannels,nchannels);
    
    mean_distance = zeros(nchannels,nchannels);
    med_distance = zeros(nchannels,nchannels);
    max_distance = zeros(nchannels,nchannels);
    min_distance = zeros(nchannels,nchannels);
    std_distance = zeros(nchannels,nchannels);

    for c1 = 1:nchannels
        for c2 = c1:nchannels
            %3D correlation
            correlations(c1,c2) = corr(Isyn{c1}(:),Isyn{c2}(:));
            
            %volume overlap
            intersection = Isyn{c1} & Isyn{c2}; %min(Isyn{c1},Isyn{c2});
            n_intersection = nnz(intersection);
            
            %scale to total volume of pixel
            overlaps(c1, c2) = (n_intersection*2)/(nnz(Isyn{c1}) + nnz(Isyn{c2}));
                    
            %union = max(Isyn{c1},Isyn{c2});
            %overlaps(c1,c2) = sum(intersection(:))./sum(union(:));

            %do for shifted control
            translated = Isyn_trans{c2}(:,1+ctrlshift_xy:end,:);
            original = Isyn{c1}(:,1+ctrlshift_xy:end,:);
            
            %3D correlation
            correlations_ctrl(c1,c2) = corr(original(:),translated(:));
            
            %volume overlap
            intersection_trans = translated & original; %min(original,translated);
            n_intersection_trans = nnz(intersection_trans);
            %union_trans = max(original,translated);
            overlaps_ctrl(c1,c2) = (n_intersection_trans*2)/(nnz(original)+nnz(translated));%sum(intersection_trans(:))./sum(union_trans(:));

            % METRICS AT THE PUNCTA LEVEL
            npuncta1 = npuncta(ss_idx,c1);
            npuncta2 = npuncta(ss_idx,c2);
            
            %loop through each pair of puncta
            frac_overlap = ones(npuncta1,npuncta2)*-1;
            distances = ones(npuncta1,npuncta2)*-1;
            for p1 = 1:npuncta1
                for p2 = 1:npuncta2
                    
                    %fraction overlap (as percent of combined puncta volume)
                    if iscell(puncta_idx_all{ss_idx,c1})
                        idx1 = puncta_idx_all{ss_idx,c1}{p1};
                    else
                        idx1 =  puncta_idx_all{ss_idx,c1}(p1);
                    end
                    
                    if iscell(puncta_idx_all{ss_idx,c2})
                        idx2 = puncta_idx_all{ss_idx,c2}{p2};
                    else
                        idx2 = puncta_idx_all{ss_idx,c2}(p2);
                    end
                    overlap = length(intersect(idx1,idx2));
                    frac_overlap(p1,p2) = (overlap*2) / (nnz(idx1) + nnz(idx2));
                    
                    %inter-puncta distance
                    cent1 = puncta_centroids_all{ss_idx,c1}(p1,:);
                    cent2 = puncta_centroids_all{ss_idx,c2}(p2,:);
                    distances(p1,p2) = norm(cent1-cent2);
                end
            end
            
            distances = distances(distances>=0);
            frac_overlap = frac_overlap(frac_overlap>=0);
            
            mean_frac_overlap(c1,c2) = mean(frac_overlap);
            med_frac_overlap(c1,c2) = median(frac_overlap);
            std_frac_overlap(c1,c2) = std(frac_overlap);%/length(frac_overlap);

            mean_distance(c1,c2) = mean(distances);
            med_distance(c1,c2) = median(distances);
            try
                max_distance(c1,c2) = max(distances);
            catch
                disp(distances);
                max_distance(c1,c2) = nan;
            end
            try
                min_distance(c1,c2) = min(distances);
            catch
                disp(distances);
                min_distance(c1,c2) = nan;
            end
            std_distance(c1,c2) = std(distances)/length(distances);
            
        end
    end

    syndata(ss_idx).correlations = correlations; % ./ correlations_ctrl;
    syndata(ss_idx).overlaps = overlaps; %./ overlaps_ctrl;
    syndata(ss_idx).correlations_ctrl = correlations_ctrl;
    syndata(ss_idx).overlaps_ctrl = overlaps_ctrl;
    %syndata(ss_idx).correlations_ctrl = correlations_ctrl;
    %syndata(ss_idx).overlaps_ctrl = overlaps_ctrl;
    
    pwdata(ss_idx).mean_frac_overlap = mean_frac_overlap + mean_frac_overlap';
    pwdata(ss_idx).med_frac_overlap = med_frac_overlap + med_frac_overlap';
    pwdata(ss_idx).std_frac_overlap = std_frac_overlap + std_frac_overlap;
    
    pwdata(ss_idx).mean_distance = mean_distance + mean_distance';
    pwdata(ss_idx).med_distance = med_distance + med_distance';
    pwdata(ss_idx).max_distance = max_distance + max_distance';
    pwdata(ss_idx).min_distance = min_distance + min_distance';
    pwdata(ss_idx).std_distance = std_distance + std_distance';
    
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

end

