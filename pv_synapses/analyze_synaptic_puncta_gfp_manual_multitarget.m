function [data] = analyze_synaptic_puncta_gfp_manual_multitarget(params)

%Extract parameters
maskfolder = params.outfolder;
targetfolder = params.targetoutfolder;
parentfolder = params.parentfolder;

%Get filenames for existing masks
gfpp_segs = dir([maskfolder filesep '*bbox_seg_GFPpos.tif']); %bounding boxes
all_masks = dir([maskfolder filesep '*affine_seg_all.tif']); %not bounding boxes

nfovs = length(gfpp_segs);

nprots = length(params.target_channels);
nsyn_estimate = sum(params.nsyn_all);

npunctas_all = zeros(nsyn_estimate,nprots);
vols_all = zeros(nsyn_estimate,nprots);
meanints_all = zeros(nsyn_estimate,nprots);
names = cell(nsyn_estimate,1);
meanints_refnorm_all = zeros(nsyn_estimate,nprots);
selected_tf = zeros(nsyn_estimate,1);
counter = 0;

%loop through each field of view in the folder
for fidx = 1:nfovs
    %load up all syanpse and GFP masks
    gfpp_seg = mat2gray(loadtiff([maskfolder gfpp_segs(fidx).name]));
    all_seg = mat2gray(loadtiff([maskfolder all_masks(fidx).name]));
    
    namesplits = strsplit(all_masks(fidx).name,"_");
    roiname = namesplits{1};

    CC_all = bwconncomp(all_seg,26);
    
    %get bounding boxes for all synapses
    bboxes_all = regionprops3(CC_all,'BoundingBox');
    nobjects_all = CC_all.NumObjects;
    nsyn_all(fidx,1) = nobjects_all;

    buffer_half = params.nbuffer/2;
    zbuffer_half = 2;

    if fidx > 1
        counter = sum(nsyn_all(1:fidx-1));
    end

    for pp = 1:nprots
        %load up background-subtracted image (NOT masked, we need this for mean
        %intensity)
        chid = params.target_channels{pp};
        chid_splits = strsplit(chid, '-');

        fname_search = [roiname '_round00' chid_splits{1} '_ch0' chid_splits{2}];
        fnames_target = dir([parentfolder filesep fname_search '*.tif']);
        ch2img = double(loadtiff([parentfolder fnames_target(1).name]));

        pname = params.target_names{pp};
        target_masks = dir([targetfolder filesep '*' pname '*.tif']);
        target_mask = double(loadtiff([targetfolder target_masks(fidx).name]));
    
        %Next, quantify signal within bounding boxes
        for ss = 1:nobjects_all
            sidx = counter + ss;
            bbox = bboxes_all(ss,1).BoundingBox;
            ulf_x = bbox(1);
            ulf_y = bbox(2);
            ulf_z = bbox(3);
            width_x = bbox(4);
            width_y = bbox(5);
            width_z = bbox(6);
    
            %Note: somehow these bounding boxes are larger than 
            % Calculate the start and end indices for each dimension
            xStart = max(1, floor(ulf_x) - buffer_half);
            xEnd = min(size(target_mask,2), floor(ulf_x + width_x) + buffer_half);
            yStart = max(1, floor(ulf_y) - buffer_half);
            yEnd = min(size(target_mask,1), floor(ulf_y + width_y) + buffer_half);
            zStart = max(1, floor(ulf_z) - zbuffer_half);
            zEnd = min(size(target_mask,3), floor(ulf_z + width_z) + zbuffer_half);
    
            %Grab the bounding box area in target channel
            ch2chunk = ch2img(yStart:yEnd, xStart:xEnd, zStart:zEnd);

            %use previously generated mask for this
            chunk_filt = target_mask(yStart:yEnd, xStart:xEnd, zStart:zEnd);
    
            if pp == 1
                %Grab the bounding box area in the reference channel, binarized
                refchunk = all_seg(yStart:yEnd, xStart:xEnd, zStart:zEnd);
    
                %grab the corresponding area in the "selected" GFP+ mask
                selected_filt = gfpp_seg(yStart:yEnd, xStart:xEnd, zStart:zEnd);
        
                if mean(selected_filt(:)) ==1
                    selected_tf(sidx,1) = 1;
                else
                    selected_tf(sidx,1)=0;
                end
            end
    
            ch2masked = ch2chunk .* chunk_filt;
            CC2 = bwconncomp(ch2masked,26);%find all connected components in the combined image
            n2 = CC2.NumObjects;
    
            if n2 > 0
                mean_int = sum(ch2masked(:))/nnz(ch2masked(:));
                mean_int_refnorm = sum(ch2masked(:))/nnz(refchunk(:));
                voltab=regionprops3(CC2,'Volume');
                volstemp=sum(voltab.Volume,1);
                vols = volstemp * params.vol_converter;
                npuncta = n2;
            else
                mean_int = 0;
                vols = 0;
                npuncta = n2;
                mean_int_refnorm = 0;
            end
    
            name = [roiname '_' pname '_syn' num2str(ss)];
            names{sidx} = name;
            meanints_all(sidx,pp) = mean_int;
            npunctas_all(sidx,pp) = npuncta;
            vols_all(sidx,pp) = vols;
            meanints_refnorm_all(sidx,pp) = mean_int_refnorm;
            if mod(sidx,10000) == 0
                disp(['Analyzed ' num2str(sidx) ' synapses!']);
            end
    
        end
    end

end
    
data.names = names;
data.total_vol = vols_all;
data.mean_int = meanints_all;
data.npuncta = npunctas_all;
data.selected_tf = selected_tf;
data.mean_int_refnorm = meanints_refnorm_all;
data.nsyn_all = nsyn_all;
