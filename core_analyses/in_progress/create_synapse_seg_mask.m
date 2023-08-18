function I_s = create_synapse_seg_mask(CC_combo,buffer,upperlim,xystep,zstep,imsize)
%Creates mask for segmenting synapses by adding bounding boxes around connected components.
%Excludes objects over the size limit

vol = regionprops3(CC_combo, 'Volume');
% sortedVol = sort([vol.Volume], 'descend');
% tokeep = sortedVol<upperlim;
% nobjects = sum(tokeep);
% CC_combo.PixelIdxList = CC_combo.PixelIdxList(tokeep);
% CC_combo.NumObjects = nobjects;
vol = vol.Volume;

%add a buffer of 100nm around each synapse
bufferx = ceil(buffer*(1/xystep));
buffery = ceil(buffer*(1/xystep));
bufferz = floor(buffer*(1/zstep));

I_s = zeros(imsize); %black background

nobjects = CC_combo.NumObjects;
indexes = regionprops3(CC_combo,'SubarrayIdx');
indexes = indexes.SubarrayIdx;

toremove = [];
for cc_idx = 1:nobjects %loop through each putative synapse and add a bounding box for segmentation
    if vol(cc_idx) < upperlim
        index_x = indexes{cc_idx,2};
        index_y = indexes{cc_idx,1};
        index_z = indexes{cc_idx,3};

        startx = index_x(1) - bufferx;
        endx = index_x(end) + bufferx;

        starty = index_y(1) - buffery;
        endy = index_y(end) + buffery;

        startz = index_z(1) - bufferz;
        endz = index_z(end) + bufferz;

        %if any of the indexes are negative
        if startx < 1
            startx = 1;
        end
        if starty < 1
            starty = 1;
        end
        if startz < 1
            startz = 1;
        end
        if endx > imsize(2)
            endx = imsize(2);
        end
        if endy > imsize(1)
            endy = imsize(1);
        end
        if endz > imsize(3)
            endz = imsize(3);
        end

        I_s(starty:endy,startx:endx,startz:endz)=1;
    else
        disp(['volume over size limit'])
        toremove=[toremove; cc_idx];
    end
end
end

