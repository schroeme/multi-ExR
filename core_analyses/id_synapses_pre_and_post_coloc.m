function Icombo = id_synapses_pre_and_post_coloc(preinds,postinds,Iall,xystep)
%Creates a mask for synapses (combined channels) based on the intersection of pre- and
%post-synaptic markers

    %intersection logic: OR within pre or post, AND between pre and post
    %separate out pre-synaptic channels and combine using OR logic
    for ii = 1:length(preinds)
        if ii == 1
            prechs = Iall{preinds(ii)};
        else
            prechs = max(prechs,Iall{preinds(ii)}); %this matrix will have a one if any presynaptic channel is located there
        end
    end

    for jj = 1:length(postinds)
        if jj == 1
            postchs = Iall{postinds(jj)};
        else
            postchs = max(postchs,Iall{postinds(jj)}); %this matrix will have a one if any presynaptic channel is located there
        end
    end
    
    CC_post = bwconncomp(postchs,26);%find all connected components in the image
    nobjects_post = CC_post.NumObjects;
    Icoloc = zeros(size(Iall{1})); %empty matrix for new mask
        
    for obj = 1:nobjects_post
        clear idx refchunk
        indexes = regionprops3(CC_post,'SubarrayIdx');
        indexes = indexes.SubarrayIdx;

        %add a small (200 nm) cushion to the bounding box to search for
        %opposing pre- and post-synaptic stains
        bbcushion_xy = ceil(.1*(1/xystep));
        bbcushion_z = 3; 
        index_x = indexes{obj,2};
        index_y = indexes{obj,1};
        index_z = indexes{obj,3};

        startx = index_x(1) - bbcushion_xy;
        endx = index_x(end) + bbcushion_xy;

        starty = index_y(1) - bbcushion_xy;
        endy = index_y(end) + bbcushion_xy;

        startz = index_z(1) - bbcushion_z;
        endz = index_z(end) + bbcushion_z;

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
        if endx > size(Iall{1},2)
            endx = size(Iall{1},2);
        end
        if endy > size(Iall{1},1)
            endy = size(Iall{1},1);
        end
        if endz > size(Iall{1},3)
            endz = size(Iall{1},3);
        end
        
        %check for close apposition of pre and post channels
        prechunk = prechs(starty:endy,startx:endx,startz:endz);
        if any(prechunk,'all')
            postchunk = postchs(starty:endy,startx:endx,startz:endz); %get the target chunk
            combochunk = max(prechunk,postchunk); %take the union
            Icoloc(starty:endy,startx:endx,startz:endz) = combochunk; %population the chunk into the segmentation
        end
    end
    Icombo = Icoloc;
end

