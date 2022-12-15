function save_chunks(parentfolder,fov,chunkfoldername,CCseg,refch,name)
%Saves segmented synapses as small chunks
mkdir([parentfolder chunkfoldername])

%extract relevant synapse properties
indexes = regionprops3(CCseg,'SubarrayIdx');
indexes = indexes.SubarrayIdx;
centroids = regionprops3(CCseg,'Centroid');
centroids = centroids.Centroid;
    
%load the raw images
fnames = dir([parentfolder '*' fov '*.tif']); %raw images
count=0;
for ch_idx = 1:length(fnames) %loop through all channels except reference
    fname = fnames(ch_idx).name;
    TF1 = contains(fname,refch);
    TF2 = contains(fname,'Seg');
    
    if ~TF1 && ~TF2
        count=count+1;
        I = loadtiff([parentfolder fname]);
        Ibin = loadtiff([parentfolder 'masks/Bin-' fname]);

        for ccidx = 1:CCseg.NumObjects %loop through the nonempty chunks (the ones that met our size criteria)
            rows = indexes{ccidx,1};
            cols = indexes{ccidx,2};
            zs = indexes{ccidx,3};

            chunk = I(rows,cols,zs);
            chunkbin = Ibin(rows,cols,zs);

            grayImage = uint16(chunk);
            binImage = uint8(chunkbin);

            centroid = (centroids(ccidx,:));
            chunktitle = [parentfolder chunkfoldername filesep name '_ch' num2str(count) '_chunk' num2str(ccidx) '_' num2str(round(centroid(1))) '_' num2str(round(centroid(2))) '_' num2str(round(centroid(3))) '.tif'];
            binchunktitle = [parentfolder chunkfoldername filesep 'Bin-' name '_ch' num2str(count) '_chunk' num2str(ccidx) '_' num2str(round(centroid(1))) '_' num2str(round(centroid(2))) '_' num2str(round(centroid(3))) '.tif'];
             
            options.overwrite = true;
            options.message = false;
            saveastiff(grayImage, chunktitle, options);
            saveastiff(binImage, binchunktitle, options);

        end
    end
end
end

