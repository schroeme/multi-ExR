function Images = load_flatten_filter(parentfolder,folders,normalization,maxes,savename,subtractmean)

for fidx = 1:length(folders)%loop through all folders
    fnames = dir([parentfolder folders{fidx} filesep '*.tif']);
    Iflat = [];
    for imidx = 1:length(fnames) %loop through all files
        fname = fnames(imidx).name;
        clear info
        info=imfinfo([parentfolder folders{fidx} filesep fname]);
        z_length = numel(info);
        %numchunks = z_length/maxes(3);

        for k = 1:z_length
            A = imread(info(1).Filename, k, 'Info', info);
            %[J,noise_out] = wiener2(A,[2 2]);
            if normalization == "global"
                Agray(:,:,k) = double(A)/65535; %normalize globally
            elseif normalization == "local"
                Agray(:,:,k) = mat2gray(A); %normalize locally
            end
        end
        Afilt = medfilt3(Agray);
        
        I = imresize3(Afilt,[maxes(1) maxes(2) maxes(3)]);
%         %chop up the image into chunks in the z direction corresponding to
%         %~34 voxels each
%         for cc = 1:numchunks
%             Im = Agray(:,:,cc:(cc+maxes(3)));
%             %im resize takes rows, columns, slices - y measures rows
%             I = imresize3(Im,[maxes(1) maxes(2) maxes(3)]);%,'method','triangle');
%             Iflat = [Iflat; I(:)'];
%         end
        clear Agray %Im
        Iflat(imidx,:) = I(:);
    end
    Iflatall{fidx}=Iflat;
    clear Iflat
end

Images = Iflatall{1};
for fidx = 2:length(folders) %loop through all folders
    Images = [Images; Iflatall{fidx}];
end

% remove all zero rows
Images = Images(any(Images,2),:);

%subtract off the mean of each pixel
if subtractmean
    Images = Images - mean(Images,1);
end


% Save them
save(savename,'Images','-v7.3')
