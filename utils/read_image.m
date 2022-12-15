function [img,nz,QC] = read_image(path,chind,nchannels,filt)

    info = imfinfo(path);
    try
        C = strsplit(info(1).ImageDescription,char(10));
        ch_str = strsplit(C{3},"=");
        z_str = strsplit(C{4},"=");
        nc = str2double(ch_str{2});
        nz = str2double(z_str{2});
    catch
        nc = nchannels;
        nz = length(info);
    end
    
    QC = 0;
    if nc ~= nchannels
        QC = 1;
        img = [];
        nz = [];
    else
        %nz = numel(info)/nchannels;
        img = zeros(info(1).Height,info(1).Width,nz);
        slices = chind:nchannels:numel(info);
        parfor j = 1:length(slices)
            img(:,:,j) = imread(path, slices(j));
        end
        if filt
            %img = img./(2e16);
            img = mat2gray(img);
            img = medfilt3(img,[5 5 3]);
        else
            img = uint16(img);
        end
    end
    
    %img = img(:,:,5:end-4);
    %nz = size(img,3);