function [Iall,Icombo] = load_and_combine_and_mask_multiplexed_channels(parentfolder,fov,fnames,refch,threshpct,medfilt)
%Loads multiple channels of the same image, e.g. from multiplexed data.
%Assumes all are the same width, height, and z-stack length
%fnames = list of filenames to load
%Iref = name of reference channel. If provided, will skip the reference
%and subtract it from the other channels

%load the reference channel for subtraction of soma/other cellular
%structures
reffnames = dir([parentfolder '*' fov '*' refch '*.tif']); %binary from test stain
reffname = reffnames(1).name;

Iref = double(loadtiff([parentfolder reffname]));
thresh = prctile(Iref(:),threshpct);
Iref_bin = imbinarize(Iref,thresh);
Iref_filt = medfilt3(Iref_bin,medfilt);
Iref = rescale(Iref_filt); %convert from 8-bit to double with zeros and ones

for ch_idx = 1:length(fnames) %loop through all images
    fname = fnames(ch_idx).name;
    TF = contains(fname,refch);
    
    if ~TF
        I = double(loadtiff([parentfolder fname]));
        thresh = prctile(I(:),threshpct);
        Imbin = imbinarize(I,thresh);
        Imfilt = medfilt3(Imbin,medfilt);
        
        I = rescale(Imfilt);
        I = I-Iref; %subtract out the reference channel
        I(I<0)=0; %set anything negative to zero
        
        if ch_idx == 1
            Icombo = I;
        else
            Icombo = Icombo + I;
        end
        
        options.overwrite = true;
        imtosave = uint8(255*I);
        saveastiff(imtosave, [parentfolder 'masks/Bin-' fname],options);
        Iall{ch_idx} = I;
        Iall = Iall(~cellfun('isempty',Iall));
    end
end
end

