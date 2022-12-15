function [Iall,Icombo] = load_and_combine_multiplexed_channels(parentfolder,fov,fnames,refch)
%Loads multiple channels of the same image, e.g. from multiplexed data.
%Assumes all are the same width, height, and z-stack length
%fnames = list of filenames to load
%Iref = name of reference channel. If provided, will skip the reference
%and subtract it from the other channels

%load the reference channel for subtraction of soma/other cellular
%structures
reffnames = dir([parentfolder '*' fov '*' refch '*.tif']); %binary from test stain
reffname = reffnames(1).name;

Iref = loadtiff([parentfolder reffname]);
Iref = rescale(Iref); %convert from 8-bit to double with zeros and ones

for ch_idx = 1:length(fnames) %loop through all images
    fname = fnames(ch_idx).name;
    TF = contains(fname,refch);
    
    if ~TF
        I = loadtiff([parentfolder fname]);
        I = rescale(I);
        I = I-Iref; %subtract out the reference channel
        I(I<0)=0; %set anything negative to zero
        
        if ch_idx == 1
            Icombo = I;
        else
            Icombo = Icombo + I;
        end
        
        Iall{ch_idx} = I;
        Iall = Iall(~cellfun('isempty',Iall));
    end
end
end

