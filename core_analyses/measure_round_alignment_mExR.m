function error_DG = measure_round_alignment_mExR(params,fov)
% Calculate's registration error (in nm, in post-expansion units) of
% registered image volumes using the method for ExSeqProcessing by Dan Goodwin
% from Alon et al., Science 2021

%extract parameters
parentfolder = params.parentfolder;
rounds = params.rounds;
nrounds = length(rounds);

imglist = dir([params.parentfolder '*' fov '*.tif']);
nrounds_measured = length(imglist)/params.nchannels;
if nrounds_measured < nrounds
    disp(['Found fewer rounds than indicated. Proceeding with lower number.'])
    nrounds = nrounds_measured;
    fnames = dir([parentfolder '*' fov '*ch01*.tif']);
    %extract the idxs of rounds that are present
    clear rounds
    for ff = 1:length(fnames)
        name = fnames(ff).name;
        splits = strsplit(name,"_");
        rounds{ff} = splits{2};
    end
end
ROIname = fov;
    
for rridx = 1:nrounds %load all reference channel images for this ROI in all rounds
    error_channel = params.error_channel{rridx};
    fname = dir([parentfolder '*' ROIname '*' rounds{rridx} '_*ch0' error_channel '*.tif']);
    
    if length(fname) > 0
        disp(fname(1).name);
        I = loadtiff([parentfolder fname(1).name]);
        dim = size(I);
        I_all{rridx,1} = I;
        
        if params.subtract_morph
            morph_channel = params.morph_channel{rridx};
            fname = dir([parentfolder '*' ROIname '*' rounds{rridx} '_*ch0' morph_channel '*.tif']);
            I_morph = loadtiff([parentfolder fname(1).name]);
            I_morph = mat2gray(I_morph);
            I_morph = medfilt3(I_morph,[5 5 3]);
            
            thresh = prctile(I_morph(:),99);
            I_morph_bin = imbinarize(I_morph,thresh);
            
            se = strel('disk',50);
            I_morph_closed = imclose(I_morph_bin,se);
            I_morph_bin = double(1 - I_morph_closed);
            
            I = mat2gray(I);
            I = double(I) .* I_morph_bin;
            I_all{rridx,1} = I;
        end
    else
        I_all{rridx,1} = [];
    end
end

rd1=1;

for rd2 = 2:nrounds
    disp(['analyzing alignment between round ' num2str(rd1) ' and round ' num2str(rd2)])
    
    img1 = mat2gray(I_all{1});
    img2 = I_all{rd2};
    if nnz(img2) >= params.nonzero_thresh %if there are enough nonzero pixels
        img2 = mat2gray(img2);
    
        %crop stack to slices with at least one nonzero pixel
        keepers = [];
        for zz = 1:dim(3)
            slice = img2(:,:,zz);
            if nnz(slice) > 0
                keepers = [keepers zz];
            end
        end
        disp([num2str(length(keepers)) ' slices kept']);
        
        if length(keepers) < 10
            error_DG{rd1,rd2} = NaN;
            error_MS(rd1,rd2)=NaN;
        else
            img2 = img2(:,:,keepers);
            img1 = img1(:,:,keepers);
    
            %find the registration error using Dan's method from ExSeq paper
            %script has been adapted by MES
            error_DG{rd1,rd2} = DG_ExSeq_reg_quality_vMS(img1,img2,params);

        end
    else
        error_DG{rd1,rd2} = NaN;
    end
end

end
