function res = run_nanocolumn_autocorr_refAlign(params,fov)
%update later with documentation

exp_factor = params.exp_factor;   %expansion factor, adjust as needed
parentdir= params.parentdir;
targetdir = params.targetdir;
nChannels = params.nChannels;
rmax = params.rmax;
pixel_size=params.pixel_size;
channels_rounds = params.channels_rounds;

res = [];       %array to keep the results
ch01_syn_files = dir([parentdir fov '*round001*ch01*.tif']);
nsynapses = length(ch01_syn_files); %number of synapses for this ROI

for sss = 1:nsynapses %loop through all synapses
    syn_name_splits = split(ch01_syn_files(sss).name,'_');
    syn_name = [syn_name_splits{5} '_' syn_name_splits{6}];
    syn_files = dir([parentdir fov '*' syn_name]);
    
    for ii = 1:nChannels%loop through all channels
        ch1str=channels_rounds{ii};
        ch1splits = split(ch1str,'-');
        ch1roundstr = ['0' ch1splits{1}];
        ch1str = ch1splits{2};
       
        ch1files = dir([parentdir fov '*round*' ch1roundstr '_ch0' ch1str '*.tif']);
        file1 = fullfile(ch1files(sss).name);

        test1 = loadtiff([parentdir file1]);
        test1 = double(test1);

        %Subdivide pixel size (physical pixel size of camera chip divided by
        %magnification of objective) and step size of stack in order to achieve
        %isometric sub-voxels. 
        im1 = expand(test1,[2,2,3]); %2,2,3 are divisors that give equal pixel values for x,y,z
        pixel = pixel_size/exp_factor; %Pixel size defined here, in nm
    
        xyzshift = [0,0,0]; %no shift because we're doing autocorrelation
        distance = [20,180]; %Min and max distance allowed for xyz shift, in nm
        flag = 0;
        cc = get_corr_3dMatrix_final(im1, im1, pixel, rmax, xyzshift, distance, flag);
        temp = [1, sss, ii, cc];
        xyzshift = cc(1:3);
        if mean(cc(10:15))~=0
            res = [res;temp];
        end
    end
end
save([targetdir 'AutoCorrelation_' fov '.txt'],'res', '-ascii', '-tabs');

%Column 1: ROI number
%Column 2: Synapse number for that ROI
%Column 3: channel 1
%Columns 5-7: xyz shift
%Column 8: Average intensity of channel 1
%Column 9: Average intensity of channel 2
%Column 10-end: Auto-correlation results
end