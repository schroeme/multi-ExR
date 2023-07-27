%% SCRIPT TO EVALUATE QUALITY OF REGISTRATION 

% Uses Dan Goodwin's method from ExSeq paper (Alon et al., Science 2021)

% Accompanies the manuscript by Kang*, Schroeder* et al., 2023
% Last modified by Margaret Schroeder on 12/14/2022

% Image processing up to this point:
% - background subtraction in Fiji and
% - registration using modified ExSeq pipeline

%% Set up and set parameters
clear all

%folder where the field of view image volumes are stored (either full or
%cropped z-stacks to mutually overlapping areas)
parentfolder = 'E:/Margaret/mExR/2023.03_synapses/cropped_z/';

fovs = {
%     'ROI1';
%     'ROI2';
%     'ROI3';
%     'ROI4';
%     'ROI5'
%     'ROI1';
%     'ROI2';
%     'ROI3';
%     'ROI4';
%     'S1ROI1';
    'S1ROI2';    
%     'S1ROI3';
%     'S1ROI4';
%     'S1ROI5';
% 
%     'S2ROI1';
%     'S2ROI2';
%     'S2ROI3';
%     'S2ROI4';    
% 
%     'S3ROI1';
%     'S3ROI2';
%     'S3ROI3';
%     'S3ROI4'; 
% 
%     'S4ROI1';
%     'S4ROI2';
%     'S4ROI3';
%     'S4ROI4'; 
    };

params.xystep = 0.1625/18; %physical pixel size divided by expansoin factor, um/voxel in x and y
params.zstep = 0.25/18; %physical z-step size divided by expansion factor, um/voxel in z
params.nchannels=3; %number of channels in each stack
params.parentfolder = parentfolder;
params.error_channel = {'1','1'};%,'4','4','2','2'}; %channel index, channel on which to calculate error
params.subtract_morph = 0; %subtract the morphology channel? no, because we are using it for reg quality analysis
params.morph_channel = nan; %morphology channel used for registration, for subtracting out
params.rounds = {'01','10'};%,'12','13','14','15','16','17','18','19'};%,'6','7','8','9','10','11','12'};%rounds to analyze for registration error

%Params for Dan Goodwin's registration quality evaluation
params.subvol_dim=100;%length of side of sub-volume, in pixels to use in Dan's registration quality measure
params.xrange = 1:2048;%in pixels, extent of volume to analyze
params.yrange = 1:2048;%in pixels, extent of volume to analyze
params.N = 1000; %number of sub-volumes for DG method
params.pct_thresh = 99;%percentage of pixel intensity values for thresholding
params.doplot = 1; %plotting flag (1 if output plots desired)
params.nonzero_thresh=.2*2048*2048*60;%number of nonzero pixels required to run registration - any
params.doplotrgb=1;

%% Running section

for fovidx = 1:length(fovs)
    fov = fovs{fovidx};
    disp(fov);
    error_DG{fovidx} = measure_round_alignment_mExR(params,fov);
end

%% Organize data for easy copy/paste
% modify as needed for convenient copy/paste into Excel or Prism for plotting
fovid = 2;
cp = [error_DG{1,fovid}{1,2}...
    error_DG{1,fovid}{1,3} error_DG{1,fovid}{1,4}...
    error_DG{1,fovid}{1,5} error_DG{1,fovid}{1,6} error_DG{1,fovid}{1,7}...
    error_DG{1,fovid}{1,8} error_DG{1,fovid}{1,9} error_DG{1,fovid}{1,10}...
    error_DG{1,fovid}{1,11} error_DG{1,fovid}{1,12}
    ];