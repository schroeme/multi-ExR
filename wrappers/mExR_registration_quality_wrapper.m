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
parentfolder = 'E:/Margaret/mExR/2022.10_5xFAD/cropped_z/';

fovs = {
    '5xFAD-Ctx-ROI1';
    '5xFAD-Ctx-ROI2';
    '5xFAD-Ctx-ROI3';
    '5xFAD-Ctx-ROI4';
    '5xFAD-Ctx-ROI5';

    'WT-Ctx-ROI1';
    'WT-Ctx-ROI2';
    'WT-Ctx-ROI3';
    'WT-Ctx-ROI4';
    'WT-Ctx-ROI5';
    };

savefolder = 'E:/Margaret/mExR/2022.10_5xFAD/registered/regQuality/';

params.xystep = 0.1625/18; %physical pixel size divided by expansoin factor, um/voxel in x and y
params.zstep = 0.25/18; %physical z-step size divided by expansion factor, um/voxel in z
params.nchannels=4; %number of channels in each stack
params.parentfolder = parentfolder;
params.error_channel = {'4','4','4','4','4','4','4'}; %channel index, channel on which to calculate error
params.subtract_morph = 0; %subtract the morphology channel? no, because we are using it for reg quality analysis
params.morph_channel = nan; %morphology channel used for registration, for subtracting out
params.rounds = {'1','2','3','4','5','6','7'};%rounds to analyze for registration error

%Params for Dan Goodwin's registration quality evaluation
params.subvol_dim=100;%length of side of sub-volume, in pixels to use in Dan's registration quality measure
params.xrange = 1:2048;%in pixels, extent of volume to analyze
params.yrange = 1:2048;%in pixels, extent of volume to analyze
params.N = 1000; %number of sub-volumes for DG method
params.pct_thresh = 99.5;%percentage of pixel intensity values for thresholding
params.doplot = 0; %plotting flag (1 if output plots desired)
params.nonzero_thresh=.2*2048*2048*60;%number of nonzero pixels required to run registration - any

%% Running section
for fovidx = 1:length(fovs)
    fov = fovs{fovidx};
    disp(fov);
    error_DG{fovidx} = measure_round_alignment_mExR(params,fov);
end

%% Organize data for easy copy/paste
% modify as needed for convenient copy/paste into Excel or Prism for
% plotting
fovid = 10;

cp = [error_DG{1,fovid}{1,2} error_DG{1,fovid}{1,3} error_DG{1,fovid}{1,4}];