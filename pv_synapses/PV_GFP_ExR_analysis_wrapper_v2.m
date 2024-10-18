%%  WRAPPER for analyzing PV-GFP (NOT multiiplexed) ExR dataset from July 2023 %% 
% Analysis version 2, using same methods as 2024 datasets (unselected vs.
% selected terminology)

% Last modified by Margaret Schroeder on 8/26/24

% Image processing up to this point:
% - background subtraction in Fiji and

% NOTE: It's best to run this wrapper script section by section and be
% mindful of the outputs

%% Set up and parameter setting - Cortex

clear all

%the directory containing all of your image files, named according to the
%following convention: [FIELD OF VIEW NAME]_[ROUND]_[CHANNEL].tif, as
%output from the ExSeqProcessing registration pipeline
params.parentfolder = 'A:/Margaret/mExR/2024-08_SSC & HPC_PV-GFP + 5 Target_50step Crop Synapse ROI/preprocessed/';

%set parameters
params.xystep = 0.1625/18; %physical pixel size divided by expansoin factor, um/voxel in x and y
params.zstep = 0.25/18; %physical z-step size divided by expansion factor, um/voxel in z
params.dofilt = 1; %1 to do median filter to eliminate high spatial frequency noise
params.filt_size = [5 5 3]; %filter size [x y z] for 3D median filter
params.filt = 'med'; %method for filtering image
params.vol_converter = ((2/3)*params.xystep + (1/3)*params.zstep)^3; %weighted average to convert voxels to um^3
params.nchannels=3;
params.gfilt=1; %whether or not to blur the synaptic reference channel for synapse segmentation
params.sigma=3;
params.sizefilt=1; %whether or not to do size filtration
params.lowerlim=450; %for bassoon synapse identification
params.target_sizefilt = 0;
params.target_lowerlim = 200; %for target channel

params.savemasks=1;
params.savetargetmasks=1;
params.maskfolder='A:/Margaret/mExR/2024-08_SSC & HPC_PV-GFP + 5 Target_50step Crop Synapse ROI/GFP_masks/'; %GFP+ PN dendrites defined by Menglong (note these are 2D for now)
params.targetoutfolder = 'A:/Margaret/mExR/2024-08_SSC & HPC_PV-GFP + 5 Target_50step Crop Synapse ROI/target_masks/';
params.outfolder = 'A:/Margaret/mExR/2024-08_SSC & HPC_PV-GFP + 5 Target_50step Crop Synapse ROI/syn_masks/';

params.thresh_method = 'zscore';
params.thresh_multiplier = 4;

params.nbuffer = 20;%number of pixels by which to buffer each bounding box (total)
params.lowerlim_inbox = 10;
params.target_names = {'GluA1','GluA4','Stg','NR1','PSD95'};

params.doplot=0;

%% Make the synapse masks

[fnames,nsyn_GFP,nsyn_all,fov_vol_all] = make_synaptic_mask_gfp_manual(params);

%% Run quantification in the target channel using previously generated masks

params.nsyn_estimate=1000; %estimated number of synapses per field of view (make generous)
data = analyze_synaptic_puncta_gfp_manual(params);

save('A:/Margaret/mExR/2024-08_SSC & HPC_PV-GFP + 5 Target_50step Crop Synapse ROI/analyzed_data_nofilt_20240830.mat','data')
