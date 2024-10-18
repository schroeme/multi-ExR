%%  WRAPPER for analyzing new PV/Thy1-GFP (NOT multiiplexed) ExR dataset from Jul / August 2024

% Last modified by Margaret Schroeder on 7/17/25

% Image processing up to this point:
% - background subtraction in Fiji and

% NOTE: It's best to run this wrapper script section by section and be
% mindful of the outputs

%% Set up and parameter setting - Cortex

clear all

%the directory containing all of your image files, named according to the
%following convention: [FIELD OF VIEW NAME]_[ROUND]_[CHANNEL].tif, as
%output from the ExSeqProcessing registration pipeline
params.parentfolder = 'A:/Margaret/mExR/2023.07_PV-GFP_ExR/preprocessed/';

%set parameters
params.xystep = 0.1625/18; %physical pixel size divided by expansoin factor, um/voxel in x and y
params.zstep = 0.25/18; %physical z-step size divided by expansion factor, um/voxel in z
params.dofilt = 1; %1 to do median filter to eliminate high spatial frequency noise
params.filt_size = [5 5 3]; %filter size [x y z] for 3D median filter
params.filt = 'med'; %method for filtering image
params.vol_converter = ((2/3)*params.xystep + (1/3)*params.zstep)^3; %weighted average to convert voxels to um^3
params.nchannels=3;
params.gfilt=1;
params.sigma=3;


%% Run quantification

%change the parameters for this analysis
params.doplot=0;
params.sizefilt=1; %in voxels
params.savemasks=1;
params.outfolder = 'A:/Margaret/mExR/2023.07_PV-GFP_ExR/synapse_masks/';
params.lowerlim=50;
params.thresh_method = 'zscore';
params.thresh_multiplier = 3;
%params.filt=0;
params.nbuffer = 20;%number of pixels by which to buffer each bounding box (total)
params.lowerlim_inbox = 10;

data = analyze_synaptic_puncta_fov(params);

save('A:/Margaret/mExR/2023.07_PV-GFP_ExR/analyzed_data_20240719.mat','data')

%% Let's repeat this only on the high-quality SSC data prepared by MZ and YDA

clear all

%the directory containing all of your image files, named according to the
%following convention: [FIELD OF VIEW NAME]_[ROUND]_[CHANNEL].tif, as
%output from the ExSeqProcessing registration pipeline
params.parentfolder = 'A:/Margaret/mExR/2023.07_PV-GFP_ExR/preprocessed_highQ_crop/SSC/';

%set parameters
params.xystep = 0.1625/18; %physical pixel size divided by expansoin factor, um/voxel in x and y
params.zstep = 0.25/18; %physical z-step size divided by expansion factor, um/voxel in z
params.dofilt = 1; %1 to do median filter to eliminate high spatial frequency noise
params.filt_size = [5 5 3]; %filter size [x y z] for 3D median filter
params.filt = 'med'; %method for filtering image
params.vol_converter = ((2/3)*params.xystep + (1/3)*params.zstep)^3; %weighted average to convert voxels to um^3
params.nchannels=3;
params.gfilt=1;
params.sigma=3;


%% Run quantification

%change the parameters for this analysis
params.doplot=0;
params.sizefilt=1; %in voxels
params.savemasks=1;
params.savetargetmasks=1;
params.targetoutfolder = 'A:/Margaret/mExR/2023.07_PV-GFP_ExR/target_masks_highQ/';
params.outfolder = 'A:/Margaret/mExR/2023.07_PV-GFP_ExR/synapse_masks_highQ/';
params.gfpfolder = 'A:/Margaret/mExR/2023.07_PV-GFP_ExR/';
params.lowerlim=250;
params.gfp_sigma=10;
params.savegfp = 0;
params.thresh_method = 'zscore';
params.thresh_multiplier = 3;
params.gfp_thresh=0.5; %0.5 standard deviations above mean for SSC, 2 standard deviations above mean for hpc
%params.filt=0;
params.nbuffer = 20;%number of pixels by which to buffer each bounding box (total)
params.lowerlim_inbox = 10;
params.nsyn_estimate=2000; %estimated number of synapses per field of view (make generous)

data = analyze_synaptic_puncta_fov(params);

%save('A:/Margaret/mExR/2023.07_PV-GFP_ExR/analyzed_data_highQ_20240728_HPC.mat','data')

