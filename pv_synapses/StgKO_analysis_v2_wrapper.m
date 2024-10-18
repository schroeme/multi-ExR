%%  WRAPPER for analyzing synapses in Ctx and TRN after Stg KO %% 

% Last modified by Margaret Schroeder on 9/4/24

% Image processing up to this point:
% - background subtraction in Fiji and
% - registration using modified ExSeq pipeline
% - manual ROI identification (for some sections) and running the Fiji
% script to automatically crop these

% NOTE: It's best to run this wrapper script section by section and be
% mindful of the outputs

% Image processing up to this point:
% - background subtraction in Fiji and

% NOTE: It's best to run this wrapper script section by section and be
% mindful of the outputs

%% Set up and parameter setting - Cortex

clear all

%the directory containing all of your image files, named according to the
%following convention: [FIELD OF VIEW NAME]_[ROUND]_[CHANNEL].tif, as
%output from the ExSeqProcessing registration pipeline

params.parentfolder = 'A:/Margaret/mExR/2023.10_StgKO_Ctx/registered/';

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
params.lowerlim=300; %for Cav2.1 synapse identification
params.upperlim=1e4; %for Lectin exclusion
params.target_sizefilt = 0;
params.target_lowerlim = 200; %for target channel

params.savemasks=1;
params.savetargetmasks=1;
params.synreffolder = 'A:/Margaret/mExR/2023.10_StgKO_Ctx/registered/round1_noSMIGFAP_ch1-3ref/'; %folder containing the Cav2.1 channel WITHOUT GFAP/SMI
params.maskfolder='A:/Margaret/mExR/2023.10_StgKO_Ctx/GFP_masks/'; %GFP+ PN dendrites defined by Menglong (note these are 2D for now)
params.targetoutfolder = 'A:/Margaret/mExR/2023.10_StgKO_Ctx/target_masks_nofilt/';
params.outfolder = 'A:/Margaret/mExR/2023.10_StgKO_Ctx/syn_masks/';
params.target_channels = {'2-2','2-3','3-2','3-3','4-2','4-3'};
params.target_names = {'GluA4','Stg','GluA1','PSD95','GluA3','NR1'};

params.thresh_method = 'zscore';
params.thresh_multiplier = 4;

params.trimz = 1; %1 to trim the stack in z, 0 to not
params.ztrim_range = 10; %number of slices to trim off each side

params.nbuffer = 20;%number of pixels by which to buffer each bounding box (total)
params.lowerlim_inbox = 10;

params.doplot=0;

%% Make the synapse masks

[fnames,nsyn_GFP,nsyn_all,fov_vol_all] = make_synaptic_mask_gfp_manual_multiplexed(params);

%% Run quantification in the target channel using previously generated masks

params.nsyn_estimate=2000; %estimated number of synapses per field of view (make generous)
params.nsyn_all = nsyn_all;
data = analyze_synaptic_puncta_gfp_manual_multitarget(params);

% Save down results
save([params.parentfolder 'analyzed_data_ctx_nofilt_ztrim_20240912.mat'],'data');


%% Set up and parameter setting - TRN

clear all

%the directory containing all of your image files, named according to the
%following convention: [FIELD OF VIEW NAME]_[ROUND]_[CHANNEL].tif, as
%output from the ExSeqProcessing registration pipeline

params.parentfolder = 'A:/Margaret/mExR/2023.10_StgKO_TRN/registered/';

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
params.lowerlim=300; %for Cav2.1 synapse identification
params.upperlim=1e4; %for Lectin exclusion
params.target_sizefilt = 0;
params.target_lowerlim = 200; %for target channel
params.trimz = 1; %1 to trim the stack in z, 0 to not
params.ztrim_range = 10; %number of slices to trim off each side

params.savemasks=1;
params.savetargetmasks=1;
params.synreffolder = 'A:/Margaret/mExR/2023.10_StgKO_TRN/registered/round1_noSMIGFAP_ch1-3ref/'; %folder containing the Cav2.1 channel WITHOUT GFAP/SMI
params.maskfolder='A:/Margaret/mExR/2023.10_StgKO_TRN/GFP_masks/'; %GFP+ PN dendrites defined by Menglong (note these are 2D for now)
params.targetoutfolder = 'A:/Margaret/mExR/2023.10_StgKO_TRN/target_masks_nofilt/';
params.outfolder = 'A:/Margaret/mExR/2023.10_StgKO_TRN/syn_masks/';
params.target_channels = {'1-2','2-2','2-3','3-2','3-3'};
params.target_names = {'GluA4','GluA3','Stg','PSD95','NR1'};

params.thresh_method = 'zscore';
params.thresh_multiplier = 4;


params.nbuffer = 20;%number of pixels by which to buffer each bounding box (total)
params.lowerlim_inbox = 10;

params.doplot=0;

%% Make the synapse masks

[fnames,nsyn_GFP,nsyn_all,fov_vol_all] = make_synaptic_mask_gfp_manual_multiplexed(params);

%% Run quantification in the target channel using previously generated masks

params.nsyn_estimate=2000; %estimated number of synapses per field of view (make generous)
params.nsyn_all = nsyn_all;
data = analyze_synaptic_puncta_gfp_manual_multitarget(params);

% Save down results
save([params.parentfolder 'analyzed_data_trn_nofilt_ztrim_20240913.mat'],'data');