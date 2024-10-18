%%  WRAPPER for analyzing Thy1-YFP vs. PV synapses %% 

% Last modified by Margaret Schroeder on 6/23/24

% Image processing up to this point:
% - background subtraction in Fiji and
% - registration using modified ExSeq pipeline
% - manual ROI identification (for some sections) and running the Fiji
% script to automatically crop these

% NOTE: It's best to run this wrapper script section by section and be
% mindful of the outputs

%% Set up and parameter setting - Cortex

clear all

%the directory containing all of your image files, named according to the
%following convention: [FIELD OF VIEW NAME]_[ROUND]_[CHANNEL].tif, as
%output from the ExSeqProcessing registration pipeline
params.parentfolder = 'A:/Margaret/mExR/2023.10_StgKO_Ctx/cropped_synaptic_rois/';

%names of fields of view within the parent folder
ROIs = {
    'mouseA-ROI1-2',
    'mouseA-ROI1-3',
    'mouseA-ROI1',
    'mouseA-ROI2',
    'mouseA-ROI3-1',
    'mouseA-ROI3-2',
    'mouseB-ROI1',
    'mouseB-ROI2',
    'mouseB-ROI3',
    'mouseB-ROI4',
    'mouseD-ROI4',
    'mouseD-ROI5',
    'mouseE-ROI1',
    'mouseE-ROI3',
    'mouseE-ROI5',
    'mouseF-ROI2',
    'mouseF-ROI4',
    'mouseG-ROI1',
    'mouseG-ROI2',
    'mouseG-ROI3',
    'mouseG-ROI4',
    'mouseG-ROI5'
    };


%set parameters
params.xystep = 0.1625/18; %physical pixel size divided by expansoin factor, um/voxel in x and y
params.zstep = 0.25/18; %physical z-step size divided by expansion factor, um/voxel in z
% params.thresh_method = 'pct'; %method for intensity thresholding; either percentile, z-score based, or absolute
% params.thresh_pctl = 99; %percentile at which to set intensity threshold
params.dofilt = 1; %1 to do median filter to eliminate high spatial frequency noise
params.filt_size = [5 5 3]; %filter size [x y z] for 3D median filter
params.filt = 'med'; %method for filtering image

% params.doplot = 0; %whether or not to plot raw and thresholded images for inspection
% params.morph_ch = '1'; %reference/morphology channel
% params.morph_round = '1'; %reference/morphology round
% params.subtract_morph=1; %whether or not to subtract out the morphology/reference channel prior to quantifying other signal in the field of view
% params.lowerlim = 100; %minimum object size, in voxels
% params.upperlim = 5000; %maximum object size, in voxels
% 
% params.syn_rad=10; %radius for morphological closing (dilation) of synapses
% params.morph_close_syn = 1; %whether or not to morphologically close synapses
params.vol_converter = ((2/3)*params.xystep + (1/3)*params.zstep)^3; %weighted average to convert voxels to um^3

%define synaptic channels using 'round-channel'
params.syn_channels = {
    '01-1';%Cav2.1
   '01-2';%Myc
   '01-3';%PV
    '02-2';%GluA4
    '02-3';%Stargazin
    '03-2';%GluA1
    '03-3';%PSD95
    '04-2';%GluA3
    '04-3';%NR1
    };


%% Run quantification of synaptic proteins in synaptic ROIs - Ctx 

%change the parameters for this analysis
params.doplot=0;
params.sizefilt=10;
params.lowerlim=0;
params.thresh_method = 'zscore';
params.thresh_multiplier = 3;
%params.filt=0;

for fidx = 1:length(ROIs)
    data(fidx).roiname = ROIs{fidx};
    [data(fidx).data] = analyze_synaptic_puncta_cropped(params,ROIs{fidx},params.syn_channels);
end

%% Compile data - synaptic proteins in Ctx
%This section puts the data in a convenient format for copy/paste into
%Excel or Prism for graphing. May require hard-coding to your preferences
topaste = [];
for ii = 1:length(ROIs)
    %mat = cell2mat(data(ii).data.AR); %change this to pull different variables of interest
    mat = data(ii).data.roi_filename;
    %topaste = [topaste; mat];
    topaste = [topaste; mat{1,1}];
end

%% Set up and parameter setting - TRN

clear all

%the directory containing all of your image files, named according to the
%following convention: [FIELD OF VIEW NAME]_[ROUND]_[CHANNEL].tif, as
%output from the ExSeqProcessing registration pipeline
params.parentfolder = 'A:/Margaret/mExR/2023.10_StgKO_TRN/cropped_synaptic_rois/';

%names of fields of view within the parent folder
ROIs = {
	'mouseA-ROI1';
	'mouseA-ROI1-2';
	'mouseA-ROI2';
	'mouseA-ROI3';
	'mouseA-ROI4';
	'mouseA-ROI5';
	
	'mouseB-ROI1';
	'mouseB-ROI2';
	'mouseB-ROI3';
	'mouseB-ROI4';
	'mouseB-ROI5';
	
	'mouseD-ROI1';
	'mouseD-ROI2';
	'mouseD-ROI4';
	'mouseD-ROI5';
	
	'mouseE-ROI2';
	'mouseE-ROI3';
	'mouseE-ROI5';
	
	'mouseG-ROI1';
	'mouseG-ROI2';
	'mouseG-ROI3';
	'mouseG-ROI4';
	'mouseG-ROI5';
    };


%set parameters
params.xystep = 0.1625/18; %physical pixel size divided by expansoin factor, um/voxel in x and y
params.zstep = 0.25/18; %physical z-step size divided by expansion factor, um/voxel in z
% params.thresh_method = 'pct'; %method for intensity thresholding; either percentile, z-score based, or absolute
% params.thresh_pctl = 99; %percentile at which to set intensity threshold
params.dofilt = 1; %1 to do median filter to eliminate high spatial frequency noise
params.filt_size = [5 5 3]; %filter size [x y z] for 3D median filter
params.filt = 'med'; %method for filtering image

% params.doplot = 0; %whether or not to plot raw and thresholded images for inspection
% params.morph_ch = '1'; %reference/morphology channel
% params.morph_round = '1'; %reference/morphology round
% params.subtract_morph=1; %whether or not to subtract out the morphology/reference channel prior to quantifying other signal in the field of view
% params.lowerlim = 100; %minimum object size, in voxels
% params.upperlim = 5000; %maximum object size, in voxels
% 
% params.syn_rad=10; %radius for morphological closing (dilation) of synapses
% params.morph_close_syn = 1; %whether or not to morphologically close synapses
params.vol_converter = ((2/3)*params.xystep + (1/3)*params.zstep)^3; %weighted average to convert voxels to um^3

%define synaptic channels using 'round-channel'
params.syn_channels = {
    '01-1';%Cav2.1
   '01-2';%GluA4
   '01-3';%PV
    '02-2';%GluA3
    '02-3';%Stargazin
    '03-2';%PSD95
    '03-3';%NR1
    };


%% Run quantification of synaptic proteins in synaptic ROIs - TRN

%change the parameters for this analysis
params.doplot=0;
params.sizefilt=10;
params.lowerlim=0;
params.thresh_method = 'zscore';
params.thresh_multiplier = 3;
%params.filt=0;

for fidx = 1:length(ROIs)
    data(fidx).roiname = ROIs{fidx};
    [data(fidx).data] = analyze_synaptic_puncta_cropped(params,ROIs{fidx},params.syn_channels);
end

%% Compile data - synaptic proteins in TRN
%This section puts the data in a convenient format for copy/paste into
%Excel or Prism for graphing. May require hard-coding to your preferences
topaste = [];
for ii = 1:length(ROIs)
    mat = cell2mat(data(ii).data.AR); %change this to pull different variables of interest
    %mat = data(ii).data.roi_filename;
    topaste = [topaste; mat];
    %topaste = [topaste; mat{1,1}];
end