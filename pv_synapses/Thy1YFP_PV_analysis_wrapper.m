%%  WRAPPER for analyzing Thy1-YFP vs. PV synapses %% 

% Last modified by Margaret Schroeder on 7/4/24

% Image processing up to this point:
% - background subtraction in Fiji and
% - registration using modified ExSeq pipeline
% - manual ROI identification (for some sections) and running the Fiji
% script to automatically crop these

% NOTE: It's best to run this wrapper script section by section and be
% mindful of the outputs

%% Set up and parameter setting

clear all

%the directory containing all of your image files, named according to the
%following convention: [FIELD OF VIEW NAME]_[ROUND]_[CHANNEL].tif, as
%output from the ExSeqProcessing registration pipeline
%params.parentfolder = 'A:/Margaret/mExR/2024.04_Thy1-YFP_PV/cropped_synaptic_rois/PN/';

params.parentfolder = 'A:/Margaret/mExR/2024.04_Thy1-YFP_PV/registered/';

%names of fields of view within the parent folder
ROIs = {
     'Sample1-CA1-fov1';
    'Sample1-CA1-fov2';
    'Sample1-CA1-fov3';
    'Sample1-CA1-fov4';
    'Sample1-CA1-fov5';
    %'Sample1-CA1-fov6'; poor registration
    'Sample1-CA1-fov7';
    'Sample1-CA1-fov8';
    'Sample1-CA1-fov9';
    'Sample1-Ctx-fov1';
    'Sample1-Ctx-fov2';
    'Sample1-Ctx-fov3';
    'Sample1-Ctx-fov4';
    'Sample1-Ctx-fov5';
   % 'Sample1-Ctx-fov6'; poor registration
   % 'Sample2-CA1-fov1'; poor registration
   % 'Sample2-CA1-fov2'; poor registration
   % 'Sample2-CA1-fov4'; poor registration
    'Sample2-CA1-fov5';
    'Sample2-CA1-fov7';
    'Sample2-CA1-fov8';
    'Sample2-CA1-fov9';
    'Sample2-Ctx-fov1';
    'Sample2-Ctx-fov2';
    'Sample2-Ctx-fov3';
    'Sample2-Ctx-fov4';
    'Sample2-Ctx-fov5';
    'Sample2-Ctx-fov6';
    'Sample2-Ctx-fov7';
    'Sample2-Ctx-fov8';
    'Sample2-Ctx-fov9'
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
   '01-2';%PV
   '01-3';%GFP
    '02-2';%GluA4
    '02-3';%NR1
    '03-2';%GluA1
    '03-3';%PSD95
    '04-2';%Stg
    '04-3';%NR2B
    '05-2';%PSD95
    '05-3';%Shank3
    '06-2';%Btbd11
    '06-3';%IRSp53
    };


%% Determine thresholds for quantification of synaptic proteins in ROIs based on whole field of view
params.thresh_method = 'absolute';
params.thresh_multiplier = 3;
threshmat = [];
for fidx = 1:length(ROIs)
    thresholds = determine_intensity_threshold(params,ROIs{fidx},params.syn_channels);
    threshmat = [threshmat; thresholds];
end

%Copied and pasted the thresholds into an Excel file and saved into the
%Dropbox folder.

%% Run quantification of synaptic proteins in synaptic ROIs - PN (pyramidal neuron)

params.parentfolder = 'A:/Margaret/mExR/2024.04_Thy1-YFP_PV/cropped_synaptic_rois/PN/';

%change the parameters for this analysis
params.doplot=0;
params.sizefilt=10;
params.lowerlim=0;
params.thresh_method = 'absolute';

for fidx = 1:length(ROIs)
    data(fidx).roiname = ROIs{fidx};
    params.thresholds=threshmat(fidx,:); %set thresholds for this field of view
    [data(fidx).data] = analyze_synaptic_puncta_cropped(params,ROIs{fidx},params.syn_channels);
end

%% Compile data - synaptic proteins in PNs
%This section puts the data in a convenient format for copy/paste into
%Excel or Prism for graphing. May require hard-coding to your preferences
topaste = [];
for ii = 1:length(ROIs)
    mat = cell2mat(data(ii).data.num_puncta); %change this to pull different variables of interest
    topaste = [topaste; mat];
end


%% Run quantification of synaptic proteins in cropped abeta nanoclusters - PV (pyramidal neuron)

params.parentfolder = 'A:/Margaret/mExR/2024.04_Thy1-YFP_PV/cropped_synaptic_rois/PV/';

%change the parameters for this analysis
params.doplot=0;
params.sizefilt=10;
params.lowerlim=0;
params.thresh_method = 'absolute';
%params.thresh_multiplier = 3;
%params.filt=0;

for fidx = 1:length(ROIs)
    data(fidx).roiname = ROIs{fidx};
    params.thresholds=threshmat(fidx,:); %set thresholds for this field of view
    [data(fidx).data] = analyze_synaptic_puncta_cropped(params,ROIs{fidx},params.syn_channels);
end

%% Compile data - synaptic proteins in pvs
%This section puts the data in a convenient format for copy/paste into
%Excel or Prism for graphing. May require hard-coding to your preferences
topaste = [];
for ii = 1:length(ROIs)
    mat = cell2mat(data(ii).data.AR); %change this to pull different variables of interest
    %mat = data(ii).data.roi_filename;
    topaste = [topaste; mat];
    %topaste = [topaste; mat{1,1}];
end

