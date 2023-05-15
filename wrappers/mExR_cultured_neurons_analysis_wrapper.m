%%  WRAPPER for analyzing multi-ExR cultured neurons %% 

% Accompanies the manuscript by Kang*, Schroeder* et al., 2023
% Last modified by Margaret Schroeder on 3/13/23

% Image processing up to this point:
% - background subtraction in Fiji and
% - registration using modified ExSeq pipeline

% NOTE: It's best to run this wrapper script section by section and be
% mindful of the outputs

%% Set parameters - 1 month cultured neuron dataset

clear all

%the directory containing all of your image files, named according to the
%following convention: [FIELD OF VIEW NAME]_[ROUND]_[CHANNEL].tif, as
%output from the ExSeqProcessing registration pipeline
params.parentfolder = 'E:/Margaret/mExR/2023.03_cultured-neurons_1mo/cropped_z/'; %path to whole fields of view, for determining intensity thresholds

%names of fields of view within the parent folder
fovs = {
    'ROI1';
    'ROI2';
    'ROI3';
    'ROI4'
    };

%set parameters
params.xystep = 0.1625/18; %physical pixel size divided by expansoin factor, um/voxel in x and y
params.zstep = 0.25/18; %physical z-step size divided by expansion factor, um/voxel in z
params.controlshift = .1; %in um, amount by which to shift the image for "control" analysis in correlations (sort of depricated as this isn't used)
params.do_binarize=1; %set to 1 if images have not been binarized yet
params.channels_rounds = { %strings defining the channels to analyze in channel-round format
    '1-2';'1-3';'1-4';
    '2-2';'2-3';'2-4';
    '3-2';'3-3';'3-4'
    };
params.subtract_morph=0;
params.dofilt=1;
params.morph_ch = '1-4';
params.thresh_method = 'zscore';%method for intensity thresholding; either percentile, z-score based, or absolute
params.thresh_multiplier = 5; %at least this many standard deviations above the mean
params.thresh_val = 98; %if thresh_method is percentile, the percentile at which to set the intensity threshold
params.do_med_filt = 1; %1 to do median filter to eliminate high spatial frequency noise
params.filt_size = [5 5 1]; %filter size [x y z] for 3D median filter
params.doreg = 0; %whether or not to perform an additional registration to round 1 based on the shared reference channel
params.pixel_size=82; %pixel size for using the nanocolumn shift method to perform this additional registration (see the associated nanocolumn code)
params.rmax = 40; %rmax for using for using the nanocolumn shift method to perform this additional registration (see the associated nanocolumn code)
params.distance = [0,500]; %minimum and maximum shift radius (in nm), based on registration error, for using the nanocolumn shift method to perform this additional registration (see the associated nanocolumn code)
params.doplot=0; %whether or not to plot output images
params.lowerlim = 300; %lower limit on connected component size (voxels), smaller than this will be filtered out of mask
params.upperlim = 1e4; %upper limit on connected component size (voxels), larger than this will be filtered out of mask
params.vol_converter = ((2/3)*params.xystep + (1/3)*params.zstep)^3; %weighted average to convert voxels to um^3

%% One month dataset - volume and number of each synaptic protein
params.doplot=1;

%define synaptic channels using 'round-channel'
params.syn_channels = {
   '1-2'; %Synapsin 1
    '1-3'; %NR1
    '2-2'; %NR2B
    '2-3'; %SynGAP
    '3-2'; %GluA1
    '3-3'; %PSD95
    '4-2'; %Bassoon
    '4-3'; %Gephyrin
    '5-2'; %RIM
    '5-3'; %CaMKIIa
    };

for fidx = 1:length(fovs)
    data(fidx).data = analyze_synapses_cultured(params,fovs{fidx},params.syn_channels);
end


%% One month dataset - volume and number of each synaptic protein, grab data
%This section puts the data in a convenient format for copy/paste into
%Excel or Prism for graphing. May require hard-coding to your preferences

protein = 10; %change this to pull data for different channels
res=[];
for ii = 1:4 %for all fields of view
    res = [res; data(ii).data.AR']; %change this to pull different variables
end

%% 2 wks - volume and number of each synaptic protein
params.doplot=1;
params.parentfolder = 'E:/Margaret/mExR/2023.03_cultured-neurons_2wks/cropped_z/';

%define synaptic channels using 'round-channel'
params.syn_channels = {
   '1-2'; %Synapsin 1
    '1-3'; %NR1
    '2-2'; %NR2B
    '2-3'; %SynGAP
    '3-2'; %GluA1
    '3-3'; %PSD95
    '4-2'; %Bassoon
    '4-3'; %Gephyrin
    '5-2'; %RIM
    '5-3'; %CaMKIIa
    };

for fidx = 1:length(fovs)
    data(fidx).data = analyze_synapses_cultured(params,fovs{fidx},params.syn_channels);
end


%% 2 wks dataset - volume and number of each synaptic protein, grab data
%This section puts the data in a convenient format for copy/paste into
%Excel or Prism for graphing. May require hard-coding to your preferences

protein = 10; %change this to pull data for different channels
res=[];
for ii = 1:4 %for all fields of view
    res = [res; data(ii).data.AR']; %change this to pull different variables
end
