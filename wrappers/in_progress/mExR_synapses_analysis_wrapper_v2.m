%%  WRAPPER for analyzing multi-ExR synapses data %% 

% Accompanies the manuscript by Kang*, Schroeder* et al., 2023
% Last modified by Margaret Schroeder on 12/14/2022

% Image processing up to this point:
% - background subtraction in Fiji and
% - registration using modified ExSeq pipeline
% - manual ROI identification (for some sections) and running the Fiji
% script to automatically crop these

% NOTE: It's best to run this wrapper script section by section and be
% mindful of the outputs

%% Set parameters

clear all

%the directory containing all of your image files, named according to the
%following convention: [FIELD OF VIEW NAME]_[ROUND]_[CHANNEL].tif, as
%output from the ExSeqProcessing registration pipeline
params.parentfolder =  'E:/Margaret/mExR/2023.03_synapses/cropped_rois/'; %path to manually-identified ROIs
%params.threshfolder = 'E:/Margaret/mExR/2022.08_synapses/cropped_z/'; %path to whole fields of view, for determining intensity thresholds

% %names of fields of view within the parent folder
% fovs = {
%     'ROI1';
%     'ROI2';
%     'ROI3';
%     'ROI4';
%     'ROI5';
%     };

fovs = {
%     'S1ROI1';
%     'S1ROI2';
%     'S1ROI3';
%     'S1ROI4';
    'S2ROI1';
    'S2ROI2';
    'S2ROI3';
    'S2ROI4';
    };
%set parameters
params.xystep = 0.1625/18; %physical pixel size divided by expansoin factor, um/voxel in x and y
params.zstep = 0.25/18; %physical z-step size divided by expansion factor, um/voxel in z
params.do_binarize=1; %set to 1 if images have not been binarized yet
params.channels_rounds = { %strings defining the channels to analyze in channel-round format
    '1-2';'1-3';
    '2-2';'2-3';
    '3-2';'3-3';
    '4-2';'4-3';
    '5-2';'5-3';
    '6-3';
    '7-3';
    '8-2';'8-3';
    '9-2';
    '10-2';'10-3';
    '11-2';'11-3';
    '12-2'
    }; %FOR BATCH 2
% params.channels_rounds = {
%     '1-1';'1-2';'1-3';
%     '2-1';'2-2';'2-3';
%     '3-1';'3-2';'3-3';
%     '4-2';'4-3';
%     '5-2';'5-3';
%     '7-2';
%     '8-2';'8-3';
%     '11-3'
%     }; %FOR BATCH 1
params.thresh_method = 'zscore';%method for intensity thresholding; either percentile, z-score based, or absolute
params.thresh_multiplier = 3; %at least this many standard deviations above the mean
params.thresh_val = 98; %if thresh_method is percentile, the percentile at which to set the intensity threshold
params.do_med_filt = 1; %1 to do median filter to eliminate high spatial frequency noise
params.filt_size = [3 3 1]; %filter size [x y z] for 3D median filter
n_synapses = [20 20 20 20 20 20 20 20]; %number of synapses that were manually identified for each field of view
params.doreg = 0; %whether or not to perform an additional registration to round 1 based on the shared reference channel
params.pixel_size=82; %pixel size for using the nanocolumn shift method to perform this additional registration (see the associated nanocolumn code)
params.rmax = 40; %rmax for using for using the nanocolumn shift method to perform this additional registration (see the associated nanocolumn code)
params.distance = [0,500]; %minimum and maximum shift radius (in nm), based on registration error, for using the nanocolumn shift method to perform this additional registration (see the associated nanocolumn code)
params.doplot=0; %whether or not to plot output images
params.minsize = 50; %minimum size (in voxels) that a puncta should be to be counted

%% Run the analysis - this extracts multiplexed synaptic properties

for fidx = 1:length(fovs)
    imdatas(fidx).filename = fovs{fidx};
    imdatas(fidx).imdata = synapse_properties_multiplexed_nopw(params,fovs{fidx},n_synapses(fidx));
end
%% Save down the data
%change the filename/date here as appropriate
save([params.parentfolder 'S1-2_multiplex_imdatas_20230609.mat'],'imdatas')

%Your next step : run mExR_synapse_data_processing.m to convert these .mat
%files to dataframe format