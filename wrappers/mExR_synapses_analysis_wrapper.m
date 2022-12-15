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
addpath('/utils');
addpath('/core_analyses');

%the directory containing all of your image files, named according to the
%following convention: [FIELD OF VIEW NAME]_[ROUND]_[CHANNEL].tif, as
%output from the ExSeqProcessing registration pipeline
params.parentfolder =  'E:/Margaret/mExR/2022.08_synapses/synapses_cropped/'; %path to manually-identified ROIs
params.threshfolder = 'E:/Margaret/mExR/2022.08_synapses/cropped_z/'; %path to whole fields of view, for determining intensity thresholds

%names of fields of view within the parent folder
fovs = {
    'ROI1';
    'ROI2';
    'ROI3';
    'ROI4';
    'ROI5';
    };

%set parameters
params.xystep = 0.1625/18; %physical pixel size divided by expansoin factor, um/voxel in x and y
params.zstep = 0.25/18; %physical z-step size divided by expansion factor, um/voxel in z
params.controlshift = .1; %in um, amount by which to shift the image for "control" analysis in correlations (sort of depricated as this isn't used)
params.do_binarize=1; %set to 1 if images have not been binarized yet
params.channels_rounds = { %strings defining the channels to analyze in channel-round format
    '1-1';'1-2';'1-3';'1-4';
    '2-1';'2-2';'2-3';'2-4';
    '3-1';'3-2';'3-3';'3-4';
    '4-2';
    '4-3';'4-4';
    '5-2';'5-3';'5-4';
    '6-2';'6-3';'6-4';
    '7-3';'7-3';'7-4';
    '8-2';'8-3';'8-4';
    '9-2';'9-3';'9-4';
    '10-2';'10-3';'10-4';
    '11-2';'11-3';'11-4'
    };
params.thresh_method = 'zscore';%method for intensity thresholding; either percentile, z-score based, or absolute
params.thresh_multiplier = 5; %at least this many standard deviations above the mean
params.thresh_val = 98; %if thresh_method is percentile, the percentile at which to set the intensity threshold
params.do_med_filt = 1; %1 to do median filter to eliminate high spatial frequency noise
params.filt_size = [5 5 1]; %filter size [x y z] for 3D median filter
n_synapses = [37 33 36 37 38]; %number of synapses that were manually identified for each field of view
params.doreg = 0; %whether or not to perform an additional registration to round 1 based on the shared reference channel
params.pixel_size=82; %pixel size for using the nanocolumn shift method to perform this additional registration (see the associated nanocolumn code)
params.rmax = 40; %rmax for using for using the nanocolumn shift method to perform this additional registration (see the associated nanocolumn code)
params.distance = [0,500]; %minimum and maximum shift radius (in nm), based on registration error, for using the nanocolumn shift method to perform this additional registration (see the associated nanocolumn code)
params.doplot=0; %whether or not to plot output images

%% Run the analysis - this extracts multiplexed synaptic properties

for fidx = 1:length(fovs)
    syndatas(fidx).fovname = fovs{fidx};
    imdatas(fidx).filename = fovs{fidx};
    pwdatas(fidx).filename = fovs{fidx};
    [imdatas(fidx).imdata,syndatas(fidx).syndata,pwdatas(fidx).pwdata] = synapse_properties_multiplexed(params,fovs{fidx},n_synapses(fidx));
end
%% Save down the data
%change the filename/date here as appropriate
save([params.parentfolder 'multiplex_pwdatas_2022119.mat'],'pwdatas')
save([params.parentfolder 'multiplex_imdatas_20221119.mat'],'imdatas')
save([params.parentfolder 'multiplex_syndatas_20221119.mat'],'syndatas')

%Your next step : run mExR_synapse_data_processing.m to convert these .mat
%files to dataframe format