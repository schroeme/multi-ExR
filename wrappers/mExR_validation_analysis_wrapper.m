%% WRAPPER CALLING SCRIPTS TO ANALYZE MULTI-EXR VALIDATION DATA

% Accompanies the manuscript by Kang*, Schroeder* et al., 2023
% Last modified by Margaret Schroeder on 12/14/2022

% Image processing up to this point:
% - background subtraction in Fiji and
% - registration using modified ExSeq pipeline
% - manual ROI identification (for some sections) and running the Fiji
% script to automatically crop these

% NOTE: It's best to run this wrapper script section by section and be
% mindful of the outputs

%% Initialize and set parameters

clear all

%folder containing registered image volumes
parentfolder = 'A:/Margaret/mExR/2022.05_validation/cropped_z/';

%fields of view within this folder to analyze. I have commented out some
%fields of view here because stripping data was not obtained for all ROIs,
%and in this instance I wanted to analyze the stripping rounds, as
%specified in params.rounds
fovs = {
    'ROI1';
    'ROI2';
%     'ROI3';
%     'ROI4';
    'ROI5';
%     'ROI6';
%     'ROI8'
    };

params.xystep = 0.1625/18; %physical pixel size divided by expansoin factor, um/voxel in x and y
params.zstep = 0.25/18; %physical z-step size divided by expansion factor, um/voxel in z
params.nchannels=4;%number of channels in each stack
params.parentfolder = parentfolder;
params.normalization = 'minmax'; %method for normalizing images
params.subtract_morph = 1; %subtract the morphology channel? no, because we are using it for reg quality analysis
params.morph_channel = '4'; %morphology channel used for registration, for subtracting out
params.nrounds = 7; %number of rounds for each field of view
params.rounds = {'1';
%     '1-strip';
    '2';
%     '2-strip';
    '3';
%     '3-strip';
    '4';'5';'6';'7'};%if analyzing strip rounds, specify here
params.thresh_method = 'pct'; %method for intensity thresholding; either percentile, z-score based, or absolute
params.doplot=1;
params.thresh_pctl = 99.5; %percentile at which to set intensity threshold
params.filt = 'med'; %method for filtering image. if 'med', is 3D median filter
params.medfilt = [9 9 3]; %size of 3D median filter in [x y z]
params.morph_close = 1; %whether or not to perform morphological closing on objects
params.morph_closingrad = 2; %radius for morphological closing disk for morphology channel, in microns
params.syn_closingrad = 0.25; %radius for morphological closing disk for synapse channel in microns
params.lowerlim = 0.05; %lower limit of size filter for object detction, in microns
params.vol_converter = ((2/3)*params.xystep + (1/3)*params.zstep)^3; %weighted average to convert voxels to um^3
params.dilationrad = 20; %radius (in pixels) for mask dilation

%% Running section for whole field of view analysis
% Will return overall SNR, SNR within/around detected objects, and number
% of objects detected

for fovidx = 1:length(fovs)
    fov = fovs{fovidx};
    disp(fov);
    [data(fovidx).SNR_overall,data(fovidx).SNR_synapses,data(fovidx).nsynapses] = analyze_mExR_validation_strip(fov,params);
end

%% Reformat data for convenient copy-paste into Prism
protein =3;
topaste = [];
for fovidx = 1:length(fovs)
    datacol = data(fovidx).nsynapses(:,protein);
    topaste = [topaste datacol];
end

%% Save down the data

%adjust the filename as needed
save([parentfolder 'validation_SNR_data_20230419.mat'],'data')

%% Running section for manually-cropped synaptic ROIs analysis

%Last run 1/24/24 with new cropped synapses
%update the folder for this analysis
params.parentfolder = 'A:/Margaret/mExR/2022.05_validation/cropped_ROIs/';

% Example of what these would look like if I were analyzing all rounds,
% non-stripping
fovs = {
    'ROI1';
    'ROI2';
    'ROI3';
    'ROI4';
    'ROI5';
    'ROI6';
    'ROI8'
    };
params.rounds = {'1';'2';'3';'4';'5';'6';'7'};%if analyzing non-strip rounds

% fovs = {
% %     'ROI1';
% %     'ROI2';
%     'ROI3';
%     'ROI4';
% %     'ROI5';
%     'ROI6';
%     'ROI8'
%     };
% 
% params.rounds = {'1';'1-strip'; '2';'2-strip';'3';'3-strip';'4';'5';'6';'7'};%if analyzing strip rounds

%params.medfilt = [5 5 1];
params.filt='none'; %don't do median filtering
params.morph_close=0; %don't do morphological closing
params.savechunks=0; %don't save down any objects
params.lowerlim = 0.02; %lower limit of size filter for object detection
params.channels = {'ch02','ch03'};
params.doplot=1;
params.normalization='none';
for fovidx = 1:length(fovs)
    fov = fovs{fovidx};
    disp(fov);
    [datasyn(fovidx).SNR,datasyn(fovidx).num_puncta,datasyn(fovidx).punctavol,datasyn(fovidx).punctaint,datasyn(fovidx).nsynapses] = analyze_mExR_validation_cropped(fov,params);
end

%% Reshape and compile data for copy/paste into Prism

%Format the data for convenient copy/paste into Excel or GraphPad Prism
%protein order = SynGAP, bassoon
%want to get - rows corresponding to rounds, columns corresponding to fovs
protein = 1; %what protein we want to extract data for
res = []; %empty for holding results

for ii = 1:length(fovs)
    for jj = 1:length(params.rounds)
        res(jj,ii) = datasyn(ii).punctaint(jj,protein);
    end
end

%% Save down the data

%change the filename here!
save([parentfolder 'validation_cropped-synapses_data_20240129.mat'],'datasyn')

%% Measure the volume occupied by each channel
% In order to quantify what the feature density should be for adequate
% registration (in our case), or (in your case) to figure out if you have
% sufficient feature density in the registration channel. Higher is usually
% always better

% directory of preprocessed images (unregistered), with channels split
params.parentfolder = 'E:/Margaret/mExR/2022.05_validation/preprocessed_split/';
params.savefolder = 'E:/Margaret/mExR/2022.05_validation/preprocessed_split/summedMasks/';

% fields of view within this folder to analyze (per the naming convention,
% this should be the first string before the underscore)
fovs = {
    'ROI1';
    'ROI2';
    'ROI3';
    'ROI4';
    'ROI5';
    'ROI6';
    'ROI8'
    };

params.morph_close=0;
params.regsummed = 1; %1 if the registration channel was created by summing all channels
params.closingrad = 0.1; %radius for morphological closing disk, in microns
params.writemask = 1;

for fovidx = 1:length(fovs)
    fov = fovs{fovidx};
    disp(fov);
    [data(fovidx).volume,data(fovidx).volume_frac,data(fovidx).min_vol,data(fovidx).max_vol] = analyze_mExR_validation_feature_density(fov,params);
end

%% Reformat feature density measurements to convenient copy/past form

topaste = [];
for fovidx = 1:length(fovs)
    topaste = [topaste;data(fovidx).max_vol];
end