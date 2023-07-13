%%  WRAPPER for analyzing multi-ExR 5xFAD vs. WT data %% 

% Accompanies the manuscript by Kang*, Schroeder* et al., 2023
% Last modified by Margaret Schroeder on 5/19/2023

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
params.parentfolder = 'E:/Margaret/mExR/2023.05_5xFAD/cropped_z/';

%names of fields of view within the parent folder
ROIs = {
    %5xFAD
    'S1ROI1';
    'S1ROI2';    
    'S1ROI3';
    'S1ROI4';
    'S1ROI5';

    'S2ROI1';
    'S2ROI2';
    'S2ROI3';
    'S2ROI4';    

    %WT
    'S3ROI1';
    'S3ROI2';
    'S3ROI3';
    'S3ROI4'; 

    'S4ROI1';
    'S4ROI2';
    'S4ROI3';
    'S4ROI4'; 
    };


%set parameters
params.xystep = 0.1625/18; %physical pixel size divided by expansoin factor, um/voxel in x and y
params.zstep = 0.25/18; %physical z-step size divided by expansion factor, um/voxel in z
params.thresh_method = 'pct'; %method for intensity thresholding; either percentile, z-score based, or absolute
params.thresh_pctl = 99; %percentile at which to set intensity threshold
params.dofilt = 1; %1 to do median filter to eliminate high spatial frequency noise
params.filt_size = [5 5 3]; %filter size [x y z] for 3D median filter
params.filt = 'med'; %method for filtering image

params.doplot = 0; %whether or not to plot raw and thresholded images for inspection
params.morph_ch = '1'; %reference/morphology channel
params.morph_round = '1'; %reference/morphology round
params.subtract_morph=1; %whether or not to subtract out the morphology/reference channel prior to quantifying other signal in the field of view
params.lowerlim = 100; %minimum object size, in voxels
params.upperlim = 5000; %maximum object size, in voxels
params.morph_rad=50; %radius for morphological closing (dilation) of reference channel
params.syn_rad=10; %radius for morphological closing (dilation) of synapses
params.morph_close_syn = 1; %whether or not to morphologically close synapses
params.vol_converter = ((2/3)*params.xystep + (1/3)*params.zstep)^3; %weighted average to convert voxels to um^3

%define synaptic channels using 'round-channel'
params.syn_channels = {
   '1-2';%RIM 
    '2-2';%GluA3
    '2-3';%NR2B
    '3-2';%RIM-BP
    '3-3';%GluA2
    '4-2';%GluA1
    '4-3';%NR1
    '5-3';%Shank3
    '6-2';%Homer1
    '8-2';%GluA4
    '9-2'; %PSD95
    '9-3';%Bassoon
    '10-2';%synGAP
    '10-3';%IRsp53

    };

%define amyloid-beta channels using 'round-channel'
params.AB_chs = {'1-3';
    '7-2';
    '8-3'};

%% Determine absolute thresholds for abeta channels based on Abeta signal

%change the parameters for this analysis
params.thresh_method = 'zscore';
params.thresh_multiplier = 5;

%run this to determine the absolute threshold
for fidx = 1:length(ROIs)
    thresh(fidx).roiname = ROIs{fidx};
    thresh(fidx).thresh = determine_intensity_threshold_5xFAD(params,ROIs{fidx},params.AB_chs);
end
%% Analyze overlap/correlation between different AB stains 

params.thresh_method = 'absolute';
params.thresholds = [388 %~rounded averages as determined above
526
189];
params.doplot=1;
params.lowerlim = 150; %minimum object size, in voxels

for fidx = 1:length(ROIs)
    data(fidx).roiname = ROIs{fidx};
    [data(fidx).data] = analyze_AB_coloc(params,ROIs{fidx},params.AB_chs);
end

%% Compile data - Abeta vs. WT abeta abundance
%This section puts the data in a convenient format for copy/paste into
%Excel or Prism for graphing. May require hard-coding to your preferences

res_5xFAD=[];
for jj = 1:9
    temp = data(jj).data.vol';%change this to pull different variables
    res_5xFAD = [res_5xFAD;temp];
end

res_WT=[];
for jj = 10:17
    temp = data(jj).data.vol';%change this to pull different variables
    res_WT = [res_WT;temp];
end

%% Quantify differences in CamKIIa signal and volume 

%change the parameters for this analysis
params.lowerlim = 50;
params.subtract_morph = 1;
params.thresh_method = 'pct';
params.thresh_pctl = 99.5;%
params.doplot=1;
camkii_chs = {'6-3'}; %round-channel format

for fidx = 1:length(ROIs)
    data(fidx).roiname = ROIs{fidx};
    [data(fidx).data] = analyze_camkii(params,ROIs{fidx},camkii_chs,params.AB_chs);
end

%% Compile data - CamKIIa in WT vs. 5xFAD
%This section puts the data in a convenient format for copy/paste into
%Excel or Prism for graphing. May require hard-coding to your preferences

res_5xFAD=[];
res_WT=[];
for ii = 10:17
    res_WT = [res_WT; data(ii).data.nobjects]; %change this to pull different variables
end

for jj = 1:9
     res_5xFAD = [res_5xFAD; data(jj).data.nobjects]; %change this to pull different variables
end

%% Determine absolute thresholds for synaptic channels based on wild-type signal

%change the parameters for this analysis
params.thresh_method = 'zscore';
params.thresh_multiplier = 5;

%run this to determine the absolute threshold
for fidx = 1:length(ROIs)
    thresh(fidx).roiname = ROIs{fidx};
    thresh(fidx).thresh = determine_intensity_threshold_5xFAD(params,ROIs{fidx},params.syn_channels);
end

%% Take averages in WT fields of view to determine thresholds

WT_inds = 10:17; %indices of WT fields of view
syn_thresholds = zeros(length(params.syn_channels),1);
for pp = 1:length(params.syn_channels)
    temp = [];
    for ff = 1:length(WT_inds)
        temp = [temp; thresh(WT_inds(ff)).thresh(pp)];
    end
    syn_thresholds(pp,1) = round(mean(temp));
end

%% Run quantification - synaptic channels, whole field of view

%change the parameters for this analysis
params.morph_rad=15;
params.subtract_morph = 1;
params.lowerlim = 100;
params.dofilt = 1;
params.mask_abeta=1;
params.doplot=1;

params.thresh_method = 'absolute';
params.thresholds = syn_thresholds;

for fidx = 1:length(ROIs)
    data(fidx).roiname = ROIs{fidx};
    [data(fidx).data] = analyze_synapses_5xFAD(params,ROIs{fidx},params.syn_channels);
end

%% Synaptic data, whole field of view - compile data by condition
%This section puts the data in a convenient format for copy/paste into
%Excel or Prism for graphing. May require hard-coding to your preferences

%1- RIM 
%2 - GluA3
%3 - NR2B
%4- RIM-BP
%5 - GluA2
%6 - GluA1
%7 - NR1
%8 Shank3
%9 - Homer1
%10 - GluA4
%11 - PSD95
%12- Bassoon
%13 - synGAP
%14- IRsp53

protein = 14; %change this to pull data for different channels
res_WT =[];
res_5xFAD = [];
for ii = 10:17
    res_WT = [res_WT; data(ii).data.vol(protein)]; %change this to pull different variables
end

for jj = 1:9
     res_5xFAD = [res_5xFAD; data(jj).data.vol(protein)];
end

% %% Run quantification of amyloid beta in cropped nanoclusters
% 
% %change the parameters for this analysis
% params.thresh_method = 'zscore';
% params.thresh_multiplier = 4;
% params.doplot=0;
% params.sizefilt=0;
% params.lowerlim=50;
% 
% %path to cropped volumes
% params.parentfolder = 'E:/Margaret/mExR/2022.10_5xFAD/cropped_abeta/';
% 
% %ROIs to quantify (in this case, only 5xFAD)
% ROIs = {
% '5xFAD-ctx-ROI1';
% '5xFAD-ctx-ROI2';
% '5xFAD-ctx-ROI3';
% '5xFAD-ctx-ROI4';
% '5xFAD-ctx-ROI5'
%     };
% 
% for fidx = 1:length(ROIs)
%     data(fidx).roiname = ROIs{fidx};
%     [data(fidx).data] = analyze_mExR_5xFAD_abeta_cropped(params,ROIs{fidx},params.AB_chs);
% end
% 
% %% Compile data - amyloid beta in cropped nanoclusters
% %This section puts the data in a convenient format for copy/paste into
% %Excel or Prism for graphing. May require hard-coding to your preferences
% 
% topaste = [];
% for ii = 1:length(ROIs)
%     mat = cell2mat(data(ii).data.num_puncta); %change this to pull different variables
%     topaste = [topaste; mat];
% end
% 
% %% Run quantification of synaptic proteins in cropped abeta nanoclusters
% 
% %change the parameters for this analysis
% params.doplot=0;
% params.sizefilt=0;
% params.lowerlim=100;
% params.thresh_method = 'absolute';
% params.thresholds = [55
% 80
% 40
% 50];
% 
% %path to cropped volumes
% params.parentfolder = 'E:/Margaret/mExR/2022.10_5xFAD/cropped_abeta/';
% 
% %ROIs to quantify (in this case, only 5xFAD)
% ROIs = {
% '5xFAD-ctx-ROI1';
% '5xFAD-ctx-ROI2';
% '5xFAD-ctx-ROI3';
% '5xFAD-ctx-ROI4';
% '5xFAD-ctx-ROI5'
%     };
% 
% for fidx = 1:length(ROIs)
%     data(fidx).roiname = ROIs{fidx};
%     [data(fidx).data] = analyze_mExR_5xFAD_abeta_cropped(params,ROIs{fidx},params.syn_channels);
% end
% 
% %% Compile data - amyloid beta in cropped nanoclusters
% %This section puts the data in a convenient format for copy/paste into
% %Excel or Prism for graphing. May require hard-coding to your preferences
% topaste = [];
% for ii = 1:length(ROIs)
%     mat = cell2mat(data(ii).data.num_puncta); %change this to pull different variables of interest
%     topaste = [topaste; mat];
% end