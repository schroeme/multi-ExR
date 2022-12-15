%% Wrapper for multi-ExR nanocolumn analysis %%

% Accompanies the manuscript by Kang*, Schroeder* et al., 2023
% Last modified by Margaret Schroeder on 12/14/2022
% Adapted from a script originally written by Tyler Tarr (Sarkar, Kang,
% Wassie et al. Nat Biomedical Engineering 2022)

% Image processing up to this point:
% - background subtraction in Fiji and
% - registration using modified ExSeq pipeline
% - manual ROI identification and running the Fiji script to automatically crop these
% - ROIs for nanocolumn analysis should contain sandwich-like views
% (lateral in the xy-plane) of synapses

% NOTE: It's best to run this wrapper script section by section and be
% mindful of the outputs

% This version of the wrapper script CAN call on a function that uses the shared 
% reference channel across all rounds to calculate x,y,z shift and align synapses
% However, this is not the way we ended up doing it in the manuscript

%% Part 1: Initialize and set parameters
clear all
params.exp_factor = 18;   %expansion factor, adjust as needed
params.parentdir = 'E:/Margaret/mExR/2022.08_synapses/synapses_cropped_nanocolumn/'; %path to cropped ROIs for this analysis
params.targetdir = params.parentdir;
params.channels_rounds = { %channels and rounds to analyze in round-channel string format
    '1-1';'1-2';'1-3';'1-4';
    '2-1';'2-2';'2-3';'2-4';
    '3-1';'3-2';'3-3';'3-4';
    '4-2';'4-3';'4-4';
    '5-2';'5-3';'5-4';
    '6-2';'6-3';'6-4';
    '7-3';'7-3';'7-4';
    '8-2';'8-3';'8-4';
    '9-2';'9-3';'9-4';
    '10-2';'10-3';'10-4';
    '11-2';'11-3';'11-4'
    };
params.prepost = [0 1 1 0 ... %pre/post identity of the synaptic protein. this is important for the xyzshift logic
    1 1 1 0 ...
    0 1 0 0 ...
    1 1 0 ...
    0 0 0 ...
    1 1 0 ...
    0 1 0 ...
    1 1 0 ...
    0 1 0 ...
    0 1 0 ...
    1 1 0
]; %0 for pre and 1 for post
params.nChannels=length(params.channels_rounds); %number of channels
params.rmax = 20;      %Maximum shift distance, in pixels (defined below)

fov_names = {'ROI1';
    'ROI2';
    'ROI3';
    'ROI4';
    'ROI5'};
params.pixel_size=82; %in nm, divide the physical pixel size of the camera chip by 2 and the z-step size by 3
%.1625/2, .250/3 in nm (not exact)

%% Part 2: Read in the saved tiffs and do the cross correlation (not done for manuscript)
% 
% for fidx = 1:length(fov_names)
%     cc_struct(fidx).res_crosscorr = run_nanocolumn_crosscor_refAlign(params,fov_names{fidx});
% end
% %Column 1: ROI number
% %Column 2: Synapse number for that ROI
% %Column 3: channel 1
% %Column 4: channel 2
% %Columns 5-7: xyz shift
% %Column 8: Average intensity of channel 1
% %Column 9: Average intensity of channel 2
% %Column 10-end: Cross-correlation results

%% Part 3: Read in the saved tiffs and do the autocorrelation
% expand the matrix to make voxel cubic
%Autocorrelation is calculated just as the cross-correlation except there
%is no xyzshift

for fidx = 1:length(fov_names)
    ac_struct(fidx).res_autocorr = run_nanocolumn_autocorr_refAlign(params,fov_names{fidx});
end

%Column 1: ROI number
%Column 2: Synapse number for that ROI
%Column 3: channel 1
%Columns 5-7: xyz shift
%Column 8: Average intensity of channel 1
%Column 9: Average intensity of channel 2
%Column 10-end: Auto-correlation results
%% Part 4- Read in the saved tiffs and do protein enrichment analysis

params.rmax=180; %Max distance over which to run the analysis, in nm
params.step=5; %Step size distance for the analysis, in nm. Ideally, this should be in
            %in increments of pixels, which are defined below. For example,
            %if the pixel size below is 5 nm, then the step could be set to
            %5, 10, 15, etc. 
params.flag=0;

params.override_shift = 1;
for fidx = 1:length(fov_names)
    enrich_struct(fidx).res_enrich = run_nanocolumn_enrichment_refAlign(params,fov_names{fidx});
end

%Column 1: ROI id
%Column 2: Synapse number in that ROI
%Column 3: Channel 1
%Column 4: Channel 2
%Column 5-end: Enrichment results
