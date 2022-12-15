%% WRAPPER SCRIPT FOR AUTOMATED SEGMENTATION OF SYNAPSES FROM MULTI-EXR IMAGES %%

% Image processing up to this point was: background subtraction, median
% filtering, and automated thresholding in ImageJ

% Last modified 6/13/22
% Note: this script has not been updated or maintained, but has been
% provided in case someone else would like to use or work on it in the meantime.

%% Initialize and set parameters
params.parentfolder =  'E:/Margaret/mExR/2022.05_synapses/registered/cropped_z/';

fovs = {
    'ROI1';
    'ROI2';
    'ROI3';
    'ROI4'
    };

%set parameters
params.xystep = 0.1625/18; %um/voxel in x and y
params.zstep = 0.40/18; % um/voxel in z
params.method = 'prepost'; %method for determining if a synapse based on colocalization
%either 'prepost' for pre/post opposition, or 'minimum' for minimum number
%of channels, regardless of identity
params.refch = 'ch04'; %name of reference channel, will get ignored
params.postinds = [4; 5; 9; 10; 12];%channel indices of post-synaptic stains (counting rounds sequentially, not counting reference)
params.preinds = [2; 8;11];%channel indices of pre-synaptic stains (counting rounds sequentially, not counting reference)
params.morphclose = 0;
params.minchannels = 3;
params.savechunks = 1;
params.chunkfoldername = 'chunks';
params.matchinginds = [4 14]; %indices of matching proteins across rounds, one row per protein (counting rounds sequentially, not counting reference)
params.controlshift = .05; %in microns
params.closingrad = 0.250; %radius for morphological closing, in microns
params.lowerlim = 0.1; %lower limit for synapse size (1 linear dimension, in microns)
params.upperlim = 0.5; %upper limit for synapse size (1 linear dimension, in microns)
params.buffersize = 0.1; %buffer to add around each synapse to create a bounding box in segmentation
params.dothresh=1; %whether images need to be thresholded
params.thresh_pct = 98; %intensity percentile for thresholding
params.dofilt = 1;
params.medfilt = [5 5 3];

%% Run the segmentation
%function segment_cc_multiplexed(parent,plotting,morphclose,minchannels,savechunks,chunkfoldername)
parfor fidx = 1:length(fovs)
    data(fidx).fov = fovs{fidx};
    [data(fidx).nsynapses,data(fidx).imdata,data(fidx).syndata] = segment_cc_multiplexed(fovs{fidx},params);
end
%%
save([params.parentfolder 'multiplex_synapse_data_20220630.mat'],'data')
