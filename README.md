# multi-EXR

Code base accompanying Kang*, Schroeder* et al., 2023.

## Broad overview of multi-ExR analysis pipeline
1. Image pre-processing (Fiji)
2. Image registration (MATLAB or Python)
3. (Optional): Z-stack shortening
4. Registration quality check (MATLAB)
5. (Optional) region of interest (ROI) identification (Fiji)
6. Image / ROI analysis (MATLAB)
7. Post-processing, graphing, statistical testing (Python, Prism, Excel)

## Additional software packages/environments needed
* Please see installation instructions on the provider websites.
1. [Fiji](https://imagej.net/software/fiji/)
2. [ExSeqProcessing MATLAB pipeline](https://github.com/dgoodwin208/ExSeqProcessing) (put this in the same working directory as this repository)
3. Python (version 3 or later), [Anaconda](https://docs.anaconda.com/free/anaconda/install/index.html) and the conda environment for linear mixed effects model (lmer) analysis in Python, located in the "environments" folder. See [these instructions](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file) on how to install a conda environment from a yaml file.
4. [MATLAB](https://www.mathworks.com/products/matlab.html), version 2023a or later, with the following tool boxes: Computer Vision Toolbox, Curve Fitting Toolbox, Image Processing Toolbox, Optimization Toolbox, Parallel Computing Toolbox, Signal Processing Toolbox, Statistics and Machine Learning Toolbox.

## Detailed multi-ExR analysis pipeline description and instructions
* Note: expected run-times are estimates based on our experience; precise runtimes will vary based on your data and your system setup.
1. Image pre-processing (Fiji)
   * Script name(s): (1) `batch_preprocess_mExR.ijm`, (2) `reformat_for_ExSeqPipeline.ijm`
   * (1) Implements background subtraction using Fiji’s rolling ball algorithm with a radius of 50 pixels, converts .nd2 files to .tif stacks, and creates a filename in proper format from the file structure in which the data was acquired. Note that depending on how you name the files, parts of the script will have to change.
   * (2) Splits channels and names them correctly for input into the ExSeqProcessing registration pipeline.
   * Expected runtime: ~1 hour or less for 4 rounds of imaging with 4 2048x2048x81 fields of view on a Windows 10 machine with 512G of RAM.
   * Expected output: preprocessed and renamed .tif stacks in the specified output directory.
2. Image registration (MATLAB or Python)
   * Script name(s): [ExSeqProcessing MATLAB pipeline](https://github.com/dgoodwin208/ExSeqProcessing) OR [ExM-Tools/mExR](https://github.com/donglaiw/ExM-Toolbox/tree/ck/mExR)
   * A script for batching registration over multiple fields of view: `mExR_ExSeqReg_wrapper.m`.
   * Starting registration parameters for ExSeqProcessing registration:
      * params.COLOR_CORRECT_CLAMP = [0,0,0];
      * params.DO_DOWNSAMPLE = true;
      * params.DOWNSAMPLE_RATE = 4.;
      * params.SCALE_PYRAMID = [1:9];
      * If a single channel (e.g. ch01) is the reference channel: regparams.REGISTERCHANNELS_SIFT = {'ch01'};
      * Otherwise, the normalized sum of all channels can be used. Run the normalization step and set that channel as the registration channel: regparams.REGISTERCHANNELS_SIFT = {‘summedNorm’};
   * Expected runtime: ~1 hour or less for 4 rounds of imaging with 4 2048x2048x81 fields of view *with GPU acceleration*, several hours to days without GPU acceleration.
   * Expected output: several; primary output is registered .tif stacks in the `4_registration/` directory.
3. (Optional) Z-stack shortening
   * Script name(s): `batch_cropstack.ijm`
   * The amount of the z-stack that is triped is hard coded (we trim the first and last 10 z-slices); adjust this if desired.
   * Expected runtime: ~1 hour or less for 4 rounds of imaging with 4 2048x2048x81 fields of view.
   * Expected output: trimmed z-stacks in the specified output directory.
4. Registration quality check (MATLAB)
   * Script name(s): mExR_registration_quality_wrapper.m
   * Expected runtime: ~3 hours for 4 rounds of imaging with 4 2048x2048x81 fields of view.
   * Expected output: a MATLAB cell array (`error_DG`) containing the estimated registration error, in nm in post-expansion units, for each field of view in a different cell.
5. (Optional) region of interest (ROI) identification (Fiji)
   * Fiji: Analyze -> Tools -> ROI Manager
   * Add ROIs manually; make sure the z-position is centered on the desired object
   * After adding all ROIs save them as a .zip file named as "[ROIname]_rois.zip". ROIs from our data are provided at the Harvard dataverse (see below).
   * Script name(s): `crop_ROIs_zrestrict.ijm`.
   * Expected output: cropped (in x,y, and z) image volumes, separate for each field of view and identified ROI and each channel, in the specified output folder. Note there will likely be thousands of small files.
6. Image / ROI analysis (MATLAB)
   * Script name(s): `mExR_5xFAD_dataset2_analysis_wrapper.m`, and `mExR_validation_analysis_wrapper.m`.
   * Each analysis script is slightly different, but the basic process is: set parameters in wrapper, call scripts and functions that do quantification. The quantification is almost always based on intensity thresholding of some sort, median and/or size filtration, identifying connected components, and pulling out properties of those connected components.
   * Expected runtime variable, usually a few hours.
   * Expected output is variable, but usually a vector/array of numbers.
7. Post-processing, graphing, statistical testing (Python, Prism, Excel)
   * Script name(s): `5xFAD_vs_WT_abeta_lmer.ipynb`, `5xFAD_vs_WT_synapses_lmer.ipynb`.
   * Runtime: minutes.

## File naming convention
Prior to pre-processing, if you wish to use our pre-processing scripts, z-stacks should be contained in a directory with the following organization: [round]/[condition]/ and named "ROI[N] [magnification].tif". Alternatively, you can name the files and organize your folders however you want, and modify the pre-processing .ijm scripts accordingly. After pre-processing (before registration), individual channel z-stacks should be named according to the following convention: [FIELD-OF-VIEW]_round[NNN]_ch[NN].tif. For example, 'WT-Ctx-ROI1_round001_ch01.tif'. After registration, the analysis scripts assume a similar naming convention, though the appendix "_warped.tif" or "_affine.tif" may be added, depending on which registration pipeline you use.

## System requirements
For image registration using the [ExSeqProcessing MATLAB pipeline](https://github.com/dgoodwin208/ExSeqProcessing), a Linux operating system is required. GPU acceleration using the `cuda` flag is recommended. The registration code was run on a machine running Linux Ubuntu 20.04 with 256G of RAM and one NVIDIA GeForce RTX 2080 Ti GPU. For other analyses, MATLAB, Fiji, and Python scripts can be run on Linux, Windows, or Mac machines. We recommend using a machine with at leaset 128G of RAM given the size of image volumes (~600MB per channel per round). The analysis code was developed and tested on a Windows 10 Pro machine with 512G RAM and a NVIDIA GeForce RTX 2080 Ti GPU.

## Example dataset and tutorial
Background-subtracted, registered images are available for download from the Harvard dataverse (link to be made public upon acceptance of publication). An end-to-end tutorial from raw data is forthcoming. In the mean time, analyses used to create the figures can be reproduced from the registered data using the scripts in this repository.

## Notes on the scripts in the folders
Note that not all scripts provided here were used for analysis provided in the manuscript. Scripts in the "in_progress" subfolders are still under development and were not used for manuscript figures. Please check the wrappers to see which scripts were actually called.

## License. 
This code is shared under an [MIT license](https://opensource.org/license/mit/), Copyright 2023.