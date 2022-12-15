# multi-EXR

Code base accompanying Kang*, Schroeder* et al., 2023 (preprint here: []).

## Broad overview of multi-ExR analysis pipeline

1. Image pre-processing (Fiji)
2. Image registration (MATLAB)
3. Registration quality check (MATLAB)
4. (Optional) region of interest (ROI) identification (Fiji)
5. Image / ROI analysis (MATLAB)
6. Post-processing, graphing, statistical testing (Python, Prism, Excel)

### Additional software packages/environments needed
1. [Fiji](https://imagej.net/software/fiji/)
2. [ExSeqProcessing MATLAB pipeline](https://github.com/dgoodwin208/ExSeqProcessing) (put this in the same working directory as this repository)
3. Python environments: [link to .env files]

### Detailed pipeline
1. Image pre-processing (Fiji)
* Script name(s): [UPDATE], split_channels
2. Image registration (MATLAB)
* Script name(s): [ExSeqProcessing MATLAB pipeline](https://github.com/dgoodwin208/ExSeqProcessing) OR [ExM-Tools/mExR](https://github.com/donglaiw/ExM-Toolbox/tree/ck/mExR)
3. Registration quality check (MATLAB)
* Script name(s): mExR_registration_quality_wrapper.m
4. (Optional) region of interest (ROI) identification (Fiji)
* Script name(s): [UPDATE]
5. Image / ROI analysis (MATLAB)
* Script name(s): mExR_5xFAD_analysis_wrapper.m, mExR_nanocolumn_analysis_wrapper.m, mExR_validation_analysis_wrapper.m, mExR_synapses_analysis_wrapper.m, mExR_synapse_data_processing.m
6. Post-processing, graphing, statistical testing (Python, Prism, Excel)
* Script name(s): [UPDATE]

### File naming convention
Prior to pre-processing, if you wish to use our pre-processing scripts, z-stacks should be contained in a directory with the following organization: [round]/[condition]/ and named "ROI[N] [magnification].tif". Alternatively, you can name the files and organize your folders however you want, and modify the pre-processing .ijm scripts accordingly. After pre-processing (before registration), individual channel z-stacks should be named according to the following convention: [FIELD-OF-VIEW]_round[NNN]_ch[NN].tif. For example, 'WT-Ctx-ROI1_round001_ch01.tif'. After registration, the analysis scripts assume a similar naming convention, though the appendix "_warped.tif" or "_affine.tif" may be added, depending on which registration pipeline you use.

### Example dataset and tutorial
[Coming soon].
