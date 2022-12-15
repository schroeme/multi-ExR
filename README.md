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
1. Fiji: https://imagej.net/software/fiji/
2. [ExSeqProcessing MATLAB pipeline](https://github.com/dgoodwin208/ExSeqProcessing) (put this in the same working directory as this repository)
3. Python environments: [link to .env files]

### Detailed pipeline
1. Image pre-processing (Fiji)
* Script name(s): [UPDATE]
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
[UPDATE]

### Example dataset and tutorial
[Coming soon].
