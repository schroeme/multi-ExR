# multi-EXR

Code base accompanying Kang*, Schroeder* et al., 2023 (preprint here: []).

## Broad overview of multi-ExR analysis pipeline

1. Image pre-processing (Fiji)
2. Image registration (MATLAB)
3. (Optional): Z-stack shortening
4. Registration quality check (MATLAB)
5. (Optional) region of interest (ROI) identification (Fiji)
6. Image / ROI analysis (MATLAB)
7. Post-processing, graphing, statistical testing (Python, Prism, Excel)

## Additional software packages/environments needed
1. [Fiji](https://imagej.net/software/fiji/)
2. [ExSeqProcessing MATLAB pipeline](https://github.com/dgoodwin208/ExSeqProcessing) (put this in the same working directory as this repository)

## Detailed multi-ExR analysis pipeline
1. Image pre-processing (Fiji)
   * Script name(s): `batch_preprocess_mExR.ijm`, `reformat_for_ExSeqPipeline.ijm`
2. Image registration (MATLAB)
   * Script name(s): [ExSeqProcessing MATLAB pipeline](https://github.com/dgoodwin208/ExSeqProcessing) OR [ExM-Tools/mExR](https://github.com/donglaiw/ExM-Toolbox/tree/ck/mExR)
3. (Optional) Z-stack shortening
   * Script name(s): `batch_cropstack.ijm`
4. Registration quality check (MATLAB)
   * Script name(s): mExR_registration_quality_wrapper.m
5. (Optional) region of interest (ROI) identification (Fiji)
   * Fiji: Analyze -> Tools -> ROI Manager
   * Add ROIs manually; make sure the z-position is centered on the desired object
   * After adding all ROIs save them as a .zip file named as "[ROIname]_rois.zip"
   * Script name(s): `crop_ROIs_zrestrict.ijm`
6. Image / ROI analysis (MATLAB)
   * Script name(s): `mExR_5xFAD_dataset2_analysis_wrapper.m`, and `mExR_validation_analysis_wrapper.m`
7. Post-processing, graphing, statistical testing (Python, Prism, Excel)
   * Script name(s): `5xFAD_vs_WT_abeta_lmer.ipynb`, `5xFAD_vs_WT_synapses_lmer.ipynb`.

## File naming convention
Prior to pre-processing, if you wish to use our pre-processing scripts, z-stacks should be contained in a directory with the following organization: [round]/[condition]/ and named "ROI[N] [magnification].tif". Alternatively, you can name the files and organize your folders however you want, and modify the pre-processing .ijm scripts accordingly. After pre-processing (before registration), individual channel z-stacks should be named according to the following convention: [FIELD-OF-VIEW]_round[NNN]_ch[NN].tif. For example, 'WT-Ctx-ROI1_round001_ch01.tif'. After registration, the analysis scripts assume a similar naming convention, though the appendix "_warped.tif" or "_affine.tif" may be added, depending on which registration pipeline you use.

## Notes on the scripts in the folders

Note that not all scripts provided here were used for analysis provided in the manuscript. Scripts in the "in_progress" subfolders are still under development and were not used for manuscript figures. Please check the wrappers to see which scripts were actually called.

## Example dataset and tutorial
[Coming soon].
