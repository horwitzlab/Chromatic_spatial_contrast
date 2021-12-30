## Thid folder contains 4 major scripts:
- `PopAnalyses_cellclassifcationSVD.m`:Script for classifying cells based on their cone weights. Cone weights are obtained by using SVD onto the STA and pulling out the 1st row and column singular vectors.
- `PopAnalysis_WN_interaction_between_subregions_GLMcrossvalidation.m`: Script for computing the spatial integration between RF subregions in response to checkerboard whitenoise stimuli
- `PopAnalysis_comparing_isoresponse_model_fits_crossvalidation.m`: Script for computing the spatial integration between RF subregions in response to isoresponse punctate flashed stimuli
- `PopAnalysis_WN_interaction_within_individual_subregions_GLMcrossvalidation.m`: Script for computing the phosphor signal integration within individual RF subregions
