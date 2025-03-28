# Code for fitting and analyzing multi-scale racing diffusion models (MS-RDMs)

### Overview
This code uses an adapted version of the Enhanced Models of Choice (EMC2) package. Code is adapted to include the multi-scale dynamic RDMs. To get started:
- `prepData.R` creates dataframes following EMC2-standards based on the online available data from four experiments, found under `datasets` 
- `fitDynamicEAM.R` fits the models
- `modelComparisons.R` does the BPIC-based model comparisons and generates the corresponding tables (Supplementary Tables)
- `makeFigures-X.R` make the figures
- `analyze-X.R` extract parameter values, analyzes the influence of initial Q-values, and analyzes the effect sizes per mechanism
- `parameter_recovery_dct.R` runs the parameter recoveries and generates the corresponding Supplementary Figures
- `dataset_descriptives.R` gets some descriptives for Supplementary Table S1.
- `extra_EMC2_functions` contains the functions that overwrite/augment EMC2 with functionality needed for dynamic models.


### Naming conventions
Compared to the paper, naming conventions here are different:
- Accuracy memory (AM) in the paper is named accuracy history (AH) in the code
- Fluency memory (FM) in the paper is named RS in the code, sometimes abbreviated with "V"
- db (threshold difference) in the paper is named z here (c.f. starting point in a DDM)
- dr (rate difference) in the paper is named v here (c.f. drift rate in DDM)
- b (average threshold) in the paper is named a here (c.f. the threshold in DDM)

Abbrevations in the code roughly combine the EAM parameter that is affected, with the mechanism. I.e.,
- zSM means that SM affects the z parameter (threshold bias)
- vSM means that SM affects the v parameter (rate bias)
- uAH means that AM affects the urgency parameter u
- aAH means that AM affects the threshold parameter a
- vAH means that AM affects the rate quality parameter
- bV means that FM affects overall threshold b (apologies for the naming inconsistencies!)

Trends/cosine models are interchangeably called 'DCT' models. They're NULL if no cosines/trends were estimated.
