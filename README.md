# Code for fitting and analyzing multi-scale racing diffusion models (MS-RDMs)

### Overview
This code uses an adapted version of the Enhanced Models of Choice (EMC2) package. Code is adapted to include the multi-scale dynamic RDMs. To get started:
- `prepData.R` creates dataframes following EMC2-standards based on the online available data from four experiments, found under `datasets` 
- `makeSamplers.R` makes the samplers objects, combining data, design and model specification.
- `fitDynamicEAM.R` fits the pre-created sampler object
- `analyzeDynamicEAMs.R` creates some standard figures for diagnostic purposes
- `makeSamplersTrends.R`, `fitDynamicEAMTrends.R`, `analyzeDynamicEAMsTrends.R` do the same but with the cosine trend models included
- `modelComparisons.R` does the BPIC-based model comparisons and generates the corresponding tables (Supplementary Tables S3-S6)
- `makeFiguresSpectra.R` makes Figure 1 (and other spectrum-related figures)
- `makeFiguresSequential.R` makes Figure 3 (and other stimulus-history related figures)
- `makeFiguresPES.R` makes Figure 4 (and other error-related figures)
- `makeFiguresRS.R` makes Figure 5
- `makeTablesParameters.R` extracts trial-by-trial parameter values from the samplers objects, generates Supplementary Table 2, and Supplementary Figures S5-S6)
- `parameter_recovery_without_dct.R` runs the parameter recoveries and generates the corresponding Supplementary Figures
- `utils_plotting.R` contains a lot of plotting utility functions
- `dataset_descriptives.R` gets some descriptives for Supplementary Table S1.


### Naming conventions
Compared to the paper, naming conventions here are different:
- Accuracy memory (AM) in the paper is named accuracy history (AH) in the code
- Fluency memory (FM) in the paper is named RS in the code
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
