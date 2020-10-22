# Diploid male analysis

R code and data used to estimate levels of genetic diversity and diploid male production from males in a UK _B. hypnorum_ population.

This analysis requires the following programs:
* R (https://www.r-project.org/)
* Cervus software (http://www.fieldgenetics.com/pages/home.jsp).

Microsatellite genotyping data for all sampled males is uploaded in both plain text format (for manipulation with R) and in an Excel format. Relevant information for each dataset is provided in Word format (files ending in _README.docx_).

The R markdown script requires the following files from this repository: _Raw_male_typing_matrix.txt_ and _Consensus_male_genotypes.txt_.

These files should then be placed in the same folder as the R markdown script (_diploid_male_analysis.rmd_) before the analysis is run, i.e.:

    local/DiploidMaleAnalysis/
    >diploid_male_analysis.rmd
    >Raw_male_typing_matrix.txt
    >Consensus_male_genotypes.txt
