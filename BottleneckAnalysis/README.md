# Bottleneck analysis

R code and data used to estimate levels of genetic diversity, recent evidence of a bottleneck (using BOTTLENECK), and calculate M-ratios for <i>B. hypnorum</i> plus five other UK <i>Bombus</i> populations (using microsatellite genotyping data from Crowther _et al._, 2019 and Dreier _et al._, 2014). This analysis requires R, Cervus (http://www.fieldgenetics.com/pages/home.jsp), and BOTTLENECK software (http://www1.montpellier.inra.fr/CBGP/software/Bottleneck/bottleneck.html).

The master script _initial popgen analysis.R_ requires all files from both the **Dreier et al data** and **data out** folders.

To run the script, the local folder setup should be set up to match this repository, i.e.:

    local/BottleneckAnalysis/
    >initial popgen analysis.R
    >Dreier et al data
    >data out

## References

Crowther, L.P., Wright, D.J., Richardson, D.S., Carvell, C. and Bourke, A.F., 2019. Spatial ecology of a range‐expanding bumble bee pollinator. _Ecology and evolution_, **9**(3), pp.986-997. https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.4722

Dreier, S., Redhead, J.W., Warren, I.A., Bourke, A.F., Heard, M.S., Jordan, W.C., Sumner, S., Wang, J. and Carvell, C., 2014. Fine‐scale spatial genetic structure of common and declining bumble bees across an agricultural landscape. _Molecular Ecology_, **23**(14), pp.3384-3395. https://onlinelibrary.wiley.com/doi/full/10.1111/mec.12823
