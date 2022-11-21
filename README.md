# GWWA_Integrated_trend

## Integrated spatial trend models for BBS and GWWA surveys

Integrated, spatially explicit trend model for Golden-winged Warblers. Using data from the BBS and a GWWA-specific survey. The species has declined in some portions of it's range to a point where it is very rarely detected during BBS surveys. As a result, there is effectively no new information coming from the BBS to help understand the species trends.

### To Do

1.  Potentially examine the change in habitat using similar methods to the [Betts et al.](https://doi.org/10.1038/s41559-022-01737-8)paper.

## Models

There are two key models being explored here. Both models account for the different observation processes in the BBS and the GWWA surveys.

1.  A spatially explicit, site-level trend model that uses recent (15-20 years) data from the BBS and the GWWA surveys to estimate trends and relative abundances at each site (BBS route, and GWWA quad). This site-level model lends itself very well to understanding and estimating the drivers of site-level trends (e.g., post-hoc analyses of factors that are associated with site-level trends, and site-level covariates of trends incorporated into the main model, like Betts et al.).

2.  A spatially explicit, 1-degree cell strata, first-difference model that uses data from the full time-series of the BBS (1966-2021) and the GWWA surveys (2009-2021) to integrate information on population change, while accounting for site and stratum relative abundance. This model does a good job of generating spatially explicit estimates of trends that also lend themselves to broader-scale summaries (e.g., BCR summaries).
