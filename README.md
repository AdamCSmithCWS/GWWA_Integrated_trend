# GWWA_Integrated_trend

## Integrated spatial trend model for BBS and GWWA surveys

Integrated, spatially explicit trend model for Golden-winged Warblers. Using data from the BBS and a GWWA-specific survey. The species has declined in some portions of it's range to a point where it is very rarely detected during BBS surveys. As a result, there is effectively no new information coming from the BBS to help understand the species trends.

### To Do

1.  See if there is a way to get the BYM model to converge - currently the data don't support for a spatial and random effect of trend.

    1.  If the BYM model will converge:

        1.  code in a leave future out (LFO) cross-validation comparing the two kinds of spatial models (simple iCAR vs BYM on the trends).

    2.  If the BYM model will not converge:

        1.  Run LFO on spatial vs random-only trend model

2.  Contact Matt to examine the change in habitat using similar methods to the [Betts et al.](https://doi.org/10.1038/s41559-022-01737-8)paper.

3.  Run the full BBS model for the species using 2021 data.
