# Propensity Matching with the `dm2200` Example 

Last Update: 2023-02-13

The entire example can be downloaded from the front page of the 500-data website.

Here is the [Quarto file](https://github.com/THOMASELOVE/500-data/blob/master/dm2200/dm2200.qmd) including all code.

I've also placed an [HTML version on RPubs](https://rpubs.com/TELOVE/dm2200-500) at https://rpubs.com/TELOVE/dm2200-500, and you can download the Quarto file from there, too (just click on the CODE box on the top right of the screen and select View Source.)

The `dm2200.csv` data file, in [raw, downloadable form](https://raw.githubusercontent.com/THOMASELOVE/500-data/master/dm2200/data/dm2200.csv).

- The Data Set is 100% fictional.

The purpose of this (simulated) example is to demonstrate a range of propensity score matching methods in R. As of the most recent update, this includes:

### Using the `Matching` package

1. 1:1 matching without replacement
2. 1:2 matching without replacement
3. 1:3 matching with replacement
4. 1:1 matching without replacement but with a caliper on the propensity score

### Using the `MatchIt` package

5. 1:1 nearest neighbor matching without replacement
6. 1:1 nearest neighbor caliper matching
7. 1:1 optimal matching without replacement
8. 1:2 optimal matching without replacement

Note that some of the outcome models developed for some of these matches do not converge, and that we're largely ignoring this in the published demonstration.

