Propensity Matching and the dm2200 Data
================
Thomas E. Love, Ph.D.
2020-01-19

  - [Setup](#setup)
  - [The `dm2200` data set](#the-dm2200-data-set)
      - [Codebook](#codebook)
      - [Comparing Exposure Groups with
        `tableone`](#comparing-exposure-groups-with-tableone)
  - [Propensity for Exposure](#propensity-for-exposure)
      - [Fitting a Propensity Model](#fitting-a-propensity-model)
          - [Storing the Propensity
            Scores](#storing-the-propensity-scores)
  - [`match1` 1:1 greedy matching without
    replacement](#match1-11-greedy-matching-without-replacement)
      - [Obtaining the Matched Sample](#obtaining-the-matched-sample)
      - [Checking Covariate Balance for our 1:1 Greedy
        Match](#checking-covariate-balance-for-our-11-greedy-match)
          - [Using `bal.tab` to obtain a balance
            table](#using-bal.tab-to-obtain-a-balance-table)
          - [Checking Rubin’s Rules 1 and
            2](#checking-rubins-rules-1-and-2)
          - [Using `bal.plot` from
            `cobalt`](#using-bal.plot-from-cobalt)
          - [Using `love.plot` to look at Standardized
            Differences](#using-love.plot-to-look-at-standardized-differences)
          - [Using `love.plot` to look at Variance
            Ratios](#using-love.plot-to-look-at-variance-ratios)
  - [Outcome Models](#outcome-models)
      - [Unadjusted Models prior to Propensity
        Matching](#unadjusted-models-prior-to-propensity-matching)
          - [Unadjusted Outcome Model for
            `bp_good`](#unadjusted-outcome-model-for-bp_good)
          - [Unadjusted Outcome Model for
            `bmi`](#unadjusted-outcome-model-for-bmi)
      - [Adjusted Outcome Models after
        `match1`](#adjusted-outcome-models-after-match1)
          - [Binary Outcome: `bp_good`](#binary-outcome-bp_good)
          - [Quantitative Outcome: `bmi`](#quantitative-outcome-bmi)
      - [Key References](#key-references)

This document will eventually demonstrate multiple matching strategies
incorporating the propensity score, including the assessment of
covariate balance before and after matching. We focus on binary and
quantitative outcomes in a (simulated) electronic health records data
setting. It uses the `cobalt` package extensively.

  - Greifer, N. (2020). cobalt: Covariate Balance Tables and Plots. R
    package version 4.0.0.

## Setup

``` r
library(skimr); library(tableone)
library(magrittr); library(janitor) 
library(broom); library(survival); library(lme4)
library(cobalt); library(Matching)
library(tidyverse)

theme_set(theme_bw())
```

``` r
dm2200 <- read_csv("data/dm2200.csv") %>% 
    type.convert() %>% # convert characters to factors
    mutate(subject = as.character(subject),
           bp_good = as.numeric(sbp < 140 & dbp < 90))

dm2200
```

    # A tibble: 2,200 x 26
       subject exposure   age race   hisp sex   insur nincome nhsgrad cleve
       <chr>   <fct>    <int> <fct> <int> <fct> <fct>   <dbl>   <int> <int>
     1 S-0001  A           67 Blac~     0 F     Medi~   24700      72     1
     2 S-0002  B           56 White     0 F     Comm~   62300      85     0
     3 S-0003  B           41 Blac~     0 M     Medi~    6500      49     1
     4 S-0004  A           56 Blac~     0 M     Medi~   45400      86     0
     5 S-0005  B           69 Blac~     0 F     Medi~   15400      72     0
     6 S-0006  B           44 Blac~     0 M     Medi~   34100      87     0
     7 S-0007  B           47 Blac~     0 F     Medi~   13900      42     1
     8 S-0008  B           60 Blac~     0 M     Medi~   30800      91     0
     9 S-0009  B           67 Other     1 F     Medi~   47100      90     1
    10 S-0010  B           70 Blac~     0 F     Medi~   20800      69     1
    # ... with 2,190 more rows, and 16 more variables: height_cm <int>,
    #   weight_kg <int>, bmi <dbl>, a1c <dbl>, sbp <int>, dbp <int>, ldl <int>,
    #   visits <int>, tobacco <fct>, statin <int>, ace_arb <int>, betab <int>,
    #   depr_dx <int>, eyeex <int>, pneumo <int>, bp_good <dbl>

# The `dm2200` data set

I’ve simulated data to match real information we’ve collected over the
years at Better Health Partnership on adults who live with diabetes.
These data mirror some of the real data colleted from electronic health
records across the region by Better Health Partnership, but individual
values have been permuted across patients, so the results are not
applicable to any population. The data I simulated from was a subset of
Better Health data that required that the subject fall into exactly one
of the two exposure groups we’ll study, that they live in Cuyahoga
County, prefer English for health-related communications, and have no
missing data on the variables we’ll study.

  - The *exposure* we’ll study is called `exposure` and consists of two
    levels: A and B. I won’t specify the details further on how the
    exposure is determined, except to say that it is uniquely
    determinable for each subject.
  - We’ll study a binary outcome, specifically whether the subject’s
    blood pressure is in control, in the sense that both their systolic
    blood pressure is below 140 mm Hg, *and* their diastolic blood
    pressure is below 90 mm Hg.
  - We’ll also study a continuous outcome, the subject’s body-mass index
    or `bmi`.

## Codebook

*Note*: I used `paste(colnames(dm2200), collapse = " | ")` to help me
make this list.

|   Variable |       Type        | Description                                                              |
| ---------: | :---------------: | ------------------------------------------------------------------------ |
|    subject |     character     | subject identifier (S-0001 to S-2200)                                    |
|   exposure | factor (2 levels) | A or B                                                                   |
|        age |      integer      | age in years                                                             |
|       race | factor (4 levels) | White, Black\_AA, Asian, Other                                           |
|       hisp |        1/0        | 1 = Hispanic or Latinx, 0 = not                                          |
|        sex |        F/M        | F = Female, M = Male                                                     |
|      insur | factor (4 levels) | Insurance: Medicare, Commercial, Medicaid or Uninsured                   |
|    nincome |      integer      | est. Neighborhood Median Income, in $                                    |
|    nhsgrad |      integer      | est. % of adults in Neighborhood who are High School graduates           |
|      cleve |        1/0        | 1 = Cleveland resident, 0 = resident of suburbs                          |
| height\_cm |      integer      | height in cm                                                             |
| weight\_kg |      integer      | weight in kg                                                             |
|        bmi |      numeric      | body mass index (kg/m<sup>2</sup>)                                       |
|        a1c |      numeric      | most recent Hemoglobin A1c (in %)                                        |
|        sbp |      numeric      | most recent systolic blood pressure (in mm Hg)                           |
|        dbp |      numeric      | most recent diastolic blood pressure (in mm Hg)                          |
|   bp\_good |        1/0        | 1 if `sbp` \< 140 and `dbp` \< 90, 0 otherwise                           |
|        ldl |      numeric      | most recent LDL cholesterol (in mg/dl)                                   |
|     visits |      integer      | primary care office visits in past year                                  |
|    tobacco | factor (3 levels) | Tobacco use: Current, Former, Never                                      |
|     statin |        1/0        | 1 if subject had a statin prescription in the past year                  |
|   ace\_arb |        1/0        | 1 if subject had an ACE inhibitor or ARB prescription in the past year   |
|      betab |        1/0        | 1 if subject had a beta-blocker prescription in the past year            |
|   depr\_dx |        1/0        | 1 if the subject has a depression diagnosis                              |
|      eyeex |        1/0        | 1 if the subject has had a retinal eye exam in the past year             |
|     pneumo |        1/0        | 1 if the subject has had a pneumococcal vaccination in the past 10 years |

## Comparing Exposure Groups with `tableone`

``` r
t1 <- CreateTableOne(
    vars = c("age", "race", "hisp", "sex", "insur", 
             "nincome", "nhsgrad", "cleve", "sbp", "dbp",
             "ldl", "visits", "tobacco", "statin", 
             "ace_arb", "betab", "depr_dx", "eyeex", 
             "pneumo", "bmi", "bp_good"), 
    factorVars = c("hisp", "cleve", "statin",
                   "ace_arb", "betab", "depr_dx", 
                   "eyeex", "pneumo", "bp_good"),
    strata = "exposure", 
    data = dm2200)

t1
```

``` 
                     Stratified by exposure
                      A                   B                   p      test
  n                        200                2000                       
  age (mean (SD))        54.90 (8.55)        58.90 (9.83)     <0.001     
  race (%)                                                    <0.001     
     Asian                   2 ( 1.0)           25 ( 1.2)                
     Black_AA              166 (83.0)         1065 (53.2)                
     Other                   1 ( 0.5)           74 ( 3.7)                
     White                  31 (15.5)          836 (41.8)                
  hisp = 1 (%)               4 ( 2.0)          104 ( 5.2)      0.068     
  sex = M (%)              128 (64.0)          842 (42.1)     <0.001     
  insur (%)                                                   <0.001     
     Commercial              7 ( 3.5)          510 (25.5)                
     Medicaid              125 (62.5)          501 (25.1)                
     Medicare               37 (18.5)          913 (45.6)                
     Uninsured              31 (15.5)           76 ( 3.8)                
  nincome (mean (SD)) 26927.50 (11642.18) 37992.55 (20854.20) <0.001     
  nhsgrad (mean (SD))    77.73 (11.50)       82.23 (11.35)    <0.001     
  cleve = 1 (%)            181 (90.5)         1109 (55.5)     <0.001     
  sbp (mean (SD))       135.88 (20.60)      132.44 (16.17)     0.005     
  dbp (mean (SD))        82.47 (9.08)        72.41 (11.94)    <0.001     
  ldl (mean (SD))        94.34 (38.62)       97.22 (36.66)     0.293     
  visits (mean (SD))      4.81 (2.62)         3.69 (1.78)     <0.001     
  tobacco (%)                                                 <0.001     
     Current               108 (54.0)          445 (22.2)                
     Former                 46 (23.0)          775 (38.8)                
     Never                  46 (23.0)          780 (39.0)                
  statin = 1 (%)           154 (77.0)         1603 (80.2)      0.334     
  ace_arb = 1 (%)          160 (80.0)         1531 (76.5)      0.310     
  betab = 1 (%)             53 (26.5)          781 (39.1)      0.001     
  depr_dx = 1 (%)           53 (26.5)          753 (37.6)      0.002     
  eyeex = 1 (%)            104 (52.0)         1225 (61.3)      0.013     
  pneumo = 1 (%)           116 (58.0)         1766 (88.3)     <0.001     
  bmi (mean (SD))        32.83 (8.22)        35.09 (8.17)     <0.001     
  bp_good = 1 (%)          117 (58.5)         1431 (71.5)     <0.001     
```

# Propensity for Exposure

We’ll fit a logistic regression model to predict propensity for exposure
`A` (as compared to `B`), on the basis of these 18 covariates:

  - age, race, hisp, sex, insur, nincome, nhsgrad, cleve,
  - a1c, ldl, visits, tobacco, statin, ace\_arb, betab,
  - depr\_dx, eyeex, pneumo

Practically, we might well fit something more complex than a simple
model with main effects, but that’s what we’ll limit ourselves to in
this setting. Note that we’re not including any direct information on
either of our outcomes, or the elements that go into them. In practical
work, we might fit different propensity scores for each outcome, but
we’re keeping things simple here.

## Fitting a Propensity Model

We’ll use the `f.build` tool from the `cobalt` package here.

``` r
dm2200 <- dm2200 %>%
    mutate(treat = as.logical(exposure == "A"))

covs_1 <- dm2200 %>%
    select(age, race, hisp, sex, insur, nincome,
           nhsgrad, cleve, a1c, ldl, visits, tobacco,
           statin, ace_arb, betab, depr_dx, eyeex, pneumo)

prop_model <- glm(f.build("treat", covs_1), data = dm2200,
                  family = binomial)

tidy(prop_model, conf.int = TRUE) %>%
    select(term, estimate, std.error, conf.low, conf.high, p.value) %>%
    knitr::kable(digits = 3)
```

| term           | estimate | std.error | conf.low | conf.high | p.value |
| :------------- | -------: | --------: | -------: | --------: | ------: |
| (Intercept)    |  \-4.285 |     1.628 |  \-7.645 |   \-1.225 |   0.008 |
| age            |    0.012 |     0.012 |  \-0.011 |     0.037 |   0.305 |
| raceBlack\_AA  |  \-0.048 |     1.011 |  \-1.787 |     2.257 |   0.962 |
| raceOther      |  \-2.176 |     1.513 |  \-5.585 |     0.789 |   0.150 |
| raceWhite      |  \-1.180 |     1.028 |  \-2.963 |     1.147 |   0.251 |
| hisp           |    0.040 |     0.633 |  \-1.347 |     1.183 |   0.949 |
| sexM           |    0.831 |     0.190 |    0.461 |     1.208 |   0.000 |
| insurMedicaid  |    2.422 |     0.418 |    1.669 |     3.330 |   0.000 |
| insurMedicare  |    0.961 |     0.452 |    0.126 |     1.922 |   0.034 |
| insurUninsured |    3.246 |     0.480 |    2.351 |     4.256 |   0.000 |
| nincome        |    0.000 |     0.000 |    0.000 |     0.000 |   0.803 |
| nhsgrad        |  \-0.004 |     0.010 |  \-0.023 |     0.015 |   0.677 |
| cleve          |    1.812 |     0.320 |    1.208 |     2.468 |   0.000 |
| a1c            |  \-0.038 |     0.047 |  \-0.132 |     0.052 |   0.412 |
| ldl            |  \-0.004 |     0.003 |  \-0.009 |     0.001 |   0.104 |
| visits         |    0.268 |     0.041 |    0.186 |     0.349 |   0.000 |
| tobaccoFormer  |  \-1.134 |     0.227 |  \-1.586 |   \-0.695 |   0.000 |
| tobaccoNever   |  \-1.023 |     0.225 |  \-1.471 |   \-0.587 |   0.000 |
| statin         |    0.006 |     0.228 |  \-0.434 |     0.461 |   0.978 |
| ace\_arb       |    0.384 |     0.229 |  \-0.055 |     0.844 |   0.093 |
| betab          |  \-0.260 |     0.204 |  \-0.665 |     0.135 |   0.202 |
| depr\_dx       |  \-0.664 |     0.208 |  \-1.078 |   \-0.263 |   0.001 |
| eyeex          |  \-0.325 |     0.191 |  \-0.699 |     0.050 |   0.089 |
| pneumo         |  \-1.612 |     0.211 |  \-2.027 |   \-1.200 |   0.000 |

``` r
glance(prop_model)
```

    # A tibble: 1 x 7
      null.deviance df.null logLik   AIC   BIC deviance df.residual
              <dbl>   <int>  <dbl> <dbl> <dbl>    <dbl>       <int>
    1         1340.    2199  -416.  880. 1016.     832.        2176

### Storing the Propensity Scores

``` r
dm2200 <- dm2200 %>%
    mutate(ps = prop_model$fitted,
           linps = prop_model$linear.predictors)

ggplot(dm2200, aes(x = exposure, y = linps)) +
    geom_violin() +
    geom_boxplot(width = 0.3)
```

![](matching_with_dm2200_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

# `match1` 1:1 greedy matching without replacement

We’re going to match on the linear propensity score, and define our
`treat` (treatment) as occurring when `exposure` is A.

``` r
match_1 <- Match(Tr = dm2200$treat, X = dm2200$linps, 
                 M = 1, replace = FALSE, ties = FALSE,
                 estimand = "ATT")

summary(match_1)
```

``` 

Estimate...  0 
SE.........  0 
T-stat.....  NaN 
p.val......  NA 

Original number of observations..............  2200 
Original number of treated obs...............  200 
Matched number of observations...............  200 
Matched number of observations  (unweighted).  200 
```

## Obtaining the Matched Sample

Now, we build a new matched sample data frame in order to do some of the
analyses to come. This will contain only the matched subjects.

``` r
match1_matches <- factor(rep(match_1$index.treated, 2))
dm2200_matched1 <- cbind(match1_matches, 
                         dm2200[c(match_1$index.control, 
                                  match_1$index.treated),])
```

Some sanity checks:

``` r
dm2200_matched1 %>% count(exposure)
```

    # A tibble: 2 x 2
      exposure     n
      <fct>    <int>
    1 A          200
    2 B          200

``` r
dm2200_matched1 %>% head()
```

``` 
  match1_matches subject exposure age     race hisp sex     insur nincome
1              1  S-0728        B  64 Black_AA    0   F  Medicaid   12600
2              4  S-1650        B  48 Black_AA    0   M  Medicaid   22100
3             17  S-0834        B  61 Black_AA    0   M  Medicaid   31500
4             26  S-1232        B  69 Black_AA    0   F  Medicare   21700
5             27  S-0678        B  36 Black_AA    0   F  Medicaid   16600
6             37  S-0765        B  39 Black_AA    0   M Uninsured   35000
  nhsgrad cleve height_cm weight_kg  bmi  a1c sbp dbp ldl visits tobacco statin
1      78     1       152        81 35.1  6.9 134  76 140      3  Former      1
2      80     1       172       134 45.3  5.1  98  66 121      3   Never      1
3      80     1       183        94 28.1  6.5 127  79 137      3 Current      1
4      77     1       168        65 23.0  5.9 101  70  76      2 Current      1
5      79     1       160        80 31.3 14.5 142  95 157      5 Current      1
6      82     0       185       125 36.5  7.3 138  83 119      4 Current      1
  ace_arb betab depr_dx eyeex pneumo bp_good treat         ps      linps
1       1     1       1     1      1       1 FALSE 0.03767221 -3.2404325
2       1     0       0     0      0       1 FALSE 0.62855155  0.5260079
3       1     0       1     0      1       1 FALSE 0.33963493 -0.6649215
4       1     0       0     1      1       1 FALSE 0.07464006 -2.5175054
5       0     1       0     1      0       0 FALSE 0.40466703 -0.3860563
6       1     0       1     1      0       1 FALSE 0.41908241 -0.3265413
```

## Checking Covariate Balance for our 1:1 Greedy Match

### Using `bal.tab` to obtain a balance table

``` r
covs_1plus <- dm2200 %>%
    select(age, race, hisp, sex, insur, nincome,
           nhsgrad, cleve, a1c, ldl, visits, tobacco,
           statin, ace_arb, betab, depr_dx, eyeex, pneumo,
           ps, linps)

bal1 <- bal.tab(M = match_1,
                treat = dm2200$exposure,
                covs = covs_1plus, quick = FALSE,
                un = TRUE, disp.v.ratio = TRUE)
bal1
```

    Balance Measures
                        Type Diff.Un V.Ratio.Un Diff.Adj V.Ratio.Adj
    age              Contin.  0.4076     1.3195  -0.0321      1.2044
    race_Asian        Binary  0.0025             -0.0050            
    race_Black_AA     Binary -0.2975              0.0050            
    race_Other        Binary  0.0320              0.0000            
    race_White        Binary  0.2630              0.0000            
    hisp              Binary  0.0320              0.0000            
    sex_M             Binary -0.2190             -0.0400            
    insur_Commercial  Binary  0.2200              0.0100            
    insur_Medicaid    Binary -0.3745              0.0100            
    insur_Medicare    Binary  0.2715              0.0050            
    insur_Uninsured   Binary -0.1170             -0.0250            
    nincome          Contin.  0.5306     3.2086  -0.0046      1.4150
    nhsgrad          Contin.  0.3958     0.9727   0.0167      0.9324
    cleve             Binary -0.3505              0.0100            
    a1c              Contin. -0.0419     0.8190   0.0008      1.0297
    ldl              Contin.  0.0783     0.9012   0.0097      0.9210
    visits           Contin. -0.6304     0.4602  -0.2529      0.8388
    tobacco_Current   Binary -0.3175             -0.0300            
    tobacco_Former    Binary  0.1575              0.0250            
    tobacco_Never     Binary  0.1600              0.0050            
    statin            Binary  0.0315             -0.0600            
    ace_arb           Binary -0.0345             -0.0450            
    betab             Binary  0.1255              0.0550            
    depr_dx           Binary  0.1115             -0.0100            
    eyeex             Binary  0.0925              0.0350            
    pneumo            Binary  0.3030              0.0900            
    ps               Contin. -2.9189     0.1576  -0.7380      0.4789
    linps            Contin. -1.7842     1.3395  -0.2125      0.5376
    
    Sample sizes
                A    B
    All       200 2000
    Matched   200  200
    Unmatched   0 1800

### Checking Rubin’s Rules 1 and 2

We’ll build a little table of the Rubin’s Rules (1 and 2) before and
after our `match_1` is applied.

``` r
covs_2plus <- dm2200 %>%
    select(linps)

rubin_m1 <- bal.tab(M = match_1,
                treat = dm2200$treat,
                covs = covs_2plus, 
                un = TRUE, disp.v.ratio = TRUE)[1]

rubin_report_m1 <- tibble(
    status = c("Rule1", "Rule2"),
    Unmatched = c(rubin_m1$Balance$Diff.Un,
                  rubin_m1$Balance$V.Ratio.Un),
    Matched = c(rubin_m1$Balance$Diff.Adj,
               rubin_m1$Balance$V.Ratio.Adj))

rubin_report_m1 %>% knitr::kable(digits = 2)
```

| status | Unmatched | Matched |
| :----- | --------: | ------: |
| Rule1  |      2.07 |    0.25 |
| Rule2  |      0.75 |    1.86 |

  - The Rule 1 results tell us about the standardized differences
    expressed as proportions, so we’d like to be certain that our
    results are as close to zero as possible, and definitely below 0.5
    in absolute value.
      - Multiply these by 100 to describe them as percentages, adjusting
        the cutoff to below 50 in absolute value.
      - Here, before matching we have a bias of 206.5005338%, and this
        is reduced to 24.5889855% after 1:1 greedy matching.
  - The Rule 2 results tell us about the variance ratio of the linear
    propensity scores. We want this to be within (0.5, 2) and ideally
    within (0.8, 1.25).
      - Here, before matching we have a variance ratio of 74.6560482%,
        and this becomes 185.9963419% after 1:1 greedy matching.

### Using `bal.plot` from `cobalt`

We can look at any particular variable with this approach, for example,
age:

``` r
bal.plot(obj = match_1,
         treat = dm2200$exposure,
         covs = covs_1plus,
         var.name = "age", 
         which = "both",
         sample.names = 
             c("Unmatched Sample", "Matched Sample"))
```

![](matching_with_dm2200_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

We could also look at the propensity scores in each group, perhaps in
mirrored histograms, with …

``` r
bal.plot(obj = match_1,
         treat = dm2200$exposure,
         covs = covs_1plus,
         var.name = "ps", 
         which = "both",
         sample.names = 
             c("Unmatched Sample", "Matched Sample"),
         type = "histogram", mirror = TRUE)
```

![](matching_with_dm2200_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

Can we look at a categorical variable this way?

``` r
bal.plot(obj = match_1,
         treat = dm2200$exposure,
         covs = covs_1plus,
         var.name = "insur", 
         which = "both",
         sample.names = 
             c("Unmatched Sample", "Matched Sample"))
```

![](matching_with_dm2200_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

### Using `love.plot` to look at Standardized Differences

``` r
love.plot(bal1, 
          threshold = .1, size = 3,
          var.order = "unadjusted",
          stats = "mean.diffs",
          stars = "raw",
          sample.names = c("Unmatched", "Matched"),
          title = "Love Plot for our 1:1 Match") +
    labs(caption = "* indicates raw mean differences (for binary variables)")
```

![](matching_with_dm2200_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
love.plot(bal1, 
          threshold = .1, size = 3,
          var.order = "unadjusted",
          stats = "mean.diffs",
          stars = "raw",
          abs = TRUE,
          sample.names = c("Unmatched", "Matched"),
          title = "Absolute Differences for 1:1 Match") +
    labs(caption = "* indicates raw mean differences (for binary variables)")
```

![](matching_with_dm2200_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

### Using `love.plot` to look at Variance Ratios

Note that this will only include the variables (and summaries like `ps`
and `linps`) that describe quantities. Categorical variables are
dropped.

``` r
love.plot(bal1, 
          threshold = .5, size = 3,
          stats = "variance.ratios",
          sample.names = c("Unmatched", "Matched"),
          title = "Variance Ratios for our 1:1 Match") 
```

![](matching_with_dm2200_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

# Outcome Models

We’ll fit two (overly simplistic) outcome models, one for `bp_good` (our
binary outcome) and another for `bmi` (our quantitative outcome.) Later,
we’ll compare the `exposure` effect estimates made here to the estimates
we obtain after propensity matching.

## Unadjusted Models prior to Propensity Matching

### Unadjusted Outcome Model for `bp_good`

``` r
unadj_mod1 <- glm(bp_good == 1 ~ exposure == "A", data = dm2200, 
                  family = binomial())

tidy(unadj_mod1, exponentiate = TRUE, 
     conf.int = TRUE, conf.level = 0.95) %>%
    select(term, estimate, std.error, 
           conf.low, conf.high) %>%
    knitr::kable(digits = 3)
```

| term                | estimate | std.error | conf.low | conf.high |
| :------------------ | -------: | --------: | -------: | --------: |
| (Intercept)         |    2.515 |     0.050 |    2.284 |     2.773 |
| exposure == “A”TRUE |    0.561 |     0.152 |    0.417 |     0.757 |

### Unadjusted Outcome Model for `bmi`

``` r
unadj_mod2 <- lm(bmi ~ exposure == "A", data = dm2200)

tidy(unadj_mod2, conf.int = TRUE, conf.level = 0.95) %>%
    select(term, estimate, std.error, 
           conf.low, conf.high) %>%
    knitr::kable(digits = 3)
```

| term                | estimate | std.error | conf.low | conf.high |
| :------------------ | -------: | --------: | -------: | --------: |
| (Intercept)         |   35.087 |     0.183 |   34.729 |    35.446 |
| exposure == “A”TRUE |  \-2.260 |     0.606 |  \-3.449 |   \-1.071 |

## Adjusted Outcome Models after `match1`

### Binary Outcome: `bp_good`

``` r
result_match1_bp <- clogit(bp_good ~ (exposure == "A") + 
                          strata(match1_matches),
                      data = dm2200_matched1)

tidy(result_match1_bp, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) %>%
    select(term, estimate, std.error, conf.low, conf.high) %>%
    knitr::kable(digits = 3)
```

| term                | estimate | std.error | conf.low | conf.high |
| :------------------ | -------: | --------: | -------: | --------: |
| exposure == “A”TRUE |     0.61 |     0.211 |    0.403 |     0.924 |

### Quantitative Outcome: `bmi`

We’ll use a mixed model to account for our 1:1 matching. The matches
here are treated as a random factor, with the exposure a fixed factor,
in the `lme4` package.

``` r
dm2200_matched1 <- dm2200_matched1 %>% 
    mutate(match1_matches_f = as.factor(match1_matches))

result_match1_bmi <- lmer(bmi ~ (exposure == "A") + 
                              (1 | match1_matches_f), 
                          data = dm2200_matched1)
```

    boundary (singular) fit: see ?isSingular

``` r
tidy(result_match1_bmi, 
     conf.int = TRUE, conf.level = 0.95) %>% 
    filter(group == "fixed") %>%
    select(term, estimate, std.error, 
           conf.low, conf.high) %>%
    knitr::kable(digits = 3)
```

    Warning in bind_rows_(x, .id): binding factor and character vector, coercing
    into character vector

    Warning in bind_rows_(x, .id): binding character and factor vector, coercing
    into character vector

| term                | estimate | std.error | conf.low | conf.high |
| :------------------ | -------: | --------: | -------: | --------: |
| (Intercept)         |   35.016 |     0.597 |   33.846 |    36.186 |
| exposure == “A”TRUE |  \-2.189 |     0.844 |  \-3.844 |   \-0.534 |

## Key References

Matching was performed using the Matching package (Sekhon, 2011), and
covariate balance was assessed using cobalt (Greifer, 2020), both in R
(R Core Team, 2019).

  - Greifer, N. (2020). cobalt: Covariate Balance Tables and Plots. R
    package version 4.0.0.
