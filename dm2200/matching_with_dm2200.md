Propensity Matching and the dm2200 Data
================
Thomas E. Love, Ph.D.
2020-03-01

  - [Setup](#setup)
  - [The `dm2200` data set](#the-dm2200-data-set)
      - [Codebook](#codebook)
      - [Comparing Exposure Groups with
        `tableone`](#comparing-exposure-groups-with-tableone)
  - [Propensity for Exposure](#propensity-for-exposure)
      - [Fitting a Propensity Model](#fitting-a-propensity-model)
          - [Storing the Propensity
            Scores](#storing-the-propensity-scores)
  - [`match_1` 1:1 greedy matching without replacement with the
    `Matching`
    package](#match_1-11-greedy-matching-without-replacement-with-the-matching-package)
      - [ATT vs. ATE vs. ATC estimates](#att-vs.-ate-vs.-atc-estimates)
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
  - [`match_2` 1:2 greedy matching without replacement with the
    `Matching`
    package](#match_2-12-greedy-matching-without-replacement-with-the-matching-package)
      - [Obtaining the Matched Sample](#obtaining-the-matched-sample-1)
      - [Checking Covariate Balance for our 1:2 Greedy
        Match](#checking-covariate-balance-for-our-12-greedy-match)
          - [Using `bal.tab` to obtain a balance
            table](#using-bal.tab-to-obtain-a-balance-table-1)
          - [Checking Rubin’s Rules 1 and
            2](#checking-rubins-rules-1-and-2-1)
          - [Using `bal.plot` from
            `cobalt`](#using-bal.plot-from-cobalt-1)
          - [Using `love.plot` to look at Standardized
            Differences](#using-love.plot-to-look-at-standardized-differences-1)
          - [Using `love.plot` to look at Variance
            Ratios](#using-love.plot-to-look-at-variance-ratios-1)
  - [`match_3` 1:3 matching, with replacement with the `Matching`
    package](#match_3-13-matching-with-replacement-with-the-matching-package)
      - [Obtaining the Matched Sample](#obtaining-the-matched-sample-2)
      - [Checking Covariate Balance for our 1:3
        Match](#checking-covariate-balance-for-our-13-match)
          - [Using `bal.tab` to obtain a balance
            table](#using-bal.tab-to-obtain-a-balance-table-2)
          - [Checking Rubin’s Rules 1 and
            2](#checking-rubins-rules-1-and-2-2)
          - [Using `bal.plot` from
            `cobalt`](#using-bal.plot-from-cobalt-2)
          - [Using `love.plot` to look at Standardized
            Differences](#using-love.plot-to-look-at-standardized-differences-2)
          - [Using `love.plot` to look at Variance
            Ratios](#using-love.plot-to-look-at-variance-ratios-2)
  - [`match_4` Caliper Matching (1:1 without replacement) with the
    `Matching`
    package](#match_4-caliper-matching-11-without-replacement-with-the-matching-package)
      - [Obtaining the Matched Sample](#obtaining-the-matched-sample-3)
      - [Checking Covariate Balance for our 1:1 Caliper
        Match](#checking-covariate-balance-for-our-11-caliper-match)
          - [Using `bal.tab` to obtain a balance
            table](#using-bal.tab-to-obtain-a-balance-table-3)
          - [Checking Rubin’s Rules 1 and
            2](#checking-rubins-rules-1-and-2-3)
          - [Using `bal.plot` from
            `cobalt`](#using-bal.plot-from-cobalt-3)
          - [Using `love.plot` to look at Standardized
            Differences](#using-love.plot-to-look-at-standardized-differences-3)
          - [Using `love.plot` to look at Variance
            Ratios](#using-love.plot-to-look-at-variance-ratios-3)
  - [`match_5` 1:1 Nearest Neighbor Matching with
    `Matchit`](#match_5-11-nearest-neighbor-matching-with-matchit)
      - [Obtaining the Matched Sample](#obtaining-the-matched-sample-4)
      - [Checking Covariate Balance for our 1:1 Nearest Neighbor
        Match](#checking-covariate-balance-for-our-11-nearest-neighbor-match)
          - [Default Numerical Balance Summary from
            `MatchIt`](#default-numerical-balance-summary-from-matchit)
          - [Using `bal.tab` to obtain a balance
            table](#using-bal.tab-to-obtain-a-balance-table-4)
          - [Checking Rubin’s Rules 1 and
            2](#checking-rubins-rules-1-and-2-4)
          - [Using `bal.plot` from
            `cobalt`](#using-bal.plot-from-cobalt-4)
          - [Using `love.plot` to look at Standardized
            Differences](#using-love.plot-to-look-at-standardized-differences-4)
          - [Using `love.plot` to look at Variance
            Ratios](#using-love.plot-to-look-at-variance-ratios-4)
  - [`match_6` 1:1 Nearest Neighbor Caliper Matching with
    `Matchit`](#match_6-11-nearest-neighbor-caliper-matching-with-matchit)
      - [Obtaining the Matched Sample](#obtaining-the-matched-sample-5)
      - [Checking Covariate Balance for our 1:1 Nearest Neighbor Caliper
        Match](#checking-covariate-balance-for-our-11-nearest-neighbor-caliper-match)
          - [Using `bal.tab` to obtain a balance
            table](#using-bal.tab-to-obtain-a-balance-table-5)
          - [Checking Rubin’s Rules 1 and
            2](#checking-rubins-rules-1-and-2-5)
          - [Using `bal.plot` from
            `cobalt`](#using-bal.plot-from-cobalt-5)
          - [Using `love.plot` to look at Standardized
            Differences](#using-love.plot-to-look-at-standardized-differences-5)
          - [Using `love.plot` to look at Variance
            Ratios](#using-love.plot-to-look-at-variance-ratios-5)
  - [`match_7` 1:1 Optimal Matching with
    `Matchit`](#match_7-11-optimal-matching-with-matchit)
      - [Obtaining the Matched Sample](#obtaining-the-matched-sample-6)
      - [Checking Covariate Balance for our 1:1 Optimal
        Match](#checking-covariate-balance-for-our-11-optimal-match)
          - [Using `bal.tab` to obtain a balance
            table](#using-bal.tab-to-obtain-a-balance-table-6)
          - [Checking Rubin’s Rules 1 and
            2](#checking-rubins-rules-1-and-2-6)
          - [Are the Optimal and Nearest Neighbor Matches the
            Same?](#are-the-optimal-and-nearest-neighbor-matches-the-same)
          - [Using `bal.plot` from
            `cobalt`](#using-bal.plot-from-cobalt-6)
          - [Using `love.plot` to look at Standardized
            Differences](#using-love.plot-to-look-at-standardized-differences-6)
          - [Using `love.plot` to look at Variance
            Ratios](#using-love.plot-to-look-at-variance-ratios-6)
  - [`match_8` 1:2 Optimal Matching with
    `Matchit`](#match_8-12-optimal-matching-with-matchit)
      - [Obtaining the Matched Sample](#obtaining-the-matched-sample-7)
      - [Checking Covariate Balance for our 1:2 Optimal
        Match](#checking-covariate-balance-for-our-12-optimal-match)
          - [Using `bal.tab` to obtain a balance
            table](#using-bal.tab-to-obtain-a-balance-table-7)
          - [Checking Rubin’s Rules 1 and
            2](#checking-rubins-rules-1-and-2-7)
          - [Using `bal.plot` from
            `cobalt`](#using-bal.plot-from-cobalt-7)
          - [Using `love.plot` to look at Standardized
            Differences](#using-love.plot-to-look-at-standardized-differences-7)
          - [Using `love.plot` to look at Variance
            Ratios](#using-love.plot-to-look-at-variance-ratios-7)
  - [Planned matches coming as soon as Dr. Love finishes
    them](#planned-matches-coming-as-soon-as-dr.-love-finishes-them)
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
      - [Adjusted Outcome Models after
        `match2`](#adjusted-outcome-models-after-match2)
          - [Binary Outcome: `bp_good`](#binary-outcome-bp_good-1)
          - [Quantitative Outcome: `bmi`](#quantitative-outcome-bmi-1)
      - [Adjusted Outcome Models after
        `match3`](#adjusted-outcome-models-after-match3)
          - [Binary Outcome: `bp_good`](#binary-outcome-bp_good-2)
          - [Quantitative Outcome: `bmi`](#quantitative-outcome-bmi-2)
      - [Adjusted Outcome Models after
        `match4`](#adjusted-outcome-models-after-match4)
          - [Binary Outcome: `bp_good`](#binary-outcome-bp_good-3)
          - [Quantitative Outcome: `bmi`](#quantitative-outcome-bmi-3)
      - [Adjusted Outcome Models after
        `match5`](#adjusted-outcome-models-after-match5)
          - [Binary Outcome: `bp_good`](#binary-outcome-bp_good-4)
          - [Quantitative Outcome: `bmi`](#quantitative-outcome-bmi-4)
      - [Adjusted Outcome Models after
        `match6`](#adjusted-outcome-models-after-match6)
          - [Binary Outcome: `bp_good`](#binary-outcome-bp_good-5)
          - [Quantitative Outcome: `bmi`](#quantitative-outcome-bmi-5)
      - [Adjusted Outcome Models after
        `match7`](#adjusted-outcome-models-after-match7)
          - [Binary Outcome: `bp_good`](#binary-outcome-bp_good-6)
          - [Quantitative Outcome: `bmi`](#quantitative-outcome-bmi-6)
      - [Adjusted Outcome Models after
        `match8`](#adjusted-outcome-models-after-match8)
          - [Binary Outcome: `bp_good`](#binary-outcome-bp_good-7)
          - [Quantitative Outcome: `bmi`](#quantitative-outcome-bmi-7)
  - [Cleanup](#cleanup)
  - [Key References](#key-references)

This document demonstrates multiple matching strategies incorporating
the propensity score, including the assessment of covariate balance
before and after matching. We focus on binary and quantitative outcomes
in a (simulated) electronic health records data setting. It uses the
`cobalt` package extensively. See the Key References section at the end
of the document.

## Setup

``` r
library(skimr); library(tableone)
library(magrittr); library(janitor) 
library(broom); library(survival); library(lme4)
library(cobalt); library(Matching); library(MatchIt)
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

# `match_1` 1:1 greedy matching without replacement with the `Matching` package

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

## ATT vs. ATE vs. ATC estimates

Note that in each of the matched samples we build, we’ll focus on ATT
estimates (average treated effect on the treated) rather than ATE
estimates. This means that in our matching we’re trying to mirror the
population represented by the “treated” sample we observed.

  - To obtain ATE estimates rather than ATT with the `Match` function
    from the `Matching` package, use `estimand = "ATE"` in the process
    of developing the matched sample.
  - To obtain ATC estimates (average treatment effect on the controls),
    use `estimand = "ATC"`.

I encourage the use of ATT estimates in your projects, where possible. I
suggest also that you define the “treated” group (the one that the
propensity score is estimating) to be the smaller of the two groups you
have, to facilitate this approach. If you estimate ATE or ATC instead of
ATT, of course, you are answering a different question than what ATT
resolves.

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
1              1  S-1640        B  56    White    0   M Uninsured   86600
2              4  S-1650        B  48 Black_AA    0   M  Medicaid   22100
3             17  S-0834        B  61 Black_AA    0   M  Medicaid   31500
4             26  S-1705        B  72 Black_AA    0   F  Medicare   19600
5             27  S-0678        B  36 Black_AA    0   F  Medicaid   16600
6             37  S-0765        B  39 Black_AA    0   M Uninsured   35000
  nhsgrad cleve height_cm weight_kg  bmi  a1c sbp dbp ldl visits tobacco statin
1      92     0       185       136 39.7  6.9 121  62  71      3  Former      1
2      80     1       172       134 45.3  5.1  98  66 121      3   Never      1
3      80     1       183        94 28.1  6.5 127  79 137      3 Current      1
4      87     1       152        74 32.0  6.8 134  69 129      4 Current      1
5      79     1       160        80 31.3 14.5 142  95 157      5 Current      1
6      82     0       185       125 36.5  7.3 138  83 119      4 Current      1
  ace_arb betab depr_dx eyeex pneumo bp_good treat         ps      linps
1       1     1       0     0      1       1 FALSE 0.03776720 -3.2378153
2       1     0       0     0      0       1 FALSE 0.62855155  0.5260079
3       1     0       1     0      1       1 FALSE 0.33963493 -0.6649215
4       1     1       0     1      1       1 FALSE 0.07539326 -2.5066506
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
    age              Contin.  0.4076     1.3195  -0.0626      1.2907
    race_Asian        Binary  0.0025             -0.0050            
    race_Black_AA     Binary -0.2975             -0.0050            
    race_Other        Binary  0.0320             -0.0050            
    race_White        Binary  0.2630              0.0150            
    hisp              Binary  0.0320             -0.0050            
    sex_M             Binary -0.2190             -0.0300            
    insur_Commercial  Binary  0.2200             -0.0100            
    insur_Medicaid    Binary -0.3745             -0.0100            
    insur_Medicare    Binary  0.2715              0.0400            
    insur_Uninsured   Binary -0.1170             -0.0200            
    nincome          Contin.  0.5306     3.2086  -0.0656      0.9947
    nhsgrad          Contin.  0.3958     0.9727   0.0044      0.8888
    cleve             Binary -0.3505              0.0200            
    a1c              Contin. -0.0419     0.8190  -0.0723      1.0127
    ldl              Contin.  0.0783     0.9012  -0.0394      0.9131
    visits           Contin. -0.6304     0.4602  -0.2979      0.7863
    tobacco_Current   Binary -0.3175              0.0050            
    tobacco_Former    Binary  0.1575              0.0300            
    tobacco_Never     Binary  0.1600             -0.0350            
    statin            Binary  0.0315             -0.0600            
    ace_arb           Binary -0.0345             -0.0500            
    betab             Binary  0.1255              0.0400            
    depr_dx           Binary  0.1115              0.0300            
    eyeex             Binary  0.0925             -0.0100            
    pneumo            Binary  0.3030              0.1100            
    ps               Contin. -2.9189     0.1576  -0.7384      0.4790
    linps            Contin. -1.7842     1.3395  -0.2126      0.5376
    
    Sample sizes
                A    B
    All       200 2000
    Matched   200  200
    Unmatched   0 1800

### Checking Rubin’s Rules 1 and 2

We’ll build a little table of the Rubin’s Rules (1 and 2) before and
after our `match_1` is applied.

``` r
covs_for_rubin <- dm2200 %>%
    select(linps)

rubin_m1 <- bal.tab(M = match_1,
                treat = dm2200$treat,
                covs = covs_for_rubin, 
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
        is reduced to 24.6044378% after 1:1 greedy matching.
  - The Rule 2 results tell us about the variance ratio of the linear
    propensity scores. We want this to be within (0.5, 2) and ideally
    within (0.8, 1.25).
      - Here, before matching we have a variance ratio of 74.6560482%,
        and this becomes 186.0251662% after 1:1 greedy matching.

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

# `match_2` 1:2 greedy matching without replacement with the `Matching` package

Again, we’ll match on the linear propensity score, and define our
`treat` (treatment) as occurring when `exposure` is A. The only
difference will be that we’ll allow each subject with exposure A to be
matched to exactly two subjects with exposure B.

``` r
match_2 <- Match(Tr = dm2200$treat, X = dm2200$linps, 
                 M = 2, replace = FALSE, ties = FALSE,
                 estimand = "ATT")

summary(match_2)
```

``` 

Estimate...  0 
SE.........  0 
T-stat.....  NaN 
p.val......  NA 

Original number of observations..............  2200 
Original number of treated obs...............  200 
Matched number of observations...............  200 
Matched number of observations  (unweighted).  400 
```

Note that we now have 400 matched exposure “B” subjects in our matched
sample.

## Obtaining the Matched Sample

As before,

``` r
match2_matches <- factor(rep(match_2$index.treated, 2))
dm2200_matched2 <- cbind(match2_matches, 
                         dm2200[c(match_2$index.control, 
                                  match_2$index.treated),])
```

How many unique subjects are in our matched sample?

``` r
dm2200_matched2 %$% n_distinct(subject)
```

    [1] 600

This match repeats each exposure A subject twice, to match up with the
400 exposure B subjects.

``` r
dm2200_matched2 %>% count(exposure)
```

    # A tibble: 2 x 2
      exposure     n
      <fct>    <int>
    1 A          400
    2 B          400

``` r
dm2200_matched2 %>% count(subject, exposure)
```

    # A tibble: 600 x 3
       subject exposure     n
       <chr>   <fct>    <int>
     1 S-0001  A            2
     2 S-0004  A            2
     3 S-0007  B            1
     4 S-0008  B            1
     5 S-0014  B            1
     6 S-0017  A            2
     7 S-0026  A            2
     8 S-0027  A            2
     9 S-0037  A            2
    10 S-0045  B            1
    # ... with 590 more rows

## Checking Covariate Balance for our 1:2 Greedy Match

### Using `bal.tab` to obtain a balance table

``` r
covs_2plus <- dm2200 %>%
    select(age, race, hisp, sex, insur, nincome,
           nhsgrad, cleve, a1c, ldl, visits, tobacco,
           statin, ace_arb, betab, depr_dx, eyeex, pneumo,
           ps, linps)

bal2 <- bal.tab(M = match_2,
                treat = dm2200$exposure,
                covs = covs_2plus, quick = FALSE,
                un = TRUE, disp.v.ratio = TRUE)
bal2
```

    Balance Measures
                        Type Diff.Un V.Ratio.Un Diff.Adj V.Ratio.Adj
    age              Contin.  0.4076     1.3195   0.0712      1.2459
    race_Asian        Binary  0.0025              0.0000            
    race_Black_AA     Binary -0.2975             -0.0225            
    race_Other        Binary  0.0320              0.0025            
    race_White        Binary  0.2630              0.0200            
    hisp              Binary  0.0320              0.0050            
    sex_M             Binary -0.2190             -0.0825            
    insur_Commercial  Binary  0.2200              0.0025            
    insur_Medicaid    Binary -0.3745             -0.0375            
    insur_Medicare    Binary  0.2715              0.0725            
    insur_Uninsured   Binary -0.1170             -0.0375            
    nincome          Contin.  0.5306     3.2086   0.0030      1.4407
    nhsgrad          Contin.  0.3958     0.9727   0.0513      1.0600
    cleve             Binary -0.3505             -0.0250            
    a1c              Contin. -0.0419     0.8190  -0.0214      0.8772
    ldl              Contin.  0.0783     0.9012   0.0037      0.9105
    visits           Contin. -0.6304     0.4602  -0.3021      0.7227
    tobacco_Current   Binary -0.3175             -0.1000            
    tobacco_Former    Binary  0.1575              0.0625            
    tobacco_Never     Binary  0.1600              0.0375            
    statin            Binary  0.0315              0.0200            
    ace_arb           Binary -0.0345             -0.0075            
    betab             Binary  0.1255              0.0350            
    depr_dx           Binary  0.1115              0.0525            
    eyeex             Binary  0.0925              0.0450            
    pneumo            Binary  0.3030              0.1500            
    ps               Contin. -2.9189     0.1576  -1.4656      0.3428
    linps            Contin. -1.7842     1.3395  -0.4322      0.4094
    
    Sample sizes
                A    B
    All       200 2000
    Matched   200  400
    Unmatched   0 1600

### Checking Rubin’s Rules 1 and 2

We’ll build a little table of the Rubin’s Rules (1 and 2) before and
after our 1:2 greedy `match_2` is applied, and compare these to the
results we found in `match_1` (the 1:1 match).

``` r
covs_for_rubin <- dm2200 %>%
    select(linps)

rubin_m2 <- bal.tab(M = match_2,
                treat = dm2200$treat,
                covs = covs_for_rubin, 
                un = TRUE, disp.v.ratio = TRUE)[1]

rubin_report_m12 <- tibble(
    status = c("Rule1", "Rule2"),
    Unmatched = c(rubin_m2$Balance$Diff.Un,
                  rubin_m2$Balance$V.Ratio.Un),
    Match1 = c(rubin_m1$Balance$Diff.Adj,
               rubin_m1$Balance$V.Ratio.Adj),
    Match2 = c(rubin_m2$Balance$Diff.Adj,
               rubin_m2$Balance$V.Ratio.Adj))

rubin_report_m12 %>% knitr::kable(digits = 2)
```

| status | Unmatched | Match1 | Match2 |
| :----- | --------: | -----: | -----: |
| Rule1  |      2.07 |   0.25 |   0.50 |
| Rule2  |      0.75 |   1.86 |   2.44 |

  - Again, we’d like to see Rule 1 as close to zero as possible, and
    definitely below 0.5 in absolute value. Unsurprisingly, when we have
    to match *two* exposure B subjects to each exposure A subject, we
    don’t get matches that are as close.
  - The Rule 2 results tell us about the variance ratio of the linear
    propensity scores. We want this to be within (0.5, 2) and ideally
    within (0.8, 1.25). Again, here the results are a bit disappointing
    in comparison to what we saw in our 1:1 match.

### Using `bal.plot` from `cobalt`

Looking at the propensity scores in each group, perhaps in mirrored
histograms, we have …

``` r
bal.plot(obj = match_2,
         treat = dm2200$exposure,
         covs = covs_2plus,
         var.name = "ps", 
         which = "both",
         sample.names = 
             c("Unmatched Sample", "match_2 Sample"),
         type = "histogram", mirror = TRUE)
```

![](matching_with_dm2200_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

### Using `love.plot` to look at Standardized Differences

``` r
love.plot(bal2, 
          threshold = .1, size = 3,
          var.order = "unadjusted",
          stats = "mean.diffs",
          stars = "raw",
          sample.names = c("Unmatched", "Matched"),
          title = "Love Plot for our 1:2 Match") +
    labs(caption = "* indicates raw mean differences (for binary variables)")
```

![](matching_with_dm2200_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

### Using `love.plot` to look at Variance Ratios

Again, the categorical variables are dropped.

``` r
love.plot(bal2, 
          threshold = .5, size = 3,
          stats = "variance.ratios",
          sample.names = c("Unmatched", "Matched"),
          title = "Variance Ratios for our 1:2 Match") 
```

![](matching_with_dm2200_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

# `match_3` 1:3 matching, with replacement with the `Matching` package

Again, we’ll match on the linear propensity score, and define our
`treat` (treatment) as occurring when `exposure` is A. But now, we’ll
match *with* replacement (which means that multiple subject with
exposure A can be matched to the same subject with exposure B) and we’ll
also match each subject with exposure A to be matched to exactly three
subjects with exposure B.

``` r
match_3 <- Match(Tr = dm2200$treat, X = dm2200$linps, 
                 M = 3, replace = TRUE, ties = FALSE,
                 estimand = "ATT")

summary(match_3)
```

``` 

Estimate...  0 
SE.........  0 
T-stat.....  NaN 
p.val......  NA 

Original number of observations..............  2200 
Original number of treated obs...............  200 
Matched number of observations...............  200 
Matched number of observations  (unweighted).  600 
```

Note that we now have 600 matched exposure “B” subjects in our matched
sample.

## Obtaining the Matched Sample

As before,

``` r
match3_matches <- factor(rep(match_3$index.treated, 2))
dm2200_matched3 <- cbind(match3_matches, 
                         dm2200[c(match_3$index.control, 
                                  match_3$index.treated),])
```

If this was being done without replacement, this would repeat each
exposure A subject three times, to match up with the 600 exposure B
subjects. But here, we have a different result.

How many unique subjects are in our matched sample?

``` r
dm2200_matched3 %$% n_distinct(subject)
```

    [1] 506

How many of those are in Exposure A?

``` r
dm2200_matched3 %>% filter(exposure == "A") %$% n_distinct(subject)
```

    [1] 200

How many of those are in Exposure B?

``` r
dm2200_matched3 %>% filter(exposure == "B") %$% n_distinct(subject)
```

    [1] 306

Among those exposure A subjects, how many times were they used in the
matches?

``` r
dm2200_matched3 %>% filter(exposure == "A") %>% 
    count(subject) %>%
    tabyl(n)
```

``` 
 n n_n percent
 3 200       1
```

Among those exposure B subjects, how many times were they used in the
matches?

``` r
dm2200_matched3 %>% filter(exposure == "B") %>% 
    count(subject) %>%
    tabyl(n)
```

``` 
  n n_n     percent
  1 203 0.663398693
  2  58 0.189542484
  3  20 0.065359477
  4   7 0.022875817
  5   4 0.013071895
  6   2 0.006535948
  7   1 0.003267974
  8   1 0.003267974
  9   2 0.006535948
 10   2 0.006535948
 13   2 0.006535948
 15   1 0.003267974
 20   1 0.003267974
 23   1 0.003267974
 24   1 0.003267974
```

## Checking Covariate Balance for our 1:3 Match

### Using `bal.tab` to obtain a balance table

``` r
covs_3plus <- dm2200 %>%
    select(age, race, hisp, sex, insur, nincome,
           nhsgrad, cleve, a1c, ldl, visits, tobacco,
           statin, ace_arb, betab, depr_dx, eyeex, pneumo,
           ps, linps)

bal3 <- bal.tab(M = match_3,
                treat = dm2200$exposure,
                covs = covs_3plus, quick = FALSE,
                un = TRUE, disp.v.ratio = TRUE)
bal3
```

    Balance Measures
                        Type Diff.Un V.Ratio.Un Diff.Adj V.Ratio.Adj
    age              Contin.  0.4076     1.3195  -0.0772      1.2937
    race_Asian        Binary  0.0025             -0.0017            
    race_Black_AA     Binary -0.2975              0.0267            
    race_Other        Binary  0.0320             -0.0033            
    race_White        Binary  0.2630             -0.0217            
    hisp              Binary  0.0320             -0.0067            
    sex_M             Binary -0.2190              0.0117            
    insur_Commercial  Binary  0.2200              0.0017            
    insur_Medicaid    Binary -0.3745             -0.0317            
    insur_Medicare    Binary  0.2715              0.0200            
    insur_Uninsured   Binary -0.1170              0.0100            
    nincome          Contin.  0.5306     3.2086   0.0325      1.3808
    nhsgrad          Contin.  0.3958     0.9727  -0.0241      0.8381
    cleve             Binary -0.3505             -0.0017            
    a1c              Contin. -0.0419     0.8190  -0.0541      1.0381
    ldl              Contin.  0.0783     0.9012   0.0820      1.0192
    visits           Contin. -0.6304     0.4602  -0.0609      0.8131
    tobacco_Current   Binary -0.3175              0.0017            
    tobacco_Former    Binary  0.1575              0.0233            
    tobacco_Never     Binary  0.1600             -0.0250            
    statin            Binary  0.0315             -0.0017            
    ace_arb           Binary -0.0345             -0.0033            
    betab             Binary  0.1255             -0.0133            
    depr_dx           Binary  0.1115              0.0100            
    eyeex             Binary  0.0925             -0.0050            
    pneumo            Binary  0.3030              0.0100            
    ps               Contin. -2.9189     0.1576  -0.0384      0.9511
    linps            Contin. -1.7842     1.3395  -0.0267      0.8893
    
    Sample sizes
                           A       B
    All                  200 2000.00
    Matched (ESS)        200  104.59
    Matched (Unweighted) 200  306.00
    Unmatched              0 1694.00

### Checking Rubin’s Rules 1 and 2

We’ll build a little table of the Rubin’s Rules (1 and 2) before and
after our 1:2 greedy `match_2` is applied, and compare these to the
results we found in `match_1` (the 1:1 match).

``` r
covs_for_rubin <- dm2200 %>%
    select(linps)

rubin_m3 <- bal.tab(M = match_3,
                treat = dm2200$treat,
                covs = covs_for_rubin, 
                un = TRUE, disp.v.ratio = TRUE)[1]

rubin_report_m123 <- tibble(
    status = c("Rule1", "Rule2"),
    Unmatched = c(rubin_m2$Balance$Diff.Un,
                  rubin_m2$Balance$V.Ratio.Un),
    Match1 = c(rubin_m1$Balance$Diff.Adj,
               rubin_m1$Balance$V.Ratio.Adj),
    Match2 = c(rubin_m2$Balance$Diff.Adj,
               rubin_m2$Balance$V.Ratio.Adj),
    Match3 = c(rubin_m3$Balance$Diff.Adj,
               rubin_m3$Balance$V.Ratio.Adj))


rubin_report_m123 %>% knitr::kable(digits = 2)
```

| status | Unmatched | Match1 | Match2 | Match3 |
| :----- | --------: | -----: | -----: | -----: |
| Rule1  |      2.07 |   0.25 |   0.50 |   0.03 |
| Rule2  |      0.75 |   1.86 |   2.44 |   1.12 |

  - Again, we’d like to see Rule 1 results as close to zero as possible,
    and definitely below 0.5 in absolute value.
  - In Rule 2, we want the variance ratio of the linear propensity
    scores to be within (0.5, 2) and ideally within (0.8, 1.25).
  - It appears that (in these data) allowing the same exposure B subject
    to be used for multiple matches (matching with replacement) more
    than makes up for the fact that matching 3 exposure B’s for each
    exposure A (1:3 matching) is a tougher job than pair (1:1) matching,
    as seen in the results for Rubin’s Rule 1 and Rule 2.

### Using `bal.plot` from `cobalt`

Looking at the propensity scores in each group, perhaps in mirrored
histograms, we have …

``` r
bal.plot(obj = match_3,
         treat = dm2200$exposure,
         covs = covs_3plus,
         var.name = "ps", 
         which = "both",
         sample.names = 
             c("Unmatched Sample", "match_3 Sample"),
         type = "histogram", mirror = TRUE)
```

![](matching_with_dm2200_files/figure-gfm/unnamed-chunk-38-1.png)<!-- -->

### Using `love.plot` to look at Standardized Differences

``` r
love.plot(bal3, 
          threshold = .1, size = 3,
          var.order = "unadjusted",
          stats = "mean.diffs",
          stars = "raw",
          abs = TRUE,
          sample.names = c("Unmatched", "Matched"),
          title = "Love Plot of |Mean Differences| for our 1:3 Match") +
    labs(caption = "* indicates raw mean differences (for binary variables)")
```

![](matching_with_dm2200_files/figure-gfm/unnamed-chunk-39-1.png)<!-- -->

### Using `love.plot` to look at Variance Ratios

Again, the categorical variables are dropped.

``` r
love.plot(bal3, 
          threshold = .5, size = 3,
          stats = "variance.ratios",
          sample.names = c("Unmatched", "Matched"),
          title = "Variance Ratios for our 1:3 Match") 
```

![](matching_with_dm2200_files/figure-gfm/unnamed-chunk-40-1.png)<!-- -->

# `match_4` Caliper Matching (1:1 without replacement) with the `Matching` package

The `Match` function in the `Matching` package allows you to specify a
caliper. From the `Matching` help file:

  - A caliper is the maximum acceptable distance (on a covariate) which
    we are willing to accept in any match. Observations for which we
    cannot find a match within the caliper are dropped.Dropping
    observations generally changes the quantity being estimated.
  - The caliper is interpreted to be in standardized units. For example,
    caliper=.25 means that all matches not equal to or within .25
    standard deviations of each covariate in X are dropped, and not
    matched.
      - If a scalar caliper is provided to the `caliper` setting in the
        `Match` function, this caliper is used for all covariates in X.
      - If a vector of calipers is provided, a caliper value should be
        provided for each covariate in X.

We’ll again perform a 1:1 match without replacement, but now we’ll do so
while only accepting matches where the linear propensity score of each
match is within 0.2 standard deviations of the linear PS.

``` r
match_4 <- Match(Tr = dm2200$treat, X = dm2200$linps, 
                 M = 1, replace = FALSE, ties = FALSE,
                 caliper = 0.2, estimand = "ATT")

summary(match_4)
```

``` 

Estimate...  0 
SE.........  0 
T-stat.....  NaN 
p.val......  NA 

Original number of observations..............  2200 
Original number of treated obs...............  200 
Matched number of observations...............  162 
Matched number of observations  (unweighted).  162 

Caliper (SDs)........................................   0.2 
Number of obs dropped by 'exact' or 'caliper'  38 
```

Note that we have now dropped 38 of the exposure “A” subjects, and
reduced our sample to the 168 remaining exposure “A” subjects, who are
paired with 162 unique matched exposure “B” subjects in our matched
sample.

## Obtaining the Matched Sample

As before,

``` r
match4_matches <- factor(rep(match_4$index.treated, 2))
dm2200_matched4 <- cbind(match4_matches, 
                         dm2200[c(match_4$index.control, 
                                  match_4$index.treated),])
```

How many unique subjects are in our matched sample?

``` r
dm2200_matched4 %$% n_distinct(subject)
```

    [1] 324

This match includes 162 pairs so 324 subjects, since we’ve done matching
without replacement.

``` r
dm2200_matched4 %>% count(exposure)
```

    # A tibble: 2 x 2
      exposure     n
      <fct>    <int>
    1 A          162
    2 B          162

## Checking Covariate Balance for our 1:1 Caliper Match

### Using `bal.tab` to obtain a balance table

``` r
covs_4plus <- dm2200 %>%
    select(age, race, hisp, sex, insur, nincome,
           nhsgrad, cleve, a1c, ldl, visits, tobacco,
           statin, ace_arb, betab, depr_dx, eyeex, pneumo,
           ps, linps)

bal4 <- bal.tab(M = match_4,
                treat = dm2200$exposure,
                covs = covs_4plus, quick = FALSE,
                un = TRUE, disp.v.ratio = TRUE)
bal4
```

    Balance Measures
                        Type Diff.Un V.Ratio.Un Diff.Adj V.Ratio.Adj
    age              Contin.  0.4076     1.3195  -0.0911      1.2192
    race_Asian        Binary  0.0025              0.0000            
    race_Black_AA     Binary -0.2975              0.0062            
    race_Other        Binary  0.0320             -0.0062            
    race_White        Binary  0.2630              0.0000            
    hisp              Binary  0.0320              0.0000            
    sex_M             Binary -0.2190             -0.0370            
    insur_Commercial  Binary  0.2200              0.0185            
    insur_Medicaid    Binary -0.3745              0.0247            
    insur_Medicare    Binary  0.2715             -0.0370            
    insur_Uninsured   Binary -0.1170             -0.0062            
    nincome          Contin.  0.5306     3.2086  -0.0572      0.8507
    nhsgrad          Contin.  0.3958     0.9727  -0.0974      1.0586
    cleve             Binary -0.3505              0.0247            
    a1c              Contin. -0.0419     0.8190   0.0115      0.9806
    ldl              Contin.  0.0783     0.9012  -0.0141      0.9002
    visits           Contin. -0.6304     0.4602  -0.0382      1.0651
    tobacco_Current   Binary -0.3175              0.0000            
    tobacco_Former    Binary  0.1575              0.0370            
    tobacco_Never     Binary  0.1600             -0.0370            
    statin            Binary  0.0315             -0.0062            
    ace_arb           Binary -0.0345             -0.0247            
    betab             Binary  0.1255              0.0432            
    depr_dx           Binary  0.1115              0.0000            
    eyeex             Binary  0.0925             -0.0185            
    pneumo            Binary  0.3030              0.0000            
    ps               Contin. -2.9189     0.1576  -0.0306      0.9598
    linps            Contin. -1.7842     1.3395  -0.0072      0.9775
    
    Sample sizes
                A    B
    All       200 2000
    Matched   162  162
    Unmatched   0 1838
    Discarded  38    0

### Checking Rubin’s Rules 1 and 2

We’ll build a little table of the Rubin’s Rules (1 and 2) before and
after our 1:2 greedy `match_4` is applied, and compare these to the
results we found in `match_1` (the 1:1 match).

``` r
covs_for_rubin <- dm2200 %>%
    select(linps)

rubin_m4 <- bal.tab(M = match_4,
                treat = dm2200$treat,
                covs = covs_for_rubin, 
                un = TRUE, disp.v.ratio = TRUE)[1]

rubin_report_m1234 <- tibble(
    status = c("Rule1", "Rule2"),
    Unmatched = c(rubin_m2$Balance$Diff.Un,
                  rubin_m2$Balance$V.Ratio.Un),
    Match1 = c(rubin_m1$Balance$Diff.Adj,
               rubin_m1$Balance$V.Ratio.Adj),
    Match2 = c(rubin_m2$Balance$Diff.Adj,
               rubin_m2$Balance$V.Ratio.Adj),
    Match3 = c(rubin_m3$Balance$Diff.Adj,
               rubin_m3$Balance$V.Ratio.Adj),
    Match4 = c(rubin_m4$Balance$Diff.Adj,
               rubin_m4$Balance$V.Ratio.Adj))

rubin_report_m1234 %>% knitr::kable(digits = 2)
```

| status | Unmatched | Match1 | Match2 | Match3 | Match4 |
| :----- | --------: | -----: | -----: | -----: | -----: |
| Rule1  |      2.07 |   0.25 |   0.50 |   0.03 |   0.01 |
| Rule2  |      0.75 |   1.86 |   2.44 |   1.12 |   1.02 |

  - This approach produces an exceptionally strong match in terms of
    balance, with Rubin’s Rule 1 being very close to 0, and Rule 2 being
    very close to 1.
  - Unfortunately, we’ve only done this by dropping the 38 “hardest to
    match” subjects receiving exposure “A”.

### Using `bal.plot` from `cobalt`

Looking at the propensity scores in each group, perhaps in mirrored
histograms, we have …

``` r
bal.plot(obj = match_4,
         treat = dm2200$exposure,
         covs = covs_4plus,
         var.name = "ps", 
         which = "both",
         sample.names = 
             c("Unmatched Sample", "match_4 Sample"),
         type = "histogram", mirror = TRUE)
```

![](matching_with_dm2200_files/figure-gfm/unnamed-chunk-47-1.png)<!-- -->

### Using `love.plot` to look at Standardized Differences

``` r
love.plot(bal4, 
          threshold = .1, size = 3,
          var.order = "unadjusted",
          stats = "mean.diffs",
          stars = "raw",
          sample.names = c("Unmatched", "Matched"),
          title = "Love Plot for our 1:1 Caliper Match") +
    labs(caption = "* indicates raw mean differences (for binary variables)")
```

![](matching_with_dm2200_files/figure-gfm/unnamed-chunk-48-1.png)<!-- -->

### Using `love.plot` to look at Variance Ratios

Again, the categorical variables are dropped.

``` r
love.plot(bal4, 
          threshold = .5, size = 3,
          stats = "variance.ratios",
          sample.names = c("Unmatched", "Matched"),
          title = "Variance Ratios for our 1:1 Caliper Match") 
```

![](matching_with_dm2200_files/figure-gfm/unnamed-chunk-49-1.png)<!-- -->

# `match_5` 1:1 Nearest Neighbor Matching with `Matchit`

The `MatchIt` package implements the suggestions of Ho, Imai, King, and
Stuart (2007) for improving parametric statistical models by
preprocessing data with nonparametric matching methods. Read more about
`MatchIt` at <https://gking.harvard.edu/matchit>.

With `MatchIt`, the matching is done using the `matchit(treat ~ X, ...)`
function, where `treat` is the vector of treatment assignments and `X`
are the covariates to be used in the matching.

We’ll start our exploration of the `MatchIt` approach to developing
matches with nearest neighbor matching, which (quoting the `MatchIt`
manual…)

> … selects the *r* (default = 1, specified by the `ratio` option) best
> control matches for each individual in the treated group, using a
> distance measure specified by the `distance` option (default = logit).
> Matches are chosen for each treated unit one at a time, with the order
> specified by the `m.order` command (default=largest to smallest). At
> each matching step we choose the control unit that is not yet matched
> but is closest to the treated unit on the distance measure.

The full syntax (with default choices indicated) is

    matchit(formula, data=NULL, discard=0, exact=FALSE, 
             replace=FALSE, ratio=1, model="logit", 
             reestimate=FALSE, nearest=TRUE, m.order=2, 
             caliper=0, calclosest=FALSE, mahvars=NULL, 
             subclass=0, sub.by="treat", counter=TRUE, 
             full=FALSE, full.options=list(), ...)

Here, we’ll match on the linear propensity score (by specifying the
`distance` to use the logistic link with the linear propensity score via
`"linear.logit"`), and we’ll perform 1:1 nearest neighbor matching
without replacement using the default ordering (largest to smallest).
Since we’ve already seen that greedy 1:1 matching without replacement
doesn’t work well in this setting, we’re not expecting a strong result
in terms of balance here, either.

``` r
dm2200 <- dm2200 %>%
    mutate(treat = as.logical(exposure == "A"))

covs_1 <- dm2200 %>%
    select(age, race, hisp, sex, insur, nincome,
           nhsgrad, cleve, a1c, ldl, visits, tobacco,
           statin, ace_arb, betab, depr_dx, eyeex, pneumo)

match_5 <- matchit(f.build("treat", covs_1), data = dm2200,
                   distance = "linear.logit", method = "nearest", 
                   ratio = 1, replace = FALSE)

match_5
```

``` 

Call: 
matchit(formula = f.build("treat", covs_1), data = dm2200, method = "nearest", 
    distance = "linear.logit", ratio = 1, replace = FALSE)

Sample sizes:
          Control Treated
All          2000     200
Matched       200     200
Unmatched    1800       0
Discarded       0       0
```

## Obtaining the Matched Sample

There is just one tricky part to doing this in `MatchIt`. The main work
can be done with a very simple command.

``` r
dm2200_matched5 <- match.data(match_5)
```

This leaves only the job of creating a matching number, for which we
have to develop some additional R code.

``` r
# Thanks to Robert McDonald at Mayo
# https://lists.gking.harvard.edu/pipermail/matchit/2012-April/000458.html

len <- dim(dm2200)[1]
matchx <- rep(NA,len)
len2 <- length(match_5$match.matrix)
count <- 1
for(i in 1:len2){
    
    match1 <- match_5$match.matrix[i]
    match2 <- row.names(match_5$match.matrix)[i]
    
    if(!is.na(match1)){
    matchx[as.numeric(match1)] <- count
    matchx[as.numeric(match2)] <- count
    count <- count+1}
    
}

dm2200_matched5 <- dm2200_matched5 %>%
    mutate(match5_matches = 
               matchx[as.numeric
                      (row.names(dm2200_matched5))]) %>%
    mutate(match5_matches = factor(match5_matches))
```

How many unique subjects are in our matched sample?

``` r
dm2200_matched5 %$% n_distinct(subject)
```

    [1] 400

This match includes 200 pairs so 400 subjects, as we’ve done matching
without replacement.

``` r
dm2200_matched5 %>% count(exposure)
```

    # A tibble: 2 x 2
      exposure     n
      <fct>    <int>
    1 A          200
    2 B          200

Do we have as many matching numbers as treated subjects?

``` r
dm2200_matched5 %$% n_distinct(match5_matches)
```

    [1] 200

## Checking Covariate Balance for our 1:1 Nearest Neighbor Match

### Default Numerical Balance Summary from `MatchIt`

``` r
summary(match_5)
```

``` 

Call:
matchit(formula = f.build("treat", covs_1), data = dm2200, method = "nearest", 
    distance = "linear.logit", ratio = 1, replace = FALSE)

Summary of balance for all data:
               Means Treated Means Control SD Control   Mean Diff eQQ Med
distance             -0.6577       -4.0979     1.9281      3.4402    3.38
age                  54.8950       58.9000     9.8257     -4.0050    4.00
raceAsian             0.0100        0.0125     0.1111     -0.0025    0.00
raceBlack_AA          0.8300        0.5325     0.4991      0.2975    0.00
raceOther             0.0050        0.0370     0.1888     -0.0320    0.00
raceWhite             0.1550        0.4180     0.4934     -0.2630    0.00
hisp                  0.0200        0.0520     0.2221     -0.0320    0.00
sexM                  0.6400        0.4210     0.4938      0.2190    0.00
insurMedicaid         0.6250        0.2505     0.4334      0.3745    0.00
insurMedicare         0.1850        0.4565     0.4982     -0.2715    0.00
insurUninsured        0.1550        0.0380     0.1912      0.1170    0.00
nincome           26927.5000    37992.5500 20854.1979 -11065.0500 9800.00
nhsgrad              77.7350       82.2260    11.3454     -4.4910    4.50
cleve                 0.9050        0.5545     0.4971      0.3505    0.00
a1c                   7.7655        7.6868     1.8808      0.0788    0.10
ldl                  94.3450       97.2160    36.6583     -2.8710    3.00
visits                4.8100        3.6885     1.7791      1.1215    1.00
tobaccoFormer         0.2300        0.3875     0.4873     -0.1575    0.00
tobaccoNever          0.2300        0.3900     0.4879     -0.1600    0.00
statin                0.7700        0.8015     0.3990     -0.0315    0.00
ace_arb               0.8000        0.7655     0.4238      0.0345    0.00
betab                 0.2650        0.3905     0.4880     -0.1255    0.00
depr_dx               0.2650        0.3765     0.4846     -0.1115    0.00
eyeex                 0.5200        0.6125     0.4873     -0.0925    0.00
pneumo                0.5800        0.8830     0.3215     -0.3030    0.00
                 eQQ Mean    eQQ Max
distance           3.4452     4.5334
age                3.9400     7.0000
raceAsian          0.0050     1.0000
raceBlack_AA       0.3000     1.0000
raceOther          0.0350     1.0000
raceWhite          0.2650     1.0000
hisp               0.0350     1.0000
sexM               0.2200     1.0000
insurMedicaid      0.3750     1.0000
insurMedicare      0.2700     1.0000
insurUninsured     0.1150     1.0000
nincome        11599.5000 74000.0000
nhsgrad            4.4100     8.0000
cleve              0.3500     1.0000
a1c                0.1530     0.8000
ldl                3.7750    27.0000
visits             1.1050     4.0000
tobaccoFormer      0.1600     1.0000
tobaccoNever       0.1600     1.0000
statin             0.0300     1.0000
ace_arb            0.0350     1.0000
betab              0.1250     1.0000
depr_dx            0.1100     1.0000
eyeex              0.0900     1.0000
pneumo             0.3000     1.0000


Summary of balance for matched data:
               Means Treated Means Control SD Control Mean Diff   eQQ Med
distance             -0.6577       -1.0679     1.2219    0.4102    0.1049
age                  54.8950       54.2500     9.5285    0.6450    1.0000
raceAsian             0.0100        0.0100     0.0997    0.0000    0.0000
raceBlack_AA          0.8300        0.8300     0.3766    0.0000    0.0000
raceOther             0.0050        0.0000     0.0000    0.0050    0.0000
raceWhite             0.1550        0.1600     0.3675   -0.0050    0.0000
hisp                  0.0200        0.0150     0.1219    0.0050    0.0000
sexM                  0.6400        0.6000     0.4911    0.0400    0.0000
insurMedicaid         0.6250        0.6450     0.4797   -0.0200    0.0000
insurMedicare         0.1850        0.1800     0.3852    0.0050    0.0000
insurUninsured        0.1550        0.1350     0.3426    0.0200    0.0000
nincome           26927.5000    25385.5000 11052.2078 1542.0000 1700.0000
nhsgrad              77.7350       77.4650    11.0726    0.2700    1.0000
cleve                 0.9050        0.9050     0.2940    0.0000    0.0000
a1c                   7.7655        7.8350     2.1865   -0.0695    0.1000
ldl                  94.3450       94.4900    37.0468   -0.1450    1.0000
visits                4.8100        4.3400     2.2870    0.4700    0.0000
tobaccoFormer         0.2300        0.2650     0.4424   -0.0350    0.0000
tobaccoNever          0.2300        0.2400     0.4282   -0.0100    0.0000
statin                0.7700        0.7150     0.4525    0.0550    0.0000
ace_arb               0.8000        0.7650     0.4251    0.0350    0.0000
betab                 0.2650        0.3100     0.4637   -0.0450    0.0000
depr_dx               0.2650        0.2850     0.4525   -0.0200    0.0000
eyeex                 0.5200        0.5300     0.5004   -0.0100    0.0000
pneumo                0.5800        0.6600     0.4749   -0.0800    0.0000
                eQQ Mean    eQQ Max
distance          0.4114     2.3376
age               0.9950     5.0000
raceAsian         0.0000     0.0000
raceBlack_AA      0.0000     0.0000
raceOther         0.0050     1.0000
raceWhite         0.0050     1.0000
hisp              0.0050     1.0000
sexM              0.0400     1.0000
insurMedicaid     0.0200     1.0000
insurMedicare     0.0050     1.0000
insurUninsured    0.0200     1.0000
nincome        2205.0000 19500.0000
nhsgrad           1.1400     8.0000
cleve             0.0000     0.0000
a1c               0.1765     0.8000
ldl               1.9650    12.0000
visits            0.4900     3.0000
tobaccoFormer     0.0350     1.0000
tobaccoNever      0.0100     1.0000
statin            0.0550     1.0000
ace_arb           0.0350     1.0000
betab             0.0450     1.0000
depr_dx           0.0200     1.0000
eyeex             0.0100     1.0000
pneumo            0.0800     1.0000

Percent Balance Improvement:
               Mean Diff.  eQQ Med eQQ Mean  eQQ Max
distance          88.0768  96.8975  88.0596  48.4367
age               83.8951  75.0000  74.7462  28.5714
raceAsian        100.0000   0.0000 100.0000 100.0000
raceBlack_AA     100.0000   0.0000 100.0000 100.0000
raceOther         84.3750   0.0000  85.7143   0.0000
raceWhite         98.0989   0.0000  98.1132   0.0000
hisp              84.3750   0.0000  85.7143   0.0000
sexM              81.7352   0.0000  81.8182   0.0000
insurMedicaid     94.6595   0.0000  94.6667   0.0000
insurMedicare     98.1584   0.0000  98.1481   0.0000
insurUninsured    82.9060   0.0000  82.6087   0.0000
nincome           86.0642  82.6531  80.9906  73.6486
nhsgrad           93.9880  77.7778  74.1497   0.0000
cleve            100.0000   0.0000 100.0000 100.0000
a1c               11.7460   0.0000 -15.3595   0.0000
ldl               94.9495  66.6667  47.9470  55.5556
visits            58.0918 100.0000  55.6561  25.0000
tobaccoFormer     77.7778   0.0000  78.1250   0.0000
tobaccoNever      93.7500   0.0000  93.7500   0.0000
statin           -74.6032   0.0000 -83.3333   0.0000
ace_arb           -1.4493   0.0000   0.0000   0.0000
betab             64.1434   0.0000  64.0000   0.0000
depr_dx           82.0628   0.0000  81.8182   0.0000
eyeex             89.1892   0.0000  88.8889   0.0000
pneumo            73.5974   0.0000  73.3333   0.0000

Sample sizes:
          Control Treated
All          2000     200
Matched       200     200
Unmatched    1800       0
Discarded       0       0
```

### Using `bal.tab` to obtain a balance table

``` r
bal5 <- bal.tab(match_5, un = TRUE, disp.v.ratio = TRUE)

bal5
```

    Call
     matchit(formula = f.build("treat", covs_1), data = dm2200, method = "nearest", 
        distance = "linear.logit", ratio = 1, replace = FALSE)
    
    Balance Measures
                         Type Diff.Un V.Ratio.Un Diff.Adj V.Ratio.Adj
    distance         Distance  2.0650     0.7466   0.2462      1.8590
    age               Contin. -0.4682     0.7579   0.0754      0.8059
    race_Asian         Binary -0.0025              0.0000            
    race_Black_AA      Binary  0.2975              0.0000            
    race_Other         Binary -0.0320              0.0050            
    race_White         Binary -0.2630             -0.0050            
    hisp               Binary -0.0320              0.0050            
    sex_M              Binary  0.2190              0.0400            
    insur_Commercial   Binary -0.2200             -0.0050            
    insur_Medicaid     Binary  0.3745             -0.0200            
    insur_Medicare     Binary -0.2715              0.0050            
    insur_Uninsured    Binary  0.1170              0.0200            
    nincome           Contin. -0.9504     0.3117   0.1324      1.1096
    nhsgrad           Contin. -0.3904     1.0281   0.0235      1.0794
    cleve              Binary  0.3505              0.0000            
    a1c               Contin.  0.0379     1.2209  -0.0334      0.9035
    ldl               Contin. -0.0743     1.1096  -0.0038      1.0865
    visits            Contin.  0.4276     2.1732   0.1792      1.3150
    tobacco_Current    Binary  0.3175              0.0450            
    tobacco_Former     Binary -0.1575             -0.0350            
    tobacco_Never      Binary -0.1600             -0.0100            
    statin             Binary -0.0315              0.0550            
    ace_arb            Binary  0.0345              0.0350            
    betab              Binary -0.1255             -0.0450            
    depr_dx            Binary -0.1115             -0.0200            
    eyeex              Binary -0.0925             -0.0100            
    pneumo             Binary -0.3030             -0.0800            
    
    Sample sizes
              FALSE TRUE
    All        2000  200
    Matched     200  200
    Unmatched  1800    0

### Checking Rubin’s Rules 1 and 2

We’ll build a little table of the Rubin’s Rules (1 and 2) results before
and after `match_5` is applied.

``` r
rubin_report_m5 <- tibble(
    status = c("Rule1", "Rule2"),
    Unmatched = c(bal5$Balance$Diff.Un[1],
                  bal5$Balance$V.Ratio.Un[1]),
    Match5 = c(bal5$Balance$Diff.Adj[1],
               bal5$Balance$V.Ratio.Adj[1]))

rubin_report_m5 %>% knitr::kable(digits = 2)
```

| status | Unmatched | Match5 |
| :----- | --------: | -----: |
| Rule1  |      2.07 |   0.25 |
| Rule2  |      0.75 |   1.86 |

  - As we’d expect based on our previous greedy 1:1 match, this isn’t
    great sufficiently strong balance.

### Using `bal.plot` from `cobalt`

Looking at the linear propensity scores in each group, in mirrored
histograms, we have …

``` r
bal.plot(obj = match_5,
         var.name = "distance", 
         which = "both",
         sample.names = 
             c("Unmatched Sample", "match_5 Sample"),
         type = "histogram", mirror = TRUE)
```

![](matching_with_dm2200_files/figure-gfm/unnamed-chunk-59-1.png)<!-- -->

### Using `love.plot` to look at Standardized Differences

``` r
love.plot(bal5, 
          threshold = .1, size = 3,
          var.order = "unadjusted",
          stats = "mean.diffs",
          stars = "raw",
          sample.names = c("Unmatched", "Matched"),
          title = "Love Plot for our 1:1 Nearest Neighbor Match") +
    labs(subtitle = "distance = linear propensity score",
         caption = "* indicates raw mean differences (for binary variables)")
```

![](matching_with_dm2200_files/figure-gfm/unnamed-chunk-60-1.png)<!-- -->

### Using `love.plot` to look at Variance Ratios

Again, the categorical variables are dropped.

``` r
love.plot(bal5, 
          threshold = .5, size = 3,
          stats = "variance.ratios",
          sample.names = c("Unmatched", "Matched"),
          title = "Variance Ratios for our 1:1 Nearest Neighbor Match") +
    labs(subtitle = "distance = linear propensity score",
         caption = "* indicates raw mean differences (for binary variables)")
```

![](matching_with_dm2200_files/figure-gfm/unnamed-chunk-61-1.png)<!-- -->

# `match_6` 1:1 Nearest Neighbor Caliper Matching with `Matchit`

As in `match_5`, we’ll match on the linear propensity score (by
specifying the `distance` to use the logistic link with the linear
propensity score via `"linear.logit"`), and we’ll perform 1:1 nearest
neighbor matching without replacement using the default ordering
(largest to smallest), but we’ll add a some arguments to build a
specific kind caliper match. Specifically, we’ll require our matches to
be within 0.25 standard deviations of each other on the linear
propensity score.

``` r
match_6 <- matchit(f.build("treat", covs_1), data = dm2200,
                   distance = "linear.logit", 
                   method = "nearest", caliper = 0.25,
                   ratio = 1, replace = FALSE)

match_6
```

``` 

Call: 
matchit(formula = f.build("treat", covs_1), data = dm2200, method = "nearest", 
    distance = "linear.logit", caliper = 0.25, ratio = 1, replace = FALSE)

Sample sizes:
          Control Treated
All          2000     200
Matched       172     172
Unmatched    1828      28
Discarded       0       0
```

## Obtaining the Matched Sample

There is just one tricky part to doing this in `MatchIt`. The main work
can be done with a very simple command.

``` r
dm2200_matched6 <- match.data(match_6)
```

This leaves only the job of creating a matching number, for which we
have to develop some additional R code.

``` r
# Thanks to Robert McDonald at Mayo
# https://lists.gking.harvard.edu/pipermail/matchit/2012-April/000458.html

len <- dim(dm2200)[1]
matchx <- rep(NA,len)
len2 <- length(match_6$match.matrix)
count <- 1
for(i in 1:len2){
    
    match1 <- match_6$match.matrix[i]
    match2 <- row.names(match_6$match.matrix)[i]
    
    if(!is.na(match1)){
    matchx[as.numeric(match1)] <- count
    matchx[as.numeric(match2)] <- count
    count <- count+1}
    
}

dm2200_matched6 <- dm2200_matched6 %>%
    mutate(match6_matches = 
               matchx[as.numeric
                      (row.names(dm2200_matched6))]) %>%
    mutate(match6_matches = factor(match6_matches))
```

How many unique subjects are in our matched sample?

``` r
dm2200_matched6 %$% n_distinct(subject)
```

    [1] 344

This match includes 172 pairs so 344 subjects, as we’ve done matching
without replacement.

``` r
dm2200_matched6 %>% count(exposure)
```

    # A tibble: 2 x 2
      exposure     n
      <fct>    <int>
    1 A          172
    2 B          172

Do we have as many matching numbers as treated subjects?

``` r
dm2200_matched6 %$% n_distinct(match6_matches)
```

    [1] 172

## Checking Covariate Balance for our 1:1 Nearest Neighbor Caliper Match

### Using `bal.tab` to obtain a balance table

``` r
bal6 <- bal.tab(match_6, un = TRUE, disp.v.ratio = TRUE)
```

    Note: s.d.denom not specified; assuming pooled.

``` r
bal6
```

    Call
     matchit(formula = f.build("treat", covs_1), data = dm2200, method = "nearest", 
        distance = "linear.logit", caliper = 0.25, ratio = 1, replace = FALSE)
    
    Balance Measures
                         Type Diff.Un V.Ratio.Un Diff.Adj V.Ratio.Adj
    distance         Distance  1.9093     0.7466   0.0963      1.2511
    age               Contin. -0.4348     0.7579   0.0309      1.0437
    race_Asian         Binary -0.0025              0.0000            
    race_Black_AA      Binary  0.2975              0.0000            
    race_Other         Binary -0.0320              0.0000            
    race_White         Binary -0.2630              0.0000            
    hisp               Binary -0.0320             -0.0233            
    sex_M              Binary  0.2190             -0.0291            
    insur_Commercial   Binary -0.2200             -0.0058            
    insur_Medicaid     Binary  0.3745             -0.0116            
    insur_Medicare     Binary -0.2715              0.0174            
    insur_Uninsured    Binary  0.1170              0.0000            
    nincome           Contin. -0.6552     0.3117  -0.0486      0.5132
    nhsgrad           Contin. -0.3931     1.0281   0.0997      0.7999
    cleve              Binary  0.3505              0.0174            
    a1c               Contin.  0.0397     1.2209   0.0282      0.9961
    ldl               Contin. -0.0763     1.1096  -0.0604      1.0657
    visits            Contin.  0.5005     2.1732   0.1920      1.2944
    tobacco_Current    Binary  0.3175              0.0174            
    tobacco_Former     Binary -0.1575             -0.0174            
    tobacco_Never      Binary -0.1600              0.0000            
    statin             Binary -0.0315              0.0233            
    ace_arb            Binary  0.0345              0.0291            
    betab              Binary -0.1255             -0.0233            
    depr_dx            Binary -0.1115              0.0465            
    eyeex              Binary -0.0925              0.0116            
    pneumo             Binary -0.3030             -0.0349            
    
    Sample sizes
              FALSE TRUE
    All        2000  200
    Matched     172  172
    Unmatched  1828   28

### Checking Rubin’s Rules 1 and 2

We’ll build a little table of the Rubin’s Rules (1 and 2) results before
and after `match_6` is applied, and compare these to the `match_5`
results.

``` r
rubin_report_m56 <- tibble(
    status = c("Rule1", "Rule2"),
    Unmatched = c(bal6$Balance$Diff.Un[1],
                  bal6$Balance$V.Ratio.Un[1]),
    Match5 = c(bal5$Balance$Diff.Adj[1],
               bal5$Balance$V.Ratio.Adj[1]),
    Match6 = c(bal6$Balance$Diff.Adj[1],
               bal6$Balance$V.Ratio.Adj[1])
    )

rubin_report_m56 %>% knitr::kable(digits = 2)
```

| status | Unmatched | Match5 | Match6 |
| :----- | --------: | -----: | -----: |
| Rule1  |      1.91 |   0.25 |   0.10 |
| Rule2  |      0.75 |   1.86 |   1.25 |

  - So the caliper appears to help quite a bit to improve the results
    that we saw in 1:1 nearest neighbor matching without replacement, at
    the cost of not including 28 of the treated subjects in the match.

### Using `bal.plot` from `cobalt`

Looking at the linear propensity scores in each group, in mirrored
histograms, we have …

``` r
bal.plot(obj = match_6,
         var.name = "distance", 
         which = "both",
         sample.names = 
             c("Unmatched Sample", "match_6 Sample"),
         type = "histogram", mirror = TRUE)
```

![](matching_with_dm2200_files/figure-gfm/unnamed-chunk-70-1.png)<!-- -->

### Using `love.plot` to look at Standardized Differences

``` r
love.plot(bal6, 
          threshold = .1, size = 3,
          var.order = "unadjusted",
          stats = "mean.diffs",
          stars = "raw",
          sample.names = c("Unmatched", "Matched"),
          title = "Love Plot for 1:1 Nearest Neighbor Caliper Match") +
    labs(subtitle = "distance = linear propensity score",
         caption = "* indicates raw mean differences (for binary variables)")
```

![](matching_with_dm2200_files/figure-gfm/unnamed-chunk-71-1.png)<!-- -->

### Using `love.plot` to look at Variance Ratios

Again, the categorical variables are dropped.

``` r
love.plot(bal6, 
          threshold = .5, size = 3,
          stats = "variance.ratios",
          sample.names = c("Unmatched", "Matched"),
          title = "Variance Ratios for 1:1 Nearest Neighbor Caliper Match") +
    labs(subtitle = "distance = linear propensity score",
         caption = "* indicates raw mean differences (for binary variables)")
```

![](matching_with_dm2200_files/figure-gfm/unnamed-chunk-72-1.png)<!-- -->

# `match_7` 1:1 Optimal Matching with `Matchit`

As in `match_5`, we’ll match on the linear propensity score (by
specifying the `distance` to use the logistic link with the linear
propensity score via `"linear.logit"`), but now we’ll use an optimal 1:1
matching without replacement.

``` r
match_7 <- matchit(f.build("treat", covs_1), data = dm2200,
                   distance = "linear.logit", 
                   method = "optimal", 
                   ratio = 1, replace = FALSE)
```

    Warning in optmatch::fullmatch(d, min.controls = ratio, max.controls = ratio, : Without 'data' argument the order of the match is not guaranteed
        to be the same as your original data.

``` r
match_7
```

``` 

Call: 
matchit(formula = f.build("treat", covs_1), data = dm2200, method = "optimal", 
    distance = "linear.logit", ratio = 1, replace = FALSE)

Sample sizes:
          Control Treated
All          2000     200
Matched       200     200
Unmatched    1800       0
Discarded       0       0
```

## Obtaining the Matched Sample

As before, much of the work can be done with `match.data`.

``` r
dm2200_matched7 <- match.data(match_7)
```

This leaves only the job of creating a matching number, for which we
have to develop some additional R code.

``` r
# Thanks to Robert McDonald at Mayo
# https://lists.gking.harvard.edu/pipermail/matchit/2012-April/000458.html

len <- dim(dm2200)[1]
matchx <- rep(NA,len)
len2 <- length(match_7$match.matrix)
count <- 1
for(i in 1:len2){
    
    match1 <- match_7$match.matrix[i]
    match2 <- row.names(match_7$match.matrix)[i]
    
    if(!is.na(match1)){
    matchx[as.numeric(match1)] <- count
    matchx[as.numeric(match2)] <- count
    count <- count+1}
    
}

dm2200_matched7 <- dm2200_matched7 %>%
    mutate(match7_matches = 
               matchx[as.numeric
                      (row.names(dm2200_matched7))]) %>%
    mutate(match7_matches = factor(match7_matches))
```

How many unique subjects are in our matched sample?

``` r
dm2200_matched7 %$% n_distinct(subject)
```

    [1] 400

This match includes 200 pairs so 400 subjects, as we’ve done matching
without replacement.

``` r
dm2200_matched7 %>% count(exposure)
```

    # A tibble: 2 x 2
      exposure     n
      <fct>    <int>
    1 A          200
    2 B          200

Do we have as many matching numbers as treated subjects?

``` r
dm2200_matched7 %$% n_distinct(match7_matches)
```

    [1] 200

## Checking Covariate Balance for our 1:1 Optimal Match

### Using `bal.tab` to obtain a balance table

``` r
bal7 <- bal.tab(match_7, un = TRUE, disp.v.ratio = TRUE)

bal7
```

    Call
     matchit(formula = f.build("treat", covs_1), data = dm2200, method = "optimal", 
        distance = "linear.logit", ratio = 1, replace = FALSE)
    
    Balance Measures
                         Type Diff.Un V.Ratio.Un Diff.Adj V.Ratio.Adj
    distance         Distance  2.0650     0.7466   0.2461      1.8593
    age               Contin. -0.4682     0.7579   0.0661      0.7978
    race_Asian         Binary -0.0025              0.0000            
    race_Black_AA      Binary  0.2975              0.0100            
    race_Other         Binary -0.0320              0.0050            
    race_White         Binary -0.2630             -0.0150            
    hisp               Binary -0.0320              0.0050            
    sex_M              Binary  0.2190              0.0350            
    insur_Commercial   Binary -0.2200              0.0050            
    insur_Medicaid     Binary  0.3745             -0.0250            
    insur_Medicare     Binary -0.2715              0.0000            
    insur_Uninsured    Binary  0.1170              0.0200            
    nincome           Contin. -0.9504     0.3117   0.1280      1.0981
    nhsgrad           Contin. -0.3904     1.0281   0.0096      1.0861
    cleve              Binary  0.3505             -0.0100            
    a1c               Contin.  0.0379     1.2209  -0.0281      0.9083
    ldl               Contin. -0.0743     1.1096   0.0287      1.1089
    visits            Contin.  0.4276     2.1732   0.1811      1.3180
    tobacco_Current    Binary  0.3175              0.0500            
    tobacco_Former     Binary -0.1575             -0.0450            
    tobacco_Never      Binary -0.1600             -0.0050            
    statin             Binary -0.0315              0.0450            
    ace_arb            Binary  0.0345              0.0450            
    betab              Binary -0.1255             -0.0500            
    depr_dx            Binary -0.1115             -0.0300            
    eyeex              Binary -0.0925             -0.0050            
    pneumo             Binary -0.3030             -0.0900            
    
    Sample sizes
              FALSE TRUE
    All        2000  200
    Matched     200  200
    Unmatched  1800    0

### Checking Rubin’s Rules 1 and 2

We’ll build a little table of the Rubin’s Rules (1 and 2) results before
and after `match_6` is applied, and compare these to the `match_5`
results.

``` r
rubin_report_m567 <- tibble(
    status = c("Rule1", "Rule2"),
    Unmatched = c(bal6$Balance$Diff.Un[1],
                  bal6$Balance$V.Ratio.Un[1]),
    Match5 = c(bal5$Balance$Diff.Adj[1],
               bal5$Balance$V.Ratio.Adj[1]),
    Match6 = c(bal6$Balance$Diff.Adj[1],
               bal6$Balance$V.Ratio.Adj[1]),
    Match7 = c(bal7$Balance$Diff.Adj[1],
               bal7$Balance$V.Ratio.Adj[1])
    )

rubin_report_m567 %>% knitr::kable(digits = 2)
```

| status | Unmatched | Match5 | Match6 | Match7 |
| :----- | --------: | -----: | -----: | -----: |
| Rule1  |      1.91 |   0.25 |   0.10 |   0.25 |
| Rule2  |      0.75 |   1.86 |   1.25 |   1.86 |

  - We note here that the optimal matching here does very little, if
    anything, to improved on the nearest neighbor match.

### Are the Optimal and Nearest Neighbor Matches the Same?

``` r
m5_subs <- dm2200_matched5 %>% select(subject)
m7_subs <- dm2200_matched7 %>% select(subject)

all_equal(m5_subs, m7_subs)
```

    [1] "Rows in x but not y: 167, 117, 352, 148. Rows in y but not x: 355, 346, 300, 137. "

Apparently not, but they are very similar.

### Using `bal.plot` from `cobalt`

Looking at the linear propensity scores in each group, in mirrored
histograms, we have …

``` r
bal.plot(obj = match_7,
         var.name = "distance", 
         which = "both",
         sample.names = 
             c("Unmatched Sample", "match_7 Sample"),
         type = "histogram", mirror = TRUE)
```

![](matching_with_dm2200_files/figure-gfm/unnamed-chunk-82-1.png)<!-- -->

### Using `love.plot` to look at Standardized Differences

``` r
love.plot(bal7, 
          threshold = .1, size = 3,
          var.order = "unadjusted",
          stats = "mean.diffs",
          stars = "raw",
          sample.names = c("Unmatched", "Matched"),
          title = "Love Plot for 1:1 Optimal Match") +
    labs(subtitle = "distance = linear propensity score",
         caption = "* indicates raw mean differences (for binary variables)")
```

![](matching_with_dm2200_files/figure-gfm/unnamed-chunk-83-1.png)<!-- -->

### Using `love.plot` to look at Variance Ratios

Again, the categorical variables are dropped.

``` r
love.plot(bal7, 
          threshold = .5, size = 3,
          stats = "variance.ratios",
          sample.names = c("Unmatched", "Matched"),
          title = "Variance Ratios for 1:1 Optimal Match") +
    labs(subtitle = "distance = linear propensity score",
         caption = "* indicates raw mean differences (for binary variables)")
```

![](matching_with_dm2200_files/figure-gfm/unnamed-chunk-84-1.png)<!-- -->

# `match_8` 1:2 Optimal Matching with `Matchit`

Again, we’ll match on the linear propensity score (by specifying the
`distance` to use the logistic link with the linear propensity score via
`"linear.logit"`), but now we’ll perform 1:2 (so `ratio = 2`) optimal
(so \`method = “optimal”) matching without replacement using the default
ordering (largest to smallest).

``` r
match_8 <- matchit(f.build("treat", covs_1), data = dm2200,
                   distance = "linear.logit", 
                   method = "optimal", ratio = 2, 
                   replace = FALSE)
```

    Warning in optmatch::fullmatch(d, min.controls = ratio, max.controls = ratio, : Without 'data' argument the order of the match is not guaranteed
        to be the same as your original data.

``` r
match_8
```

``` 

Call: 
matchit(formula = f.build("treat", covs_1), data = dm2200, method = "optimal", 
    distance = "linear.logit", ratio = 2, replace = FALSE)

Sample sizes:
          Control Treated
All          2000     200
Matched       400     200
Unmatched    1600       0
Discarded       0       0
```

## Obtaining the Matched Sample

As before, much of the work can be done with a very simple command.

``` r
dm2200_matched8 <- match.data(match_8)
```

This leaves only the job of creating a matching number, for which we
have to develop some additional R code, and this needs some modification
because we now have two matched controls for each treated subject.

``` r
# Thanks to Robert McDonald at Mayo
# https://lists.gking.harvard.edu/pipermail/matchit/2012-April/000458.html

mm1 <- match_8$match.matrix[,1]
mm2 <- match_8$match.matrix[,2]

len <- dim(dm2200)[1]
matchx <- rep(NA,len)
len2 <- length(match_8$match.matrix)
count <- 1
for(i in 1:len2){
    
    match_tr <- row.names(match_8$match.matrix)[i]
    match_con1 <- mm1[i]
    match_con2 <- mm2[i]
    
    if(!is.na(match_tr) & !is.na(match_con1) & !is.na(match_con2)){
    matchx[as.numeric(match_con1)] <- count
    matchx[as.numeric(match_con2)] <- count
    matchx[as.numeric(match_tr)] <- count
    count <- count+1}
    
}

dm2200_matched8 <- dm2200_matched8 %>%
     mutate(match8_matches = 
            matchx[as.numeric(row.names(dm2200_matched8))]) 
```

How many unique subjects are in our matched sample?

``` r
dm2200_matched8 %$% n_distinct(subject)
```

    [1] 600

This match includes 200 triplets (1 treated and 2 control) so 600
subjects, and again we’ve done matching without replacement.

``` r
dm2200_matched8 %>% count(exposure)
```

    # A tibble: 2 x 2
      exposure     n
      <fct>    <int>
    1 A          200
    2 B          400

Do we have as many matching numbers as treated subjects?

``` r
dm2200_matched8 %$% n_distinct(match8_matches)
```

    [1] 200

## Checking Covariate Balance for our 1:2 Optimal Match

### Using `bal.tab` to obtain a balance table

``` r
bal8 <- bal.tab(match_8, un = TRUE, disp.v.ratio = TRUE)

bal8
```

    Call
     matchit(formula = f.build("treat", covs_1), data = dm2200, method = "optimal", 
        distance = "linear.logit", ratio = 2, replace = FALSE)
    
    Balance Measures
                         Type Diff.Un V.Ratio.Un Diff.Adj V.Ratio.Adj
    distance         Distance  2.0650     0.7466   0.5002      2.4426
    age               Contin. -0.4682     0.7579  -0.0722      0.8303
    race_Asian         Binary -0.0025              0.0000            
    race_Black_AA      Binary  0.2975              0.0175            
    race_Other         Binary -0.0320              0.0050            
    race_White         Binary -0.2630             -0.0225            
    hisp               Binary -0.0320              0.0000            
    sex_M              Binary  0.2190              0.0750            
    insur_Commercial   Binary -0.2200              0.0000            
    insur_Medicaid     Binary  0.3745              0.0375            
    insur_Medicare     Binary -0.2715             -0.0850            
    insur_Uninsured    Binary  0.1170              0.0475            
    nincome           Contin. -0.9504     0.3117   0.0171      0.7121
    nhsgrad           Contin. -0.3904     1.0281  -0.0306      0.9398
    cleve              Binary  0.3505              0.0175            
    a1c               Contin.  0.0379     1.2209   0.0485      1.1913
    ldl               Contin. -0.0743     1.1096  -0.0269      1.0734
    visits            Contin.  0.4276     2.1732   0.2135      1.3368
    tobacco_Current    Binary  0.3175              0.1075            
    tobacco_Former     Binary -0.1575             -0.0575            
    tobacco_Never      Binary -0.1600             -0.0500            
    statin             Binary -0.0315             -0.0050            
    ace_arb            Binary  0.0345              0.0100            
    betab              Binary -0.1255             -0.0375            
    depr_dx            Binary -0.1115             -0.0625            
    eyeex              Binary -0.0925             -0.0225            
    pneumo             Binary -0.3030             -0.1475            
    
    Sample sizes
              FALSE TRUE
    All        2000  200
    Matched     400  200
    Unmatched  1600    0

### Checking Rubin’s Rules 1 and 2

We’ll build a little table of the Rubin’s Rules (1 and 2) results before
and after `match_8` is applied, and compare these to the other `MatchIt`
results we’ve developed.

``` r
rubin_report_m5678 <- tibble(
    status = c("Rule1", "Rule2"),
    Unmatched = c(bal8$Balance$Diff.Un[1],
                  bal8$Balance$V.Ratio.Un[1]),
    Match5 = c(bal5$Balance$Diff.Adj[1],
               bal5$Balance$V.Ratio.Adj[1]),
    Match6 = c(bal6$Balance$Diff.Adj[1],
               bal6$Balance$V.Ratio.Adj[1]),
    Match7 = c(bal7$Balance$Diff.Adj[1],
               bal7$Balance$V.Ratio.Adj[1]),
    Match8 = c(bal8$Balance$Diff.Adj[1],
               bal8$Balance$V.Ratio.Adj[1])
    )

rubin_report_m5678 %>% knitr::kable(digits = 2)
```

| status | Unmatched | Match5 | Match6 | Match7 | Match8 |
| :----- | --------: | -----: | -----: | -----: | -----: |
| Rule1  |      2.07 |   0.25 |   0.10 |   0.25 |   0.50 |
| Rule2  |      0.75 |   1.86 |   1.25 |   1.86 |   2.44 |

  - So the optimal matching doesn’t look very strong here. Of course,
    we’re matching 1:2 in this situation.

### Using `bal.plot` from `cobalt`

Looking at the linear propensity scores in each group, in mirrored
histograms, we have …

``` r
bal.plot(obj = match_8,
         var.name = "distance", 
         which = "both",
         sample.names = 
             c("Unmatched Sample", "match_8 Sample"),
         type = "histogram", mirror = TRUE)
```

![](matching_with_dm2200_files/figure-gfm/unnamed-chunk-93-1.png)<!-- -->

### Using `love.plot` to look at Standardized Differences

``` r
love.plot(bal8, 
          threshold = .1, size = 3,
          var.order = "unadjusted",
          stats = "mean.diffs",
          stars = "raw",
          sample.names = c("Unmatched", "Matched"),
          title = "Love Plot for 1:2 Optimal Match") +
    labs(subtitle = "distance = linear propensity score",
         caption = "* indicates raw mean differences (for binary variables)")
```

![](matching_with_dm2200_files/figure-gfm/unnamed-chunk-94-1.png)<!-- -->

### Using `love.plot` to look at Variance Ratios

Again, the categorical variables are dropped.

``` r
love.plot(bal8, 
          threshold = .5, size = 3,
          stats = "variance.ratios",
          sample.names = c("Unmatched", "Matched"),
          title = "Variance Ratios for 1:2 Optimal Match") +
    labs(subtitle = "distance = linear propensity score",
         caption = "* indicates raw mean differences (for binary variables)")
```

![](matching_with_dm2200_files/figure-gfm/unnamed-chunk-95-1.png)<!-- -->

# Planned matches coming as soon as Dr. Love finishes them

  - Optimal Matching using the `MatchIt` package
  - Full Matching using the `MatchIt` package
  - Genetic Matching using the `MatchIt` package
  - Coarsened Exact Matching using the `MatchIt` package

# Outcome Models

We’ll fit two (overly simplistic) outcome models, one for `bp_good` (our
binary outcome) and another for `bmi` (our quantitative outcome.) Later,
we’ll compare the `exposure` effect estimates made here to the estimates
we obtain after propensity matching. In each case, we’ll focus on ATT
estimates (average treated effect on the treated) rather than ATE
estimates.

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
| exposure == “A”TRUE |    0.636 |     0.216 |    0.417 |     0.972 |

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
| (Intercept)         |   35.691 |     0.619 |   34.478 |    36.904 |
| exposure == “A”TRUE |  \-2.864 |     0.867 |  \-4.563 |   \-1.165 |

## Adjusted Outcome Models after `match2`

### Binary Outcome: `bp_good`

``` r
result_match2_bp <- clogit(bp_good ~ (exposure == "A") + 
                          strata(match2_matches),
                      data = dm2200_matched2)

tidy(result_match2_bp, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) %>%
    select(term, estimate, std.error, conf.low, conf.high) %>%
    knitr::kable(digits = 3)
```

| term                | estimate | std.error | conf.low | conf.high |
| :------------------ | -------: | --------: | -------: | --------: |
| exposure == “A”TRUE |    0.451 |     0.174 |    0.321 |     0.634 |

### Quantitative Outcome: `bmi`

We’ll use a mixed model to account for our 1:1 matching. The matches
here are treated as a random factor, with the exposure a fixed factor,
in the `lme4` package.

``` r
dm2200_matched2 <- dm2200_matched2 %>% 
    mutate(match2_matches_f = as.factor(match2_matches))

result_match2_bmi <- lmer(bmi ~ (exposure == "A") + 
                              (1 | match2_matches_f), 
                          data = dm2200_matched2)

tidy(result_match2_bmi, 
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
| (Intercept)         |   34.719 |     0.476 |   33.787 |    35.651 |
| exposure == “A”TRUE |  \-1.892 |     0.552 |  \-2.973 |   \-0.811 |

## Adjusted Outcome Models after `match3`

### Binary Outcome: `bp_good`

``` r
result_match3_bp <- clogit(bp_good ~ (exposure == "A") + 
                          strata(match3_matches),
                      data = dm2200_matched3)

tidy(result_match3_bp, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) %>%
    select(term, estimate, std.error, conf.low, conf.high) %>%
    knitr::kable(digits = 3)
```

| term                | estimate | std.error | conf.low | conf.high |
| :------------------ | -------: | --------: | -------: | --------: |
| exposure == “A”TRUE |    0.609 |     0.137 |    0.466 |     0.797 |

### Quantitative Outcome: `bmi`

We’ll use a mixed model to account for our 1:1 matching. The matches
here are treated as a random factor, with the exposure a fixed factor,
in the `lme4` package.

``` r
dm2200_matched3 <- dm2200_matched3 %>% 
    mutate(match3_matches_f = as.factor(match3_matches))

result_match3_bmi <- lmer(bmi ~ (exposure == "A") + 
                              (1 | match3_matches_f), 
                          data = dm2200_matched3)

tidy(result_match3_bmi, 
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
| (Intercept)         |   34.637 |     0.416 |   33.821 |    35.452 |
| exposure == “A”TRUE |  \-1.810 |     0.458 |  \-2.708 |   \-0.912 |

## Adjusted Outcome Models after `match4`

### Binary Outcome: `bp_good`

``` r
result_match4_bp <- clogit(bp_good ~ (exposure == "A") + 
                          strata(match4_matches),
                      data = dm2200_matched4)

tidy(result_match4_bp, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) %>%
    select(term, estimate, std.error, conf.low, conf.high) %>%
    knitr::kable(digits = 3)
```

| term                | estimate | std.error | conf.low | conf.high |
| :------------------ | -------: | --------: | -------: | --------: |
| exposure == “A”TRUE |    0.659 |     0.239 |    0.412 |     1.053 |

### Quantitative Outcome: `bmi`

We’ll use a mixed model to account for our 1:1 matching. The matches
here are treated as a random factor, with the exposure a fixed factor,
in the `lme4` package.

``` r
dm2200_matched4 <- dm2200_matched4 %>% 
    mutate(match4_matches_f = as.factor(match4_matches))

result_match4_bmi <- lmer(bmi ~ (exposure == "A") + 
                              (1 | match4_matches_f), 
                          data = dm2200_matched4)

tidy(result_match4_bmi, 
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
| (Intercept)         |   35.799 |     0.700 |   34.427 |    37.171 |
| exposure == “A”TRUE |  \-2.548 |     0.953 |  \-4.416 |   \-0.679 |

## Adjusted Outcome Models after `match5`

### Binary Outcome: `bp_good`

``` r
result_match5_bp <- clogit(bp_good ~ (exposure == "A") + 
                          strata(match5_matches),
                      data = dm2200_matched5)

tidy(result_match5_bp, exponentiate = TRUE, 
     conf.int = TRUE, conf.level = 0.95) %>%
    select(term, estimate, std.error, conf.low, conf.high) %>%
    knitr::kable(digits = 3)
```

| term                | estimate | std.error | conf.low | conf.high |
| :------------------ | -------: | --------: | -------: | --------: |
| exposure == “A”TRUE |     0.63 |     0.219 |     0.41 |     0.967 |

### Quantitative Outcome: `bmi`

We’ll use a mixed model to account for our 1:1 matching. The matches
here are treated as a random factor, with the exposure a fixed factor,
in the `lme4` package.

``` r
result_match5_bmi <- lmer(bmi ~ (exposure == "A") + 
                              (1 | match5_matches), 
                          data = dm2200_matched5)
```

    boundary (singular) fit: see ?isSingular

``` r
tidy(result_match5_bmi, 
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
| (Intercept)         |   35.408 |     0.618 |   34.197 |    36.619 |
| exposure == “A”TRUE |  \-2.581 |     0.874 |  \-4.294 |   \-0.868 |

## Adjusted Outcome Models after `match6`

### Binary Outcome: `bp_good`

``` r
result_match6_bp <- clogit(bp_good ~ (exposure == "A") + 
                          strata(match6_matches),
                      data = dm2200_matched6)

tidy(result_match6_bp, exponentiate = TRUE, 
     conf.int = TRUE, conf.level = 0.95) %>%
    select(term, estimate, std.error, conf.low, conf.high) %>%
    knitr::kable(digits = 3)
```

| term                | estimate | std.error | conf.low | conf.high |
| :------------------ | -------: | --------: | -------: | --------: |
| exposure == “A”TRUE |    0.533 |     0.253 |    0.325 |     0.875 |

### Quantitative Outcome: `bmi`

We’ll use a mixed model to account for our 1:1 matching. The matches
here are treated as a random factor, with the exposure a fixed factor,
in the `lme4` package.

``` r
result_match6_bmi <- lmer(bmi ~ (exposure == "A") + 
                              (1 | match6_matches), 
                          data = dm2200_matched6)

tidy(result_match6_bmi, 
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
| (Intercept)         |   34.694 |     0.622 |   33.474 |    35.914 |
| exposure == “A”TRUE |  \-1.470 |     0.844 |  \-3.124 |     0.184 |

## Adjusted Outcome Models after `match7`

### Binary Outcome: `bp_good`

``` r
result_match7_bp <- clogit(bp_good ~ (exposure == "A") + 
                          strata(match7_matches),
                      data = dm2200_matched7)

tidy(result_match7_bp, exponentiate = TRUE, 
     conf.int = TRUE, conf.level = 0.95) %>%
    select(term, estimate, std.error, conf.low, conf.high) %>%
    knitr::kable(digits = 3)
```

| term                | estimate | std.error | conf.low | conf.high |
| :------------------ | -------: | --------: | -------: | --------: |
| exposure == “A”TRUE |    0.558 |     0.232 |    0.354 |     0.878 |

### Quantitative Outcome: `bmi`

We’ll use a mixed model to account for our 1:1 matching. The matches
here are treated as a random factor, with the exposure a fixed factor,
in the `lme4` package.

``` r
result_match7_bmi <- lmer(bmi ~ (exposure == "A") + 
                              (1 | match7_matches), 
                          data = dm2200_matched7)

tidy(result_match7_bmi, 
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
| (Intercept)         |   35.362 |     0.618 |   34.150 |    36.573 |
| exposure == “A”TRUE |  \-2.535 |     0.871 |  \-4.242 |   \-0.827 |

## Adjusted Outcome Models after `match8`

### Binary Outcome: `bp_good`

``` r
result_match8_bp <- clogit(bp_good ~ (exposure == "A") + 
                          strata(match8_matches),
                      data = dm2200_matched8)

tidy(result_match8_bp, exponentiate = TRUE, 
     conf.int = TRUE, conf.level = 0.95) %>%
    select(term, estimate, std.error, conf.low, conf.high) %>%
    knitr::kable(digits = 3)
```

| term                | estimate | std.error | conf.low | conf.high |
| :------------------ | -------: | --------: | -------: | --------: |
| exposure == “A”TRUE |     0.55 |     0.191 |    0.378 |       0.8 |

### Quantitative Outcome: `bmi`

We’ll use a mixed model to account for our 1:1 matching. The matches
here are treated as a random factor, with the exposure a fixed factor,
in the `lme4` package.

``` r
result_match8_bmi <- lmer(bmi ~ (exposure == "A") + 
                              (1 | match8_matches), 
                          data = dm2200_matched8)
```

    boundary (singular) fit: see ?isSingular

``` r
tidy(result_match8_bmi, 
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
| (Intercept)         |   34.738 |     0.433 |   33.888 |    35.587 |
| exposure == “A”TRUE |  \-1.911 |     0.750 |  \-3.381 |   \-0.440 |

# Cleanup

We’ve created a lot of variables here that we don’t actually need going
forward. So I’ll remove them here:

``` r
rm(bal1, bal2, bal3, bal4,
   covs_1, covs_1plus, covs_2plus, covs_3plus, covs_4plus,
   covs_for_rubin, dm2200, 
   dm2200_matched1, dm2200_matched2, dm2200_matched3,
   dm2200_matched4,
   match_1, match_2, match_3, match_4,
   prop_model, 
   result_match1_bmi, result_match2_bmi, result_match3_bmi,
   result_match4_bmi,
   result_match1_bp, result_match2_bp, result_match3_bp,
   result_match4_bp,
   rubin_m1, rubin_m2, rubin_m3, rubin_m4,
   rubin_report_m1, rubin_report_m12, rubin_report_m123,
   rubin_report_m1234, rubin_report_m5, 
   rubin_report_m56, rubin_report_m567, rubin_report_m5678, 
   t1, unadj_mod1, unadj_mod2,
   match1_matches, match2_matches, match3_matches,
   count, i, len, len2, match_con1, match_con2, match_tr,
   match1, match2, match4_matches, matchx, mm1, mm2, 
   bal5, bal6, bal7, bal8, dm2200_matched5, 
   dm2200_matched6, dm2200_matched7, dm2200_matched8,
   m5_subs, m7_subs, match_5, match_6, match_7, match_8,
   result_match5_bmi, result_match6_bmi, result_match7_bmi,
   result_match8_bmi,
   result_match5_bp, result_match6_bp, result_match7_bp,
   result_match8_bp)
```

# Key References

Matching in these examples was performed using the Matching package
(Sekhon, 2011), and covariate balance was assessed using cobalt
(Greifer, 2020), both in R (R Core Team, 2019).

  - Greifer, N. (2020). cobalt: Covariate Balance Tables and Plots. R
    package version 4.0.0.
  - Sekhon, J.S. (2011) Multivariate and Propensity Score Matching
    Software with Automated Balance Optimization: The Matching Package
    for R, *J of Statistical Software*, 2011, 42: 7,
    <http://www.jstatsoft.org/>. R package version 4.9-6.
  - R Core Team (2019). R: A language and environment for statistical
    computing. R Foundation for Statistical Computing, Vienna, Austria.
    URL <https://www.R-project.org/>.
