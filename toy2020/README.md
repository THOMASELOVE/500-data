# The Toy Example (2020 version)

The entire example can be downloaded from the front page of the 500-data website.

Here is the [R Markdown file](https://github.com/THOMASELOVE/500-data/blob/master/toy2020/toy2020_analysis.Rmd) including all code.

Here is the [Github Markdown file](https://github.com/THOMASELOVE/500-data/blob/master/toy2020/toy2020_analysis.md) which displays the code and its results.

I've also placed an [HTML version on RPubs](https://rpubs.com/TELOVE/toy2020-500) at https://rpubs.com/TELOVE/toy2020-500.

The `toy.csv` data file, in [raw, downloadable form](https://raw.githubusercontent.com/THOMASELOVE/500-data/master/toy2020/data/toy2020.csv).

- The Data Set is 100% fictional.
    - It contains data on 400 subjects (140 treated and 260 controls) on treatment status, six covariates, and three outcomes, with no     missing observations anywhere.
    - We assume that a logical argument suggests that the square of `covA`, as well as the interactions of `covB` with `covC` and with   `covD` should be related to treatment assignment, and thus should be included in our propensity model.
    - Our objective is to estimate the average causal effect of treatment (as compared to control) on each of the three outcomes, without   propensity adjustment, and then with propensity matching, subclassification, weighting and regression adjustment using the propensity score.
    
# The 13 Tasks We’ll Tackle in this Example

1.  Ignoring the covariate information, what is the unadjusted point estimate (and 95% confidence interval) for the effect of the     treatment on each of the three outcomes (`out1.cost`, `out2.event`, and `out3.time`)?
2.  Assume that theory suggests that the square of `covA`, as well as the interactions of `covB` with `covC` and `covB` with `covD` should be related to treatment assignment. Fit a propensity score model to the data, using the six covariates (A-F) and the three     transformations (A<sup>2</sup>, and the B-C and B-D interactions.) Plot the resulting propensity scores, by treatment group, in an     attractive and useful way. 
3.  Use Rubin’s Rules to assess the overlap of the propensity scores and the individual covariates prior to the use of any propensity score adjustments.
4.  Use 1:1 greedy matching to match all 140 treated subjects to control subjects without replacement on the basis of the linear propensity for treatment. Evaluate the degree of covariate imbalance before and after propensity matching for each of the six covariates, and present the pre- and post-match standardized differences and variance ratios for the covariates, as well as the square term and interactions, as well as both the raw and linear propensity score in appropriate plots. Now, build a new data frame containing the propensity-matched sample, and use it to first check Rubin’s Rules after matching.
5.  Now, use the matched sample data set to evaluate the treatment’s average causal effect on each of the three outcomes. In each case, specify a point estimate (and associated 95% confidence interval) for the effect of being treated (as compared to being a control     subject) on the outcome. Compare your results to the automatic versions reported by the Matching package when you include the     outcome in the matching process.
6.  Now, instead of matching, instead subclassify the subjects into quintiles by the raw propensity score. Display the balance in terms of standardized differences by quintile for the covariates, their transformations, and the propensity score in an appropriate table or plot(s). Are you satisfied?
7.  Regardless of your answer to the previous question, use the propensity score quintile subclassification approach to find a point estimate (and 95% confidence interval) for the effect of the treatment on each outcome.
8.  Now using a reasonable propensity score weighting strategy, assess the balance of each covariate, the transformations and the linear propensity score prior to and after propensity weighting. Is the balance after weighting satisfactory?
9.  Using propensity score weighting to evaluate the treatment’s effect, developing a point estimate and 95% CI for the average causal effect of treatment on each outcome.
10. Finally, use direct adjustment for the linear propensity score on the entire sample to evaluate the treatment’s effect, developing a point estimate and 95% CI for each outcome.
11. Now, try a double robust approach. Weight, then adjust for linear propensity score.
12. Compare your conclusions about the average causal effect obtained in the following six ways to each other. What happens and why? Which of these methods seems most appropriate given the available information?
    - without propensity adjustment,
    - after propensity matching,
    - after propensity score subclassification,
    - after propensity score weighting,
    - after adjusting for the propensity score directly, and
    - after weighting then adjusting for the PS, to each other.
13. Perform a sensitivity analysis for your matched samples analysis and the first outcome (`out1.cost`) if it turns out to show a     statistically detectable treatment effect.

