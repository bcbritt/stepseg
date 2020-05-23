# stepseg: Stepwise Segmented Regression Analysis
Performs stepwise segmented regression analysis as described in Britt (2015). This analysis is designed to identify and evaluate the statistical significance of any number of continuous and discontinuous breakpoints. It is especially useful for analyzing data containing substantial error variance, such as volatile longitudinal observations of human behavior.

If you use this package, please cite the original manuscript:

Britt, B. C. (2015). Stepwise segmented regression analysis: An iterative statistical algorithm to detect and quantify evolutionary and revolutionary transformations in longitudinal data. In S. A. Matei, M. G. Russell, & E. Bertino (Eds.), _Transparency in social media: Tools, methods, and algorithms for mediating online interactions_ (pp. 125-144). Heidelberg, Germany: Springer.



## Installation

To use this package, first install and load it in R with

```r
install.packages("devtools")
library(devtools)
install_github("bcbritt/stepseg")
library("stepseg")
```



## Usage

The `stepseg` function requires a vector or 1-column data.frame containing values for the dependent variable (`dv`). You may optionally specify a data.frame whose columns contain independent variables to be used in the analysis (`ivs`). If the `ivs` argument is omitted, `stepseg` will treat the observations in `dv` as being sequentially ordered and will use the row number as the sole independent variable, which is most appropriate for analyzing univariate time series data.

Several other arguments may be provided as well:

- `start_formula`: The formula to be used at the beginning of the stepwise model selection procedure. Variables inputted into the model using this argument will never be removed from the model. This argument should only include linear coefficients for independent variables, not indicator functions representing breakpoints or interactions among multiple independent variables. If interactions between multiple independent variables are desired, you should create another column in the `ivs` data.frame to represent that interaction. If this argument is not specified, `stepseg` will add an intercept and a linear coefficient for each independent variable to the regression model. This is sufficient for most use cases, so the default value of this argument should usually be left unchanged.

- `add` and `remove` indicate the _p_-value thresholds for adding variables to the regression model (forward selection) and removing variables from the regression model (backward selection), respectively. Both values must be between `0` and `1`, and `add` cannot be greater than `remove`.

- `add_mode`: When using stepwise segmented regression analysis, coefficients are added to and removed from the model in "blocks," with a given block consisting of the indicator function representing a given breakpoint location as well as all interaction terms between that indicator function and independent variables that are eligible to be added to the model. `add_mode` indicates the mode that should be used to identify the "best" block of coefficients to be considered for addition to the model. If `add_mode = "mse"` (the default), then the candidate block whose coefficients would jointly reduce the model mean squared error by the greatest amount is treated as the "best" block, and all of its coefficients are added to the model if the resulting `p`-value for at least one of those coefficients would be less than `add`. If `add_mode = "p"`, then all possible candidate blocks that could be added to the model are evaluated, and whichever block has a `p`-value less than or equal to all `p`-values in all other blocks is treated as the "best" block. As long as that `p`-value is less than the threshold specified by the \code{add} argument, all coefficients included in the block are added to the model. `add_mode = "mse"` is generally recommended, especially when `order > 0`, as `add_mode = "p"` is prone to false positives in the breakpoint detection procedure and can sometimes fall into an infinite loop, failing to ever converge.

- `family`: The distribution family to be used in the analysis. Any family specifications that are valid inputs to [**stats**::glm](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/glm) may be used here, including link functions.

- `order`: The maximum exponent that will be applied to independent variables in regression coefficients. All interaction terms that are considered valid, with exponents ranging from `0` to `order`, are included in a given block. `order = 0` will only add indicator functions alone to the model (essentially representing changes in the intercept), while `order = 1` will also add interactions between the indicator functions and one or more independent variables (representing changes in the slope at different levels of the specified independent variable). Setting `order` to a value greater than `1` is generally not recommended, as higher-order terms may compete with lower-level terms in the regression model and lead to unpredictable or difficult-to-interpret results; see Britt (2015) for more information.

- `allow_interactions`: Whether or not interactions may be created between indicator functions representing values of a given independent variable and other independent variables in the regression model. This is `FALSE` by default, as such interactions may be difficult to interpret and are likely to be spurious in many use cases.

- `update` and `verbose` allow you to specify how many iterations should elapse between progress reports from the algorithm and how much information should be printed in each report. By default, `update = 0` and `verbose = FALSE`, which does not print any output. Setting `update` to a positive integer will print a one-line notation every `update` iterations. Setting `verbose = TRUE` will add the regression model, as it appears during the current iteration, to that progress report. However, the regression model at each iteration of the algorithm is also provided in the `model_iterations` output element, so there is rarely a need to set `update = TRUE`; this should generally only be done for debugging purposes, such as to diagnose the causes of unexpected results or to detect infinite loops if they emerge when `add_mode = "p"`.

The output for the `stepseg` function contains the final regression model as well as many other elements:

- `output`: The summary of the final regression model that resulted from the stepwise segmented regression analysis, constructed using
[**stats**::glm](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/glm) and passed through [**stats**::drop1](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/add1) to obtain results using Type II sums of squares (since Type I sums of squares are especially prone to infinite loops when performing stepwise model selection)

- `raw_glm`: The raw output for the final regression model using [**stats**::glm](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/glm), without passing it through [**stats**::drop1](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/add1), thereby using Type I sums of squares rather than Type II

- `cleaned_output`: The summary of the final regression model that resulted from the stepwise segmented regression analysis, removing any coefficients that had \code{NA} \emph{p}-values due to singularities or other issues, constructed using [**stats**::glm](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/glm) and passed through [**stats**::drop1](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/add1) to obtain results using Type II sums of squares

- `cleaned_raw_glm`: The raw output for the final regression model using [**stats**::glm](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/glm), removing any coefficients that had `NA` _p_-values due to singularities or other issues, but without passing it through [**stats**::drop1](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/add1), thereby using Type I sums of squares rather than Type II

- `start_formula`: The formula for the regression model prior to the first iteration of the stepwise segmented regression analysis

- `final_formula`: The final formula for the regression model that resulted from the stepwise segmented regression analysis

- `cleaned_formula`: The final formula for the regression model that resulted from the stepwise segmented regression analysis, removing any coefficients that had `NA` _p_-values due to singularities or other issues

- `start_factors`: The coefficients used in the regression model prior to the first iteration of the stepwise segmented regression analysis

- `final_factors`: The coefficients used in the regression model that resulted from the stepwise segmented regression analysis

- `cleaned_factors`: The coefficients used in the regression model that resulted from the stepwise segmented regression analysis, removing any coefficients that had `NA` _p_-values due to singularities or other issues

- `add_iterations`: The number of forward selection iterations elapsed

- `remove_iterations`: The number of backward selection iterations elapsed

- `model_iterations`: A list showing how the regression model appeared at the end of each forward selection and backward selection iteration

- `dv`: The value set for the `dv` argument

- `ivs`: The value set for the `ivs` argument

- `add`: The value set for the `add` argument

- `remove`: The value set for the `remove` argument

- `add_mode`: The value set for the `add_mode` argument

- `family`: The value set for the `family` argument

- `order`: The value set for the `order` argument

- `allow_interactions`: The value set for the `allow_interactions` argument

- `update`: The value set for the `update` argument

- `verbose`: The value set for the `verbose` argument



## Examples

```r
set.seed(123)
dv <- c(11:30, seq(70, 12, by = -2)) + (rnorm(50)/10)
results1 <- stepseg(dv)
print(results1)
# Single term deletions
# 
# Model:
# dv ~ 1 + I(ivs[, 1]^1) + I(ivs[, 1] > 20) + I(ivs[, 1] > 20):I(ivs[, 
#     1]^1)
#                                Df Deviance    AIC F value    Pr(>F)    
# <none>                                 0.4 -89.11                      
# I(ivs[, 1]^1)                   1    663.3 279.15   75604 < 2.2e-16 ***
# I(ivs[, 1] > 20)                1  12818.4 427.22 1461966 < 2.2e-16 ***
# I(ivs[, 1]^1):I(ivs[, 1] > 20)  1   4607.1 376.06  525415 < 2.2e-16 ***
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

ivs <- cbind(c(1:10, 21:30, 1:30), c(1:40,50:41))
results2 <- stepseg(dv, ivs, order = 1, allow_interactions = TRUE, update = 2)
# "Beginning the stepwise procedure..."
# "Completed iteration #2"
# "Completed iteration #4"

print(results2)
# Single term deletions
# 
# Model:
# dv ~ 1 + I(ivs[, 1]^1) + I(ivs[, 2]^1) + I(ivs[, 2] > 20) + I(ivs[, 
#     2] > 20):I(ivs[, 1]^1) + I(ivs[, 2] > 20):I(ivs[, 2]^1) + 
#     I(ivs[, 1] > 25) + I(ivs[, 1] > 25):I(ivs[, 1]^1) + I(ivs[, 
#     1] > 25):I(ivs[, 2]^1)
#                                Df Deviance     AIC    F value    Pr(>F)    
# <none>                                0.31 -92.709                         
# I(ivs[, 1]^1)                   1     0.32 -92.585 1.7790e+00  0.189631    
# I(ivs[, 2]^1)                   1    28.33 131.492 3.7394e+03 < 2.2e-16 ***
# I(ivs[, 2] > 20)                1  1104.78 314.663 1.4737e+05 < 2.2e-16 ***
# I(ivs[, 1] > 25)                1     0.32 -92.518 1.8364e+00  0.182794    
# I(ivs[, 1]^1):I(ivs[, 2] > 20)  1   182.43 224.610 2.4301e+04 < 2.2e-16 ***
# I(ivs[, 2]^1):I(ivs[, 2] > 20)  1    21.79 118.370 2.8668e+03 < 2.2e-16 ***
# I(ivs[, 1]^1):I(ivs[, 1] > 25)  1     0.31 -93.556 9.5640e-01  0.333834    
# I(ivs[, 2]^1):I(ivs[, 1] > 25)  1     0.37 -85.378 8.4115e+00  0.005969 ** 
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

results3 <- stepseg(dv, ivs, order = 2, allow_interactions = FALSE,
  family = "poisson", update = 10) #Generates warnings due to non-Poisson DV
# "Beginning the stepwise procedure..."
# There were 50 or more warnings (use warnings() to see the first 50)

print(results3)
# Single term deletions

# Model:
# dv ~ 1 + I(ivs[, 1]^1) + I(ivs[, 1]^2) + I(ivs[, 2]^1) + I(ivs[, 
#     2]^2) + I(ivs[, 2] > 20) + I(ivs[, 2] > 20):I(ivs[, 2]^1) + 
#     I(ivs[, 2] > 20):I(ivs[, 2]^2)
#                                Df Deviance AIC     LRT  Pr(>Chi)    
# <none>                              0.8352 Inf                      
# I(ivs[, 1]^1)                   1   6.3978 Inf  5.5627 0.0183474 *  
# I(ivs[, 1]^2)                   1  13.0428 Inf 12.2076 0.0004759 ***
# I(ivs[, 2]^1)                   1   1.7744 Inf  0.9392 0.3324749    
# I(ivs[, 2]^2)                   1   5.8674 Inf  5.0322 0.0248799 *  
# I(ivs[, 2] > 20)                1  16.5465 Inf 15.7114 7.378e-05 ***
# I(ivs[, 2]^1):I(ivs[, 2] > 20)  1  12.8583 Inf 12.0231 0.0005255 ***
# I(ivs[, 2]^2):I(ivs[, 2] > 20)  1   3.7525 Inf  2.9174 0.0876314 .  
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```r
