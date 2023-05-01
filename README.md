# stepseg: Stepwise Segmented Regression Analysis
Performs stepwise segmented regression analysis as described in Britt (2015). This analysis is designed to identify and evaluate the statistical significance of any number of continuous and discontinuous breakpoints. It is especially useful for analyzing data containing substantial error variance, such as volatile longitudinal observations of human behavior.

If you use this package, please cite the original publication:

Britt, B. C. (2015). Stepwise segmented regression analysis: An iterative statistical algorithm to detect and quantify evolutionary and revolutionary transformations in longitudinal data. In S. A. Matei, M. G. Russell, & E. Bertino (Eds.), _Transparency in social media: Tools, methods, and algorithms for mediating online interactions_ (pp. 125-144). Heidelberg, Germany: Springer.



## Installation

To use this package, first install and load it in R with

```r
install.packages("devtools")
library(devtools)
install_github("bcbritt/stepseg@master")
library("stepseg")
```



## Usage

To analyze univariate data (normal, Poisson, etc.), call the `stepseg` function. If your data are Dirichlet-distributed, then use `stepsegdir` instead. Otherwise, if your data are compositional, then `stepsegcomp` is the most appropriate choice.

### stepseg: Univariate Data

The `stepseg` function requires a vector or 1-column data.frame containing values for the dependent variable (`dv`). You may optionally specify a data.frame whose columns contain independent variables to be used in the analysis (`ivs`). If the `ivs` argument is omitted, `stepseg` will treat the observations in `dv` as being sequentially ordered and will use the row number as the sole independent variable, which is most appropriate for analyzing univariate time series data.

Several other arguments may be provided as well:

- `start_formula`: The formula to be used at the beginning of the stepwise model selection procedure. Variables inputted into the model using this argument will never be removed from the model. This argument should only include linear coefficients for independent variables, not indicator functions representing breakpoints or interactions among multiple independent variables. If interactions between multiple independent variables are desired, you should create another column in the `ivs` data.frame to represent that interaction. If this argument is not specified, `stepseg` will add an intercept and a linear coefficient for each independent variable to the regression model. This is sufficient for most use cases, so the default value of this argument should usually be left unchanged.

- `add` and `remove` indicate the _p_-value thresholds for adding variables to the regression model (forward selection) and removing variables from the regression model (backward selection), respectively. Both values must be between `0` and `1`, and `add` cannot be greater than `remove`.

- `add_mode`: When using stepwise segmented regression analysis, coefficients are added to and removed from the model in "blocks," with a given block consisting of the indicator function representing a given breakpoint location as well as all interaction terms between that indicator function and independent variables that are eligible to be added to the model. `add_mode` indicates the mode that should be used to identify the "best" block of coefficients to be considered for addition to the model. If `add_mode = "mse"` (the default), then the candidate block whose coefficients would jointly reduce the model mean squared error by the greatest amount is treated as the "best" block, and all of its coefficients are added to the model if the resulting _p_-value for at least one of those coefficients would be less than `add`. If `add_mode = "p"`, then all possible candidate blocks that could be added to the model are evaluated, and whichever block has a _p_-value less than or equal to all _p_-values in all other blocks is treated as the "best" block. As long as that _p_-value is less than the threshold specified by the \code{add} argument, all coefficients included in the block are added to the model. `add_mode = "mse"` is generally recommended, especially when `order > 0`, as `add_mode = "p"` is prone to false positives in the breakpoint detection procedure and can sometimes fall into an infinite loop, failing to ever converge.

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

### stepsegdir and stepsegcomp: Dirichlet and Compositional Data

The `stepsegdir` and `stepsegcomp` functions each requires a data.frame containing values for the multivariate dependent variable (`dv`). These functions may be used, for instance, to assess longitudinal changes in the relative representation of topics obtained via topic modeling. `stepsegdir` is most suitable for data that are known to be Dirichlet-distributed, such as results obtained via latent Dirichlet allocation. `stepsegcomp` is more versatile, as it does not require data to be Dirichlet-distributed, so it may be used for the output of other topic modeling procedures such as Dirichlet multinomial mixture models (e.g., LF-DMM, GS-DMM).

As with `stepseg`, you may optionally specify a data.frame whose columns contain independent variables to be used in the analysis (`ivs`). If the `ivs` argument is omitted, `stepseg` will treat the observations in `dv` as being sequentially ordered and will use the row number as the sole independent variable, which is most appropriate for analyzing univariate time series data.

Several other arguments may be provided as well:

- `start_breakpoints`: A vector of breakpoints to be added to the model before the stepwise procedure begins. This can be used to resume a stopped analysis (in which case `retain_start_breakpoints` should be set to `FALSE`) or to ensure that a set of breakpoints established _a priori_ persist in the final model (in which case `retain_start_breakpoints` should be set to `TRUE`).

- `retain_start_breakpoints`: If `TRUE`, then the breakpoints added to the model from `start_breakpoints` are never removed.

- `add` and `remove` indicate the _p_-value thresholds for adding variables to the regression model (forward selection) and removing variables from the regression model (backward selection), respectively. Both values must be between `0` and `1`, and `add` cannot be greater than `remove`.

- `bonferroni`: If `TRUE`, then the values of `add` and `remove` are adjusted using a Bonferroni correction based on the number of terms added to the that correspond to each breakpoint

- `add_mode`: The mode that should be used to identify the "best" breakpoint to be considered for addition to the model in each iteration of the stepwise procedure. If `add_mode = "mse"` (the default), then the candidate block whose coefficients would jointly reduce the model mean squared error (MSE) by the greatest amount is treated as the "best" block, and all of its coefficients are added to the model if the resulting _p_-value for at least one of those coefficients would be less than `add`. If `add_mode = "p"`, then all possible candidate blocks that could be added to the model are evaluated, and whichever block has a _p_-value less than or equal to all _p_-values in all other blocks is treated as the "best" block. As long as that _p_-value is less than the threshold specified by the \code{add} argument, all coefficients included in the block are added to the model. As with the `stepseg` procedure, `add_mode = "mse"` is generally recommended, especially when `order > 0`, as `add_mode = "p"` is prone to false positives in the breakpoint detection procedure and can sometimes fall into an infinite loop, failing to ever converge. However, in rare cases, the breakpoint that would be selected based on MSE will not correspond to a statistically significant change in any single category of a multivariate dependent variable. This can cause the stepwise procedure to prematurely terminate, in which case it may be worth considering `add_mode = "p"` as a fallback option.

- `ord`: The maximum exponent that will be applied to independent variables in regression coefficients. All interaction terms that are considered valid, with exponents ranging from `0` to `ord`, are included in a given block. `ord = 0` will only add indicator functions alone to the model (essentially representing changes in the intercept), while `ord = 1` will also add interactions between the indicator functions and one or more independent variables (representing changes in the slope at different levels of the specified independent variable). Setting `ord` to a value greater than `1` is generally not recommended, as higher-order terms may compete with lower-level terms in the regression model and lead to unpredictable or difficult-to-interpret results; see Britt (2015) for more information. Additionally, in `stepsegdir`, `ord` must be greater than `0`.

- `parameterization` (`stepsegdir` only): The parameterization used to perform Dirichlet regression. Permitted values are `"common"` and `"alternative"`.

- `type` (`stepsegcomp` only): The analysis type used for compositional regression. [**Compositional**::comp.reg](https://www.rdocumentation.org/packages/Compositional/versions/1.1/topics/comp.reg) does not return standard errors for `"lmfit"` or `"spatial"`, so at present, this must be set to `"classical"` to perform a classical multivariate regression.

- `allow_interactions`: Whether or not interactions may be created between indicator functions representing values of a given independent variable and other independent variables in the regression model. This is `FALSE` by default, as such interactions may be difficult to interpret and are likely to be spurious in many use cases.

- `iterlim`, `tol1`, and `tol2` (`stepsegdir` only) indicate the number of iterations that may be used for convergence when building each Dirichlet regression model as well as the convergence criteria for BFGS and NR optimization, respectively.

- `update` and `verbose` allow you to specify how many loop iterations should elapse between progress reports from the algorithm and how much information should be printed at the end of each stepwise model selection step. By default, `update = 0` and `verbose = FALSE`, which does not print any output. Setting `update` to a positive integer will print a one-line notation every `update` iterations. Setting `verbose = TRUE` will add the regression model, as it appears during the current iteration, to that progress report. However, the regression model at each iteration of the algorithm is also provided in the `model_iterations` output element, so there is rarely a need to set `update = TRUE`; this should generally only be done for debugging purposes, such as to diagnose the causes of unexpected results or to detect infinite loops if they emerge when `add_mode = "p"`.

Both `stepsegdir` and `stepsegcomp` output the final regression model obtained via the procedure.



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
#
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

set.seed(2023)
dv <- rbind(DirichletReg::rdirichlet(12, c(8,8,32)),
            DirichletReg::rdirichlet(18, c(40,24,4)))
model <- stepsegdir(dv, add_mode="mse", update=10, verbose=TRUE)
# "Beginning the stepwise procedure..."
# "Testing observation 10"
# "Testing observation 20"
# "End of step 1 ~ Current Breakpoints: 10"
# "Testing observation 10"
# "Testing observation 20"
# "No additional breakpoints can be added to the model."
# Warning message:
# In min(abs(summary(model)$coef.mat[indices, 4]), na.rm = T) :
#   no non-missing arguments to min; returning Inf

summary(model)
#
# Call:
# DirichReg(formula = as.formula(paste0("dv ~ ", paste(current_terms, collapse
# = " + "))), data = df, model = parameterization, control = list(iterlim =
# iterlim, tol1 = tol1, tol2 = tol2))
#
# Standardized Residuals:
#         Min       1Q   Median      3Q     Max
# V1  -2.1452  -0.6651  -0.2261  0.6674  2.1566
# V2  -1.7372  -1.1093   0.3604  0.8497  1.6465
# V3  -1.7572  -0.6547  -0.2742  0.6269  2.2728
#
# ------------------------------------------------------------------
# Beta-Coefficients for variable no. 1: V1
#                                    Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                         4.57845    0.68780   6.657  2.8e-11 ***
# I(ivs[, 1]^1)                      -0.22604    0.09559  -2.365   0.0180 *  
# I(ivs[, 1] > 12)TRUE               -0.54752    1.24725  -0.439   0.6607    
# I(ivs[, 1]^1):I(ivs[, 1] > 12)TRUE  0.21705    0.10660   2.036   0.0417 *  
# ------------------------------------------------------------------
# Beta-Coefficients for variable no. 2: V2
#                                    Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                         4.35351    0.69899   6.228 4.72e-10 ***
# I(ivs[, 1]^1)                      -0.18497    0.09834  -1.881    0.060 .  
# I(ivs[, 1] > 12)TRUE               -0.84702    1.24069  -0.683    0.495    
# I(ivs[, 1]^1):I(ivs[, 1] > 12)TRUE  0.17734    0.10875   1.631    0.103    
# ------------------------------------------------------------------
# Beta-Coefficients for variable no. 3: V3
#                                    Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                         5.97326    0.68569   8.711  < 2e-16 ***
# I(ivs[, 1]^1)                      -0.21326    0.09505  -2.244  0.02485 *  
# I(ivs[, 1] > 12)TRUE               -3.70315    1.21792  -3.041  0.00236 ** 
# I(ivs[, 1]^1):I(ivs[, 1] > 12)TRUE  0.18075    0.10538   1.715  0.08630 .  
# ------------------------------------------------------------------
# Significance codes: 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Log-likelihood: 118.8 on 12 df (83 BFGS + 2 NR Iterations)
# AIC: -213.5, BIC: -196.7
# Number of Observations: 30
# Link: Log
# Parametrization: common

dv <- rbind(c(0.058, 0.050, 0.010, 0.221, 0.186, 0.068, 0.407),
            c(0.054, 0.065, 0.009, 0.262, 0.196, 0.069, 0.345),
            c(0.040, 0.114, 0.004, 0.235, 0.191, 0.084, 0.332),
            c(0.044, 0.070, 0.011, 0.356, 0.172, 0.049, 0.298),
            c(0.048, 0.123, 0.007, 0.277, 0.175, 0.053, 0.317),
            c(0.046, 0.101, 0.013, 0.240, 0.212, 0.105, 0.283),
            c(0.058, 0.125, 0.010, 0.257, 0.191, 0.050, 0.309),
            c(0.048, 0.065, 0.010, 0.264, 0.247, 0.083, 0.283),
            c(0.072, 0.067, 0.006, 0.237, 0.197, 0.097, 0.324),
            c(0.040, 0.192, 0.002, 0.243, 0.174, 0.069, 0.280),
            c(0.041, 0.175, 0.005, 0.208, 0.174, 0.075, 0.322),
            c(0.046, 0.098, 0.004, 0.274, 0.206, 0.058, 0.314),
            c(0.037, 0.087, 0.004, 0.237, 0.222, 0.092, 0.321),
            c(0.035, 0.084, 0.004, 0.287, 0.196, 0.067, 0.327),
            c(0.040, 0.061, 0.004, 0.271, 0.226, 0.059, 0.339),
            c(0.038, 0.239, 0.006, 0.228, 0.148, 0.052, 0.289),
            c(0.051, 0.141, 0.004, 0.243, 0.180, 0.068, 0.313),
            c(0.035, 0.156, 0.007, 0.227, 0.199, 0.063, 0.313),
            c(0.037, 0.188, 0.008, 0.258, 0.197, 0.048, 0.264),
            c(0.037, 0.126, 0.002, 0.309, 0.180, 0.032, 0.314),
            c(0.035, 0.092, 0.005, 0.237, 0.242, 0.059, 0.330),
            c(0.043, 0.114, 0.005, 0.225, 0.245, 0.057, 0.311),
            c(0.038, 0.139, 0.002, 0.270, 0.184, 0.053, 0.314),
            c(0.049, 0.080, 0.010, 0.243, 0.215, 0.068, 0.335),
            c(0.041, 0.081, 0.005, 0.219, 0.193, 0.053, 0.408),
            c(0.045, 0.115, 0.001, 0.266, 0.181, 0.052, 0.340),
            c(0.047, 0.140, 0.003, 0.253, 0.170, 0.066, 0.321),
            c(0.049, 0.226, 0.003, 0.239, 0.170, 0.047, 0.266),
            c(0.041, 0.119, 0.002, 0.358, 0.153, 0.035, 0.292),
            c(0.030, 0.165, 0.010, 0.391, 0.090, 0.030, 0.284))
ivs <- c(1:30)
model <- stepsegcomp(dv, ivs, verbose=TRUE)
# "Beginning the stepwise procedure..."
# "End of step 1 ~ Current Breakpoints: 23"
# "End of step 2 ~ Current Breakpoints: 4, 23"
# "No additional breakpoints can be added to the model."

model
# $runtime
#    user  system elapsed 
#       0       0       0 
#
# $be
#                                  [,1]        [,2]       [,3]        [,4]
# (Intercept)                -0.2878712 -1.96059047  1.0820294  1.12821655
# xivs[,1]^1                  0.2700093  0.06038650  0.2450410  0.08682669
# xI(ivs[,1]>4)*(ivs[,1]^0)   0.7953717  0.19964276  0.4145045  0.09793102
# xI(ivs[,1]>4)*(ivs[,1]^1)  -0.2344747 -0.08772856 -0.2252737 -0.06469739
# xI(ivs[,1]>23)*(ivs[,1]^0) -4.2608171 -1.32045504 -3.1891707  1.80206300
# xI(ivs[,1]>23)*(ivs[,1]^1)  0.1431590  0.05065426  0.1120019 -0.08473563
#                                   [,5]         [,6]
# (Intercept)                 0.22781059  1.919198434
# xivs[,1]^1                  0.03425128  0.015529618
# xI(ivs[,1]>4)*(ivs[,1]^0)   0.23577716 -0.241304627
# xI(ivs[,1]>4)*(ivs[,1]^1)  -0.03957347  0.005926647
# xI(ivs[,1]>23)*(ivs[,1]^0)  1.58830152  0.353969914
# xI(ivs[,1]>23)*(ivs[,1]^1) -0.06607795 -0.022286050
#
# $seb
#                                  [,1]      [,2]       [,3]       [,4]
# (Intercept)                0.50524000 0.7699158 0.24702365 0.22129975
# xivs[,1]^1                 0.18448756 0.2811335 0.09020028 0.08080724
# xI(ivs[,1]>4)*(ivs[,1]^0)  0.56810356 0.8657112 0.27775912 0.24883456
# xI(ivs[,1]>4)*(ivs[,1]^1)  0.18529495 0.2823639 0.09059504 0.08116089
# xI(ivs[,1]>23)*(ivs[,1]^0) 2.12661681 3.2406697 1.03975270 0.93147764
# xI(ivs[,1]>23)*(ivs[,1]^1) 0.07985208 0.1216835 0.03904155 0.03497594
#                                  [,5]       [,6]
# (Intercept)                0.32977623 0.21663586
# xivs[,1]^1                 0.12041725 0.07910423
# xI(ivs[,1]>4)*(ivs[,1]^0)  0.37080803 0.24359039
# xI(ivs[,1]>4)*(ivs[,1]^1)  0.12094425 0.07945042
# xI(ivs[,1]>23)*(ivs[,1]^0) 1.38806839 0.91184680
# xI(ivs[,1]>23)*(ivs[,1]^1) 0.05212042 0.03423883
#
# $est
# NULL
# 
# $t
#                                  [,1]       [,2]      [,3]       [,4]
# (Intercept)                -0.5697712 -2.5464998  4.380266  5.0981376
# xivs[,1]^1                  1.4635635  0.2147965  2.716632  1.0744914
# xI(ivs[,1]>4)*(ivs[,1]^0)   1.4000470  0.2306113  1.492316  0.3935588
# xI(ivs[,1]>4)*(ivs[,1]^1)  -1.2654136 -0.3106933 -2.486601 -0.7971499
# xI(ivs[,1]>23)*(ivs[,1]^0) -2.0035660 -0.4074636 -3.067240  1.9346283
# xI(ivs[,1]>23)*(ivs[,1]^1)  1.7928026  0.4162787  2.868788 -2.4226831
#                                  [,5]        [,6]
# (Intercept)                 0.6908036  8.85909836
# xivs[,1]^1                  0.2844383  0.19631842
# xI(ivs[,1]>4)*(ivs[,1]^0)   0.6358470 -0.99061637
# xI(ivs[,1]>4)*(ivs[,1]^1)  -0.3272042  0.07459553
# xI(ivs[,1]>23)*(ivs[,1]^0)  1.1442531  0.38819012
# xI(ivs[,1]>23)*(ivs[,1]^1) -1.2677940 -0.65089993
#
# $p
#                                  [,1]       [,2]         [,3]         [,4]
# (Intercept)                0.57412808 0.01771475 0.0002008412 3.242277e-05
# xivs[,1]^1                 0.15628846 0.83174067 0.0120382387 2.932869e-01
# xI(ivs[,1]>4)*(ivs[,1]^0)  0.17429429 0.81957130 0.1486457861 6.973815e-01
# xI(ivs[,1]>4)*(ivs[,1]^1)  0.21787098 0.75871609 0.0202485998 4.331800e-01
# xI(ivs[,1]>23)*(ivs[,1]^0) 0.05653118 0.68727920 0.0052865773 6.490810e-02
# xI(ivs[,1]>23)*(ivs[,1]^1) 0.08561900 0.68090526 0.0084573677 2.332019e-02
#                                 [,5]         [,6]
# (Intercept)                0.4963166 4.953975e-09
# xivs[,1]^1                 0.7785143 8.460138e-01
# xI(ivs[,1]>4)*(ivs[,1]^0)  0.5308933 3.317559e-01
# xI(ivs[,1]>4)*(ivs[,1]^1)  0.7463503 9.411548e-01
# xI(ivs[,1]>23)*(ivs[,1]^0) 0.2638031 7.012974e-01
# xI(ivs[,1]>23)*(ivs[,1]^1) 0.2170338 5.212948e-01


```
