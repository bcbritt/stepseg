#' Stepwise Segmented Regression Analysis
#'
#' @description This function performs stepwise segmented regression analysis,
#' as described by Britt (2015).
#'
#' @details When using stepwise segmented regression analysis, coefficients are
#' added to and removed from the model in "blocks," with a given block
#' consisting of the indicator function representing a given breakpoint
#' location as well as all interaction terms between that indicator function
#' and independent variables that are eligible to be added to the model.
#'
#' The \code{order} and \code{allow_interactions} arguments jointly indicate
#' what interaction terms are valid. \code{order} indicates the maximum
#' exponent that can be applied to the independent variables, while
#' \code{allow_interactions} indicates whether interaction terms can be created
#' between a given independent variable and an indicator function that is
#' defined by values of a different independent variable.
#'
#' All interaction terms that are considered valid, with exponents ranging from
#' \code{0} to \code{order}, are included in a given block. Consider, for
#' instance, an analysis in which \code{ivs}, the data.frame representing the
#' independent variables, has three columns signifying three variables. If
#' \code{order = 1} and \code{allow_interactions = FALSE}, which are their
#' default values in \code{stepseg}, then the block corresponding to a
#' potential breakpoint at \code{ivs[, 1] = 10} would include two coefficients:
#' \itemize{
#'   \item{\code{I(ivs[, 1] = 10)}}
#'   \item{\code{I(ivs[, 1] = 10):(ivs[, 1]^1)}}
#' }
#' Setting \code{order = 2} would retain both of the aforementioned
#' coefficients and also add
#' \itemize{
#'  \item{\code{I(ivs[, 1] = 10):(ivs[, 1]^2)}}
#' }
#' to the block. Maintaining \code{order = 2} and setting
#' \code{allow_interactions = TRUE} would result in the block containing a
#' total of seven coefficients:
#' \itemize{
#'   \item{\code{ivs[, 1] = 10}}
#'   \item{\code{I(ivs[, 1] = 10):(ivs[, 1]^1)}}
#'   \item{\code{I(ivs[, 1] = 10):(ivs[, 1]^2)}}
#'   \item{\code{I(ivs[, 1] = 10):(ivs[, 2]^1)}}
#'   \item{\code{I(ivs[, 1] = 10):(ivs[, 2]^2)}}
#'   \item{\code{I(ivs[, 1] = 10):(ivs[, 3]^1)}}
#'   \item{\code{I(ivs[, 1] = 10):(ivs[, 3]^2)}}
#' }
#' Note that if \code{order = 0}, then each block will only include the
#' indicator function itself (in this example, \code{I(ivs[, 1] = 10)})
#' regardless of the value of \code{allow_interactions}.
#' 
#' During each forward selection iteration, the algorithm selects the "best"
#' block of coefficients that can be added. If \code{add_mode = "p"}, then all
#' possible candidate blocks that could be added to the model are evaluated,
#' and whichever block has a \emph{p}-value less than or equal to all
#' \emph{p}-values in all other blocks is treated as the "best" block. As long
#' as that \emph{p}-value is less than the threshold specified by the
#' \code{add} argument, all coefficients included in the block are added to the
#' model. If \code{add_mode = "mse"}, then the candidate block whose
#' coefficients would jointly reduce the model mean squared error by the
#' greatest amount is treated as the `best` block, and all of its coefficients
#' are added to the model if the resulting \emph{p}-value for at least one of
#' those coefficients would be less than \code{add}.
#' 
#' During each backward selection iteration, all coefficients listed in
#' \code{start_formula} are retained in the regression model. Among all
#' remaining coefficients, any blocks whose \emph{p}-values are all greater
#' than or equal to \code{remove} are removed from the model. If at least one
#' \emph{p}-value contained in a given block is less than \code{remove}, then
#' none of the coefficients in the block are removed from the model.
#' 
#' When \code{order = 0}, \code{add_mode = "p"} is generally acceptable to
#' follow common conventions of stepwise model selection. When
#' \code{order > 0}, however, \code{add_mode = "p"} is more likely to result in
#' spurious breakpoints being added to the model due to intercept and
#' higher-order terms competing with one another, and in rare cases the
#' algorithm may entirely fail to converge. \code{add_mode = "mse"} is more
#' robust against these issues and is strongly recommended whenever
#' \code{order > 0}.
#'
#' Note that this function uses Type II sums of squares (via
#' \code{\link[stats]{drop1}}) for parameter estimation, as Type I sums of
#' squares (the default in \code{\link[stats]{glm}}) can result in an infinite
#' loop. If Type I sums of squares for the final model are desired, the
#' corresponding \code{glm} object can be retrieved via the \code{raw_glm}
#' element of the output from \code{stepseg}.
#'
#' The \code{start_formula} argument can be used to provide a regression model
#' whose coefficients will be inputted prior to the first stepwise model
#' selection iteration. However, only first-order independent variables should
#' be listed in this formula. Higher-order terms are not outputted in a
#' predictable manner by some of the external functions used called by
#' \code{stepseg} and may lead to invalid results. If an interaction term is
#' desired, you should create a new column in the data.frame submitted as the
#' \code{ivs} argument to represent that interaction, and that column can
#' subsequently be used in the analysis. Similarly, listing indicator functions
#' representing breakpoints in this argument may lead to unexpected behavior
#' during backward selection iterations of the \code{stepseg} function, so
#' forcibly adding them to the model via the \code{start_formula} argument is
#' not advised. For most use cases, this argument can and should be omitted.
#' 
#' The larger the value of \code{order}, the more likely spurious
#' breakpoints are to emerge. As with other regression contexts, you should
#' only increase the complexity of your model (such as by increasing the value
#' of \code{order}) when you have a clear reason to do so. Moreover, whenever
#' \code{order > 0}, singularities may occur that render some coefficients
#' inestimable. For instance, if \code{order = 3} and a pair of breakpoints is
#' identified along the same independent variable, with those two breakpoints
#' occurring two data points apart, there will be insufficient data between the
#' breakpoints to estimate the standard error of the interaction between the
#' indicator function and a cubic term, so it will be reported as \code{NA}.
#' This is normal, expected behavior in stepwise segmented regression analysis.
#' In such cases, those \code{NA} coefficients can be treated as though they
#' were absent from the model. The \code{cleaned_output},
#' \code{cleaned_formula}, and \code{cleaned_factors} elements of the output
#' from \code{stepseg} remove any \code{NA} coefficients accordingly, while
#' \code{cleaned_raw_glm} uses \code{cleaned_formula} in a
#' \code{\link[stats]{glm}} function call that, like \code{raw_glm}, uses Type
#' I sums of squares.
#' 
#' However, if \code{family = "binomial"} or \code{family = "poisson"}, any
#' singularities will generate an error. This is due to a known limitation in
#' \code{\link[stats]{drop1}}, which \code{stepseg} calls to generate
#' statistics using Type II sums of squares for general linear models. A
#' warning is generated when either of these two distribution families are used
#' as input arguments.
#' 
#' Finally, use caution when setting \code{order > 1}, as increasingly
#' high-order terms may compete with lower-order terms in the model and yield
#' unpredictable results. Refer to Britt (2015) for more information.
#' 
#' @param dv A vector or data.frame representing values of the dependent
#'   variable, which is coerced to a data.frame if not provided as one
#' @param ivs A data.frame representing values of the independent variables,
#'   which is coerced to a data.frame if not provided as one (default =
#'   \code{c(1:nrow(dv))}, which implicitly treats \code{dv} as sequentially
#'   ordered data and attempts to detect breakpoints in that sequence)
#' @param start_formula A \code{\link[stats]{formula}} to be used for the
#'   regression model prior to the first stepwise model selection iteration,
#'   with all listed coefficients permanently retained in the model (default =
#'   \code{NULL}, which adds an intercept and a linear coefficient for each
#'   independent variable to the regression model)
#' @param add A numeric atomic vector between \code{0} and \code{1} indicating
#'   the \emph{p}-value threshold to add coefficients to the model during each
#'   forward selection iteration, which must be less than or equal to the
#'   \code{remove} argument in order to avoid infinite loops (default =
#'   \code{.15}, as recommended by Britt, 2015)
#' @param remove A numeric atomic vector between \code{0} and \code{1}
#'   indicating the \emph{p}-value threshold to remove coefficients from the
#'   model during each backward selection iteration, which must be greater than
#'   or equal to the \code{add} argument in order to avoid infinite loops
#'   (default = \code{.20}, as recommended by Britt, 2015)
#' @param add_mode A character atomic vector (either \code{"p"} or \code{"mse"})
#'   indicating what criterion should be used to select the best candidate
#'   block during each forward selection iteration (default = \code{"mse"})
#' @param family A \code{\link[stats]{family}} object specifying the
#'   distribution and link function to be used in the general linear model
#'   (default = \code{stats::gaussian(link = "identity")}, which is appropriate
#'   for normally-distributed data and which facilitates a simple or multiple
#'   linear regression)
#' @param order A non-negative numeric atomic vector indicating the maximum
#'   exponent that will be applied to coefficients added to the regression
#'   model (default = \code{1}; it is generally recommended to use either
#'   \code{order = 0} or \code{order = 1})
#' @param allow_interactions A boolean value indicating whether interactions
#'   may be created between indicator functions representing breakpoints along
#'   one independent variable with other independent variables in the model
#'   (default = \code{FALSE}, which restricts each indicator function to only
#'   interact with the independent variable along which its breakpoint was
#'   defined)
#' @param update A numeric value indicating how many loop iterations should
#'   elapse between progress updates (default = \code{0}, which suppresses
#'   output)
#' @param verbose A boolean value indicating whether the current regression
#'   model and results should be outputted as part of the progress report every
#'   \code{update} iterations (default = \code{FALSE}; note that if
#'   \code{update = 0}, no updates will be printed even if
#'   \code{verbose = TRUE})
#' @return An object with the \code{stepseg_output} class that contains the
#'   following list elements:
#'
#'   \code{output}: The summary of the final regression model that resulted
#'     from the stepwise segmented regression analysis, constructed using
#'     \code{\link[stats]{glm}} and passed through \code{\link[stats]{drop1}})
#'     to obtain results using Type II sums of squares
#'
#'   \code{raw_glm}: The raw output for the final regression model using
#'     \code{\link[stats]{glm}}, without passing it through
#'     \code{\link[stats]{drop1}}, thereby using Type I sums of squares rather
#'     than Type II
#'
#'   \code{cleaned_output}: The summary of the final regression model that
#'     resulted from the stepwise segmented regression analysis, removing any
#'     coefficients that had \code{NA} \emph{p}-values due to singularities or
#'     other issues, constructed using \code{\link[stats]{glm}} and passed
#'     through \code{\link[stats]{drop1}}) to obtain results using Type II sums
#'     of squares
#'
#'   \code{cleaned_raw_glm}: The raw output for the final regression model
#'     using \code{\link[stats]{glm}}, removing any coefficients that had
#'     \code{NA} \emph{p}-values due to singularities or other issues, but
#'     without passing it through \code{\link[stats]{drop1}}, thereby using
#'     Type I sums of squares rather than Type II
#'
#'   \code{start_formula}: The formula for the regression model prior to the
#'     first iteration of the stepwise segmented regression analysis
#'
#'   \code{final_formula}: The final formula for the regression model that
#'     resulted from the stepwise segmented regression analysis
#'
#'   \code{cleaned_formula}: The final formula for the regression model that
#'     resulted from the stepwise segmented regression analysis, removing any
#'     coefficients that had \code{NA} \emph{p}-values due to singularities or
#'     other issues
#'
#'   \code{start_factors}: The coefficients used in the regression model prior
#'     to the first iteration of the stepwise segmented regression analysis
#'
#'   \code{final_factors}: The coefficients used in the regression model that
#'     resulted from the stepwise segmented regression analysis
#'
#'   \code{cleaned_factors}: The coefficients used in the regression model that
#'     resulted from the stepwise segmented regression analysis, removing any
#'     coefficients that had \code{NA} \emph{p}-values due to singularities or
#'     other issues
#'
#'   \code{add_iterations}: The number of forward selection iterations elapsed
#'
#'   \code{remove_iterations}: The number of backward selection iterations
#'     elapsed
#'
#'   \code{model_iterations}: A list showing how the regression model appeared
#'     at the end of each forward selection and backward selection iteration
#'
#'   \code{dv}: The value set for the \code{dv} argument
#'
#'   \code{ivs}: The value set for the \code{ivs} argument
#'
#'   \code{add}: The value set for the \code{add} argument
#'
#'   \code{remove}: The value set for the \code{remove} argument
#'
#'   \code{add_mode}: The value set for the \code{add_mode} argument
#'
#'   \code{family}: The value set for the \code{family} argument
#'
#'   \code{order}: The value set for the \code{order} argument
#'
#'   \code{allow_interactions}: The value set for the \code{allow_interactions}
#'     argument
#'
#'   \code{update}: The value set for the \code{update} argument
#'
#'   \code{verbose}: The value set for the \code{verbose} argument
#' @importFrom "stats" "glm"
#' @importFrom "stats" "drop1"
#' @importFrom "stats" "terms"
#' @importFrom "stats" "as.formula"
#' @importFrom "stats" "gaussian"
#' @section References:
#'   Britt, B. C. (2015). Stepwise segmented regression analysis: An iterative
#'   statistical algorithm to detect and quantify evolutionary and
#'   revolutionary transformations in longitudinal data. In S. A. Matei, M. G.
#'   Russell, & E. Bertino (Eds.), \emph{Transparency in social media: Tools,
#'   methods, and algorithms for mediating online interactions} (pp. 125-144).
#'   Heidelberg, Germany: Springer.
#' @examples
#' set.seed(123)
#' dv <- c(11:30, seq(70, 12, by = -2)) + (rnorm(50)/10)
#' results1 <- stepseg(dv)
#' print(results1)
#' # Single term deletions
#' # 
#' # Model:
#' # dv ~ 1 + I(ivs[, 1]^1) + I(ivs[, 1] > 20) + I(ivs[, 1] > 20):I(ivs[, 
#' #     1]^1)
#' #                                Df Deviance    AIC F value    Pr(>F)    
#' # <none>                                 0.4 -89.11                      
#' # I(ivs[, 1]^1)                   1    663.3 279.15   75604 < 2.2e-16 ***
#' # I(ivs[, 1] > 20)                1  12818.4 427.22 1461966 < 2.2e-16 ***
#' # I(ivs[, 1]^1):I(ivs[, 1] > 20)  1   4607.1 376.06  525415 < 2.2e-16 ***
#' # ---
#' # Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#'
#' ivs <- cbind(c(1:10, 21:30, 1:30), c(1:40,50:41))
#' results2 <- stepseg(dv, ivs, order = 1, allow_interactions = TRUE, update = 2)
#' # "Beginning the stepwise procedure..."
#' # "Completed iteration #2"
#' # "Completed iteration #4"
#'
#' print(results2)
#' # Single term deletions
#' # 
#' # Model:
#' # dv ~ 1 + I(ivs[, 1]^1) + I(ivs[, 2]^1) + I(ivs[, 2] > 20) + I(ivs[, 
#' #     2] > 20):I(ivs[, 1]^1) + I(ivs[, 2] > 20):I(ivs[, 2]^1) + 
#' #     I(ivs[, 1] > 25) + I(ivs[, 1] > 25):I(ivs[, 1]^1) + I(ivs[, 
#' #     1] > 25):I(ivs[, 2]^1)
#' #                                Df Deviance     AIC    F value    Pr(>F)    
#' # <none>                                0.31 -92.709                         
#' # I(ivs[, 1]^1)                   1     0.32 -92.585 1.7790e+00  0.189631    
#' # I(ivs[, 2]^1)                   1    28.33 131.492 3.7394e+03 < 2.2e-16 ***
#' # I(ivs[, 2] > 20)                1  1104.78 314.663 1.4737e+05 < 2.2e-16 ***
#' # I(ivs[, 1] > 25)                1     0.32 -92.518 1.8364e+00  0.182794    
#' # I(ivs[, 1]^1):I(ivs[, 2] > 20)  1   182.43 224.610 2.4301e+04 < 2.2e-16 ***
#' # I(ivs[, 2]^1):I(ivs[, 2] > 20)  1    21.79 118.370 2.8668e+03 < 2.2e-16 ***
#' # I(ivs[, 1]^1):I(ivs[, 1] > 25)  1     0.31 -93.556 9.5640e-01  0.333834    
#' # I(ivs[, 2]^1):I(ivs[, 1] > 25)  1     0.37 -85.378 8.4115e+00  0.005969 ** 
#' # ---
#' # Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#'
#' results3 <- stepseg(dv, ivs, order = 2, allow_interactions = FALSE,
#'   family = "poisson", update = 10) #Generates warnings due to non-Poisson DV
#' # "Beginning the stepwise procedure..."
#' # There were 50 or more warnings (use warnings() to see the first 50)
#'
#' print(results3)
#' # Single term deletions
#' # 
#' # Model:
#' # dv ~ 1 + I(ivs[, 1]^1) + I(ivs[, 1]^2) + I(ivs[, 2]^1) + I(ivs[, 
#' #     2]^2) + I(ivs[, 2] > 20) + I(ivs[, 2] > 20):I(ivs[, 2]^1) + 
#' #     I(ivs[, 2] > 20):I(ivs[, 2]^2)
#' #                                Df Deviance AIC     LRT  Pr(>Chi)    
#' # <none>                              0.8352 Inf                      
#' # I(ivs[, 1]^1)                   1   6.3978 Inf  5.5627 0.0183474 *  
#' # I(ivs[, 1]^2)                   1  13.0428 Inf 12.2076 0.0004759 ***
#' # I(ivs[, 2]^1)                   1   1.7744 Inf  0.9392 0.3324749    
#' # I(ivs[, 2]^2)                   1   5.8674 Inf  5.0322 0.0248799 *  
#' # I(ivs[, 2] > 20)                1  16.5465 Inf 15.7114 7.378e-05 ***
#' # I(ivs[, 2]^1):I(ivs[, 2] > 20)  1  12.8583 Inf 12.0231 0.0005255 ***
#' # I(ivs[, 2]^2):I(ivs[, 2] > 20)  1   3.7525 Inf  2.9174 0.0876314 .  
#' # ---
#' # Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#' @export
#' stepseg

stepseg <- function(dv, ivs = data.frame(1:nrow(data.frame(dv))), start_formula = NULL, add = .15, remove = .20, add_mode = "mse", family = stats::gaussian(link = "identity"), order = 1, allow_interactions = FALSE, update = 0, verbose = FALSE) {
   #Initial Setup
   if(add > remove) {
      stop("'add' cannot be set to a value greater than 'remove,' as doing so may result in an infinite loop.")
   }
   if(order > 1) {
      warning("Stepwise segmented regression may produce unpredictable results when order > 1.")
   }
   if((order > 0) & (add_mode == "p")) {
      warning("When order > 0, setting add_mode = 'p' may yield spurious results. In some cases, the algorithm may entirely fail to converge. Consider setting add_mode = 'mse' instead.")
   }
   if(!is.atomic(family)) { #If family is an atomic vector, that means it only indicates the family for the glm and nothing else; if it is not, then (if properly constructed) the specific family element needs to be retrieved
      family_type = family$family
   } else {
      family_type = family
   }
   if(family_type %in% c("binomial", "poisson")) { #For models with known deviance, a Chi-square test is most appropriate
      test = "Chi"
      warning("When 'family_type' is set to 'binomial' or 'poisson,' this function generates an error if a singularity occurs at any stage of the analysis. This is due to a known limitation in the stats::drop1() function used to obtain p-values with Type II sums of squares (rather than Type I, which is the default in the stats::glm() function).")
   } else if(family_type %in% c("gaussian", "Gamma", "inverse.gaussian", "quasi", "quasibinomial", "quasipoisson")) { #For models with estimated deviance, an F-test is most appropriate
      test = "F"
   } else if(family_type %in% c("Dirichlet", "dirichlet")) {
      stop("To analyze Dirichlet-distributed data, use stepsegdir().")
   } else if(family_type %in% c("compositional")) {
      stop("To analyze compositional data, use stepsegcomp().")
   } else {
      stop("The distribution family must be one of the following: 'gaussian', 'binomial', 'Gamma', 'inverse.gaussian', 'poisson', 'quasi', 'quasibinomial', 'quasipoisson'")
   }
   utils::flush.console()

   add_iterations <- 0
   remove_iterations <- 0
   ivs <- data.frame(ivs)
   if(is.null(start_formula)) {
      factors = c("1")
      for(i in 1:ncol(ivs)) {  #If no formula was provided, then by default, set the starting formula to include an intercept and all IVs (with no indicator functions)
         for(j in 1:order) {
            factors <- append(factors, paste("I(ivs[,", i, "]^", j, ")", sep=""))
         }
      }
      start_formula <- stats::as.formula(paste("dv ~ ", paste(factors, collapse="+")))
   } else { #If a formula was provided, then use that set of IVs as the complete set of factors for the starting function
      factors <- all.vars(formula)[2:length(all.vars(formula))]
   }
   start_factors <- factors
   formula <- stats::as.formula(paste("dv ~ ", paste(factors, collapse="+")))
   permanent_factors <- length(factors) #The number of factors that were manually inputted into the model from the start; these should never be removed from the model
   if(allow_interactions) { #This block of code sets block_size, which indicates how many IVs are added to the model together (corresponding to a single breakpoint) and must be removed together
      block_size <- 1 + (order*ncol(ivs))
   } else {
      block_size <- 1 + order
   }

   added <- TRUE #Indicates whether any variables were added to the model in the last addition step
   removed <- TRUE #Indicates whether any variables were removed from the model in the last removal step
   model_iterations <- list()
   model_iterations[[1]] <- stats::drop1(stats::glm(stats::terms(formula, specials = "s", keep.order = TRUE), family = family), scope = formula, test=test) #Records the model at each step, starting with how the model appeared at the beginning of the process

   if(update > 0) {
      print("Beginning the stepwise procedure...")
      utils::flush.console()
   }

   repeat {
      #Forward selection
      add_iterations <- add_iterations + 1
      results <- data.frame(rep(NA,nrow(ivs)))
      if(ncol(ivs) > 1) {
         for(i in 2:ncol(ivs)) {
            results <- cbind(results, rep(NA,nrow(ivs)))
         }
      }

      for(i in 1:ncol(ivs)) {
         for(j in 1:nrow(ivs)) {
            temp_factors <- factors
            temp_factors <- append(temp_factors, paste("I(ivs[,", i, "] > ivs[", j, ",", i, "])", sep=""))
            if(order > 0) {
               if(allow_interactions) {  #If allow_interactions is TRUE, then we test how the model would look if we added all possible terms, including (linear, quadratic, etc.) interactions between this indicator function and unrelated IVs
                  for(m in 1:ncol(ivs)) {
                     for(k in 1:order) {
                        temp_factors <- append(temp_factors, paste("I(ivs[,", i, "] > ivs[", j, ",", i, "]):I(ivs[,", m, "]^", k, ")", sep=""))
                     }
                  }
               } else {
                  for(k in 1:order) {
                     temp_factors <- append(temp_factors, paste("I(ivs[,", i, "] > ivs[", j, ",", i, "]):I(ivs[,", i, "]^", k, ")", sep=""))
                  }
               }
            }
            temp_formula <- stats::as.formula(paste("dv ~ ", paste(temp_factors, collapse="+")))
            if(add_mode == "p") { #Record the lowest p-value that this block of coefficients would yield
               temp_pvalues <- rep(NA,length(temp_factors))
               #NOTE: When test = "chi" (which is true when family = "binomial" or "poisson"), the below line trigger an error if any coefficients are not defined due to singularities. There is not a good workaround in R, as all other functions that I have found to perform GLM either do not provide p-values or use Type I sums of squares.
               suppressWarnings(temp_pvalues[1:permanent_factors] <- stats::drop1(stats::glm(stats::terms(temp_formula, specials = "s", keep.order = TRUE), family = family), scope = temp_formula, test=test)[,ncol(stats::drop1(stats::glm(stats::terms(temp_formula, specials = "s", keep.order = TRUE), family = family), scope = temp_formula, test=test))][1:permanent_factors]) #The stats::drop1() function provides Type II sums of squares
               k <- permanent_factors + 1
               m <- k + ((length(temp_pvalues) - permanent_factors) / block_size)
               n <- k
               temp_results <- stats::drop1(stats::glm(stats::terms(temp_formula, specials = "s", keep.order = TRUE), family = family), scope = temp_formula, test=test)
               while(n <= length(temp_pvalues)) {
                  suppressWarnings(temp_pvalues[n] <- temp_results[k,ncol(temp_results)])
                  suppressWarnings(temp_pvalues[(n+1):(n+block_size-1)] <- temp_results[m:(m+block_size-2),ncol(temp_results)])
                  k <- k + 1
                  m <- m + block_size - 1
                  n <- n + block_size
               }
               suppressWarnings(results[j,i] <- min(temp_pvalues[(length(temp_pvalues)-block_size+1):length(temp_pvalues)], na.rm=TRUE))
            } else {
               if(add_mode == "mse") { #Record the MSE that would result from adding this block of coefficients
                  old_mse <- mean(stats::glm(stats::terms(formula, specials = "s", keep.order = TRUE), family = family)$residuals^2)
                  new_mse <- mean(stats::glm(stats::terms(temp_formula, specials = "s", keep.order = TRUE), family = family)$residuals^2)
                  if(new_mse < old_mse) { #We only want to consider adding this block of coefficients if it would reduce the MSE, not increase MSE or leave it unchanged
                     results[j,i] <- mean(stats::glm(stats::terms(temp_formula, specials = "s", keep.order = TRUE), family = family)$residuals^2)
                  } else {
                     results[j,i] <- Inf
                  }
               }
            }
         }
      }

      candidate_row <- which.min(apply(results,MARGIN=1,min)) #Finds the cell of results that has the lowest value and returns its row number...
      candidate_column <- which.min(apply(results,MARGIN=2,min)) #...and column number

      temp_factors <- factors
      temp_factors <- append(temp_factors, paste("I(ivs[,", candidate_column, "] > ivs[", candidate_row, ",", candidate_column, "])", sep=""))
      if(order > 0) {
         if(allow_interactions) { #If allow_interactions is TRUE, then we test how the model would look if we added all possible terms, including (linear, quadratic, etc.) interactions between this indicator function and unrelated IVs
            for(m in 1:ncol(ivs)) {
               for(k in 1:order) {
                  temp_factors <- append(temp_factors, paste("I(ivs[,", candidate_column, "] > ivs[", candidate_row, ",", candidate_column, "]):I(ivs[,", m, "]^", k, ")", sep=""))
               }
            }
         } else {
            for(k in 1:order) {
               temp_factors <- append(temp_factors, paste("I(ivs[,", candidate_column, "] > ivs[", candidate_row, ",", candidate_column, "]):I(ivs[,", candidate_column, "]^", k, ")", sep=""))
            }
         }
      }
      temp_formula <- stats::as.formula(paste("dv ~ ", paste(temp_factors, collapse="+")))
      pvalues <- rep(NA,length(temp_factors))
      suppressWarnings(pvalues[1:permanent_factors] <- stats::drop1(stats::glm(stats::terms(temp_formula, specials = "s", keep.order = TRUE), family = family), scope = temp_formula, test=test)[,ncol(stats::drop1(stats::glm(stats::terms(temp_formula, specials = "s", keep.order = TRUE), family = family), scope = temp_formula, test=test))][1:permanent_factors]) #The stats::drop1() function provides Type II sums of squares
      k <- permanent_factors + 1
      m <- k + ((length(pvalues) - permanent_factors) / block_size)
      n <- k
      temp_results <- stats::drop1(stats::glm(stats::terms(temp_formula, specials = "s", keep.order = TRUE), family = family), scope = temp_formula, test=test)
      while(n <= length(pvalues)) {
         suppressWarnings(pvalues[n] <- temp_results[k,ncol(temp_results)])
         suppressWarnings(pvalues[(n+1):(n+block_size-1)] <- temp_results[m:(m+block_size-2),ncol(temp_results)])
         k <- k + 1
         m <- m + block_size - 1
         n <- n + block_size
      }
      if(sum(pvalues[(length(pvalues)-block_size+1):length(pvalues)] < remove, na.rm=TRUE) > 0) { #If at least one p-value in this block of coefficients to be added to the model would have a p-value less than the "remove" threshold, then we add the block of coefficients to the model
         apply1 <- which.min(apply(results,MARGIN=1,min)) #Row number with the lowest p-value
         apply2 <- which.min(apply(results,MARGIN=2,min)) #Column number with the lowest p-value
         factors <- append(factors, paste("I(ivs[,", apply2, "] > ", ivs[apply1,apply2], ")", sep=""))
         if(allow_interactions) { #If allow_interactions is TRUE, then we add all possible terms, including (linear, quadratic, etc.) interactions between this indicator function and unrelated IVs
            if(order > 0) {
               for(j in 1:ncol(ivs)) {
                  for(k in 1:order) {
                     factors <- append(factors, paste("I(ivs[,", apply2, "] > ", ivs[apply1,apply2], "):I(ivs[,", j, "]^", k, ")", sep=""))
                  }
               }
            }
         } else { #If allow_interactions is FALSE, then we add this indicator function and its (linear, quadratic, etc.) interactions with the relevant IV alone
            if(order > 0) {
               for(k in 1:order) {
                  factors <- append(factors, paste("I(ivs[,", apply2, "] > ", ivs[apply1,apply2], "):I(ivs[,", apply2, "]^", k, ")", sep=""))
               }
            }
         }
         formula <- stats::as.formula(paste("dv ~ ", paste(factors, collapse="+")))
         added <- TRUE
      } else { #Otherwise, note that no variables were added during this step
         added <- FALSE
      }

      model_iterations[[1+add_iterations+remove_iterations]] <- stats::drop1(stats::glm(stats::terms(formula, specials = "s", keep.order = TRUE), family = family), scope = formula, test=test)

      if((update > 0) & ((add_iterations + remove_iterations) %% update == 0)) {
         print(paste("Completed iteration #", add_iterations + remove_iterations, sep=""))
         if(verbose) {
            print(model_iterations[[1+add_iterations+remove_iterations]])
         }
         utils::flush.console()
      }

      if(!added & !removed) { #If no variables can be added or removed from the model, terminate the function
         break
      }

      #Backward selection
      remove_iterations <- remove_iterations + 1

      keep_factors <- rep(NA,length(factors))
      pvalues <- rep(NA,length(factors))
      suppressWarnings(pvalues[1:permanent_factors] <- stats::drop1(stats::glm(stats::terms(formula, specials = "s", keep.order = TRUE), family = family), scope = formula, test=test)[,ncol(stats::drop1(stats::glm(stats::terms(formula, specials = "s", keep.order = TRUE), family = family), scope = formula, test=test))][1:permanent_factors]) #The stats::drop1() function provides Type II sums of squares
      k <- permanent_factors + 1
      m <- k + ((length(pvalues) - permanent_factors) / block_size)
      n <- k

      temp_results <- stats::drop1(stats::glm(stats::terms(formula, specials = "s", keep.order = TRUE), family = family), scope = formula, test=test)
      while(n <= length(pvalues)) {
         suppressWarnings(pvalues[n] <- temp_results[k,ncol(temp_results)])
         suppressWarnings(pvalues[(n+1):(n+block_size-1)] <- temp_results[m:(m+block_size-2),ncol(temp_results)])
         k <- k + 1
         m <- m + block_size - 1
         n <- n + block_size
      }

      keep_factors[1:permanent_factors] <- TRUE
      j <- permanent_factors + 1
      while(j < length(pvalues)) {
         if(sum(pvalues[j:(j+block_size-1)] < remove, na.rm=TRUE) > 0) { #If at least one p-value in this block is less than the "remove" threshold, then none of the factors in this block can be removed
            keep_factors[j:(j+block_size-1)] <- TRUE
         } else { #Otherwise, all of the factors in this block should be removed
            keep_factors[j:(j+block_size-1)] <- FALSE
         }
         j <- j + block_size
      }

      if(length(pvalues) > sum(keep_factors)) { #If there are any factors that can be removed, then do so
         factors <- factors[keep_factors]
         formula <- stats::as.formula(paste("dv ~ ", paste(factors, collapse="+")))
         removed <- TRUE
      } else { #Otherwise, note that no factors were removed during this step
         removed <- FALSE
      }

      if(!added & !removed) { #If no variables can be added or removed from the model, terminate the function
         break
      }

      model_iterations[[1+add_iterations+remove_iterations]] <- stats::drop1(stats::glm(stats::terms(formula, specials = "s", keep.order = TRUE), family = family), scope = formula, test=test)
      if((update > 0) & ((add_iterations + remove_iterations) %% update == 0)) {
         print(paste("Completed iteration #", add_iterations + remove_iterations, sep=""))
         if(verbose) {
            print(model_iterations[[1+add_iterations+remove_iterations]])
         }
         utils::flush.console()
      }
   }
   output <- stats::drop1(stats::glm(stats::terms(formula, specials = "s", keep.order = TRUE), family = family), scope = formula, test=test)
   pvalues <- rep(NA,length(factors))
   suppressWarnings(pvalues[1:permanent_factors] <- output[1:permanent_factors,ncol(output)])
   k <- permanent_factors + 1
   m <- k + ((length(pvalues) - permanent_factors) / block_size)
   n <- k
   while(n <= length(pvalues)) {
      suppressWarnings(pvalues[n] <- output[,ncol(output)][k])
      suppressWarnings(pvalues[(n+1):(n+block_size-1)] <- output[,ncol(output)][m:(m+block_size-2)])
      k <- k + 1
      m <- m + block_size - 1
      n <- n + block_size
   }
   keep_factors <- rep(NA,length(factors))
   keep_factors[1:permanent_factors] <- TRUE
   keep_factors[(permanent_factors+1):length(keep_factors)] <- !is.na(pvalues[(permanent_factors+1):length(pvalues)])
   cleaned_factors <- factors[which(keep_factors)]
   cleaned_formula <- stats::as.formula(paste("dv ~ ", paste(cleaned_factors, collapse="+")))
   suppressWarnings(output <- list(output = output,
                                   raw_glm = stats::glm(stats::terms(formula, specials = "s", keep.order = TRUE), family = family),
                                   cleaned_output = stats::drop1(stats::glm(stats::terms(cleaned_formula, specials = "s", keep.order = TRUE), family = family), scope = cleaned_formula, test=test),
                                   cleaned_raw_glm = stats::glm(stats::terms(cleaned_formula, specials = "s", keep.order = TRUE), family = family),
                                   start_formula = start_formula,
                                   final_formula = formula,
                                   cleaned_formula = cleaned_formula,
                                   start_factors = start_factors,
                                   final_factors = factors,
                                   cleaned_factors = cleaned_factors,
                                   add_iterations = add_iterations,
                                   remove_iterations = remove_iterations,
                                   model_iterations = model_iterations,
                                   dv = dv,
                                   ivs = ivs,
                                   add = add,
                                   remove = remove,
                                   add_mode = add_mode,
                                   family = family,
                                   order = order,
                                   allow_interactions = allow_interactions,
                                   update = update,
                                   verbose = verbose
                   ))
   class(output) <- c("stepseg_output", class(output))
   return(output) #A list containing a glm object, the number of iterations, the model at each iteration, the final formula, and the inputs for the function call
}

#' Stepwise Segmented Compositional Regression Analysis
#'
#' @description This function performs stepwise segmented compositional
#' regression analysis, a modified form of the procedure described by Britt
#' (2015).
#'
#' @details When using any form of stepwise segmented regression analysis,
#' coefficients are added to and removed from the model in "blocks," with a
#' given block consisting of the indicator function representing a given
#' breakpoint location as well as all interaction terms between that indicator
#' function and independent variables that are eligible to be added to the
#' model.
#'
#' The \code{ord} and \code{allow_interactions} arguments jointly indicate what
#' interaction terms are valid. \code{ord} indicates the maximum exponent that
#' can be applied to the independent variables, while \code{allow_interactions}
#' indicates whether interaction terms can be created between a given
#' independent variable and an indicator function that is defined by values of a
#' different independent variable.
#'
#' All interaction terms that are considered valid, with exponents ranging from
#' \code{0} to \code{ord}, are included in a given block. Consider, for
#' instance, an analysis in which \code{ivs}, the data.frame representing the
#' independent variables, has three columns signifying three variables. If
#' \code{ord = 1} and \code{allow_interactions = FALSE}, which are their default
#' values in \code{stepseg}, then the block corresponding to a potential
#' breakpoint at \code{ivs[, 1] = 10} would include two coefficients:
#' \itemize{
#'   \item{\code{I(ivs[, 1] = 10)}}
#'   \item{\code{I(ivs[, 1] = 10) * (ivs[, 1]^1)}}
#' }
#' Setting \code{ord = 2} would retain both of the aforementioned
#' coefficients and also add
#' \itemize{
#'  \item{\code{I(ivs[, 1] = 10) * (ivs[, 1]^2)}}
#' }
#' to the block. Maintaining \code{ord = 2} and setting
#' \code{allow_interactions = TRUE} would result in the block containing a
#' total of seven coefficients:
#' \itemize{
#'   \item{\code{ivs[, 1] = 10}}
#'   \item{\code{I(ivs[, 1] = 10) * (ivs[, 1]^1)}}
#'   \item{\code{I(ivs[, 1] = 10) * (ivs[, 1]^2)}}
#'   \item{\code{I(ivs[, 1] = 10) * (ivs[, 2]^1)}}
#'   \item{\code{I(ivs[, 1] = 10) * (ivs[, 2]^2)}}
#'   \item{\code{I(ivs[, 1] = 10) * (ivs[, 3]^1)}}
#'   \item{\code{I(ivs[, 1] = 10) * (ivs[, 3]^2)}}
#' }
#' Note that if \code{ord = 0}, then each block will only include the
#' indicator function itself (in this example, \code{I(ivs[, 1] = 10)})
#' regardless of the value of \code{allow_interactions}.
#' 
#' During each forward selection iteration, the algorithm selects the "best"
#' block of coefficients that can be added. If \code{add_mode = "p"}, then all
#' possible candidate blocks that could be added to the model are evaluated,
#' and whichever block has a \emph{p}-value less than or equal to all
#' \emph{p}-values in all other blocks is treated as the "best" block. As long
#' as that \emph{p}-value is less than the threshold specified by the
#' \code{add} argument, all coefficients included in the block are added to the
#' model. If \code{add_mode = "mse"}, then the candidate block whose
#' coefficients would jointly reduce the model mean squared error by the
#' greatest amount is treated as the `best` block, and all of its coefficients
#' are added to the model if the resulting \emph{p}-value for at least one of
#' those coefficients would be less than \code{add}.
#' 
#' During each backward selection iteration, all coefficients listed in
#' \code{start_formula} are retained in the regression model. Among all
#' remaining coefficients, any blocks whose \emph{p}-values are all greater
#' than or equal to \code{remove} are removed from the model. If at least one
#' \emph{p}-value contained in a given block is less than \code{remove}, then
#' none of the coefficients in the block are removed from the model.
#' 
#' When \code{ord = 0}, \code{add_mode = "p"} is generally acceptable to follow
#' common conventions of stepwise model selection. When \code{ord > 0}, however,
#' \code{add_mode = "p"} is more likely to result in spurious breakpoints being
#' added to the model due to intercept and higher-order terms competing with one
#' another, and in rare cases the algorithm may entirely fail to converge.
#' \code{add_mode = "mse"} is more robust against these issues and is generally
#' recommended as whenever \code{ord > 0}.
#'
#' As a caveat, since compositional data are multidimensional in nature, there
#' may be rare cases in which the most suitable breakpoint selected based on MSE
#' does not correspond to a statistically significant change in any individual
#' category of the dependent variable. This can cause the model to prematurely
#' terminate. If an inspection of the final model suggests that such an issue
#' has occurred, then setting \code{add_mode = "p"} may allow that problem to be
#' overcome, regardless of other weaknesses in that option.
#'
#' The larger the value of \code{ord}, the more likely spurious breakpoints are
#' to emerge. As with other regression contexts, you should only increase the
#' complexity of your model (such as by increasing the value of \code{ord}) when
#' you have a clear reason to do so. Moreover, whenever \code{ord > 0},
#' singularities may occur that render some coefficients inestimable. For
#' instance, if \code{ord = 3} and a pair of breakpoints is identified along the
#' same independent variable, with those two breakpoints occurring two data
#' points apart, there will be insufficient data between the breakpoints to
#' estimate the standard error of the interaction between the indicator function
#' and a cubic term, so it will be reported as \code{NA}. This is normal,
#' expected behavior in stepwise segmented regression analysis. In such cases,
#' those \code{NA} coefficients can be treated as though they were absent from
#' the model.
#' 
#' Finally, use caution when setting \code{ord > 1}, as increasingly high-order
#' terms may compete with lower-order terms in the model and yield unpredictable
#' results. Refer to Britt (2015) for more information.
#' 
#' @param dv A data.frame representing values of the dependent variable
#' @param ivs A data.frame representing values of the independent variables,
#'   coerced to a data.frame if not provided as one (default =
#'   \code{data.frame(1:nrow(data.frame(dv)))}, which implicitly treats
#'   \code{dv} as sequentially ordered data and attempts to detect breakpoints
#'   in that sequence)
#' @param start_breakpoints A numeric vector of breakpoints to be added to the
#'   regression model before the first iteration of the stepwise procedure
#' @param retain_start_breakpoints If \code{TRUE}, then the breakpoints
#'   specified in \code{start_breakpoints} can never be removed from the model
#' @param add A numeric atomic vector between \code{0} and \code{1} indicating
#'   the \emph{p}-value threshold to add coefficients to the model during each
#'   forward selection iteration, which must be less than or equal to the
#'   \code{remove} argument in order to avoid infinite loops (default =
#'   \code{.15}, as recommended by Britt, 2015)
#' @param remove A numeric atomic vector between \code{0} and \code{1}
#'   indicating the \emph{p}-value threshold to remove coefficients from the
#'   model during each backward selection iteration, which must be greater than
#'   or equal to the \code{add} argument in order to avoid infinite loops
#'   (default = \code{.20}, as recommended by Britt, 2015)
#' @param bonferroni If \code{TRUE}, then the \emph{p}-values specified in
#'   \code{add} and \code{remove} are divided by the number of terms added to
#'   the model for each breakpoint, including all categories of the dependent
#'   variable, all independent variables (including any interactions), and all
#'   orders of those independent variables (\code{0:ord})
#' @param add_mode A character atomic vector (either \code{"p"} or \code{"mse"})
#'   indicating what criterion should be used to select the best candidate
#'   block during each forward selection iteration (default = \code{"mse"})
#' @param ord A non-negative numeric atomic vector indicating the maximum
#'   exponent that will be applied to coefficients added to the regression
#'   model (default = \code{1}; it is generally recommended to use either
#'   \code{order = 0} or \code{order = 1})
#' @param type The method used to conduct compositional regression; since
#'   \code{\link[Compositional]{comp.reg}} does not return standard errors for
#'   \code{"lmfit"} or \code{"spatial"}, at present, this must be set to
#'   \code{"classical"}
#' @param allow_interactions A boolean value indicating whether interactions
#'   may be created between indicator functions representing breakpoints along
#'   one independent variable with other independent variables in the model
#'   (default = \code{FALSE}, which restricts each indicator function to only
#'   interact with the independent variable along which its breakpoint was
#'   defined)
#' @param update A numeric value indicating how many loop iterations should
#'   elapse between progress updates (default = \code{0}, which suppresses
#'   output)
#' @param verbose A boolean value indicating whether the current regression
#'   model and results should be outputted as part of the progress report after
#'   every step of the stepwise procedure (default = \code{FALSE}
#' @return The final compositional regression model outputted from
#'   \code{\link[Compositional]{comp.reg}}
#' @section References:
#'   Britt, B. C. (2015). Stepwise segmented regression analysis: An iterative
#'   statistical algorithm to detect and quantify evolutionary and
#'   revolutionary transformations in longitudinal data. In S. A. Matei, M. G.
#'   Russell, & E. Bertino (Eds.), \emph{Transparency in social media: Tools,
#'   methods, and algorithms for mediating online interactions} (pp. 125-144).
#'   Heidelberg, Germany: Springer.
#' @examples
#' dv <- rbind(c(0.058, 0.050, 0.010, 0.221, 0.186, 0.068, 0.407),
#'             c(0.054, 0.065, 0.009, 0.262, 0.196, 0.069, 0.345),
#'             c(0.040, 0.114, 0.004, 0.235, 0.191, 0.084, 0.332),
#'             c(0.044, 0.070, 0.011, 0.356, 0.172, 0.049, 0.298),
#'             c(0.048, 0.123, 0.007, 0.277, 0.175, 0.053, 0.317),
#'             c(0.046, 0.101, 0.013, 0.240, 0.212, 0.105, 0.283),
#'             c(0.058, 0.125, 0.010, 0.257, 0.191, 0.050, 0.309),
#'             c(0.048, 0.065, 0.010, 0.264, 0.247, 0.083, 0.283),
#'             c(0.072, 0.067, 0.006, 0.237, 0.197, 0.097, 0.324),
#'             c(0.040, 0.192, 0.002, 0.243, 0.174, 0.069, 0.280),
#'             c(0.041, 0.175, 0.005, 0.208, 0.174, 0.075, 0.322),
#'             c(0.046, 0.098, 0.004, 0.274, 0.206, 0.058, 0.314),
#'             c(0.037, 0.087, 0.004, 0.237, 0.222, 0.092, 0.321),
#'             c(0.035, 0.084, 0.004, 0.287, 0.196, 0.067, 0.327),
#'             c(0.040, 0.061, 0.004, 0.271, 0.226, 0.059, 0.339),
#'             c(0.038, 0.239, 0.006, 0.228, 0.148, 0.052, 0.289),
#'             c(0.051, 0.141, 0.004, 0.243, 0.180, 0.068, 0.313),
#'             c(0.035, 0.156, 0.007, 0.227, 0.199, 0.063, 0.313),
#'             c(0.037, 0.188, 0.008, 0.258, 0.197, 0.048, 0.264),
#'             c(0.037, 0.126, 0.002, 0.309, 0.180, 0.032, 0.314),
#'             c(0.035, 0.092, 0.005, 0.237, 0.242, 0.059, 0.330),
#'             c(0.043, 0.114, 0.005, 0.225, 0.245, 0.057, 0.311),
#'             c(0.038, 0.139, 0.002, 0.270, 0.184, 0.053, 0.314),
#'             c(0.049, 0.080, 0.010, 0.243, 0.215, 0.068, 0.335),
#'             c(0.041, 0.081, 0.005, 0.219, 0.193, 0.053, 0.408),
#'             c(0.045, 0.115, 0.001, 0.266, 0.181, 0.052, 0.340),
#'             c(0.047, 0.140, 0.003, 0.253, 0.170, 0.066, 0.321),
#'             c(0.049, 0.226, 0.003, 0.239, 0.170, 0.047, 0.266),
#'             c(0.041, 0.119, 0.002, 0.358, 0.153, 0.035, 0.292),
#'             c(0.030, 0.165, 0.010, 0.391, 0.090, 0.030, 0.284))
#' ivs <- c(1:30)
#' model <- stepsegcomp(dv, ivs, verbose=TRUE)
#' # "Beginning the stepwise procedure..."
#' # "End of step 1 ~ Current Breakpoints: 23"
#' # "End of step 2 ~ Current Breakpoints: 4, 23"
#' # "No additional breakpoints can be added to the model."
#' model
#' # $runtime
#' #    user  system elapsed 
#' #       0       0       0 
#' #
#' # $be
#' #                                  [,1]        [,2]       [,3]        [,4]
#' # (Intercept)                -0.2878712 -1.96059047  1.0820294  1.12821655
#' # xivs[,1]^1                  0.2700093  0.06038650  0.2450410  0.08682669
#' # xI(ivs[,1]>4)*(ivs[,1]^0)   0.7953717  0.19964276  0.4145045  0.09793102
#' # xI(ivs[,1]>4)*(ivs[,1]^1)  -0.2344747 -0.08772856 -0.2252737 -0.06469739
#' # xI(ivs[,1]>23)*(ivs[,1]^0) -4.2608171 -1.32045504 -3.1891707  1.80206300
#' # xI(ivs[,1]>23)*(ivs[,1]^1)  0.1431590  0.05065426  0.1120019 -0.08473563
#' #                                   [,5]         [,6]
#' # (Intercept)                 0.22781059  1.919198434
#' # xivs[,1]^1                  0.03425128  0.015529618
#' # xI(ivs[,1]>4)*(ivs[,1]^0)   0.23577716 -0.241304627
#' # xI(ivs[,1]>4)*(ivs[,1]^1)  -0.03957347  0.005926647
#' # xI(ivs[,1]>23)*(ivs[,1]^0)  1.58830152  0.353969914
#' # xI(ivs[,1]>23)*(ivs[,1]^1) -0.06607795 -0.022286050
#' #
#' # $seb
#' #                                  [,1]      [,2]       [,3]       [,4]
#' # (Intercept)                0.50524000 0.7699158 0.24702365 0.22129975
#' # xivs[,1]^1                 0.18448756 0.2811335 0.09020028 0.08080724
#' # xI(ivs[,1]>4)*(ivs[,1]^0)  0.56810356 0.8657112 0.27775912 0.24883456
#' # xI(ivs[,1]>4)*(ivs[,1]^1)  0.18529495 0.2823639 0.09059504 0.08116089
#' # xI(ivs[,1]>23)*(ivs[,1]^0) 2.12661681 3.2406697 1.03975270 0.93147764
#' # xI(ivs[,1]>23)*(ivs[,1]^1) 0.07985208 0.1216835 0.03904155 0.03497594
#' #                                  [,5]       [,6]
#' # (Intercept)                0.32977623 0.21663586
#' # xivs[,1]^1                 0.12041725 0.07910423
#' # xI(ivs[,1]>4)*(ivs[,1]^0)  0.37080803 0.24359039
#' # xI(ivs[,1]>4)*(ivs[,1]^1)  0.12094425 0.07945042
#' # xI(ivs[,1]>23)*(ivs[,1]^0) 1.38806839 0.91184680
#' # xI(ivs[,1]>23)*(ivs[,1]^1) 0.05212042 0.03423883
#' #
#' # $est
#' # NULL
#' # 
#' # $t
#' #                                  [,1]       [,2]      [,3]       [,4]
#' # (Intercept)                -0.5697712 -2.5464998  4.380266  5.0981376
#' # xivs[,1]^1                  1.4635635  0.2147965  2.716632  1.0744914
#' # xI(ivs[,1]>4)*(ivs[,1]^0)   1.4000470  0.2306113  1.492316  0.3935588
#' # xI(ivs[,1]>4)*(ivs[,1]^1)  -1.2654136 -0.3106933 -2.486601 -0.7971499
#' # xI(ivs[,1]>23)*(ivs[,1]^0) -2.0035660 -0.4074636 -3.067240  1.9346283
#' # xI(ivs[,1]>23)*(ivs[,1]^1)  1.7928026  0.4162787  2.868788 -2.4226831
#' #                                  [,5]        [,6]
#' # (Intercept)                 0.6908036  8.85909836
#' # xivs[,1]^1                  0.2844383  0.19631842
#' # xI(ivs[,1]>4)*(ivs[,1]^0)   0.6358470 -0.99061637
#' # xI(ivs[,1]>4)*(ivs[,1]^1)  -0.3272042  0.07459553
#' # xI(ivs[,1]>23)*(ivs[,1]^0)  1.1442531  0.38819012
#' # xI(ivs[,1]>23)*(ivs[,1]^1) -1.2677940 -0.65089993
#' #
#' # $p
#' #                                  [,1]       [,2]         [,3]         [,4]
#' # (Intercept)                0.57412808 0.01771475 0.0002008412 3.242277e-05
#' # xivs[,1]^1                 0.15628846 0.83174067 0.0120382387 2.932869e-01
#' # xI(ivs[,1]>4)*(ivs[,1]^0)  0.17429429 0.81957130 0.1486457861 6.973815e-01
#' # xI(ivs[,1]>4)*(ivs[,1]^1)  0.21787098 0.75871609 0.0202485998 4.331800e-01
#' # xI(ivs[,1]>23)*(ivs[,1]^0) 0.05653118 0.68727920 0.0052865773 6.490810e-02
#' # xI(ivs[,1]>23)*(ivs[,1]^1) 0.08561900 0.68090526 0.0084573677 2.332019e-02
#' #                                 [,5]         [,6]
#' # (Intercept)                0.4963166 4.953975e-09
#' # xivs[,1]^1                 0.7785143 8.460138e-01
#' # xI(ivs[,1]>4)*(ivs[,1]^0)  0.5308933 3.317559e-01
#' # xI(ivs[,1]>4)*(ivs[,1]^1)  0.7463503 9.411548e-01
#' # xI(ivs[,1]>23)*(ivs[,1]^0) 0.2638031 7.012974e-01
#' # xI(ivs[,1]>23)*(ivs[,1]^1) 0.2170338 5.212948e-01
#' @export
#' stepsegcomp

stepsegcomp <- function(dv, ivs = data.frame(1:nrow(data.frame(dv))),
                        start_breakpoints = numeric(0),
                        retain_start_breakpoints = TRUE, add = .15,
                        remove = .20, bonferroni = TRUE, add_mode = "mse",
                        ord = 1, type = "classical", allow_interactions = FALSE,
                        update = 0, verbose = FALSE) {

   if(add > remove) {
      stop("'add' cannot be set to a value greater than 'remove,' as doing so may result in an infinite loop.")
   }
   if(ord > 1) {
      warning("Stepwise segmented regression may produce unpredictable results when ord > 1.")
   }
   if((ord > 0) & (add_mode == "p")) {
      warning("When ord > 0, setting add_mode = 'p' may yield spurious results. In some cases, the algorithm may entirely fail to converge. In most cases, add_mode = 'mse' is a better choice.")
   }
   ivs <- as.data.frame(ivs)
   if(ncol(ivs) > 1) {
      warning("The use of multiple independent variables is currently experimental and may result in errors in some cases.")
   }
   if(allow_interactions) {
      warning("The allow_interactions parameter is currently experimental.")
   }

   if(verbose | (update > 0)) {
      print("Beginning the stepwise procedure...")
      utils::flush.console()
   }

   #Initial setup
   dv <- as.matrix(dv)
   current_step <- 0
   current_breakpoints <- start_breakpoints
   start_terms <- data.frame(matrix(NA, nrow = nrow(as.data.frame(dv)), ncol = 0))
   if(ord > 0) { #We don't add an intercept here, as that's automatically incorporated into the model
      for(i in 1:ord) {
         for(j in 1:ncol(ivs)) {
            start_terms <- cbind(start_terms, as.numeric(ivs[,j]^i))
            colnames(start_terms)[ncol(start_terms)] <- paste0("ivs[,", j, "]^", i)
         }
      }
   }

   #Conduct stepwise procedure
   repeat {
      current_step <- current_step + 1
      current_terms <- start_terms

      #Set up data set based on current breakpoints
      results_best_addition_tvalues <- rep(0, nrow(dv))
      results_best_addition_sse <- rep(Inf, nrow(dv))

      #Test all candidate breakpoints for forward selection
      for(i in (ord+1):(nrow(dv)-(ord+1))) {
         if(update > 0) {
            if(i %% update == 0) {
               print(paste0("Testing observation ", i))
               utils::flush.console()
            }
         }
         if(suppressWarnings({min(abs(current_breakpoints-i))}) > ord) {
            temp_terms <- start_terms
            temp_breakpoints <- c(current_breakpoints, i)
            temp_breakpoints <- temp_breakpoints[order(temp_breakpoints)]
            for(j in 1:length(temp_breakpoints)) {
               for(k in 0:ord) {
                  for(l in 1:ncol(ivs)) {
                     if(allow_interactions) {
                        for(m in 1:ncol(ivs)) {
                           temp_terms <- cbind(temp_terms, (as.numeric(I(ivs[,l] > ivs[temp_breakpoints[j],l]))*(ivs[,m]^k)))
                           colnames(temp_terms)[ncol(temp_terms)] <- paste0("I(ivs[,", l, "]>", ivs[temp_breakpoints[j],l], ")*(ivs[,", m, "]^", k, ")")                           
                        }
                     } else {
                        temp_terms <- cbind(temp_terms, (as.numeric(I(ivs[,l] > ivs[temp_breakpoints[j],l]))*(ivs[,l]^k)))
                        colnames(temp_terms)[ncol(temp_terms)] <- paste0("I(ivs[,", l, "]>", ivs[temp_breakpoints[j],l], ")*(ivs[,", l, "]^", k, ")")
                     }
                  }
               }
            }
            model <- Compositional::comp.reg(y=dv, x=temp_terms, type=type, xnew=as.matrix(temp_terms))
            results_best_addition_tvalues[i] <- get_t_from_model(model, ord, which(temp_breakpoints==i), allow_interactions, ivs)
            results_best_addition_sse[i] <- sum((model$est - dv)^2)
         }
      }

      #Identify the best candidate breakpoint
      if(add_mode=="mse") {
         if(min(results_best_addition_sse) < Inf) {
            best_index <- which(results_best_addition_sse == min(results_best_addition_sse, na.rm=T))[1]
         } else {
            print("No additional breakpoints can be added to the model.")
            utils::flush.console()
            break
         }
      } else if(add_mode=="p") {
         if(max(abs(results_best_addition_tvalues)) > 0) {
            best_index <- which(abs(results_best_addition_tvalues) == max(abs(results_best_addition_tvalues), na.rm=T))[1]
         } else {
            print("No additional breakpoints can be added to the model.")
            utils::flush.console()
            break
         }
      }

      #Prepare a model that includes the best candidate breakpoint
      temp_terms <- start_terms
      temp_breakpoints <- c(current_breakpoints, best_index)
      temp_breakpoints <- temp_breakpoints[order(temp_breakpoints)]
      for(i in 1:length(temp_breakpoints)) {
         if(update > 0) {
            if(i %% update == 0) {
               print(paste0("Testing observation ", i))
               utils::flush.console()
            }
         }
         for(j in 0:ord) {
            for(k in 1:ncol(ivs)) {
               if(allow_interactions) {
                  for(l in 1:ncol(ivs)) {
                     temp_terms <- cbind(temp_terms, (as.numeric(I(ivs[,k] > ivs[temp_breakpoints[i],k]))*(ivs[,l]^j)))
                     colnames(temp_terms)[ncol(temp_terms)] <- paste0("I(ivs[,", k, "]>", ivs[temp_breakpoints[i],k], ")*(ivs[,", l, "]^", j, ")")
                  }
               } else {
                  temp_terms <- cbind(temp_terms, (as.numeric(I(ivs[,k] > ivs[temp_breakpoints[i],k]))*(ivs[,k]^j)))
                  colnames(temp_terms)[ncol(temp_terms)] <- paste0("I(ivs[,", k, "]>", ivs[temp_breakpoints[i],k], ")*(ivs[,", k, "]^", j, ")")
               }
            }
         }
      }
      model <- Compositional::comp.reg(y=dv, x=temp_terms, type=type, xnew=as.matrix(temp_terms))

      #If possible, add the best candidate breakpoint to the model
      if(allow_interactions) {
         if(bonferroni) {
            threshold <- stats::qt(1-(add/(ncol(dv)*(ncol(ivs)^2)*(ord+1))), nrow(dv)-(nrow(model$be)))
         } else {
            threshold <- stats::qt(1-add, nrow(dv)-(nrow(model$be)))
         }
      } else {
         if(bonferroni) {
            threshold <- stats::qt(1-(add/(ncol(dv)*ncol(ivs)*(ord+1))), nrow(dv)-(nrow(model$be)))
         } else {
            threshold <- stats::qt(1-add, nrow(dv)-(nrow(model$be)))
         }
      }
      if(abs(results_best_addition_tvalues[best_index]) > threshold) {
         current_breakpoints <- c(current_breakpoints, best_index) #Add the breakpoint
         current_breakpoints <- current_breakpoints[order(current_breakpoints)]
      } else {
         if(verbose) {
            print("No additional breakpoints can be added to the model.")
            utils::flush.console()
         }
         break
      }

      #Backward selection
      current_terms <- start_terms
      for(i in 1:length(current_breakpoints)) {
         for(j in 0:ord) {
            for(k in 1:ncol(ivs)) {
               if(allow_interactions) {
                  for(l in 1:ncol(ivs)) {
                     current_terms <- cbind(current_terms, (as.numeric(I(ivs[,k] > ivs[current_breakpoints[i],k]))*(ivs[,l]^j)))
                     colnames(current_terms)[ncol(current_terms)] <- paste0("I(ivs[,", k, "]>", ivs[current_breakpoints[i],k], ")*(ivs[,", l, "]^", j, ")")
                  }
               } else {
                  current_terms <- cbind(current_terms, (as.numeric(I(ivs[,k] > ivs[current_breakpoints[i],k]))*(ivs[,k]^j)))
                  colnames(current_terms)[ncol(current_terms)] <- paste0("I(ivs[,", k, "]>", ivs[current_breakpoints[i],k], ")*(ivs[,", k, "]^", j, ")")
               }
            }
         }
      }
      model <- Compositional::comp.reg(dv, current_terms, type)
      tstatistics <- numeric(0)
      for(i in 1:length(current_breakpoints)) {
         tstatistics[i] <- get_t_from_model(model, ord, i, allow_interactions, ivs)
      }
      if(allow_interactions) {
         if(bonferroni) {
            threshold <- stats::qt(1-(remove/(ncol(dv)*(ncol(ivs)^2)*(ord+1))), nrow(dv)-(nrow(model$be)))
         } else {
            threshold <- stats::qt(1-remove, nrow(dv)-(nrow(model$be)))
         }
      } else {
         if(bonferroni) {
            threshold <- stats::qt(1-(remove/(ncol(dv)*ncol(ivs)*(ord+1))), nrow(dv)-(nrow(model$be)))
         } else {
            threshold <- stats::qt(1-remove, nrow(dv)-(nrow(model$be)))
         }
      }
      current_breakpoints <- current_breakpoints[abs(tstatistics) > threshold]
      if(retain_start_breakpoints && (length(start_breakpoints) > 0)) {
         for(i in 1:length(start_breakpoints)) {
            if(!(start_breakpoints[i] %in% current_breakpoints)) { #If any starting breakpoint has been erroneously removed, we must immediately add it back into the model
               current_breakpoints <- c(current_breakpoints, start_breakpoints[i])
            }
         }
         current_breakpoints <- current_breakpoints[order(current_breakpoints)]
      }

      #Status Update
      if(verbose) {
         print(paste0("End of step ", current_step, " ~ Current Breakpoints: ", paste(current_breakpoints, collapse=", ")))
         utils::flush.console()
      }
   }

   #Final model
   current_terms <- start_terms
   if(length(current_breakpoints) > 0) {
      for(i in 1:length(current_breakpoints)) {
         for(j in 0:ord) {
            for(k in 1:ncol(ivs)) {
               current_terms <- cbind(current_terms, (as.numeric(I(ivs[,k] > ivs[current_breakpoints[i],k]))*(ivs[,k]^j)))
               colnames(current_terms)[ncol(current_terms)] <- paste0("I(ivs[,", k, "]>", ivs[current_breakpoints[i],k], ")*(ivs[,", k, "]^", j, ")")
            }
         }
      }
   }
   model <- Compositional::comp.reg(y=dv, x=current_terms, type=type)
   model$t <- model$be / model$seb
   model$p <- (1-stats::pt(abs(model$t), nrow(dv)-(nrow(model$t))))*2
   return(model)
}

#' Stepwise Segmented Dirichlet Regression Analysis
#'
#' @description This function performs stepwise segmented Dirichlet regression
#' analysis, a modified form of the procedure described by Britt (2015).
#'
#' @details When using any form of stepwise segmented regression analysis,
#' coefficients are added to and removed from the model in "blocks," with a
#' given block consisting of the indicator function representing a given
#' breakpoint location as well as all interaction terms between that indicator
#' function and independent variables that are eligible to be added to the
#' model.
#'
#' The \code{ord} and \code{allow_interactions} arguments jointly indicate what
#' interaction terms are valid. \code{ord} indicates the maximum exponent that
#' can be applied to the independent variables, while \code{allow_interactions}
#' indicates whether interaction terms can be created between a given
#' independent variable and an indicator function that is defined by values of a
#' different independent variable.
#'
#' All interaction terms that are considered valid, with exponents ranging from
#' \code{0} to \code{ord}, are included in a given block. Consider, for
#' instance, an analysis in which \code{ivs}, the data.frame representing the
#' independent variables, has three columns signifying three variables. If
#' \code{ord = 1} and \code{allow_interactions = FALSE}, which are their default
#' values in \code{stepseg}, then the block corresponding to a potential
#' breakpoint at \code{ivs[, 1] = 10} would include two coefficients:
#' \itemize{
#'   \item{\code{I(ivs[, 1] = 10)}}
#'   \item{\code{I(ivs[, 1] = 10) * (ivs[, 1]^1)}}
#' }
#' Setting \code{ord = 2} would retain both of the aforementioned
#' coefficients and also add
#' \itemize{
#'  \item{\code{I(ivs[, 1] = 10) * (ivs[, 1]^2)}}
#' }
#' to the block. Maintaining \code{ord = 2} and setting
#' \code{allow_interactions = TRUE} would result in the block containing a
#' total of seven coefficients:
#' \itemize{
#'   \item{\code{ivs[, 1] = 10}}
#'   \item{\code{I(ivs[, 1] = 10) * (ivs[, 1]^1)}}
#'   \item{\code{I(ivs[, 1] = 10) * (ivs[, 1]^2)}}
#'   \item{\code{I(ivs[, 1] = 10) * (ivs[, 2]^1)}}
#'   \item{\code{I(ivs[, 1] = 10) * (ivs[, 2]^2)}}
#'   \item{\code{I(ivs[, 1] = 10) * (ivs[, 3]^1)}}
#'   \item{\code{I(ivs[, 1] = 10) * (ivs[, 3]^2)}}
#' }
#' Note that if \code{ord = 0}, then each block will only include the
#' indicator function itself (in this example, \code{I(ivs[, 1] = 10)})
#' regardless of the value of \code{allow_interactions}.
#' 
#' During each forward selection iteration, the algorithm selects the "best"
#' block of coefficients that can be added. If \code{add_mode = "p"}, then all
#' possible candidate blocks that could be added to the model are evaluated,
#' and whichever block has a \emph{p}-value less than or equal to all
#' \emph{p}-values in all other blocks is treated as the "best" block. As long
#' as that \emph{p}-value is less than the threshold specified by the
#' \code{add} argument, all coefficients included in the block are added to the
#' model. If \code{add_mode = "mse"}, then the candidate block whose
#' coefficients would jointly reduce the model mean squared error by the
#' greatest amount is treated as the `best` block, and all of its coefficients
#' are added to the model if the resulting \emph{p}-value for at least one of
#' those coefficients would be less than \code{add}.
#' 
#' During each backward selection iteration, all coefficients listed in
#' \code{start_formula} are retained in the regression model. Among all
#' remaining coefficients, any blocks whose \emph{p}-values are all greater
#' than or equal to \code{remove} are removed from the model. If at least one
#' \emph{p}-value contained in a given block is less than \code{remove}, then
#' none of the coefficients in the block are removed from the model.
#' 
#' When \code{ord = 0}, \code{add_mode = "p"} is generally acceptable to follow
#' common conventions of stepwise model selection. When \code{ord > 0}, however,
#' \code{add_mode = "p"} is more likely to result in spurious breakpoints being
#' added to the model due to intercept and higher-order terms competing with one
#' another, and in rare cases the algorithm may entirely fail to converge.
#' \code{add_mode = "mse"} is more robust against these issues and is generally
#' recommended as whenever \code{ord > 0}.
#'
#' As a caveat, since the Dirichlet distribution is multidimensional in nature,
#' there may be rare cases in which the most suitable breakpoint selected based
#' on MSE does not correspond to a statistically significant change in any
#' individual category of the dependent variable. This can cause the model to
#' prematurely terminate. If an inspection of the final model suggests that such
#' an issue has occurred, then setting \code{add_mode = "p"} may allow that
#' problem to be overcome, regardless of other weaknesses in that option.
#'
#' The larger the value of \code{ord}, the more likely spurious breakpoints are
#' to emerge. As with other regression contexts, you should only increase the
#' complexity of your model (such as by increasing the value of \code{ord}) when
#' you have a clear reason to do so. Moreover, whenever \code{ord > 0},
#' singularities may occur that render some coefficients inestimable. For
#' instance, if \code{ord = 3} and a pair of breakpoints is identified along the
#' same independent variable, with those two breakpoints occurring two data
#' points apart, there will be insufficient data between the breakpoints to
#' estimate the standard error of the interaction between the indicator function
#' and a cubic term, so it will be reported as \code{NA}. This is normal,
#' expected behavior in stepwise segmented regression analysis. In such cases,
#' those \code{NA} coefficients can be treated as though they were absent from
#' the model.
#' 
#' Finally, use caution when setting \code{ord > 1}, as increasingly high-order
#' terms may compete with lower-order terms in the model and yield unpredictable
#' results. Refer to Britt (2015) for more information.
#' 
#' @param dv A data.frame representing values of the dependent variable
#' @param ivs A data.frame representing values of the independent variables,
#'   coerced to a data.frame if not provided as one (default =
#'   \code{data.frame(1:nrow(data.frame(dv)))}, which implicitly treats
#'   \code{dv} as sequentially ordered data and attempts to detect breakpoints
#'   in that sequence)
#' @param start_breakpoints A numeric vector of breakpoints to be added to the
#'   regression model before the first iteration of the stepwise procedure
#' @param retain_start_breakpoints If \code{TRUE}, then the breakpoints
#'   specified in \code{start_breakpoints} can never be removed from the model
#' @param add A numeric atomic vector between \code{0} and \code{1} indicating
#'   the \emph{p}-value threshold to add coefficients to the model during each
#'   forward selection iteration, which must be less than or equal to the
#'   \code{remove} argument in order to avoid infinite loops (default =
#'   \code{.15}, as recommended by Britt, 2015)
#' @param remove A numeric atomic vector between \code{0} and \code{1}
#'   indicating the \emph{p}-value threshold to remove coefficients from the
#'   model during each backward selection iteration, which must be greater than
#'   or equal to the \code{add} argument in order to avoid infinite loops
#'   (default = \code{.20}, as recommended by Britt, 2015)
#' @param bonferroni If \code{TRUE}, then the \emph{p}-values specified in
#'   \code{add} and \code{remove} are divided by the number of terms added to
#'   the model for each breakpoint, including all categories of the dependent
#'   variable, all independent variables (including any interactions), and all
#'   orders of those independent variables (\code{0:ord})
#' @param add_mode A character atomic vector (either \code{"p"} or \code{"mse"})
#'   indicating what criterion should be used to select the best candidate
#'   block during each forward selection iteration (default = \code{"mse"})
#' @param ord A non-negative numeric atomic vector indicating the maximum
#'   exponent that will be applied to coefficients added to the regression
#'   model (default = \code{1}; it is generally recommended to use either
#'   \code{order = 0} or \code{order = 1})
#' @param parameterization The parameterization used to perform Dirichlet
#'   regression; since \link{DirichReg} does not return standard errors for
#'   \code{"lmfit"} or \code{"spatial"}, at present, this
#'   must be set to \code{"classical"}
#' @param allow_interactions A boolean value indicating whether interactions
#'   may be created between indicator functions representing breakpoints along
#'   one independent variable with other independent variables in the model
#'   (default = \code{FALSE}, which restricts each indicator function to only
#'   interact with the independent variable along which its breakpoint was
#'   defined)
#' @param iterlim The maximum number of iterations permitted for convergence in
#'   each regression, as conducted by \link{DirichReg}
#' @param tol1 The convergence criterion for BFGS optimization in each
#'   regression, as conducted by \link{DirichReg}
#' @param tol2 The convergence criterion for NR optimization in each regression,
#'   as conducted by \link{DirichReg}
#' @param update A numeric value indicating how many loop iterations should
#'   elapse between progress updates (default = \code{0}, which suppresses
#'   output)
#' @param verbose A boolean value indicating whether the current regression
#'   model and results should be outputted as part of the progress report after
#'   every step of the stepwise procedure (default = \code{FALSE}
#' @return The final Dirichlet regression model outputted from \link{DirichReg}
#' @section References:
#'   Britt, B. C. (2015). Stepwise segmented regression analysis: An iterative
#'   statistical algorithm to detect and quantify evolutionary and
#'   revolutionary transformations in longitudinal data. In S. A. Matei, M. G.
#'   Russell, & E. Bertino (Eds.), \emph{Transparency in social media: Tools,
#'   methods, and algorithms for mediating online interactions} (pp. 125-144).
#'   Heidelberg, Germany: Springer.
#' @examples
#' set.seed(2023)
#' dv <- rbind(DirichletReg::rdirichlet(12, c(8,8,32)),
#'             DirichletReg::rdirichlet(18, c(40,24,4)))
#' model <- stepsegdir(dv, add_mode="mse", update=10, verbose=TRUE)
#' # "Beginning the stepwise procedure..."
#' # "Testing observation 10"
#' # "Testing observation 20"
#' # "End of step 1 ~ Current Breakpoints: 10"
#' # "Testing observation 10"
#' # "Testing observation 20"
#' # "No additional breakpoints can be added to the model."
#' # Warning message:
#' # In min(abs(summary(model)$coef.mat[indices, 4]), na.rm = T) :
#' #   no non-missing arguments to min; returning Inf
#' summary(model)
#' #
#' # Call:
#' # DirichReg(formula = as.formula(paste0("dv ~ ", paste(current_terms, collapse
#' # = " + "))), data = df, model = parameterization, control = list(iterlim =
#' # iterlim, tol1 = tol1, tol2 = tol2))
#' #
#' # Standardized Residuals:
#' #         Min       1Q   Median      3Q     Max
#' # V1  -2.1452  -0.6651  -0.2261  0.6674  2.1566
#' # V2  -1.7372  -1.1093   0.3604  0.8497  1.6465
#' # V3  -1.7572  -0.6547  -0.2742  0.6269  2.2728
#' #
#' # ------------------------------------------------------------------
#' # Beta-Coefficients for variable no. 1: V1
#' #                                    Estimate Std. Error z value Pr(>|z|)    
#' # (Intercept)                         4.57845    0.68780   6.657  2.8e-11 ***
#' # I(ivs[, 1]^1)                      -0.22604    0.09559  -2.365   0.0180 *  
#' # I(ivs[, 1] > 12)TRUE               -0.54752    1.24725  -0.439   0.6607    
#' # I(ivs[, 1]^1):I(ivs[, 1] > 12)TRUE  0.21705    0.10660   2.036   0.0417 *  
#' # ------------------------------------------------------------------
#' # Beta-Coefficients for variable no. 2: V2
#' #                                    Estimate Std. Error z value Pr(>|z|)    
#' # (Intercept)                         4.35351    0.69899   6.228 4.72e-10 ***
#' # I(ivs[, 1]^1)                      -0.18497    0.09834  -1.881    0.060 .  
#' # I(ivs[, 1] > 12)TRUE               -0.84702    1.24069  -0.683    0.495    
#' # I(ivs[, 1]^1):I(ivs[, 1] > 12)TRUE  0.17734    0.10875   1.631    0.103    
#' # ------------------------------------------------------------------
#' # Beta-Coefficients for variable no. 3: V3
#' #                                    Estimate Std. Error z value Pr(>|z|)    
#' # (Intercept)                         5.97326    0.68569   8.711  < 2e-16 ***
#' # I(ivs[, 1]^1)                      -0.21326    0.09505  -2.244  0.02485 *  
#' # I(ivs[, 1] > 12)TRUE               -3.70315    1.21792  -3.041  0.00236 ** 
#' # I(ivs[, 1]^1):I(ivs[, 1] > 12)TRUE  0.18075    0.10538   1.715  0.08630 .  
#' # ------------------------------------------------------------------
#' # Significance codes: 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#' #
#' # Log-likelihood: 118.8 on 12 df (83 BFGS + 2 NR Iterations)
#' # AIC: -213.5, BIC: -196.7
#' # Number of Observations: 30
#' # Link: Log
#' # Parametrization: common
#' @export
#' stepsegdir

stepsegdir <- function(dv, ivs = data.frame(1:nrow(data.frame(dv))),
                       start_breakpoints = numeric(0),
                       retain_start_breakpoints = TRUE, add = .15,
                       remove = .20, bonferroni = TRUE, add_mode = "mse",
                       ord = 1, parameterization = "common",
                       allow_interactions = FALSE, iterlim = 1000, tol1 = 1e-5,
                       tol2 = 1e-10, update = 0, verbose = FALSE) {

   if(ord == 0) {
      stop("Due to a limitation in the DirichletReg package, intercept terms alone cannot be modeled for individual breakpoints.")
   }
   if(add > remove) {
      stop("'add' cannot be set to a value greater than 'remove,' as doing so may result in an infinite loop.")
   }
   if(ord > 1) {
      warning("Stepwise segmented regression may produce unpredictable results when ord > 1.")
   }
   if((ord > 0) & (add_mode == "p")) {
      warning("When ord > 0, setting add_mode = 'p' may yield spurious results. In some cases, the algorithm may entirely fail to converge. In most cases, add_mode = 'mse' is a better choice.")
   }
   ivs <- as.data.frame(ivs)
   if(ncol(ivs) > 1) {
      warning("The use of multiple independent variables is currently experimental and may result in errors in some cases.")
   }
   if(allow_interactions) {
      if(ncol(ivs) > 2) {
         warning("The allow_interactions parameter is currently experimental.")
      } else {
         warning("The ivs object contains only one column, so there are not multiple independent variables from which to construct interactions. Therefore, the allow_interactions argument has been ignored.")
         allow_interactions = FALSE
      }
   }

   if(verbose | (update > 0)) {
      print("Beginning the stepwise procedure...")
      utils::flush.console()
   }

   #Initial setup
   df <- data.frame(matrix(NA, nrow = nrow(as.data.frame(dv)), ncol = 0))
   df$dv <- DirichletReg::DR_data(as.data.frame(dv))
   df$ivs <- ivs
   current_step <- 0
   current_breakpoints <- start_breakpoints
   start_terms <- character(0)
   if(ord > 0) { #We don't add an intercept here, as that's automatically incorporated into the model
      for(i in 1:ord) {
         for(j in 1:ncol(df$ivs)) {
            start_terms <- c(start_terms, paste0("I(ivs[,", j, "]^", i, ")"))
         }
      }
   }

   #Conduct stepwise procedure
   repeat {
      current_step <- current_step + 1
      current_terms <- start_terms

      #Set up data set based on current breakpoints
      results_best_addition_pvalues <- rep(Inf, nrow(df$dv))
      results_best_addition_sse <- rep(Inf, nrow(df$dv))

      #Test all candidate breakpoints for forward selection
      for(i in (ord+1):(nrow(df$dv)-(ord+1))) {
         if(update > 0) {
            if(i %% update == 0) {
               print(paste0("Testing observation ", i))
               utils::flush.console()
            }
         }
         if(suppressWarnings({min(abs(current_breakpoints-i))}) > ord) {
            temp_terms <- start_terms
            temp_breakpoints <- c(current_breakpoints, i)
            temp_breakpoints <- temp_breakpoints[order(temp_breakpoints)]
            for(j in 1:length(temp_breakpoints)) {
               for(k in 1:ord) {
                  for(l in 1:ncol(df$ivs)) {
                     if(allow_interactions) {
                        for(m in 1:ncol(df$ivs)) {
                           temp_terms <- c(temp_terms, paste0("I(ivs[,", l, "]>", ivs[temp_breakpoints[j],l], ")*(I(ivs[,", m, "]^", k, "))"))
                        }
                     } else {
                        temp_terms <- c(temp_terms, paste0("I(ivs[,", l, "]>", ivs[temp_breakpoints[j],l], ")*(I(ivs[,", l, "]^", k, "))"))
                     }
                  }
               }
            }
            tryCatch({
               model <- DirichReg(as.formula(paste0("dv ~ ", paste(temp_terms, collapse=" + "))), df, model=parameterization, control=list(iterlim=iterlim, tol1=tol1, tol2=tol2))
               results_best_addition_pvalues[i] <- get_p_from_model(model, ord, temp_breakpoints, which(temp_breakpoints==i), allow_interactions, df$ivs, parameterization)
               results_best_addition_sse[i] <- sum((model$fitted.values$mu-df$dv)^2)
            }, error = function(e) {
               results_best_addition_pvalues[i] <- Inf
               results_best_addition_sse[i] <- Inf
            })
         }
      }

      #Identify the best candidate breakpoint
      if(add_mode=="mse") {
         if(min(results_best_addition_sse) < Inf) {
            best_index <- which(results_best_addition_sse == min(results_best_addition_sse, na.rm=T))[1]
         } else {
            print("No additional breakpoints can be added to the model.")
            utils::flush.console()
            break
         }
      } else if(add_mode=="p") {
         if(min(abs(results_best_addition_pvalues)) < Inf) {
            best_index <- which(abs(results_best_addition_pvalues) == min(abs(results_best_addition_pvalues), na.rm=T))[1]
         } else {
            print("No additional breakpoints can be added to the model.")
            utils::flush.console()
            break
         }
      }

      #Prepare a model that includes the best candidate breakpoint
      temp_terms <- start_terms
      temp_breakpoints <- c(current_breakpoints, best_index)
      temp_breakpoints <- temp_breakpoints[order(temp_breakpoints)]
      for(i in 1:length(temp_breakpoints)) {
         for(j in 1:ord) {
            for(k in 1:ncol(df$ivs)) {
               if(allow_interactions) {
                  for(l in 1:ncol(df$ivs)) {
                     temp_terms <- c(temp_terms, paste0("I(ivs[,", k, "]>", ivs[temp_breakpoints[i],k], ")*(I(ivs[,", l, "]^", j, "))"))
                  }
               } else {
                  temp_terms <- c(temp_terms, paste0("I(ivs[,", k, "]>", ivs[temp_breakpoints[i],k], ")*(I(ivs[,", k, "]^", j, "))"))
               }
            }
         }
      }
      tryCatch({
         model <- DirichReg(as.formula(paste0("dv ~ ", paste(temp_terms, collapse=" + "))), df, model=parameterization, control=list(iterlim=iterlim, tol1=tol1, tol2=tol2))
      }, error = function(e) {
      })

      #If possible, add the best candidate breakpoint to the model
      if(allow_interactions) {
         if(bonferroni) {
            threshold <- add/((ncol(df$dv)-1)*(ncol(df$ivs) + ((ncol(df$ivs)^2)*ord)))
         } else {
            threshold <- add
         }
      } else {
         if(bonferroni) {
            threshold <- add/(ncol(df$dv)*(ncol(df$ivs)*(ord+1)))
         } else {
            threshold <- add
         }
      }
      if(abs(results_best_addition_pvalues[best_index]) < threshold) {
         current_breakpoints <- c(current_breakpoints, best_index) #Add the breakpoint
         current_breakpoints <- current_breakpoints[order(current_breakpoints)]
      } else {
         if(verbose) {
            print("No additional breakpoints can be added to the model.")
            utils::flush.console()
         }
         break
      }

      #Backward selection
      current_terms <- start_terms
      for(i in 1:length(current_breakpoints)) {
         for(j in 1:ord) {
            for(k in 1:ncol(df$ivs)) {
               if(allow_interactions) {
                  for(l in 1:ncol(df$ivs)) {
                     current_terms <- c(current_terms, paste0("I(ivs[,", k, "]>", ivs[current_breakpoints[i],k], ")*(I(ivs[,", l, "]^", j, "))"))
                  }
               } else {
                  current_terms <- c(current_terms, paste0("I(ivs[,", k, "]>", ivs[current_breakpoints[i],k], ")*(I(ivs[,", k, "]^", j, "))"))
               }
            }
         }
      }
      tryCatch({
         model <- DirichReg(as.formula(paste0("dv ~ ", paste(current_terms, collapse=" + "))), df, model=parameterization, control=list(iterlim=iterlim, tol1=tol1, tol2=tol2))
      }, error = function(e) {
      })
      pvalues <- numeric(0)
      for(i in 1:length(current_breakpoints)) {
         pvalues[i] <- get_p_from_model(model, ord, current_breakpoints, i, allow_interactions, df$ivs, parameterization)
      }
      if(allow_interactions) {
         if(bonferroni) {
            threshold <- remove/((ncol(df$dv)-1)*(ncol(df$ivs) + ((ncol(df$ivs)^2)*ord)))
         } else {
            threshold <- remove
         }
      } else {
         if(bonferroni) {
            threshold <- remove/(ncol(df$dv)*(ncol(df$ivs)*(ord+1)))
         } else {
            threshold <- remove
         }
      }
      current_breakpoints <- current_breakpoints[pvalues < threshold]
      if(retain_start_breakpoints && (length(start_breakpoints) > 0)) {
         for(i in 1:length(start_breakpoints)) {
            if(!(start_breakpoints[i] %in% current_breakpoints)) { #If any starting breakpoint has been erroneously removed, we must immediately add it back into the model
               current_breakpoints <- c(current_breakpoints, start_breakpoints[i])
            }
         }
         current_breakpoints <- current_breakpoints[order(current_breakpoints)]
      }

      #Status Update
      if(verbose) {
         print(paste0("End of step ", current_step, " ~ Current Breakpoints: ", paste(current_breakpoints, collapse=", ")))
         utils::flush.console()
      }
   }

   #Final model
   current_terms <- start_terms
   if(length(current_breakpoints) > 0) {
      for(i in 1:length(current_breakpoints)) {
         for(j in 1:ord) {
            for(k in 1:ncol(df$ivs)) {
               if(allow_interactions) {
                  for(l in 1:ncol(df$ivs)) {
                     current_terms <- c(current_terms, paste0("I(ivs[,", k, "]>", ivs[current_breakpoints[i],k], ")*(I(ivs[,", l, "]^", j, "))"))
                  }
               } else {
                  current_terms <- c(current_terms, paste0("I(ivs[,", k, "]>", ivs[current_breakpoints[i],k], ")*(I(ivs[,", k, "]^", j, "))"))
               }
            }
         }
      }
   }
   model <- DirichReg(as.formula(paste0("dv ~ ", paste(current_terms, collapse=" + "))), df, model=parameterization, control=list(iterlim=iterlim, tol1=tol1, tol2=tol2))
   return(model)
}

#' Retrieve Largest t-statistic
#'
#' @description This function obtains the largest absolute value among all
#' \emph{t}-statistics for a given breakpoint in a regression model summary
#' prepared by \code{\link[Compositional]{comp.reg}}.
#'
#' @param model The regression model outputted from
#'   \code{\link[Compositional]{comp.reg}}
#' @param ord A non-negative numeric atomic vector indicating the maximum
#'   exponent that was applied to coefficients added to the regression model
#' @param i The numeric index of the breakpoint for which to obtain the
#'   \emph{t}-statistic
#' @param allow_interactions A boolean value indicating whether interactions
#'   were created between independent variables
#' @param ivs A data.frame representing values of the independent variables
#' @return A numeric value representing the largest absolute value among all
#'   \emph{t}-statistics corresponding to the specified breakpoint
#' @export
#' get_t_from_model

get_t_from_model <- function(model, ord, i, allow_interactions, ivs) {
   if(allow_interactions) {   
      return(max(abs(model$be[(ord+2+(ord*ncol(ivs)*(i-1))):(ord+1+(ord*ncol(ivs)*(i))),] / model$seb[(ord+2+(ord*ncol(ivs)*(i-1))):(ord+1+(ord*ncol(ivs)*(i))),])))
   } else {
      return(max(abs(model$be[(1+((ord+1)*i)):((ord+1)+((ord+1)*i)),] / model$seb[(1+((ord+1)*i)):((ord+1)+((ord+1)*i)),])))
   }
}

#' Retrieve Smallest p-value
#'
#' @description This function obtains the smallest \emph{p}-value among all
#' \emph{p}-value for a given breakpoint in a regression model prepared by
#' prepared by \code{\link[Compositional]{comp.reg}}.
#'
#' @param model The regression model outputted from \link{DirichReg}
#' @param ord A non-negative numeric atomic vector indicating the maximum
#'   exponent that was applied to coefficients added to the regression model
#' @param breakpoints A vector of all breakpoints used in \code{model}
#' @param i The numeric index of the breakpoint for which to obtain the
#'   \emph{t}-statistic
#' @param allow_interactions A boolean value indicating whether interactions
#'   were created between independent variables
#' @param ivs A data.frame representing values of the independent variables
#' @param parameterization The specific parameterization used to create
#'   \code{model}
#' @return A numeric value representing the smallest \emph{p}-value
#'   corresponding to the specified breakpoint
#' @export
#' get_p_from_model

get_p_from_model <- function(model, ord, breakpoints, i, allow_interactions, ivs, parameterization) {
   base_indices <- c()
   indices <- c()
   if(parameterization == "common") {
      vars <- length(model$varnames)
   } else if (parameterization == "alternative") {
      vars <- length(model$varnames)-1
   }
   if(allow_interactions) {
      base_indices <- c(base_indices, 1+(i*ncol(ivs))+(1:ncol(ivs)))
      if(ord > 0) {
         for(j in 1:ord) { #Get indices corresponding to this breakpoint for the first DV
            base_indices <- c(base_indices, 1+((length(breakpoints)+1)*ncol(ivs))+((i-1)*(ncol(ivs)^2))+(1:(ncol(ivs)^2)))
         }
      }
      for(j in 1:vars) { #Extend the indices to cover all DVs in the model
         indices <- c(indices, base_indices+(j-1)*(1+ncol(ivs)+(ncol(ivs)*length(breakpoints))+((ncol(ivs)^2)*length(breakpoints)*(ord))))
      }
   } else {
      for(j in 0:ord) { #Get indices corresponding to this breakpoint for the first DV
         base_indices <- c(base_indices, 1+(i*ncol(ivs))+(1:ncol(ivs))+(j*ncol(ivs)*length(breakpoints)))
      }
      for(j in 1:vars) { #Extend the indices to cover all DVs in the model
         indices <- c(indices, base_indices+(j-1)*(1+ncol(ivs)+(ncol(ivs)*length(breakpoints)*(ord+1))))
      }
   }
   return(min(abs(summary(model)$coef.mat[indices,4]), na.rm=T))
}

#' Fitting a Dirichlet Regression
#'
#' @description This is a modified form of the
#'   \code{\link[DirichletReg]{DirichReg}} function, adjusted to avoid "model
#'   frame and formula mismatch in model.matrix()" errors that frequently result
#'   when the formula is long and the original function is called from within
#'   another function.
#'
#' @param formula The formula specifying the regression model
#' @param data The data set from which the model is constructed
#' @param model The parameterization of the model
#' @param subset Not used
#' @param sub.comp Not used
#' @param base Not used
#' @param weights Not used
#' @param control A list of variables to control the convergence process
#' @param verbosity Not used
#' @return A Dirichlet regression model as specified by
#'   \code{\link[DirichletReg]{DirichReg}}
#' @export
#' DirichReg

DirichReg <- function (formula, data, model = c("common", "alternative"), 
    subset, sub.comp, base, weights, control, verbosity = 0) 
{
    this.call <- match.call()
    if (!(verbosity %in% 0:4)) {
        verbosity <- 0L
        warning("invalid value for verbosity.")
    }
    storage.mode(verbosity) <- "integer"
    if (verbosity > 0) {
        cat("- PREPARING DATA\n")
        if (interactive()) 
            utils::flush.console()
    }
    oformula <- formula
    if (missing(control)) {
        control <- list(sv = NULL, iterlim = 10000L, tol1 = .Machine$double.eps^(1/2), 
            tol2 = .Machine$double.eps^(3/4))
    } else {
        if (is.null(control$sv)) 
            control$sv <- NULL
        if (is.null(control$iterlim)) 
            control$iterlim <- 10000L
        if (is.null(control$tol1)) 
            control$tol1 <- .Machine$double.eps^(1/2)
        if (is.null(control$tol2)) 
            control$tol2 <- .Machine$double.eps^(3/4)
    }
    resp_lang <- oformula[[2L]]
    #resp_char <- deparse_nocutoff(resp_lang)
    resp_char <- deparse(resp_lang)
    has_data <- !missing(data)
    Y_in_data <- ifelse(has_data, resp_char %in% names(data), 
        FALSE)
    has_DR_call <- grepl("DR_data", resp_char, fixed = TRUE)
    if (Y_in_data) {
        Y_full <- data[[resp_char]]
    } else if (has_DR_call) {
        Y_full <- eval(resp_lang)
        warning(paste0(strwrap("The response was transformed by DR_data() on the fly. This is not recommended, consider adapting your code.", 
            width = getOption("width") - 9L, exdent = 9L), collapse = "\n"), 
            call. = FALSE, immediate. = TRUE)
        oformula[[2L]] <- as.symbol("Y_full")
    } else {
        Y_full <- get(resp_char, environment(oformula))
    }
    formula <- Formula::as.Formula(oformula)
    if (has_data) {
        assign(resp_char, Y_full)
    } else {
        data[[resp_char]] <- Y_full
    }
    repar <- ifelse(model == "common", FALSE, TRUE)
    mf <- match.call(expand.dots = FALSE)
    if (has_DR_call) {
        mf[["formula"]][[2L]] <- as.symbol("Y_full")
    }
    mf <- mf[c(1L, match(c("formula", "data", "subset", "weights"), 
        names(mf), 0L))]
    mf[["formula"]] <- Formula::as.Formula(formula)
    mf[["drop.unused.levels"]] <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf_formula <- mf
    d <- mf <- eval(mf, parent.frame())
    if ("(weights)" %in% names(mf)) 
        weights <- mf[["(weights)"]]
    else weights <- rep(1, nrow(mf))
    storage.mode(weights) <- "double"
    Y <- stats::model.response(mf, "numeric")
    if (missing(sub.comp)) {
        sub.comp <- seq_len(ncol(Y))
    } else {
        if (length(sub.comp) == ncol(Y)) 
            warning("no subcomposition made, because all variables were selected")
        if (length(sub.comp) == (ncol(Y) - 1)) 
            stop("no subcomposition made, because all variables except one were selected")
        if (any((sub.comp < 1) | (sub.comp > ncol(Y)))) 
            stop("subcompositions must contain indices of variables of the Dirichlet data object")
        y_in <- (seq_len(ncol(Y)))[sub.comp]
        y_out <- (seq_len(ncol(Y)))[-sub.comp]
        y_in_labels <- colnames(Y)[y_in]
        y_out_labels <- paste(colnames(Y)[y_out], sep = "", collapse = " + ")
        Y <- cbind(rowSums(Y[, y_out]), Y[, y_in])
        colnames(Y) <- c(y_out_labels, y_in_labels)
    }
    base <- ifelse(missing(base), attr(Y_full, "base"), base)
    if (!(base %in% seq_len(ncol(Y)))) 
        stop("the base variable lies outside the number of variables")
    n.dim <- ncol(Y)
    if (length(formula)[1] != 1) 
        stop("the left hand side of the model must contain one object prepared by DR_data()")
    if (!repar) {
        if (length(formula)[2] == 1) 
            for (i in 2:ncol(Y)) attr(formula, "rhs") <- lapply(seq_len(ncol(Y)), 
                function(i) attr(formula, "rhs")[[1]])
        if (length(formula)[2] > ncol(Y)) 
            stop("the right hand side must contain specifications for either one or all variables")
    } else {
        if (length(formula)[2] == 1) 
            formula <- Formula::as.Formula(formula(formula), ~1)
        if (length(formula)[2] > 2) 
            stop("the right hand side can only contain one or two specifications in the alternative parametrization")
    }
    if (!repar) {
        X.mats <- lapply(seq_len(ncol(Y)), function(i) {
            stats::model.matrix(terms(formula, data = data, rhs = i), 
                mf)
        })
        Z.mat <- NULL
        n.vars <- unlist(lapply(X.mats, ncol))
    } else {
        X.mats <- lapply(seq_len(ncol(Y)), function(i) stats::model.matrix(terms(formula, 
            data = data, rhs = 1), mf))
        Z.mat <- stats::model.matrix(terms(formula, data = data, rhs = 2), 
            mf)
        n.vars <- c(unlist(lapply(X.mats, ncol))[-1], ncol(Z.mat))
    }
    Y_fit <- unclass(Y)
    attributes(Y_fit) <- NULL
    dim(Y_fit) <- dim(Y)
    storage.mode(Y_fit) <- "double"
    X_fit <- lapply(X.mats, function(this_mat) {
        attr(this_mat, "dimnames") <- NULL
        attr(this_mat, "assign") <- NULL
        return(this_mat)
    })
    for (i in seq_along(X_fit)) storage.mode(X_fit[[i]]) <- "double"
    if (!is.null(Z.mat)) 
        storage.mode(Z.mat) <- "double"
    storage.mode(n.dim) <- "integer"
    storage.mode(n.vars) <- "integer"
    storage.mode(base) <- "integer"
    if (verbosity > 0) {
        cat("- COMPUTING STARTING VALUES\n")
        if (interactive()) 
            utils::flush.console()
    }
    if (is.null(control$sv)) {
        starting.vals <- get_starting_values(Y = Y_fit, X.mats = X_fit, 
            Z.mat = {
                if (repar) 
                  as.matrix(Z.mat)
                else Z.mat
            }, repar = repar, base = base, weights = weights) * 
            if (repar) {
                1
            } else {
                1/n.dim
            }
    } else {
        if (length(control$sv) != n.vars) 
            stop("wrong number of starting values supplied.")
        starting.vals <- control$sv
    }
    parametrization <- ifelse(repar, "alternative", "common")
    if (verbosity > 0) {
        cat("- ESTIMATING PARAMETERS\n")
        if (interactive()) 
            utils::flush.console()
    }
    fit.res <- DirichReg_fit(Y = Y_fit, X = X_fit, Z = as.matrix(Z.mat), 
        sv = starting.vals, d = n.dim, k = n.vars, w = as.vector(weights), 
        ctls = control, repar = repar, base = base, vrb = verbosity)
    varnames <- colnames(Y)
    coefs <- fit.res$estimate
    if (repar) {
        names(coefs) <- unlist(as.vector(c(rep(colnames(X.mats[[1]]), 
            n.dim - 1), colnames(Z.mat))))
    } else {
        names(coefs) <- unlist(lapply(X.mats, colnames))
    }
    if (repar) {
        B <- matrix(0, nrow = n.vars[1L], ncol = n.dim)
        B[cbind(rep(seq_len(n.vars[1L]), (n.dim - 1L)), rep(seq_len(n.dim)[-base], 
            each = n.vars[1]))] <- coefs[1:((n.dim - 1) * n.vars[1])]
        g <- matrix(coefs[((n.dim - 1) * n.vars[1] + 1):length(coefs)], 
            ncol = 1)
        XB <- exp(apply(B, 2L, function(b) {
            as.matrix(X.mats[[1L]]) %*% b
        }))
        MU <- apply(XB, 2L, function(x) {
            x/rowSums(XB)
        })
        PHI <- exp(as.matrix(Z.mat) %*% g)
        ALPHA <- apply(MU, 2L, "*", PHI)
    } else {
        B <- sapply(seq_len(n.dim), function(i) {
            coefs[(cumsum(c(0, n.vars))[i] + 1):cumsum(n.vars)[i]]
        }, simplify = FALSE)
        ALPHA <- sapply(seq_len(n.dim), function(i) {
            exp(as.matrix(X.mats[[i]]) %*% matrix(B[[i]], ncol = 1))
        })
        PHI <- rowSums(ALPHA)
        MU <- apply(ALPHA, 2L, "/", PHI)
    }
    colnames(ALPHA) <- varnames
    colnames(MU) <- varnames
    hessian <- fit.res$hessian
    vcov <- tryCatch(solve(-fit.res$hessian), error = function(x) {
        return(matrix(NA, nrow = nrow(hessian), ncol = ncol(hessian)))
    }, silent = TRUE)
    if (!repar) {
        coefnames <- apply(cbind(rep(varnames, n.vars), unlist(lapply(X.mats, 
            colnames))), 1, paste, collapse = ":")
    } else {
        coefnames <- apply(cbind(rep(c(varnames[-base], "(phi)"), 
            n.vars), c(unlist(lapply(X.mats, colnames)[-base]), 
            colnames(Z.mat))), 1, paste, collapse = ":")
    }
    dimnames(hessian) <- list(coefnames, coefnames)
    dimnames(vcov) <- list(coefnames, coefnames)
    shortnames <- names(coefs)
    names(coefs) <- coefnames
    se <- if (!any(is.na(vcov))) 
        sqrt(diag(vcov))
    else rep(NA, length(coefs))
    res <- structure(list(call = this.call, parametrization = parametrization, 
        varnames = varnames, n.vars = n.vars, dims = length(varnames), 
        Y = Y, X = X.mats, Z = Z.mat, sub.comp = sub.comp, base = base, 
        weights = weights, orig.resp = Y_full, data = data, d = d, 
        formula = formula, mf_formula = mf_formula, npar = length(coefs), 
        coefficients = coefs, coefnames = shortnames, fitted.values = list(mu = MU, 
            phi = PHI, alpha = ALPHA), logLik = fit.res$maximum, 
        vcov = vcov, hessian = hessian, se = se, optimization = list(convergence = fit.res$code, 
            iterations = fit.res$iterations, bfgs.it = fit.res$bfgs.it, 
            message = fit.res$message)), class = "DirichletRegModel")
    for (maxLik_ob in c("lastFuncGrad", "lastFuncParam")) {
        if (exists(maxLik_ob, envir = parent.frame(), inherits = FALSE)) 
            rm(list = maxLik_ob, envir = parent.frame(), inherits = FALSE)
    }
    used_objects <- ls(all.names = TRUE)
    rm(list = c("used_objects", used_objects[used_objects != 
        "res"]))
    on.exit(gc(verbose = FALSE, reset = TRUE), add = TRUE)
    return(res)
}

#' Get Starting Values
#'
#' @description This function is copied from
#' \url{https://rdrr.io/rforge/DirichletReg/src/R/get_starting_values.R}. This
#' was necessary in order to enable \link{DirichReg} to be run, as this function
#' is not an exported namespace from the \code{DirichletReg} package.

get_starting_values <- function(Y, X.mats, Z.mat, repar, base, weights){

  ops <- options(warn = -1L)
  on.exit(options(ops))

  if(!repar){###################################################### COMMON MODEL

    ### collinearity check begin
    exclude_par <- lapply(X.mats, function(list_el){
      lin_coef <- stats::lm.fit(x = list_el, y = stats::runif(nrow(list_el)))[["coefficients"]]
      if(any(na_pos <- is.na(lin_coef))){
        temp_names <- names(lin_coef)
        lin_coef <- seq_along(lin_coef)
        names(lin_coef) <- temp_names
        return(lin_coef[na_pos])
      } else {
        return(NULL)
      }
    })
    ### collinearity check end

    beta.LL <- function(x, y, X, w){
      b <- matrix(x, ncol = 2L)
      if(ncol(X) > 1L){
        LL <- w * stats::dbeta(y, exp(X%*%b[,1]), exp(X%*%b[,2]), log=TRUE)
      } else {
        LL <- w * stats::dbeta(y, unlist(exp(X*x[1])), unlist(exp(X*x[2])), log=TRUE)
      }
      return(LL)
    }

    beta.LL.deriv <- function(x, y, X, w){
      b <- matrix(x, ncol=2L)
      grad <- matrix(0.0, nrow=nrow(X), ncol=prod(dim(b)))
      element <- 1L
      if(ncol(X) > 1L){
        for(cc in seq_len(ncol(b)))for(rr in seq_len(nrow(b))){
          grad[,element] <- w * X[,rr]*(psigamma(exp(X%*%b[,1L]+X%*%b[,2L]))-psigamma(exp(X%*%b[,cc]))+log(y))
          element <- element + 1L
        }
      } else {
        for(cc in seq_len(ncol(b))){
          grad[,element] <- w * X*(psigamma(exp(X%*%x[1L]+X%*%x[2L]))-psigamma(exp(X%*%x[cc]))+log(y))
          element <- element + 1L
        }
      }
      return(grad)
    }

    unidim_fit <- lapply(seq_len(ncol(Y)), function(i){
      if(is.null(exclude_par[[i]])){
        correctX <- X.mats[[i]]
      } else {
        correctX <- X.mats[[i]][ , -exclude_par[[i]], drop = FALSE]
      }
      #suppressWarnings(
        maxLik::maxBFGS(beta.LL, beta.LL.deriv,
          start        = rep(0, 2*ncol(correctX)),
          tol          = 1e-05,
          finalHessian = FALSE,
          X            = correctX,
          y            = Y[,i],
          w            = weights)$estimate[seq_len(ncol(correctX))]
      #)
    })

    for(cmp in seq_len(ncol(Y))){
      if(is.null(exclude_par[[cmp]])){
        break
      } else {
        for(NA_vars in seq_along(exclude_par[[cmp]])) unidim_fit[[cmp]] <- append(unidim_fit[[cmp]], NA, exclude_par[[cmp]][NA_vars] - 1L)
      }
    }

    unidim_fit <- unlist(unidim_fit)

  } else {#################################################### ALTERNATIVE MODEL

    Y_logr <- log(Y[,-base,drop=FALSE]/(Y[,base,drop=TRUE]))
    unidim_fit <- as.numeric(stats::lm(Y_logr ~ X.mats[[1L]] - 1, weights = weights)[["coefficients"]])

    epsilon <- matrix(1.0, nrow = nrow(Y), ncol = ncol(Y))

    start_var <- 1L
    n_mean_par <- ncol(X.mats[[1L]])

    for(i in seq_len(ncol(Y))){
      if(i == base){
        next
      } else {
        epsilon[, i] <- as.numeric(exp(X.mats[[1L]] %*% unidim_fit[seq.int(start_var, start_var + n_mean_par - 1L)]))
        start_var <- start_var + n_mean_par
      }
    }

    MU <- epsilon / rowSums(epsilon)

    log_phi <- stats::optimize(function(x){ sum(weights*DirichletReg::ddirichlet(Y, MU * exp(x), log = TRUE)) }, c(-20, 20), maximum = TRUE)[["maximum"]]
    gammas <- as.numeric(stats::lm(I(rep(log_phi, nrow(Y))) ~ Z.mat - 1, weights = weights)[["coefficients"]])

    unidim_fit <- c(unidim_fit, gammas)

  }

  return(as.numeric(unlist(unidim_fit)))

}

#' DirichReg Fit
#'
#' @description This function is copied from
#' \url{https://rdrr.io/rforge/DirichletReg/src/R/DirichReg_fit.R}. This was
#' necessary in order to enable \link{DirichReg} to be run, as this function is
#' not an exported namespace from the \code{DirichletReg} package.

DirichReg_fit <- function(Y, X, Z, sv, d, k, w, ctls, repar, base, vrb){

  n <- nrow(Y)
  npar <- length(sv)

  ops <- options(warn = -1L)
  on.exit(options(ops))

  ##############################################################################
  ############################################## alternative parametrization ###
  if(repar){
    k <- k[1L]
    beta_ind <- as.integer(matrix(seq_len(d*k), nrow = k, ncol = d)[,-base])
    beta_x_ind <- seq_len((d-1L)*k)
    gamma_ind <- seq.int((d-1L)*k + 1L, npar)
    ncolX <- ncol(X[[1L]])
    ncolZ <- ncol(Z)
    hessian.ind <- rbind(as.matrix(expand.grid(seq_len(k), seq_len(d)[-base])[,c(2L, 1L)]), cbind(-1L, seq_len(ncolZ)))

#    ## first a couple of iterations only for precisions to accelerate
#      bfgs1 <- maxLik::maxBFGS(fn=DReg.repar,
#        start=sv, fixed=1:((d-1)*ncol(X[[1]])),
#        finalHessian=FALSE, iterlim=5L, tol=1e-2, reltol=1e-2,
#        logY=log(Y), X=X[[1L]], ncolX=ncolX, Z=Z, ncolZ=ncolZ, n=n, d=d, k=k, w=w, base=base, npar=npar, bi=beta_ind, bx=beta_x_ind, gi=gamma_ind, NR=FALSE)

      bfgs <- maxLik::maxBFGS(fn=DReg.repar,
        start=sv,#bfgs1$estimate,
        finalHessian=FALSE, iterlim=ctls$iterlim, tol=ctls$tol1, reltol=ctls$tol1, print.level=ifelse(vrb == 0, 0, vrb - 1),
        logY=log(Y), X=X[[1L]], ncolX=ncolX, Z=Z, ncolZ=ncolZ, n=n, d=d, k=k, w=w, base=base, npar=npar, bi=beta_ind, bx=beta_x_ind, gi=gamma_ind, NR=FALSE)

      res <- maxLik::maxNR(fn=DReg.repar,
        start=bfgs$estimate,
        iterlim=ctls$iterlim, tol=ctls$tol2, reltol=ctls$tol2, print.level=ifelse(vrb == 0, 0, vrb - 1),
        logY=log(Y), X=X[[1L]], ncolX=ncolX, Z=Z, ncolZ=ncolZ, n=n, d=d, k=k, w=w, base=base, npar=npar, bi=beta_ind, bx=beta_x_ind, gi=gamma_ind, NR=TRUE, h_dims=hessian.ind[,1L], h_vars=hessian.ind[,2L])

  ##############################################################################
  ################################################### common parametrization ###
  } else {
    seq_along_d <- seq_len(d)
    ncolX <- unlist(lapply(X, ncol))
    beta_x_ind <- lapply(seq_along_d, function(i){ seq.int(cumsum(c(0L, k))[i] + 1L, cumsum(k)[i]) })
    hessian.ind <- cbind(rep(seq_along_d, k), unlist(lapply(k, function(i){ seq_len(i) })))
    
      bfgs <- maxLik::maxBFGS(fn=DReg,
        start=sv,
        finalHessian=FALSE, iterlim=ctls$iterlim, tol=ctls$tol1, reltol=ctls$tol1, print.level=ifelse(vrb == 0, 0, vrb - 1),
        logY=log(Y), X=X, ncolX=ncolX, n=n, d=d, k=k, w=w, npar=npar, seq_along_d=seq_along_d, bx=beta_x_ind, NR=FALSE)

      res <- maxLik::maxNR(fn=DReg,
        start=bfgs$estimate,
        iterlim=ctls$iterlim, tol=ctls$tol2, reltol=ctls$tol2, print.level=ifelse(vrb == 0, 0, vrb - 1),
        logY=log(Y), X=X, ncolX=ncolX, n=n, d=d, k=k, w=w, npar=npar, seq_along_d=seq_along_d, bx=beta_x_ind, NR=TRUE, h_dims=hessian.ind[,1L], h_vars=hessian.ind[,2L])
  }

  res$bfgs.it <- bfgs$iterations

  return(res)

}

#' DirichReg Repar
#'
#' @description This function is copied from
#' \url{https://rdrr.io/cran/DirichletReg/src/R/DR_LL_alt.R}. This was necessary
#' in order to enable \link{DirichReg} to be run, as this function is not an
#' exported namespace from the \code{DirichletReg} package.

DReg.repar <- function(x, logY, X, ncolX, Z, ncolZ, n, d, k, w, base, npar, bi, bx, gi, NR, h_dims, h_vars){
################################################################################
### PREPARATION ################################################################
################################################################################

  B <- matrix(0.0, nrow = k, ncol = d)
  B[bi] <- x[bx]

  eps  <- apply(B, 2L, function(b){ exp(X %*% b) })

  mu <- eps / .rowSums(eps, n, d, FALSE)

  phi <- as.numeric(exp( Z %*% matrix(x[gi], ncol=1L) ))

  A <- mu * phi

  digamma_A <- digamma(A)
  trigamma_A <- trigamma(A)
  digamma_phi <- digamma(phi)



################################################################################
### LOG-LIKELIHOOD & GRADIENT ##################################################
################################################################################

  LL <- .Call(DirichletReg:::wght_LL_grad_alternative, logY, A, mu, phi, digamma_A, digamma_phi, X, ncolX, Z, ncolZ, n, d, base, npar, w)



################################################################################
### HESSIAN ####################################################################
################################################################################

  if(NR){

  trigamma_phi <- trigamma(phi)

  hessian <- matrix(NA_real_, nrow=npar, ncol=npar)

  for(hess.j in seq_len(npar)){
    for(hess.i in seq_len(npar)){
      if(hess.i < hess.j){ next }

      v1 <- h_vars[hess.i]
      v2 <- h_vars[hess.j]

      derv <- h_dims[c(hess.i, hess.j)]
      d1 <- derv[1L]
      d2 <- derv[2L]

      ##########################################################################
      ############################################### BETAs - SAME RESPONSES ###
      if((derv[1L] == derv[2L]) & all(derv != -1L)) {
        derv <- derv[1L]

        hessian[hess.i, hess.j] <- hessian[hess.j, hess.i] <-
          sum(w*(
            X[,v1] * X[,v2] * A[,derv] * (
              (2.0*mu[,derv] - 1.0) * (
                rowSums(mu[,-derv,drop=FALSE]*(logY[,-derv,drop=FALSE]-digamma_A[,-derv,drop=FALSE]))
              - (1.0 - mu[,derv]) * (logY[,derv]-digamma_A[,derv])
              )
              -
              A[,derv]*(
                (1.0 - mu[,derv])^2 * trigamma_A[,derv]
              + rowSums(mu[,-derv,drop=FALSE]^2*trigamma_A[,-derv,drop=FALSE])
              )
            )
          ))
      ##########################################################################
      ########################################## BETAs - DIFFERENT RESPONSES ###
      } else if((derv[1L] != derv[2L]) & all(derv != -1L)) {
        hessian[hess.i, hess.j] <- hessian[hess.j, hess.i] <-
          sum(w*(
            X[,v1] * X[,v2] * mu[,d1] * mu[,d2] * phi * (
              rowSums(
                mu[,-derv,drop=FALSE] * (
                  2.0 * (logY[,-derv,drop=FALSE] - digamma_A[,-derv,drop=FALSE])
                - A[,-derv,drop=FALSE] * trigamma_A[,-derv,drop=FALSE]
                )
              )
            +
              rowSums(
                (2*mu[,derv,drop=FALSE] - 1.0) * (logY[,derv,drop=FALSE] - digamma_A[,derv,drop=FALSE])
              - A[,derv,drop=FALSE] * (mu[,derv,drop=FALSE] - 1.0) * trigamma_A[,derv,drop=FALSE]
              )
            )
          ))
      ##########################################################################
      ######################################################### BETA / GAMMA ###
      } else if(any(derv != -1L) & any(derv == -1L)) {
        derv <- derv[which(derv != -1L)]

        hessian[hess.i, hess.j] <- hessian[hess.j, hess.i] <-
          sum(w*(
            Z[,v1] * X[,v2] * A[,derv] * (
              rowSums(mu[,-derv,drop=FALSE] * (
                digamma_A[,-derv,drop=FALSE] + A[,-derv,drop=FALSE]*trigamma_A[,-derv,drop=FALSE] - logY[,-derv,drop=FALSE]
              ))
            +
              (mu[,derv] - 1.0) * (
                digamma_A[,derv] + A[,derv]*trigamma_A[,derv] - logY[,derv]
              )
            )
          ))
      ##########################################################################
      ############################################################### GAMMAs ###
      } else if(all(derv == -1)){
        hessian[hess.i, hess.j] <- hessian[hess.j, hess.i] <-
          sum(w*(
            Z[,v1] * Z[,v2] * phi * ( digamma_phi + phi * trigamma_phi + rowSums(mu * (logY - digamma_A - A * trigamma_A)) )
          ))
      }
    }
  }

  attr(LL, "hessian") <- hessian

  }

  return(LL)

}

#' DReg
#'
#' @description This function is copied from
#' \url{https://rdrr.io/cran/DirichletReg/src/R/DR_LL_common.R}. This was
#' necessary in order to enable \link{DirichReg} to be run, as this function is
#' not an exported namespace from the \code{DirichletReg} package.

DReg <- function(x, logY, X, ncolX, n, d, k, w, npar, seq_along_d, bx, NR, h_dims, h_vars){

  B <- lapply(bx, function(b_ind){ x[b_ind] })

  A <- matrix(unlist(lapply(seq_along_d, function(i){ exp(X[[i]] %*% B[[i]]) })), nrow = n, ncol = d)
  Aplus <- .rowSums(A, n, d, FALSE)

  digamma_A  <- digamma(A)
  trigamma_A <- trigamma(A)
  digamma_Aplus  <- digamma(Aplus)
  trigamma_Aplus <- trigamma(Aplus)



################################################################################
### LOG-LIKELIHOOD & GRADIENT ##################################################
################################################################################

  LL <- .Call(DirichletReg:::wght_LL_grad_common, logY, A, Aplus, digamma_A, digamma_Aplus, X, ncolX, c(n, d), npar, w)



################################################################################
### HESSIAN ####################################################################
################################################################################

  if(NR){

    hessian <- matrix(NA_real_, nrow = npar, ncol = npar)

    for(hess.j in seq_len(npar)){
      for(hess.i in seq_len(npar)){
        if(hess.i < hess.j) next

        derv <- h_dims[c(hess.i, hess.j)]

        vars <- h_vars[c(hess.i, hess.j)]

        ########################################################################
        ##################################################### SAME RESPONSES ###
        if(derv[1L] == derv[2L]) {
          derv <- derv[1L]

          hessian[hess.i, hess.j] <- hessian[hess.j, hess.i] <-
          sum(w*(
            X[[derv]][,vars[1L]] * X[[derv]][,vars[2L]] * A[,derv] * (
            logY[,derv] + digamma_Aplus - digamma_A[,derv] + A[,derv] * (
              trigamma_Aplus - trigamma_A[,derv]
              )
            )
          ))
        ########################################################################
        ################################################ DIFFERENT RESPONSES ###
        } else {
          hessian[hess.i, hess.j] <- hessian[hess.j, hess.i] <-
          sum(w*(
            X[[derv[1L]]][,vars[1L]]*X[[derv[2L]]][,vars[2L]]*A[,derv[1L]]*A[,derv[2L]]*trigamma_Aplus
          ))
        }
      }
    }

    attr(LL, "hessian") <- hessian

  }

  return(LL)

}

#' Print
#'
#' This method prints the \code{output} element of an object with the
#' \code{stepseg_output} class, such as output obtained from the
#' \code{\link{stepseg}} function.
#'
#' @param x An object with the \code{stepseg_output} class
#' @param ... Other arguments inherited from the generic \code{print} function
#' @method print stepseg_output
#' @exportS3Method print stepseg_output
#' @export
#' print.stepseg_output

print.stepseg_output <- function(x, ...) {
    print(x$output)
}

#' @evalNamespace "S3method(print,stepseg_output)"
