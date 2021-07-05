#' funmediation:  Fit funmediation model
#'
#' Calculate indirect effect of a binary treatment on a scalar
#' response as mediated by a longitudinal functional trajectory
#' (see Baron & Kenny, 1986; Lindquist, 2012; Coffman et al., 2019).
#'
#' @note This function calls the tvem function in the tvem package.
#' It also calls the pfr function in the refund package (see
#' Goldsmith et al., 2011) to perform penalized functional regression.
#' Some suggestions on interpreting the output from penalized functional
#' regression are given by Dziak et al. (2019).
#'
#' @references
#' Baron, R.M., & Kenny, D.A. (1986). The moderator-mediator variable
#' distinction in social psychological research: Conceptual, strategic,
#' and statistical considerations. Journal of Personality & Social
#' Psychology, 51: 1173-1182.
#' @references
#' Coffman, D. L., Dziak, J. J., Li, R., & Piper, M. (2019). Functional
#' regression mediation analysis with application to a smoking
#' cessation intervention. Joint Statistical Meetings presentation,
#' August 2019.
#' @references
#' Dziak, J. J., Coffman, D. L., Reimherr, M., Petrovich, J., Li, R.,
#' Shiffman, S., & Shiyko, M. P. (2019). Scalar-on-function regression
#' for predicting distal outcomes from intensively gathered longitudinal
#' data: interpretability for applied scientists.
#' Statistics Surveys, 13, 150-180. <doi:10.1214/19-SS126>
#' @references
#' Goldsmith, J., Bobb, J., Crainiceanu, C., Caffo, B., & Reich, D.
#' (2011). Penalized functional regression. Journal of Computational
#' and Graphical Statistics, 20(4), 830-851.
#' @references
#' Lindquist, M. A. (2012). Functional Causal Mediation Analysis
#' With an Application to Brain Connectivity. Journal of the American
#' Statistical Association, 107: 1297-1309. <doi:10.1080/01621459.2012.695640>
#'
#' @param data The dataset containing the data to be analyzed, in long
#'  format (one row per observation, multiple per individual).
#' @param treatment  The name of the variable containing the treatment
#'  assignment, assumed to be unidimensional (either binary
#'  or else numeric).  We recommend a binary (dichotomous) treatment with
#'  0 for control and 1 for experimental).  The values of this variable
#'  should be the same for each row for a given subject.  If there are
#'  more than one treatment variables, such as a dummy-coded exposure with
#'  more than two levels, specify them as a formula such as ~x1+x2.
#' @param mediator The name of the mediator variable. The values of
#' this variable can (and should) vary within each subject.
#' @param outcome The name of the outcome variable. The values of this
#' variable should be the same for each row for a given subject.
#' @param id The name of the variable identifying each subject.
#' @param time The name of the time variable.
#' @param tve_covariates_on_mediator The covariates with time-varying-effects,
#' if any, to be included in the model predicting the mediator from
#' the treatment.
#' @param tie_covariates_on_mediator The covariates with time-invariant
#' effects, if any, to be included in the
#' model predicting the mediator from the treatment.
#' @param covariates_on_outcome The covariates, if any, to be included
#' in the model predicting the outcome from the treatment.
#' They are assumed to be subject-level (time-invariant both in value
#' and in effect).
#' @param binary_mediator Whether the mediator should be modeled as
#' dichotomous with a logistic model (TRUE),
#' or numerical with a normal model (FALSE).
#' @param binary_outcome Whether the outcome should be modeled as
#' dichotomous with a logistic model (TRUE),
#' or numerical with a normal model (FALSE).
#' @param interpolate What kind of presmoothing to use in the
#' penalized functional regression -- specifically, whether
#' to interpolate each subject's trajectory on the mediator
#' (TRUE) or fit a spline to each subject's trajectory on the
#' mediator (FALSE).  This will be counted as TRUE if
#' binary_mediator is TRUE because it does
#' not make as much sense to interpolate a binary outcome.
#' @param nboot Number of bootstrap samples for bootstrap significance
#' test of the overall effect. This test is done using the boot
#' function from the boot package by Angelo Canty and Brian Ripley.
#' It differs somewhat from the bootstrap approach used in a similar
#' context by Lindquist (2012).  We recommend using at least 200 bootstrap
#' samples and preferably 500 or more if time permits.
#' @param boot_level One minus the nominal coverage for
#' the bootstrap confidence interval estimates.
#' @param tvem_spline_order Input to be passed on to the tvem function
#' @param tvem_penalty_order Input to be passed on to the tvem function
#' @param tvem_penalize Input to be passed on to the tvem function
#' @param tvem_do_loop Whether to use a loop to select the number of knots
#' with a pseudo-AIC or pseudo-BIC, passed on to the tvem function
#' @param tvem_num_knots If tvem_do_loop is FALSE, then tvem_num_knots
#' is passed on to the tvem function as num_knots, an integer representing
#' the number of interior knots per B-spline. If tvem_do_loop is
#' TRUE then tvem_num_knots is reinterpreted as the highest number of
#' interior knots to try.
#' @param tvem_use_bic This parameter only matters if tvem_do_loop is TRUE.
#' If tvem_do_loop is TRUE
#' and tvem_use_bic is TRUE, then the information criterion used will be
#' a pseudolikelihood version of BIC.
#' If tvem_do_loop is TRUE and tvem_use_bic is FALSE, then the information
#' criterion used will be
#' a pseudolikelihood version of AIC instead. If tvem_do_loop is FALSE then
#' tvem_use_bic is ignored.
#'
#' @return An object of type funmediation. The components of an object of
#' type funmediation are as follows:
#' \describe{
#' \item{original_results}{The estimates from the fitted models for
#' predicting the mediator from the treatment, predicting the outcome
#' from the mediator and treatment, and predicting the outcome from the
#' treatment alone.}
#' \item{bootstrap_results}{The estimate and confidence interval of the
#' indirect effect using a bootstrap approach.}
#' }
#'
#' The original_results component has these components within it:
#' \describe{
#' \item{time_grid}{Grid of time points on which the functional
#' coefficients are estimated.}
#' \item{alpha_int_estimate}{Estimated intercept function (as a vector
#'  of estimates) from the TVEM regression of the mediator, M, on treatment, X.}
#' \item{alpha_int_se}{Estimated pointwise standard errors associated
#' with the above.}
#' \item{alpha_X_estimate}{Estimated time-varying treatment effect
#' from the TVEM regression of the mediator, M, on the treatment, X.}
#' \item{alpha_X_se}{Estimated pointwise standard errors associated
#' with the above.}
#' \item{beta_int_estimate}{Estimated scalar intercept from the scalar-on-
#' function regression of the outcome, Y, on the mediator, M, and treatment, X.}
#' \item{beta_int_se}{Estimated standard error for the above.}
#' \item{beta_X_estimate}{Estimated scalar coefficient for the treatment, X,
#' from the scalar-on-function regression of the outcome, Y, on the mediator,
#' M, and treatment, X.}
#' \item{beta_X_se}{Estimated standard error for the above.}
#' \item{beta_M_estimate}{Estimated functional coefficient for the mediator,
#'  M, from the scalar-on-function regression of the outcome, Y, on the mediator,
#'  M, and treatment, X.}
#' \item{beta_M_se}{Estimated pointwise standard errors associated with the above}
#' \item{beta_M_pvalue}{The p-value for significance of the mediator, M, in
#' predicting outcome, Y, after adjusting for treatment, X.}
#' \item{tau_int_estimate}{Intercept from simple model predicting outcome, Y,
#' directly from treatment, X.}
#' \item{tau_int_se}{Estimated standard error for the above.}
#' \item{tau_X_estimate}{Coefficient for treatment in model predicting outcome,
#' Y, directly from treatment, X.}
#' \item{tau_X_se}{Estimated standard error for the above.}
#' \item{indirect_effect_estimate}{Estimated indirect effect, calculated as
#' the dot product of the effect of treatment on mediator and the treatment-
#' adjusted effect of mediator on outcome.  It is a scalar, even though the
#' two component effects are functions of time.}
#' \item{tvem_XM_details}{Detailed output from the tvem function for the time-
#' varying-effect model predicting the mediator, M, from the treatment, X.}
#' \item{funreg_MY_details}{Detailed output from the refund::pfr function for
#' the scalar-on-function functional regression predicting the outcome, Y, from
#' the treatment, X, and mediator, M.}
#' \item{total_effect_details}{Detailed output from the linear or generalized
#' linear model predicting the outcome from the treatment alone, ignoring the
#' mediator (i.e., total effect)}
#' }
#'
#' The bootstrap_results component has these components within it:
#' \describe{
#' \item{indirect_effect_boot_estimate}{Bootstrap point estimate of the
#' indirect effect (average of bootstrap sample estimates).}
#' \item{indirect_effect_boot_se}{Bootstrap standard error for the
#' indirect effect (standard deviation of bootstrap sample estimates).}
#' \item{indirect_effect_boot_norm_lower}{Lower end of the bootstrap
#' confidence interval using the normal method in boot.ci in the boot package.}
#' \item{indirect_effect_boot_norm_upper}{Upper end of the bootstrap
#' confidence interval using the normal method.}
#' \item{indirect_effect_boot_basic_lower}{Lower end of the bootstrap
#' confidence interval using the basic method in boot.ci in the boot package.}
#' \item{indirect_effect_boot_basic_upper}{Upper end of the bootstrap
#' confidence interval using the basic method.}
#' \item{indirect_effect_boot_perc_lower}{Lower end of the bootstrap
#' confidence interval using the percentile method in boot.ci in the boot package.}
#' \item{indirect_effect_boot_perc_upper}{Upper end of the bootstrap
#' confidence interval using the percentile method.}
#' \item{boot_level}{The alpha level used for the bootstrap confidence interval.}
#' \item{boot1}{The output returned from the boot function.}
#' \item{time.required}{The amount of time spent doing the bootstrap test,
#' including generating and analyzing all samples.}
#' }
#'
#' @importFrom tvem tvem
#' @importFrom boot boot boot.ci
#' @importFrom refund pfr
#'
#' @importFrom stats as.formula binomial coef gaussian
#'  glm rbinom rnorm sd terms update var
#' @export

funmediation <- function(data,
                         treatment,
                         mediator,
                         outcome,
                         id,
                         time,
                         tve_covariates_on_mediator=NULL,
                         tie_covariates_on_mediator=NULL,
                         covariates_on_outcome=NULL,
                         interpolate=TRUE,
                         tvem_penalize=TRUE,
                         tvem_penalty_order=1,
                         tvem_spline_order=3,
                         tvem_num_knots=3,
                         tvem_do_loop=FALSE,
                         tvem_use_bic=FALSE,
                         binary_mediator=FALSE, # FALSE for numerical mediator, TRUE for dichotomous 0/1;
                         binary_outcome=FALSE, # FALSE for numerical outcome, TRUE for dichotomous 0/1;
                         nboot=200,
                         boot_level=.05) {
  #-------------------------------------------;
  #--- PROCESSING OF INPUT -------------------;
  #-------------------------------------------;
  m <- match.call(expand.dots = FALSE);
  m$treatment <- NULL;
  m$mediator <- NULL;
  m$outcome <- NULL;
  m$tve_covariates_on_mediator <- NULL;
  m$tie_covariates_on_mediator <- NULL;
  m$covariates_on_outcome <- NULL;
  m$binary_mediator <- NULL;
  m$binary_outcome <- NULL;
  m$tvem_num_knots <- NULL;
  m$tvem_penalty_order <- NULL;
  m$interpolate <- NULL;
  m$tvem_spline_order <- NULL;
  m$tvem_penalize <- NULL;
  m$tvem_do_loop <- NULL;
  m$grid <- NULL;
  m$nboot <- NULL;
  if (is.matrix(eval.parent(m$data))) {
    m$data <- as.data.frame(data);
  }
  m[[1]] <- quote(stats::model.frame);
  m <- eval.parent(m);
  id_variable_name <- as.character(substitute(id));
  time_variable_name <- as.character(substitute(time));
  if(class(substitute(treatment))=="call") {
    # Exposure(s) were specified as a formula with one or more variables.
    treatment_variable_names <- attr(terms(as.formula(treatment)),"term.labels");
  } else {
    if(class(substitute(treatment))=="name") {
      # Exposure was specified as a single variable.
      treatment_variable_names <- as.character(substitute(treatment));
    } else {
      stop(paste("Please specify the_predictors as a name or formula",
                 "i.e., in the form x1 or ~x1+x2."));
    }
  }
  mediator_variable_name <- as.character(substitute(mediator));
  outcome_variable_name <- as.character(substitute(outcome));
  long_data_for_analysis <- as.data.frame(eval.parent(data));
  id_variable <- long_data_for_analysis[,id_variable_name];
  if (min(table(id_variable))<2) {
    message <- "At least one subject has less than two measurement occasions.";
    if (interpolate==TRUE) {
      message <- paste(message, "We suggest using interpolate=FALSE.");
    }
    warning(message);
  }
  time_variable <- long_data_for_analysis[,time_variable_name];
  treatment_variables <- long_data_for_analysis[,treatment_variable_names,drop=FALSE];
  num_treatment_variables <- length(treatment_variable_names);
  mediator_variable <- long_data_for_analysis[,mediator_variable_name];
  outcome_variable <- long_data_for_analysis[,outcome_variable_name];
  if (sum(is.na(id_variable>0))) {
    stop("Please remove data rows which have missing data on the id variable.");
  }
  #-------------------------------------------;
  # ----- Process covariates -----------------;
  if (is.null(covariates_on_outcome)) {
    covariates_on_outcome_names <- NULL;
  } else {
    covariates_on_outcome_names <- attr(terms(covariates_on_outcome),"term.labels");
  }
  num_covariates_on_outcome <- length(covariates_on_outcome_names);
  covariates_on_outcome_data <- long_data_for_analysis[,covariates_on_outcome_names,drop=FALSE];
  if (is.null(tve_covariates_on_mediator)) {
    tve_covariates_on_mediator_names <- NULL;
  } else {
    tve_covariates_on_mediator_names <- attr(terms(tve_covariates_on_mediator),"term.labels");
  }
  num_tve_covariates_on_mediator <- length(tve_covariates_on_mediator_names);
  tve_covariates_on_mediator_data <- long_data_for_analysis[,tve_covariates_on_mediator_names,drop=FALSE];
  if (is.null(tie_covariates_on_mediator)) {
    tie_covariates_on_mediator_names <- NULL;
  } else {
    tie_covariates_on_mediator_names <- attr(terms(tie_covariates_on_mediator),"term.labels");
  }
  num_tie_covariates_on_mediator <- length(tie_covariates_on_mediator_names);
  tie_covariates_on_mediator_data <- long_data_for_analysis[,tie_covariates_on_mediator_names,drop=FALSE];
  if (length(intersect(tve_covariates_on_mediator_names, tie_covariates_on_mediator_names))>0) {
    stop(paste("Please make sure that no covariate is specified as having",
               "both time-varying and time-invariant effects."));
  }
  #-------------------------------------------;
  #--- CONVERT MEDIATOR M TO WIDE FORM -------;
  #-------------------------------------------;
  for (this_treatment_variable_index in 1:num_treatment_variables) {
    this_treatment_variable <- treatment_variables[,this_treatment_variable_index];
    variance_check_vector_1 <- unlist(lapply(split(this_treatment_variable,f=id_variable),var));
    if (max(variance_check_vector_1, na.rm=TRUE)>1e-10) {
      stop("Please make sure that the subject-level treatment is constant within subject.")
    }
  }
  variance_check_vector_2 <- unlist(lapply(split(outcome_variable,f=id_variable),var));
  if (max(variance_check_vector_2, na.rm=TRUE)>1e-10) {
    stop("Please make sure that the subject-level outcome is constant within subject.")
  }
  observed_time_grid <- sort(unique(time_variable));
  wide_id <- sort(unique(id_variable[which(!is.na(id_variable))]));
  # Should the user provide the data as long or wide?
  nobs <- length(observed_time_grid);
  nsub <- length(wide_id);
  MEDIATOR <- matrix(NA,nsub,nobs);
  OUTCOME <- rep(NA,nsub);
  if (num_covariates_on_outcome>0) {
    wide_covariates_on_outcome <- matrix(NA,nsub,num_covariates_on_outcome);
  }
  if (num_tve_covariates_on_mediator>0) {
    wide_tve_covariates_on_mediator <- matrix(NA,nsub,num_tve_covariates_on_mediator);
  }
  if (num_tie_covariates_on_mediator>0) {
    wide_tie_covariates_on_mediator <- matrix(NA,nsub,num_tie_covariates_on_mediator);
  }
  stopifnot(nsub>0);
  if (!is.null(covariates_on_outcome_data)) {
    stopifnot(length(id_variable)==nrow(covariates_on_outcome_data));
  }
  if (!is.null(tie_covariates_on_mediator_data)) {
    stopifnot(length(id_variable)==nrow(tie_covariates_on_mediator_data));
  }
  if (!is.null(tve_covariates_on_mediator_data)) {
    stopifnot(length(id_variable)==nrow(tve_covariates_on_mediator_data));
  }
  temp_treatment_variables_matrix <- matrix(NA,nsub,num_treatment_variables);
  for (this_id in 1:length(wide_id)) {
    these_rows <- which(id_variable==wide_id[this_id]);
    if (length(these_rows)>0) {
      for (j in 1:num_treatment_variables) {
        treatment_check_vector <- treatment_variables[these_rows,j,drop=FALSE];
        if (sum(!is.na(treatment_check_vector))>0) {
          if (min(treatment_check_vector, na.rm=TRUE)!=
              max(treatment_check_vector, na.rm=TRUE)) {
            print(treatment_check_vector);
            stop("Please make sure the treatment is the same for each observation within subject.")
          }
        }
        temp_treatment_variables_matrix[this_id,j] <- as.numeric(treatment_variables[min(these_rows),j,drop=FALSE]);
      }
      if (sum(!is.na(outcome_variable[these_rows]))>0) {
        if (min(outcome_variable[these_rows], na.rm=TRUE)!=
            max(outcome_variable[these_rows], na.rm=TRUE)) {
          stop("Please make sure the outcome is the same for each observation within subject.")
        }
      }
      OUTCOME[this_id] <- outcome_variable[min(these_rows)];
      if (num_covariates_on_outcome>0) {
        if (max(apply(covariates_on_outcome_data[these_rows,,drop=FALSE],2,var, na.rm = TRUE), na.rm = TRUE)>0) {
          print(paste("Possible problem in covariates_on_outcome_data for participant",wide_id[this_id]));
          print(covariates_on_outcome_data[these_rows,,drop=FALSE]);
          stop("Please make sure that the covariates on the outcome do not vary within subject for this model.")
        }
      }
      if (num_tve_covariates_on_mediator>0) {
        if (max(apply(tve_covariates_on_mediator_data[these_rows,,drop=FALSE],2,var, na.rm = TRUE), na.rm = TRUE)>0) {
          print(paste("Possible problem in tve_covariates_on_mediator_data for participant",wide_id[this_id]));
          print(tve_covariates_on_mediator_data[these_rows,,drop=FALSE]);
          stop("Please make sure that the covariates on the mediator do not vary within subject for this model.")
        }
      }
      if (num_tie_covariates_on_mediator>0) {
        if (max(apply(tie_covariates_on_mediator_data[these_rows,,drop=FALSE],2,var, na.rm = TRUE), na.rm = TRUE)>0) {
          print(paste("Possible problem in tie_covariates_on_mediator_data for participant",wide_id[this_id]));
          print(tie_covariates_on_mediator_data[these_rows,,drop=FALSE]);
          stop("Please make sure that the covariates on the mediator do not vary within subject for this model.")
        }
      }
      stopifnot(length(observed_time_grid)>0);
      for (this_time in 1:length(observed_time_grid)) {
        this_data <- which(id_variable==wide_id[this_id] &
                             time_variable==observed_time_grid[this_time]);
        if (length(this_data)==1) {
          MEDIATOR[this_id,this_time] <- data[this_data,mediator_variable_name];
          if (num_covariates_on_outcome>0) {
            wide_covariates_on_outcome[this_id,] <- as.matrix(covariates_on_outcome_data[this_data,,drop=FALSE]);
          }
          if (num_tve_covariates_on_mediator>0) {
            wide_tve_covariates_on_mediator[this_id,] <- as.matrix(tve_covariates_on_mediator_data[this_data,,drop=FALSE]);
          }
          if (num_tie_covariates_on_mediator>0) {
            wide_tie_covariates_on_mediator[this_id,] <- as.matrix(tie_covariates_on_mediator_data[this_data,,drop=FALSE]);
          }
        }
        if (length(this_data)==2) {
          stop("There seems to be more than one measurement on the same subject and time.")
        }
      }
    }
  }
  wide_data <- data.frame(wide_id=wide_id);
  for (j in 1:num_treatment_variables) {
    assign(paste("TREATMENT",j,sep=""),
           temp_treatment_variables_matrix);
    wide_data <- cbind(wide_data,
                       temp_treatment_variables_matrix[,j]);
    colnames(wide_data)[ncol(wide_data)] <- paste("TREATMENT",j,sep="");
  }
  wide_data <- cbind(wide_data,
                     OUTCOME=OUTCOME,
                     MEDIATOR=MEDIATOR);
  wide_id_column <- 1;
  treatment_columns <- 2:(1+num_treatment_variables);
  outcome_column <- 2+num_treatment_variables;
  mediator_columns <- (3+num_treatment_variables):ncol(wide_data);
  if (num_covariates_on_outcome>0) {
    wide_data <- data.frame(cbind(wide_data,
                                  wide_covariates_on_outcome));
    wide_covariates_on_outcome_columns <- (ncol(wide_data)-num_covariates_on_outcome+1):ncol(wide_data);
  }
  if (num_tve_covariates_on_mediator>0) {
    wide_data <- cbind(wide_data,
                       wide_tve_covariates_on_mediator);
    wide_tve_covariates_on_mediator_columns <- (ncol(wide_data)-num_tve_covariates_on_mediator+1):ncol(wide_data);
  }
  if (num_tie_covariates_on_mediator>0) {
    wide_data <- data.frame(cbind(wide_data,
                                  wide_tie_covariates_on_mediator))
    wide_tie_covariates_on_mediator_columns <- (ncol(wide_data)-num_tie_covariates_on_mediator+1):ncol(wide_data);
  }
  #-------------------------------------------;
  #--- DEFINE MAIN ANALYSIS ------------------;
  #-------------------------------------------;
  analyze_data_for_mediation <- function(wide_data,
                                         indices,
                                         get_details=FALSE) {
    local_wide_data <- wide_data[indices,];
    usable <- which((apply(!is.na(local_wide_data[,mediator_columns]),1,sum)>=1) & 
                      !is.na(local_wide_data[,outcome_column])); 
    local_wide_data <- local_wide_data[usable,];
    #--- Take data frame apart into pieces to use with pfr function
    MEDIATOR <- as.matrix(local_wide_data[,mediator_columns]);
    wide_id <- unlist(local_wide_data[,wide_id_column]);
    if (length(treatment_columns)>1) {
      for (j in 1:length(treatment_columns)) {
        assign(paste("TREATMENT",j,sep=""), unlist(local_wide_data[,treatment_columns[j]]));
      }
    } else {
      TREATMENT <- unlist(local_wide_data[,treatment_columns]);
    }
    OUTCOME <- unlist(local_wide_data[,outcome_column]);
    if (num_covariates_on_outcome>0) {
      wide_covariates_on_outcome <- local_wide_data[,wide_covariates_on_outcome_columns,drop=FALSE];
    }
    if (num_tve_covariates_on_mediator>0) {
      wide_tve_covariates_on_mediator <- local_wide_data[,wide_tve_covariates_on_mediator_columns,drop=FALSE];
    }
    if (num_tie_covariates_on_mediator>0) {
      wide_tie_covariates_on_mediator <- local_wide_data[,wide_tie_covariates_on_mediator_columns,drop=FALSE];
    }
    nobs <- length(mediator_columns);
    nsub <- length(wide_id);
    #--- EFFECT OF MEDIATOR M AND TREATMENT X ON OUTCOME Y ---;
    if (binary_mediator==TRUE | interpolate==FALSE) {
      pfr_formula <- OUTCOME ~ lf(MEDIATOR,
                                  argvals=observed_time_grid,
                                  presmooth="bspline",
                                  presmooth.opts=list(nbasis=4));
    } else {
      pfr_formula <- OUTCOME ~ lf(MEDIATOR,
                                  argvals=observed_time_grid,
                                  presmooth="interpolate");
    }
    if (length(treatment_columns)==1) {
      new_one <- as.formula("~.+TREATMENT");
      pfr_formula <- update(pfr_formula,new_one);
    } else {
      for (j in 1:length(treatment_columns)) {
        new_one <- as.formula(paste("~.+TREATMENT",j,sep=""));
        pfr_formula <- update(pfr_formula,new_one);
      }
    }
    if (num_covariates_on_outcome>0) {
      for (this_one in 1:num_covariates_on_outcome) {
        assign(covariates_on_outcome_names[this_one],
               wide_covariates_on_outcome[,this_one]);
        new_one <- as.formula(paste("~.+",covariates_on_outcome_names[this_one],sep=""));
        pfr_formula <- update(pfr_formula,new_one);
      }
    }
    if (binary_outcome) {
      funreg_MY <- try(refund::pfr(pfr_formula,
                                   scale=1,
                                   family=binomial()));
      # I wish I could let the user send the data and family in from outside the function,
      # but the pfr function does not allow this due to its unusual implementation
      # as a wrap-around for a hidden call to gam.;
    } else {
      funreg_MY <- try(refund::pfr(pfr_formula,
                                   family=gaussian()));
    }
    if (any(class(funreg_MY)=="try-error")) {
      print(pfr_formula);
      print(funreg_MY);
      stop("Error in running pfr.")
    }
    beta_int_estimate <- as.numeric(funreg_MY$coefficient["(Intercept)"]);
    beta_int_se <- as.numeric(summary(funreg_MY)$se["(Intercept)"]);
    if (length(treatment_columns)==1) {
      beta_X_estimate <-  as.numeric(funreg_MY$coefficient["TREATMENT"]);
      beta_X_se <- as.numeric(summary(funreg_MY)$se["TREATMENT"]);
    } else {
      beta_X_estimate <- rep(NA,length(treatment_columns));
      beta_X_se <- rep(NA,length(treatment_columns));
      for (j in 1:length(treatment_columns)) {
        beta_X_estimate[j] <-  as.numeric(funreg_MY$coefficient[paste("TREATMENT",j,sep="")]);
        beta_X_se[j] <-  as.numeric(summary(funreg_MY)$se[paste("TREATMENT",j,sep="")]);
      }
    }
    temp_coefs <- coef(funreg_MY, coords=list(observed_time_grid));
    beta_M_estimate <- as.numeric(temp_coefs[,"value"]);
    beta_M_se <- as.numeric(temp_coefs[,"se"]);
    time_grid_for_fitting <- observed_time_grid;
    beta_M_pvalue <- summary(funreg_MY)$s.table[1,"p-value"];
    if (binary_mediator) {
      tvem_family <- binomial(link="logit");
    } else {
      tvem_family <- gaussian();
    }
    #--- TOTAL EFFECT OF TREATMENT X ON OUTCOME Y ---;
    if (length(treatment_columns)==1) {
      glm_formula <- OUTCOME ~ TREATMENT;
    } else {
      temp_string <- "OUTCOME ~ ";
      for (j in 1:length(treatment_columns)) {
        temp_string <- paste(temp_string,
                             " TREATMENT",
                             j,
                             ifelse(j<length(treatment_columns),"+",""),
                             sep="")
      }
      glm_formula <- as.formula(temp_string);
    }
    if (num_covariates_on_outcome>0) {
      for (this_one in 1:num_covariates_on_outcome) {
        new_one <- as.formula(paste("~.+",covariates_on_outcome_names[this_one],sep=""));
        glm_formula <- update(glm_formula,new_one);
      }
    }
    if (binary_outcome) {
      model_for_total_effect_XY <- glm(glm_formula,
                                       family=binomial);
    } else {
      model_for_total_effect_XY <- glm(glm_formula);
    }
    tau_int_estimate <- as.numeric(model_for_total_effect_XY$coefficients["(Intercept)"]);
    tau_int_se <- summary(model_for_total_effect_XY)$coefficients["(Intercept)","Std. Error"];
    if (length(treatment_columns)==1) {
      tau_X_estimate <- as.numeric(model_for_total_effect_XY$coefficients["TREATMENT"]);
      tau_X_se <- summary(model_for_total_effect_XY)$coefficients["TREATMENT","Std. Error"];
      if (binary_outcome) {
        tau_X_pvalue <- summary(model_for_total_effect_XY)$coefficients["TREATMENT","Pr(>|z|)"];
      } else {
        tau_X_pvalue <- summary(model_for_total_effect_XY)$coefficients["TREATMENT","Pr(>|t|)"];
      }
    } else {
      tau_X_estimate <- rep(NA,length(treatment_columns));
      tau_X_se <- rep(NA,length(treatment_columns));
      tau_X_pvalue <- rep(NA,length(treatment_columns));
      for (j in 1:length(treatment_columns)) {
        this_treatment_variable_name <- paste("TREATMENT",j,sep="");
        tau_X_estimate[j] <- as.numeric(model_for_total_effect_XY$coefficients[this_treatment_variable_name]);
        tau_X_se[j] <- summary(model_for_total_effect_XY)$coefficients[this_treatment_variable_name,"Std. Error"];
        if (binary_outcome) {
          tau_X_pvalue[j] <- summary(model_for_total_effect_XY)$coefficients[this_treatment_variable_name,"Pr(>|z|)"];
        } else {
          tau_X_pvalue[j] <- summary(model_for_total_effect_XY)$coefficients[this_treatment_variable_name,"Pr(>|t|)"];
        }
      }
    }
    #--- EFFECT OF TREATMENT X ON MEDIATOR M ---;
    if (length(treatment_columns)==1) {
      local_long_data <- data.frame(id=rep(wide_id, each=nobs),
                                    time=rep(observed_time_grid, times=nsub),
                                    outcome=rep(OUTCOME, each=nobs),
                                    TREATMENT=rep(TREATMENT, each=nobs),
                                    MEDIATOR=as.vector(t(MEDIATOR)));
      tvem_formula1 <- MEDIATOR ~ TREATMENT;
    } else {
      local_long_data <- data.frame(id=rep(wide_id, each=nobs),
                                    time=rep(observed_time_grid, times=nsub),
                                    outcome=rep(OUTCOME, each=nobs));
      temp_string <- "MEDIATOR ~ ";
      for (j in 1:length(treatment_columns)) {
        temp <- rep(unlist(local_wide_data[,treatment_columns[j]]),each=nobs);
        local_long_data <- cbind(local_long_data,
                                 temp);
        colnames(local_long_data)[ncol(local_long_data)] <- paste("TREATMENT",j,sep="");
        temp_string <- paste(temp_string,
                             ifelse(j>1,"+",""),
                             paste("TREATMENT",j,sep=""));
      }
      tvem_formula1 <- as.formula(temp_string);
      local_long_data <- cbind(local_long_data,
                               MEDIATOR=as.vector(t(MEDIATOR)));
    }
    tvem_formula2 <- tie_covariates_on_mediator;
    if (num_tve_covariates_on_mediator>0) {
      for (this_one in 1:num_tve_covariates_on_mediator) {
        new_name <- tve_covariates_on_mediator_names[this_one];
        temp_data_frame <- data.frame(temp=rep(wide_tve_covariates_on_mediator[,this_one],each=nobs));
        colnames(temp_data_frame) <- new_name;
        local_long_data <- cbind(local_long_data, temp_data_frame);
        new_term <- as.formula(paste("~.+",new_name,sep=""));
        tvem_formula1 <- update(tvem_formula1,new_term);
      }
    }
    if (num_tie_covariates_on_mediator>0) {
      for (this_one in 1:num_tie_covariates_on_mediator) {
        new_name <- tie_covariates_on_mediator_names[this_one];
        temp_data_frame <- data.frame(temp=rep(wide_tie_covariates_on_mediator[,this_one],each=nobs));
        colnames(temp_data_frame) <- new_name;
        local_long_data <- cbind(local_long_data, temp_data_frame);
      }
    }
    local_long_data <- local_long_data[which(!is.na(local_long_data$MEDIATOR)),];
    # listwise deletion to remove empty observations;



    if (min(table(local_long_data$id))<2) {
      message <- "At least one subject has less than two non-missing measurement occasions.";
      if (interpolate==TRUE) {
        message <- paste(message, "We suggest using interpolate=FALSE.");
      }
      warning(message);
    }


    if (tvem_do_loop) {
      tvem_results_list <- list();
      max_knots <- tvem_num_knots;
      IC_values <- rep(Inf,max_knots+1);
      for (this_num_knots in 0:max_knots) {
        tvem_XM <- suppressWarnings(tvem::tvem(data=local_long_data,
                                               formula=tvem_formula1,
                                               time=time,
                                               id=id,
                                               invar_effects=tvem_formula2,
                                               spline_order=tvem_spline_order,
                                               family=tvem_family,
                                               penalty_function_order=tvem_penalty_order,
                                               penalize=tvem_penalize,
                                               num_knots=this_num_knots,
                                               grid=time_grid_for_fitting));
        IC_values[1+this_num_knots] <- ifelse(tvem_use_bic,
                                              tvem_XM$model_information$pseudo_bic,
                                              tvem_XM$model_information$pseudo_aic);
        tvem_results_list[[1+this_num_knots]] <- tvem_XM;
      }
      IC_table <- data.frame(0:max_knots, IC_values);
      colnames(IC_table) <- c("Number_Of_Interior_Knots",ifelse(tvem_use_bic,
                                                                "Pseudo_BIC",
                                                                "Pseudo_AIC"));
      tvem_XM <- tvem_results_list[[which.min(IC_values)]];
    } else {
      if (get_details) {
        tvem_XM <- tvem::tvem(data=local_long_data,
                              formula=tvem_formula1,
                              time=time,
                              id=id,
                              invar_effects=tvem_formula2,
                              spline_order=tvem_spline_order,
                              family=tvem_family,
                              penalty_function_order=tvem_penalty_order,
                              penalize=tvem_penalize,
                              num_knots=tvem_num_knots,
                              grid=time_grid_for_fitting);
      } else {
        tvem_XM <- suppressWarnings(tvem::tvem(data=local_long_data,
                                               formula=tvem_formula1,
                                               time=time,
                                               id=id,
                                               invar_effects=tvem_formula2,
                                               spline_order=tvem_spline_order,
                                               penalty_function_order=tvem_penalty_order,
                                               penalize=tvem_penalize,
                                               num_knots=tvem_num_knots,
                                               grid=time_grid_for_fitting));
      }
    }
    #--- MEDIATED EFFECT OF TREATMENT X THROUGH MEDIATOR M ON OUTCOME Y ---;
    alpha_int_estimate <- tvem_XM$grid_fitted_coefficients[[1]]$estimate;
    alpha_int_se <- tvem_XM$grid_fitted_coefficients[[1]]$standard_error;
    if (length(treatment_columns)==1) {
      alpha_X_estimate <- tvem_XM$grid_fitted_coefficients[[2]]$estimate;
      alpha_X_se <- tvem_XM$grid_fitted_coefficients[[2]]$standard_error;
      if (length(alpha_X_estimate)!=length(beta_M_estimate)) {
        stop("Dimension error in functional mediation function;")
      }
      indirect_effect_estimate <- mean(alpha_X_estimate*beta_M_estimate);
    } else {
      alpha_X_estimate <- list();
      alpha_X_se <- list();
      indirect_effect_estimate <- rep(NA,length(treatment_columns));
      for (j in 1:length(treatment_columns)) {
        alpha_X_estimate[[j]] <- tvem_XM$grid_fitted_coefficients[[1+j]]$estimate;
        alpha_X_se[[j]] <- tvem_XM$grid_fitted_coefficients[[1+j]]$standard_error;
        indirect_effect_estimate[j] <- mean(alpha_X_estimate[[j]]*beta_M_estimate[[j]]);
      }
    }
    if (get_details) {
      answer_list <- list(time_grid=time_grid_for_fitting,
                          alpha_int_estimate=alpha_int_estimate,
                          alpha_int_se=alpha_int_se,
                          alpha_X_estimate=alpha_X_estimate,
                          alpha_X_se=alpha_X_se,
                          beta_int_estimate=beta_int_estimate,
                          beta_int_se=beta_int_se,
                          beta_X_estimate=beta_X_estimate,
                          beta_X_se=beta_X_se,
                          beta_M_estimate=beta_M_estimate,
                          beta_M_se=beta_M_se,
                          beta_M_pvalue=beta_M_pvalue,
                          tau_int_estimate=tau_int_estimate,
                          tau_int_se=tau_int_se,
                          tau_X_estimate=tau_X_estimate,
                          tau_X_se=tau_X_se,
                          tau_X_pvalue=tau_X_pvalue,
                          indirect_effect_estimate=indirect_effect_estimate,
                          tvem_XM_details=tvem_XM,
                          funreg_MY_details=funreg_MY,
                          total_effect_details=model_for_total_effect_XY,
                          binary_mediator=binary_mediator,
                          binary_outcome=binary_outcome);
      if (tvem_do_loop) {
        answer_list$tvem_IC_table <- IC_table;
      }
      return(answer_list);
    } else {
      if (this_bootstrap>0) {cat(this_bootstrap);cat(".");}
      this_bootstrap <<- this_bootstrap + 1;
      if (tvem_do_loop) {
        ICs_table_from_bootstraps <<- rbind(ICs_table_from_bootstraps,IC_table[,2]);
      }
      return(indirect_effect_estimate);
    }
  }
  #------------------------------------------;
  #---------------- Call function -----------;
  #------------------------------------------;
  ############debug(analyze_data_for_mediation);
  original_results <- analyze_data_for_mediation(wide_data,
                                                 indices=1:nrow(wide_data),
                                                 get_details=TRUE);
  before_boot <- Sys.time();
  ICs_table_from_bootstraps <- NULL;
  cat("Ran original results.\n");
  cat("Working on bootstrap results:\n");
  this_bootstrap <- 0;
  boot1 <- boot::boot(data=wide_data,
                      statistic=analyze_data_for_mediation,
                      R=nboot);
  cat("Done bootstrapping.\n");
  if (nboot<50) {
    warning("The number of bootstrap samples was very small and results will be dubious.");
  }
  boot_norm <- list();
  boot_basic <- list();
  boot_perc <- list();
  indirect_effect_boot_estimate <- rep(NA, num_treatment_variables);
  indirect_effect_boot_se <- rep(NA, num_treatment_variables);
  for (j in 1:num_treatment_variables) {
    boot_norm[[j]] <- boot::boot.ci(boot1,conf=1-boot_level,index=j,type="norm");
    boot_basic[[j]] <- boot::boot.ci(boot1,conf=1-boot_level,index=j,type="basic");
    boot_perc[[j]] <- boot::boot.ci(boot1,conf=1-boot_level,index=j,type="perc");
    indirect_effect_boot_estimate[j] <- boot::norm.ci(boot1,conf=.0001,index=j)[2];
    indirect_effect_boot_se[j] <- sd(boot1$t[,j]);
  }
  after_boot <- Sys.time();
  bootstrap_results <- list(indirect_effect_boot_estimate=indirect_effect_boot_estimate,
                            indirect_effect_boot_se=indirect_effect_boot_se,
                            bootstrap_output=boot1,
                            boot_norm=boot_norm,
                            boot_basic=boot_basic,
                            boot_perc=boot_perc,
                            boot_level=boot_level,
                            nboot=nboot,
                            time_required=difftime(after_boot,before_boot));
  if (tvem_do_loop) {
    colnames(ICs_table_from_bootstraps) <- paste(0:(ncol(ICs_table_from_bootstraps)-1),"InteriorKnots",sep="");
    bootstrap_results$ICs_table_from_bootstraps <- ICs_table_from_bootstraps;
    bootstrap_results$num_knots_from_bootstraps <- table(paste(apply(ICs_table_from_bootstraps,1,which.min)+1,"knots"));
  }
  important_variable_names <- list(time = time_variable_name,
                                   treatment = treatment_variable_names,
                                   mediator = mediator_variable_name,
                                   outcome = outcome_variable_name);
  answer <- list(observed_time_grid_for_debug=observed_time_grid,
                 wide_data_for_debug=wide_data,
                 original_results=original_results,
                 bootstrap_results=bootstrap_results,
                 important_variable_names=important_variable_names);
  class(answer) <- "funmediation";
  return(answer);
}
