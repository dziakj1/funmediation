#' simulate_funmediation_example function
#'
#' Simulates a dataset for demonstrating the funmediation function.
#'
#' @param nsub Number of subjects
#' @param ntimes Number of potential times that could be observed on each subject
#' @param nlevels Number of treatment groups or levels on the treatment variable X.
#' Subjects are assumed to be randomly assigned to each level with equal
#' probability (i.e., the probability per level is 1/nlevel).  Default is 2 for
#' a randomized controlled trial with a control group X=0 and an experimental
#' group X=1.  There should not be less than 2 or more than 5 groups for purposes
#' of this function.
#' @param observe_rate Proportion of potential times on which there are actually
#' observations. Not all times are observed; this is assumed to be completely random and
#' to be done by design to reduce participant burden.
#' @param alpha_int  Function representing the time-varying mean of mediator variable
#' for the level of treatment with all treatment dummy codes X set to 0 (e.g., the
#' control group).
#' @param alpha_X Function representing the time-varying effect of X on the
#' mediator (if there are two treatment levels) or a list of nlevels-1
#' functions representing the effect of receiving each nonzero level of X rather
#' than control (if there are more than two treatment levels).
#' @param beta_M Function representing the functional coefficient for cumulative
#' (scalar-on-function) effect of the mediator M on the treatment Y adjusting
#' for the treatment X
#' @param beta_int Mean of Y if the X is zero and M is the 0 function
#' @param beta_X Numeric value representing the direct effect of X on Y after adjusting
#' for M (if there are two treatment levels) or a vector of nlevels-1 numeric values
#'  (if there are more than two treatment levels)
#' @param sigma_Y Error standard deviation of the outcome Y (conditional on
#' treatment and mediator trajectory)
#' @param sigma_M_error Error standard deviation of the mediator M (conditional
#' on treatment and time)
#' @param rho_M_error Autoregressive correlation coefficient of the error
#'  in the mediator M, from one observation to the next
#' @param simulate_binary_Y  Whether Y should be generated from a binary
#' logistic (TRUE) or Gaussian (FALSE) model
#'
#' @return A list with the following components:
#' \describe{
#' \item{time_grid}{The time grid for interpreting functional coefficients.}
#' \item{true_alpha_int}{True value of the time-varying alpha_int parameter,
#' representing the time-specific mean of the mediator M when the
#' treatment value X is 0. }
#' \item{true_alpha_X}{True value of the time-varying alpha_X parameter,
#' representing the effect of X on M.  This is a single number if nlevels=2,
#' or a vector of effects if nlevels>2.}
#' \item{true_beta_int}{True value of the beta_M parameter, representing
#' the mean of the outcome Y when X=0 and M=0.}
#' \item{true_beta_M}{True value of the beta_M parameter, representing the
#'  functional effect of treatment on the outcome Y.}
#' \item{true_beta_X}{True value of the beta_X parameter, representing the
#'  effect of treatment on the outcome Y adjusting for the mediator. This
#'  is a single function if nlevels=2, or a vector of functions if nlevels>2.}
#' \item{true_indirect}{True value of the indirect parameter, representing
#'  the indirect (mediated) effect of treatment on the outcome Y.  This is
#'  a single number if nlevels=2, or a vector of effects if nlevels>2.}
#' \item{dataset}{The simulated longitudinal dataset in long form.}
#' }
#'
#' @examples
#' set.seed(123)
#' # Simplest way to call the function:
#' simulation_all_defaults <- simulate_funmediation_example()
#' summary(simulation_all_defaults)
#' head(simulation_all_defaults)
#' # Changing the sample size to be larger:
#' simulation_larger <- simulate_funmediation_example(nsub=10000)
#' summary(simulation_larger)
#' # Changing the effect of the mediator to be null:
#' simulation_null <- simulate_funmediation_example(beta_M=function(t) {return(0*t)})
#' summary(simulation_null)
#' # Simulating a exposure variable with three levels (two dichotomous dummy codes)
#' simulation_three_group <- simulate_funmediation_example(nlevels=3,
#'                               alpha_X = list(function(t) {return(.1*t)},
#'                                              function(t) {return(-(t/2)^.5)}),
#'                               beta_X = c(-.2,.2))
#' print(summary(simulation_three_group));
#'
#' @export

simulate_funmediation_example <- function(
  nsub = 500,
  nlevels = 2,
  ntimes = 100,
  observe_rate = .4,
  alpha_int = function(t) {return(t^.5)}, # time-varying mean of mediator variable for the X=0 group;
  alpha_X = function(t) {return(-(t/2)^.5)}, # time-varying effect of X on the mediator;
  beta_M = function(t) {(1/2)*(exp(t)-1)}, # functional (funreg) coefficient for cumulative effect of M on Y;
  beta_int = 0,  # mean of Y if X=0 and M(t)=0;
  beta_X = .2,  # direct effect of X on Y after adjusting for M;
  sigma_Y = 1,
  sigma_M_error = 2,
  rho_M_error = .8,
  simulate_binary_Y=FALSE,
  make_covariate_S=FALSE)
{
  #############################################
  # Check that everything has the right length:
  if (nlevels<2) {stop("The treatment variable needs to have at least two levels.")};
  if (nlevels>5) {stop("Please choose 5 or fewer levels for the treatment variable.")};
  if (nlevels==2) {
    if (!is.list(alpha_X)) {alpha_X <- list(alpha_X);}
    # For convenience, we'll always treat alpha_X as a list of functions, even
    # if it is a list with only one item on the list (which still counts as a
    # list in R)
  }
  if (!is.list(alpha_X)) {stop("Please specify alpha_X as a list of functions.")};
  if (length(alpha_X)<(nlevels-1)) {
    stop(paste("Please specify a function in alpha_X for the",
               "effect of each level of X other than the control level."));
  }
  if (length(alpha_X)>(nlevels-1)) {
    stop(paste("There are too many functions in the list alpha_X."));
  }
  if (length(beta_X)<(nlevels-1)) {
    stop(paste("Please specify a value in beta_X for the",
               "effect of each level of X other than the control level."));
  }
  #############################################
  # Generate time grid:
  time_grid <- (1:ntimes)/ntimes;  # vector of all possible times, scaled within 0 to 1;
  # Determine true value of indirect effects for each non-control level of X:
  true_indirect <- rep(NA,nlevels-1);
  for (this_nonreference_level in 1:(nlevels-1)) {
    this_alpha_X_function <- alpha_X[[this_nonreference_level]];
    true_indirect[this_nonreference_level] <- mean(beta_M(time_grid)*
                                                     this_alpha_X_function(time_grid));
  }
  if (nlevels==2) {
  short_X <- rbinom(nsub,size=1,prob=.5);
  } else {
    short_X_as_multinomial <- sample.int(n=nlevels,
                                         size=nsub,
                                         replace=TRUE)-1;
    short_X <- matrix(NA,
                      nrow=length(short_X_as_multinomial),
                      ncol=nlevels-1);
    for (j in 1:(nlevels-1)) {
      short_X[,j] <- 1*(short_X_as_multinomial==j);
    }
  }
  if (make_covariate_S) {
    short_S <- rbinom(nsub,size=1,prob=.5);
  }
  #############################################
  # Simulate M from X...
  autoreg_error <- matrix(0,nsub,ntimes);
  autoreg_error[,1] <- rnorm(n=nsub,
                             mean=0,
                             sd=sigma_M_error);
  for (j in 2:ntimes) {
    autoreg_error[,j] <- rho_M_error*autoreg_error[,j-1] +
      sqrt(1-rho_M_error^2)*rnorm(n=nsub,
                                  mean=0,
                                  sd=sigma_M_error);
  }
  all_M <- matrix(0,nsub,ntimes);  # time-varying mediator;
  if (nlevels==2) {
    for (i in 1:nsub) {
      all_M[i,] <- alpha_int(time_grid) + short_X[i]*alpha_X[[1]](time_grid);
    }
  } else {
    stopifnot(length(alpha_X)==ncol(short_X));
    for (i in 1:nsub) {
      all_M[i,] <- alpha_int(time_grid);
      for (j in 1:ncol(short_X)) {
        all_M[i,] <- all_M[i,] + short_X[i,j]*alpha_X[[j]](time_grid);
      }
    }
  }
  all_M <- all_M + autoreg_error;
  eta <- rep(NA,nsub); # = E(Y|X,M);
  for (i in 1:nsub) {
    eta[i] <- beta_int + mean(beta_M(time_grid) * all_M[i,])
    if (nlevels==2) {
      eta[i] <- eta[i] + beta_X*short_X[i];
    } else {
      for (j in 1:ncol(short_X)) {
        eta[i] <- eta[i] + beta_X[j]*short_X[i,j];
      }
    }
  }
  if (simulate_binary_Y) {
    # Simulate Y from M and X...
    mu <- exp(eta)/(1+exp(eta));
    short_Y <- unlist(lapply(X=mu,FUN=rbinom,size=1,n=1));
  } else {
    mu <- eta;
    short_Y <- round(mu + rnorm(n=nsub,mean=0,sd=sigma_Y),5);
  }

  # Assemble simulated data into a long-form dataset:
  M <- all_M;
  for (i in 1:nsub) {
    which.missing.for.this.person <- which(rbinom(ntimes,1,1-observe_rate)==1);
    M[i,which.missing.for.this.person] <- NA;
  }
  if (nlevels==2) {
    temp <- data.frame(
      subject_id=rep(1:nsub,each=ntimes),
      t=rep(time_grid,times=nsub),
      X=rep(short_X,each=ntimes),
      M=as.vector(t(M)),
      Y=rep(short_Y,each=ntimes));
  } else {
    temp <- data.frame(subject_id=rep(1:nsub,each=ntimes),
                       t=rep(time_grid,times=nsub));
    for (j in 1:ncol(short_X)) {
      this_long_X <- rep(short_X[,j],each=ntimes);
      temp <- cbind(temp,this_long_X);
      colnames(temp)[ncol(temp)] <- paste("X",j,sep="");
    }
    temp <- cbind(temp,
                  data.frame(M=as.vector(t(M)),
                             Y=rep(short_Y,each=ntimes)));
  }
  if (make_covariate_S) {temp$S <- rep(short_S,each=ntimes);}
  long_simulated_data <- temp[which(!is.na(temp$M)),];
  rownames(long_simulated_data) <- NULL;
  true_alpha_X <- lapply(1:length(alpha_X),function(j) {return(alpha_X[[j]](time_grid))});
  if (nlevels==2) {
    true_alpha_X <- unlist(true_alpha_X);
  }
  return(list(time_grid=time_grid,
              true_alpha_int=alpha_int(time_grid),
              true_alpha_X=true_alpha_X,
              true_beta_int=beta_int,
              true_beta_M=beta_M(time_grid),
              true_beta_X=beta_X,
              true_indirect=true_indirect,
              dataset=long_simulated_data));
}
