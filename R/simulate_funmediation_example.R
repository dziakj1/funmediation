#' simulate_funmediation_example function
#'
#' Simulates a dataset for demonstrating the funmediation function.
#'
#' @param nsub Number of subjects
#' @param ntimes Number of potential times that could be observed on each subject
#' @param observe_rate Proportion of potential times on which there are actually
#' observations. Not all times are observed; this is assumed to be completely random and
#' to be done by design to reduce participant burden.
#' @param alpha_int  Function representing the time-varying mean of mediator variable
#'  for the X=0 group
#' @param alpha_X Function representing the time-varying effect of X on the mediator
#' @param beta_M Function representing the functional coefficient for cumulative
#' (scalar-on-function) effect of M on Y adjusting for X
#' @param beta_int  Mean of Y if X is zero and M is the 0 function
#' @param beta_X  Direct effect of X on Y after adjusting for M
#' @param sigma_Y Error standard deviation of the outcome Y
#' @param sigma_M_error Error standard deviation of the mediator M
#' @param rho_M_error Autoregressive correlation coefficient of the mediator M
#' from one observation to the next
#' @param simulate_binary_Y  Whether Y should be generated from a binary
#' logistic (TRUE) or Gaussian (FALSE) model
#' @param make_covariate_S Whether to generate a random binary covariate S  
#'  at the subject (i.e., time-invariant) level. It will be generated to have 
#'  zero population-level relationship to the outcome. 
#'
#' @return A list with the following components:
#' \describe{
#' \item{time_grid}{The time grid for interpreting functional coefficients.}
#' \item{true_alpha_int}{True value of the time-varying alpha_int parameter,
#' representing the time-specific mean of the mediator M when X is 0.}
#' \item{true_alpha_X}{True value of the time-varying alpha_X parameter,
#' representing the effect of X on M.}
#' \item{true_beta_int}{True value of the beta_M parameter, representing
#' the mean of the outcome, Y, when X=0 and M=0.}
#' \item{true_beta_M}{True value of the beta_M parameter, representing the
#' functional effect of the mediator, M, on the outcome, Y.}
#' \item{true_beta_X}{True value of the beta_X parameter, representing the
#' effect of X on the outcome, Y, adjusting for the mediator, M.}
#' \item{true_indirect}{True value of the indirect parameter, representing
#' the indirect (mediated) effect of X on the outcome, Y.}
#' \item{dataset}{The simulated longitudinal dataset in long form.}
#' }
#'
#' @importFrom stats rbinom rnorm
#'
#' @export

simulate_funmediation_example <- function(
  nsub = 500,
  ntimes = 100,
  observe_rate = .4,
  alpha_int = function(t) {return(t^.5)}, # time-varying mean of M for the X=0 group;
  alpha_X = function(t) {return(-(t/2)^.5)}, # time-varying effect of X on M;
  beta_M = function(t) {(1/2)*(exp(t)-1)}, # functional coefficient for cumulative effect of M on Y;
  beta_int = 0,  # mean of Y if the X = 0 and M is the 0 function;
  beta_X = .2,  # direct effect of X on Y after adjusting for M;
  sigma_Y = 1,
  sigma_M_error = 2,
  rho_M_error = .8,
  simulate_binary_Y=FALSE,
  make_covariate_S=FALSE)
{
  time_grid <- (1:ntimes)/ntimes;  # vector of all possible times, scaled within 0 to 1;
  true_indirect <- mean(beta_M(time_grid)*alpha_X(time_grid));
  short_X <- rbinom(nsub,size=1,prob=.5);
  if (make_covariate_S) {
    short_S <- rbinom(nsub,size=1,prob=.5);
  }
  # Simulate M from X...
  autoreg_error <- matrix(0,nsub,ntimes);
  autoreg_error[,1] <- rnorm(n=nsub,mean=0,sd=sigma_M_error);
  for (j in 2:ntimes) {
    autoreg_error[,j] <- rho_M_error*autoreg_error[,j-1] +
      sqrt(1-rho_M_error^2)*rnorm(n=nsub,mean=0,sd=sigma_M_error);
  }
  all_M <- matrix(0,nsub,ntimes);  # time-varying mediator;
  for (i in 1:nsub) {
    all_M[i,] <- alpha_int(time_grid) + short_X[i]*alpha_X(time_grid);
  }
  all_M <- all_M + autoreg_error;
  if (simulate_binary_Y) {
    # Simulate Y from M and X...
    eta <- rep(NA,nsub); # = E(Y|X,M);
    for (i in 1:nsub) {
      eta[i] <- beta_int + mean(beta_M(time_grid) * all_M[i,]) + beta_X*short_X[i];
    }
    mu <- exp(eta)/(1+exp(eta));
    short_Y <- unlist(lapply(X=mu,FUN=rbinom,size=1,n=1));
  } else {
    mu <- rep(NA,nsub); # = E(Y|X,M);
    for (i in 1:nsub) {
      mu[i] <- beta_int + mean(beta_M(time_grid) * all_M[i,]) + beta_X*short_X[i];
    }
    short_Y <- round(mu + rnorm(n=nsub,mean=0,sd=sigma_Y),5);
  }
  # Assemble simulated data into a long-form dataset:
  M <- all_M;
  for (i in 1:nsub) {
    which.missing.for.this.person <- which(rbinom(ntimes,1,1-observe_rate)==1);
    M[i,which.missing.for.this.person] <- NA;
  }
  temp <- data.frame(
    subject_id=rep(1:nsub,each=ntimes),
    t=rep(time_grid,times=nsub),
    X=rep(short_X,each=ntimes),
    M=as.vector(t(M)),
    Y=rep(short_Y,each=ntimes)
  );
  if (make_covariate_S) {temp$S <- rep(short_S,each=ntimes);}
  long_simulated_data <- temp[which(!is.na(temp$M)),];
  rownames(long_simulated_data) <- NULL;
  return(list(time_grid=time_grid,
              true_alpha_int=alpha_int(time_grid),
              true_alpha_X=alpha_X(time_grid),
              true_beta_int=beta_int,
              true_beta_M=beta_M(time_grid),
              true_beta_X=beta_X,
              true_indirect=true_indirect,
              dataset=long_simulated_data));
}
