#' plot.funmediation: Produces plots for a funmediation model.
#'
#' Produces plots from a funmediation object produced by
#' the funmediation function. These plots will be shown on the default
#' output device (likely the screen); they can of course be
#' written to a file instead, by preceding the call to plot
#' with a call to png(), pdf(), or other R graphic file output functions.
#'
#' @param x The funmediation object to be plotted.
#' @param use_panes Whether to plot multiple coefficient functions in a single image.
#' @param what_plot One of "pfr","coefs", or "tvem."  These options are as follows:
#' \describe{
#' \item{pfr}{For a "pfr" plot, the functional coefficient for predicting the
#' outcome, Y, from the mediator, M (conditional on X), is shown, by calling
#' the plot method for the penalized functional regression results. See the
#' documentation for plot.gam() in the mgcv function for more information,
#' because the pfr() function uses gam() as a back end.}
#' \item{coefs}{For a "coefs" plot, the three important functional coefficients
#' in the model (intercept for predicting M, effect of X on M,
#' and the effect of M on Y adjusting for X) are plotted one after
#' another. That is, the plots are shown for the alpha_int_estimate, alpha_X_estimate,
#' and beta_M_estimate, each as a function of time_grid.  Approximate pointwise
#' 95% confidence intervals are also shown if possible.  If there is only
#' one dichotomous treatment variable and panes are being used, the lower right
#' pane will be free, so the indirect effect will be printed there even
#' though it is a scalar.}
#' \item{tvem}{For a "tvem" plot, the functional coefficients in the TVEM model
#' predicting M from X are displayed.}
#' }
#' @param ... Further arguments currently not supported
#'
#' @import tvem
#' @importFrom graphics par plot text lines
#' @export
#' @method plot funmediation
#'
plot.funmediation <- function(x,
                              use_panes=TRUE,
                              what_plot=c("pfr",
                                          "coefs",
                                          "tvem"), ...) {
  num_treatment_variables <- length(x$important_variable_names$treatment);
  if (use_panes==FALSE) {par(mfrow=c(1,1));}
  what <- match.arg(what_plot)
  if (what=="pfr") {
    if (use_panes) {par(mfrow=c(1,1));}
    plot(x$original_results$funreg_MY_details,
         xlab="Time",
         ylab="Coefficient function",
         main="Functional effect of mediator on outcome"); # plot.gam function;
  }
  if (what=="coefs") {
    if (use_panes) {
      if (num_treatment_variables==1 |
          num_treatment_variables==2) {
        par(mfrow=c(2,2),mar=c(2,2,2,2));
      }
      if (num_treatment_variables==3 |
          num_treatment_variables==4) {
        par(mfrow=c(3,2),mar=c(2,2,2,2));
      }
      if (num_treatment_variables>4) {
        stop("There are too many treatment variables to draw this plot legibly.")
      }
    }
    plot(x$original_results$time_grid,
         x$original_results$alpha_int_estimate,
         xlab=expression(t),
         ylab=expression(alpha[int](t)),
         main="Intercept for Mediator");
    lines(x$original_results$time_grid,
          x$original_results$alpha_int_estimate-1.96*x$original_results$alpha_int_se);
    lines(x$original_results$time_grid,
          x$original_results$alpha_int_estimate+1.96*x$original_results$alpha_int_se);
    if (num_treatment_variables==1) {
      plot(x$original_results$time_grid,
           x$original_results$alpha_X_estimate,
           xlab=expression(t),
           ylab=expression(alpha[X](t)),
           main="Treatment on Mediator");
      lines(x$original_results$time_grid,
            x$original_results$alpha_X_estimate-1.96*x$original_results$alpha_X_se);
      lines(x$original_results$time_grid,
            x$original_results$alpha_X_estimate+1.96*x$original_results$alpha_X_se);
    }
    if (num_treatment_variables==2 |
        num_treatment_variables==3 |
        num_treatment_variables==4 ) {
      for (j in 1:num_treatment_variables) {
        plot(x$original_results$time_grid,
             x$original_results$alpha_X_estimate[[j]],
             xlab=expression(t),
             ylab=expression(alpha[X](t)),
             main=paste(x$important_variable_names$treatment[j],
                        "on Mediator"));
        lines(x$original_results$time_grid,
              x$original_results$alpha_X_estimate[[j]]-1.96*x$original_results$alpha_X_se[[j]]);
        lines(x$original_results$time_grid,
              x$original_results$alpha_X_estimate[[j]]+1.96*x$original_results$alpha_X_se[[j]]);
      }
    }
    plot(x$original_results$time_grid,
         x$original_results$beta_M_estimate,
         xlab=expression(t),
         ylab=expression(beta[M](t)),
         main="Mediator on Outcome");
    lines(x$original_results$time_grid,
          x$original_results$beta_M_estimate-1.96*x$original_results$beta_M_se);
    lines(x$original_results$time_grid,
          x$original_results$beta_M_estimate+1.96*x$original_results$beta_M_se);
    if (num_treatment_variables==1 & use_panes==TRUE) {
      # There is an extra pane, so we can write the indirect effect
      # estimate in the white space.
      plot(x=0:5,y=0:5,type="n",bty="n",xaxt="n",yaxt="n",xlab="",ylab="");
      text(x=1.1,y=2,pos=4,labels=paste("Indirect effect:"));
      est2 <- round(x$bootstrap_results$indirect_effect_boot_estimate,3);
      text(x=1.2,y=1,pos=4,labels=est2);
    }
  }
  if (what=="tvem") {
    par(mar=rep(2,4),cex.main=.7);
    plot(x$original_results$tvem_XM_details);
  }
}
