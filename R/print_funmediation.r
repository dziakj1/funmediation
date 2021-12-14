#' print.funmediation:  Print output from a model that was fit by the funmediation function.
#'
#' @param x The funmediation object (output of the funmediation function)
#' @param ... Further arguments currently not supported
#'
#' @import tvem
#' @importFrom utils getS3method
#' @return This function does not return an object, but is called for its side effect of printing information.
#'
#' @export
#' @method print funmediation

print.funmediation <- function(x, ...) {
  cat("======================================================= \n");
  cat("Functional Regression Mediation Function Output \n");
  cat("======================================================= \n");
  treatment_variable_names <- x$important_variable_names[["treatment"]];
  num_treatment_variables <- length(treatment_variable_names);
  cat(paste("TREATMENT:",
            paste(treatment_variable_names,collapse=" "),
            "\n"));
  cat(paste("MEDIATOR:",
            x$important_variable_names[["mediator"]],
            ifelse(x$original_results$binary_mediator,
                   "(Assumed Binary)",
                   "(Assumed Continuous)"),
            "\n"));
  cat(paste("OUTCOME:",
            x$important_variable_names[["outcome"]],
            ifelse(x$original_results$binary_outcome,
                   "(Assumed Binary)",
                   "(Assumed Continuous)"),
            "\n"));
  cat("=======================================================");
  for (j in 1:num_treatment_variables) {
    cat(paste(" \nDIRECT EFFECT OF",treatment_variable_names[j],"...\n"))
    cat(" Direct effect estimate:\n" );
    cat(x$original_results$beta_X_estimate[j]);
    cat("\n Direct effect standard error:\n" );
    cat(x$original_results$beta_X_se[j]);
    cat("\n Direct effect bootstrap estimate:\n" );
    cat(x$bootstrap_results$direct_effect_boot_estimate[j]);
    cat("\n Direct effect bootstrap confidence interval: ");
    cat("\n  ... by normal method:\n ");
    cat(c(round(x$bootstrap_results$direct_effect_boot_norm[[j]]$normal[,2],4), ", ",
          round(x$bootstrap_results$direct_effect_boot_norm[[j]]$normal[,3],4)));
    if (x$bootstrap_results$nboot>=50) {
      cat("\n  ... by percentile method:\n");
      cat(c(round(x$bootstrap_results$direct_effect_boot_perc[[j]]$percent[,4],4), ", ",
            round(x$bootstrap_results$direct_effect_boot_perc[[j]]$percent[,5],4)));
    }
  }
  cat("\n======================================================= ");
  for (j in 1:num_treatment_variables) {
    cat(paste("\nINDIRECT EFFECT OF",treatment_variable_names[j],"...\n"))
    cat(" Indirect effect parametric estimate:\n ");
    cat(x$original_results$indirect_effect_estimate[j]);
    cat("\n Indirect effect bootstrap estimate:\n ");
    cat(x$bootstrap_results$indirect_effect_boot_estimate[j]);
    cat("\n Indirect effect bootstrap confidence interval: ");
    cat("\n  ... by normal method:\n ");
    cat(c(round(x$bootstrap_results$indirect_effect_boot_norm[[j]]$normal[,2],4), ", ",
          round(x$bootstrap_results$indirect_effect_boot_norm[[j]]$normal[,3],4)));
    if (x$bootstrap_results$nboot>=50) {
      cat("\n  ... by percentile method:\n");
      cat(c(round(x$bootstrap_results$indirect_effect_boot_perc[[j]]$percent[,4],4), ", ",
            round(x$bootstrap_results$indirect_effect_boot_perc[[j]]$percent[,5],4)));
    }
  }
  cat("\n======================================================= \n");
  cat("Computation time:\n");
  print(x$bootstrap_results$time_required);
  cat("======================================================= \n");
  cat("Time-varying Effects Model Predicting MEDIATOR from TREATMENT:\n");
  workaround_print <- getS3method("print", "tvem");
  workaround_print(x$original_results$tvem_XM_details,ornate=FALSE);
  cat("======================================================= \n");
  cat("Functional Regression Model Predicting OUTCOME from MEDIATOR \n");
  cat("  and TREATMENT: \n");
  print(x$original_results$funreg_MY_details);
  cat("Scalar terms in functional regression model:\n")
  temp <- x$original_results$funreg_MY_details$coefficients;
  print(temp[-grep(x=names(temp),pattern=":")]);
  cat("======================================================= \n");
  cat("Parametric model for Predicting OUTCOME from TREATMENT \n");
  cat("  ignoring MEDIATOR: \n");
  print(summary(x$original_results$total_effect_details));
  cat("======================================================= \n");
  if(!is.null(x$original_results$tvem_IC_table)) {
    cat("ICs table for selecting number of interior knots in TVEM:\n");
    print(x$original_results$tvem_IC_table);
    cat("======================================================= \n");
  }
  if(!is.null(x$bootstrap_results$num_knots_from_bootstraps)) {
    cat("Numbers of interior knots selected in bootstrap samples:\n");
    print(x$bootstrap_results$num_knots_from_bootstraps);
    cat("======================================================= \n");
  }
}
