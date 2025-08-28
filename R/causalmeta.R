#' Causal Meta Analysis
#'
#' @description
#' Function to perform causal meta-analysis for binary treatment and binary outcomes.
#'
#' @param measure a character string to specify which effect size or outcome measure should be calculated
#' @param ai vector with the 2x2 table frequencies (upper left cell)
#' @param bi vector with the 2x2 table frequencies (upper right cell)
#' @param ci vector with the 2x2 table frequencies (lower left cell)
#' @param di vector with the 2x2 table frequencies (lower right cell)
#' @param n1i vector with the group sizes or row totals (first group/row)
#' @param n2i vector with the group sizes or row totals (second group/row)
#' @param slab optional vector with labels for the studies
#' @param data optional data of type data frame, matrix, or dgCMatrix. If NULL, vectors must be provided directly
#' @param weights optional vector of study weights, a value of NULL (the default) corresponds to study size / total size
#' @param plot set to TRUE to print forestplot, default is FALSE
#' @param log.scale set to TRUE to have log scale of the measure (only for RR, OR and SR), default is FALSE
#' @param ...
#'
#' @details The function can be used by either inputting "ai,bi,ci,di" or "ai,ci,n1i,n2i". Accepted measures are: "RD" = Risk Difference, "RR" = Risk Ratio, "OR" = Odds Ratio, "SR" = Survival Ratio
#'
#' @return A list with the following elements:
#' \describe{
#' \item{study_results}{list of effect sizes or outcome measures for each study and associated standard errors}
#' \item{final_result}{aggregated effect size or outcome measure and associated standard error}
#' }
#'
#' @examples
#' ## Example 1: With data frame
#' data <- data.frame(
#'   study = paste("Study", 1:5),
#'   treated_events = c(10, 15, 8, 12, 20),
#'   treated_total = c(100, 120, 80, 110, 150),
#'   control_events = c(15, 20, 12, 18, 25),
#'   control_total = c(100, 115, 85, 105, 145)
#' )
#'
#' # Risk difference
#' result <- causalmeta(measure = "RD", ai = treated_events, n1i = treated_total,
#'                     ci = control_events, n2i = control_total,
#'                     data = data, slab = study)
#'
#' ## Example 2: With vectors
#' treated_positives <- c(13, 5, 14, 18, 9)
#' treated_negatives <- c(80, 63, 72, 130, 100)
#' control_positives <- c(25, 18, 34, 23, 16)
#' control_negatives <- c(125, 98, 165, 117, 85)
#'
#'
#' # Risk ratio
#' result <- causalmeta(measure = "RR", ai = treated_positives, bi = treated_negatives,
#'                     ci = control_positives, di = control_negatives)
#'
#'
#'
#' @references
#' Berenfeld, C., Boughdiri, A., Colnet, B., van Amsterdam, W. A. C., Bellet, A.,
#' Khellaf, R., Scornet, E., & Josse, J. (2025). Causal Meta-Analysis: Rethinking
#' the Foundations of Evidence-Based Medicine. arXiv:2505.20168.
#' \url{https://arxiv.org/abs/2505.20168}
#'
#'
#' @importFrom purrr pmap_dfr
#' @importFrom metafor rma escalc forest addpoly
#' @importFrom tibble tibble
#'
#' @export
#'
#'
causalmeta <- function(measure, ai, bi, ci, di, n1i, n2i, slab, data = NULL, weights = NULL, plot = FALSE, log.scale = FALSE, ...){

  # Input validation on required parameters (measure, ai, ci)

  if(missing(measure)){stop("Input 'measure' is required", call. = FALSE)}

  if(!is.element(measure,c("RD","RR","OR","SR"))){stop("Unsupported measure:", measure, call. = FALSE)}

  if(missing(ai)){stop("Input 'ai' is required", call. = FALSE)}
  if(missing(ci)){stop("Input 'ci' is required", call. = FALSE)}

  # Create logical for type of parameters in input
  has_n_totals <- !(missing(n1i) || missing(n2i))
  has_table_entries <- !(missing(bi) || missing(di))

  if (!(has_n_totals || has_table_entries)) {
    stop("Must provide either ('bi', 'di') or ('n1i', 'n2i')", call. = FALSE)
  }

  # Stop if measure is Risk Differene and log scale is TRUE
  if(measure == "RD" && log.scale == TRUE){stop("Risk Difference cannot be computed in log scale", call. = FALSE)}

  # Validate values of inputs and create list needed for next part of function
  if(is.null(data)){

    # Validate inputs
    dat <- .validate_vector_inputs(ai = ai, bi = bi, ci = ci, di = di, n1i = n1i, n2i = n2i, has_n = has_n_totals, has_t = has_table_entries)
    if(!missing(slab)) {slab_vals <- slab}

  } else if(is.data.frame(data)){

    # Check: is variable in dataset?
    .validate_variables(names = names(data), ai = deparse(substitute(ai)), bi = deparse(substitute(bi)),
                        ci = deparse(substitute(ci)), di = deparse(substitute(di)), n1i = deparse(substitute(n1i)),
                        n2i = deparse(substitute(n2i)), has_n = has_n_totals, has_t = has_table_entries)

    # Validate inputs
    dat <- .validate_vector_inputs(ai = data[[deparse(substitute(ai))]], bi = data[[deparse(substitute(bi))]],
                                   ci = data[[deparse(substitute(ci))]], di = data[[deparse(substitute(di))]],
                                   n1i = data[[deparse(substitute(n1i))]], n2i = data[[deparse(substitute(n2i))]],
                                   has_n = has_n_totals, has_t = has_table_entries)
    if(!missing(slab)) {slab_vals <- data[[deparse(substitute(slab))]]}

  } else if(is.matrix(data) || inherits(data, "Matrix")){

    # Check: is variable in matrix?
    .validate_variables(names = colnames(data), ai = deparse(substitute(ai)), bi = deparse(substitute(bi)),
                        ci = deparse(substitute(ci)), di = deparse(substitute(di)), n1i = deparse(substitute(n1i)),
                        n2i = deparse(substitute(n2i)), has_n = has_n_totals, has_t = has_table_entries)

    # Validate inputs
    dat <- .validate_vector_inputs(ai = data[, deparse(substitute(ai))], bi = data[, deparse(substitute(bi))],
                                   ci = data[, deparse(substitute(ci))], di = data[, deparse(substitute(di))],
                                   n1i = data[, deparse(substitute(n1i))], n2i = data[, deparse(substitute(n2i))],
                                   has_n = has_n_totals, has_t = has_table_entries)
    if(!missing(slab)) {slab_vals <- data[, deparse(substitute(slab))]}

  } else {
    stop("Input must be vectors, data frame, matrix, or sparse matrix", call. = FALSE)
  }

  # Calculate the total sample size across all studies
  n_total <- sum(dat$n1i_vals + dat$n2i_vals)

  # For each study (row in 'data'), calculate required statistics
  results <- purrr::pmap_dfr(dat, function(ai_vals, ci_vals, n1i_vals, n2i_vals, ...) {
    # Skip studies with zero sample size in either arm
    if (n1i_vals == 0 || n2i_vals == 0) return(NULL)

    # Calculate control and experimental event rates for each study
    mu_0_k <- ci_vals / n2i_vals   # Control arm event proportion
    mu_1_k <- ai_vals / n1i_vals  # Experimental arm event proportion

    # Weight (study size / total size)

    proba_k <- if (is.null(weights)) (n1i_vals + n2i_vals) / n_total else weights / sum(weights)

    # if(plot){
      if(measure == "RD"){
        yi <- mu_1_k - mu_0_k
        vi <- (mu_0_k*(1- mu_0_k) / n1i_vals) + (mu_1_k*(1 - mu_1_k) / n2i_vals)
        ci.lb <- yi - 1.96*sqrt(vi)
        ci.ub <- yi + 1.96*sqrt(vi)
      }
      else if(measure == "RR"){
        ratio_i <- mu_1_k / mu_0_k
        vi_log <- (n1i_vals - ai_vals) / (n1i_vals * ai_vals) + (n2i_vals - ci_vals) / (n2i_vals * ci_vals)
        ci_log_i <- 1.96 * sqrt(vi_log)

        if(log.scale){ yi <- log(ratio_i); vi <- vi_log; ci.lb <- yi - ci_log_i; ci.ub <- yi + ci_log_i}
        else if(!log.scale){ yi <- ratio_i; vi <- vi_log*ratio_i^2; ef <- exp(ci_log_i); ci.lb <- yi / ef; ci.ub <- yi * ef}

      }
      else if(measure == "OR"){
        ratio_i <- (ai_vals * (n2i_vals - ci_vals)) / (ci_vals * (n1i_vals - ai_vals))
        vi_log <- (1/ai_vals) + (1/(n1i_vals - ai_vals)) + (1/ci_vals) + (1/(n2i_vals - ci_vals))
        ci_log_i <- 1.96 * sqrt(vi_log)

        if(log.scale){ yi <- log(ratio_i); vi <- vi_log; ci.lb <- yi - ci_log_i; ci.ub <- yi + ci_log_i}
        else if(!log.scale){ yi <- ratio_i; vi <- vi_log*ratio_i^2; ef <- exp(ci_log_i); ci.lb <- yi / ef; ci.ub <- yi * ef}

      }
      else if(measure == "SR"){
        ratio_i <- (1 - mu_1_k) / (1 - mu_0_k)
        vi_log <- (ai_vals / (n1i_vals * (n1i_vals - ai_vals))) + (ci_vals / (n2i_vals * (n2i_vals - ci_vals)))
        ci_log_i <- 1.96 * sqrt(vi_log)

        if(log.scale){ yi <- log(ratio_i); vi <- vi_log; ci.lb <- yi - ci_log_i; ci.ub <- yi + ci_log_i}
        else if(!log.scale){ yi <- ratio_i; vi <- vi_log*ratio_i^2; ef <- exp(ci_log_i); ci.lb <- yi / ef; ci.ub <- yi * ef}

      }

      tibble::tibble(
        proba_k = proba_k,
        mu_0_k = mu_0_k,
        mu_1_k = mu_1_k,
        yi = yi,
        sei = sqrt(vi),
        ci.lb = ci.lb,
        ci.ub = ci.ub,
        n1i = n1i_vals,
        n2i = n2i_vals,
        n.k = n1i_vals + n2i_vals
      )

    # }
    #
    # # Store computed statistics in a tibble (row)
    # else{tibble::tibble(
    #   proba_k = proba_k,
    #   mu_0_k = mu_0_k,
    #   mu_1_k = mu_1_k,
    #   n1i = n1i_vals,
    #   n2i = n2i_vals,
    #   n.k = n1i_vals + n2i_vals
    # )
    # }
  })

  if(!missing(slab)) {results$study <- slab_vals}
  else {results$study <- paste0("Study ", seq_along(results$yi))}

  # Number of studies with non-zero sample size
  n_study <- nrow(results)

  if (n_study == 0) {
    # Return NA if no studies are left
    return(tibble::tibble(estimate = NA, variance = NA, n_study = 0))
  }

  # Calculate pooled (weighted average) event rates across studies
  EY0 <- sum(results$mu_0_k * results$proba_k) # Control
  EY1 <- sum(results$mu_1_k * results$proba_k) # Experimental

  # Compute the plug-in estimator
  final_result <- switch(measure,
                         "RD" = EY1 - EY0,
                         "RR" = if (EY0 != 0) EY1 / EY0 else NA_real_,
                         "OR" = {
                           if (EY1 != 1 && EY0 != 1 && EY0 != 0) {
                             (EY1 / (1 - EY1)) / (EY0 / (1 - EY0))
                           } else NA_real_
                         },
                         "SR" = if (EY0 != 1) (1 - EY1) / (1 - EY0) else NA_real_
  )

  # Variance estimation (depends on effect measure)
  variance <- switch(measure,
                     "RD" = {
                       theta_k <- results$mu_1_k - results$mu_0_k
                       theta <- sum(results$n.k * theta_k / n_total)
                       xi1_sq <- sum((results$n.k^2) / (results$n1i * n_total) * results$mu_1_k * (1 - results$mu_1_k))
                       xi0_sq <- sum((results$n.k^2) / (results$n2i * n_total) * results$mu_0_k * (1 - results$mu_0_k))
                       part2 <- sum(results$n.k / n_total * theta_k^2) - theta^2
                       (xi1_sq + xi0_sq + part2) / n_total
                     },
                     "RR" = {
                       if (!is.na(final_result)) {
                         zeta1_sq <- sum((results$n.k^2) / (results$n1i * n_total) * results$mu_1_k) / (EY1^2)
                         zeta0_sq <- sum((results$n.k^2) / (results$n2i * n_total) * results$mu_0_k) / (EY0^2)
                         var <- (zeta1_sq + zeta0_sq) / n_total
                         if(!log.scale){var <- var*final_result^2} # delta method on log(RR): divide by ratio^2 variance of RR
                         var
                       } else NA_real_
                     },
                     "OR" = {
                       if(!is.na(final_result)) {
                         theta_k <- (results$mu_1_k * (1 -  results$mu_0_k)) / (results$mu_0_k * (1 -  results$mu_1_k))
                         theta <- sum(results$n.k * theta_k / n_total)
                         xi_k <- results$n.k * (((2*results$mu_0_k - 1)^2 / (results$mu_0_k * (1 - results$mu_0_k) * results$n2i)) + ((2*results$mu_1_k - 1)^2 / (results$mu_1_k * (1 - results$mu_1_k) * results$n1i)) )
                         var <- (sum(results$n.k * xi_k / n_total) + sum(results$n.k * log(theta_k)^2 / n_total) - log(theta)^2) / n_total
                         if(!log.scale){var <- var * final_result^2}
                         # if(log.scale){
                         #   xi_k <- results$n.k * ((1 / (results$mu_0_k * (1 - results$mu_0_k) * results$n2i)) + (1 / (results$mu_1_k * (1 - results$mu_1_k) * results$n1i)) )
                         #   var <- (sum(results$n.k * xi_k / n_total) + sum(results$n.k * log(theta_k)^2 / n_total) - log(theta)^2) / n_total
                         # } else {
                         #   xi_k <- (results$n.k^2 / n_total) * ( ((1 - results$mu_0_k)^3 / (results$n2i * results$mu_0_k * (1 - results$mu_1_k)^4)) + ( results$mu_1_k^3 / (results$n1i * (1 - results$mu_1_k) * results$mu_0_k^4)) )
                         #   var <- (sum(results$n.k * xi_k / n_total) + sum(results$n.k * theta_k^2 / n_total) - theta^2) / n_total
                         # }
                         var
                       } else NA_real_
                     },
                     "SR" = {
                       if (!is.na(final_result)) {
                        zeta1_sq <- sum((results$n.k^2) / (results$n1i * n_total) * (1 - results$mu_1_k)) / (EY1^2)
                        zeta0_sq <- sum((results$n.k^2) / (results$n2i * n_total) * (1 - results$mu_0_k)) / (EY0^2)
                        var <- (zeta1_sq + zeta0_sq) / n_total
                        if(!log.scale){var <- var*final_result^2} # delta method on log(SR): divide by ratio^2 variance of SR
                        var
                       } else NA_real_
                     }
  )

  if(log.scale){final_result <- log(final_result)} # Change it after computing the variance since we need ratios to compute variance

  result <- tibble::tibble(estimate = final_result, se = sqrt(variance))

  if(plot){

    res <- metafor::rma(yi = result$estimate, sei = result$se, method="FE")

    if(measure!="SR"){
      # Random effects model
      dat <- metafor::escalc(measure=measure, ai=ai_vals, n1i=n1i_vals, ci=ci_vals, n2i=n2i_vals, data=dat)
      res_random_effects <- metafor::rma(dat$yi, dat$vi, method="DL")
      if(is.element(measure,c("OR","RR")) && log.scale == FALSE){
        res_random_effects_beta <- exp(res_random_effects$beta) # estimate of the measure with random effects model
        ef_random_effects <- exp(1.96 * res_random_effects$se) # error factor
        # upper and lower 95% confidence bounds for the estimate
        res_random_effects_ci.lb <- res_random_effects_beta / ef_random_effects
        res_random_effects_ci.ub <- res_random_effects_beta * ef_random_effects
      }
    }

    measure_name <- if(log.scale) paste("log(",measure,")", sep = "") else measure
    refline <- if(is.element(measure,c("RD")) || log.scale == TRUE) 0 else 1

    par(mar=c(2,2,2,2))
    if(!missing(slab)){metafor::forest(x = results$yi, ci.lb = results$ci.lb, ci.ub = results$ci.ub, header = c("Study",paste(measure_name,"[95% CI]")), top=3,ylim=c(-3, nrow(results)+3), slab = slab_vals, refline = refline)}
    else {metafor::forest(x = results$yi, ci.lb = results$ci.lb, ci.ub = results$ci.ub, header = c("Study",paste(measure_name,"[95% CI]")), top=3,ylim=c(-3, nrow(results)+3), refline = refline)}
    metafor::addpoly(res, row = -1, mlab="Causal meta-analysis")
    if(measure!="SR"){
    if(is.element(measure,c("OR","RR")) && log.scale == FALSE){metafor::addpoly(x = res_random_effects_beta, ci.lb = res_random_effects_ci.lb, ci.ub = res_random_effects_ci.ub, row = -2, mlab="Random-effects model")}
    else {metafor::addpoly(res_random_effects, row = -2, mlab="Random-effects model")}
    }

  }

  return(list(study_results = select(results,c(study,yi,sei)), final_result = result))

}



.validate_variables <- function(names,ai,bi,ci,di,n1i,n2i,has_n,has_t){

  if(has_n){ vector_variables <- list(ai = ai, ci = ci, n1i= n1i, n2i = n2i) }
  else if(has_t){ vector_variables <- list(ai = ai, bi = bi, ci = ci, di = di) }

  # Check if all variables are present in dataset/matrix
  not_in_names <- setdiff(vector_variables,names)
  if(length(not_in_names) > 0L){
    stop("Variable(s) not present in 'data': ",paste(not_in_names, collapse = ", "), call. = FALSE)
  }

  invisible(TRUE)

}


.validate_vector_inputs <- function(ai,bi,ci,di,n1i,n2i,has_n,has_t){

  if(has_n){ vector_inputs <- list(ai = ai, ci = ci,n1i= n1i, n2i = n2i) }
  else if(has_t){ vector_inputs <- list(ai = ai, bi = bi, ci = ci, di = di) }

  # Check if numeric
  non_numeric <- names(vector_inputs)[!vapply(vector_inputs, is.numeric, logical(1L))]
  if (length(non_numeric) > 0L) {
    stop("Non-numeric inputs: ", paste(non_numeric, collapse = ", "), call. = FALSE)
  }

  # Check if integer
  non_integer <- names(vector_inputs)[!vapply(vector_inputs, function(x) all(x%%1==0), logical(1L))]
  if (length(non_integer) > 0L) {
    stop("Non-integer inputs: ", paste(non_integer, collapse = ", "), call. = FALSE)
  }

  # Check if positive or null
  non_positive <- names(vector_inputs)[!vapply(vector_inputs, function(x) all(x>=0), logical(1L))]
  if (length(non_positive) > 0L) {
    stop("Negative inputs: ", paste(non_positive, collapse = ", "), call. = FALSE)
  }

  # Check equal length vectors
  lengths <- vapply(vector_inputs, length, integer(1L))
  if (length(unique(lengths)) > 1L) {
    stop("All input vectors must have the same length", call. = FALSE)
  }

  # Check missing values
  missing <- names(vector_inputs)[vapply(vector_inputs, anyNA, logical(1L))]
  if (length(missing) > 0L) {
    stop("Missing data in inputs: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  if(has_n){ return(list(ai_vals = ai, ci_vals = ci, n1i_vals = n1i, n2i_vals = n2i)) }
  else if(has_t){
    n1i <- ai + bi
    n2i <- ci + di
    return(list(ai_vals = ai, ci_vals = ci, n1i_vals = n1i, n2i_vals = n2i))
  }
}

