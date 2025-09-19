#' Causal Meta Analysis for Aggregated Data
#'
#' @description
#' Function to perform causal meta-analysis on aggregated data of trials with binary treatment and binary outcome.
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
#' @param random.effects set to TRUE to have random effects meta-analysis estimate (along with causal meta-analysis estimate) both in the output and in the plot, default is TRUE (possible with all measures except SR)
#' @param ...
#'
#' @details The function can be used by either inputting "ai,bi,ci,di" or "ai,ci,n1i,n2i". Accepted measures are: "RD" = Risk Difference, "RR" = Risk Ratio, "OR" = Odds Ratio, "SR" = Survival Ratio
#'
#' @return A list with the following elements:
#' \describe{
#' \item{study_results}{list of effect sizes or outcome measures for each study and associated standard errors}
#' \item{final_result}{aggregated effect size or outcome measure and associated standard error}
#' \item{random_effects_model}{random-effect model estimate and standard error, if "random.effects = TRUE"}
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
#' result <- camea(measure = "RD", ai = treated_events, n1i = treated_total,
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
#' result <- camea(measure = "RR", ai = treated_positives, bi = treated_negatives,
#'                     ci = control_positives, di = control_negatives)
#'
#' @author Clement Berenfeld, Ahmed Boughdiri, Julie Josse, Charif El Gataa.
#'
#' \strong{Maintainer}: Clement Berenfeld <clement.berenfeld@inria.fr>
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
#' @importFrom dplyr select
#'
#' @export
#'
#'
camea <- function(measure, ai, bi, ci, di, n1i, n2i, slab, data = NULL, weights = NULL, plot = FALSE, log.scale = FALSE, random.effects = TRUE, ...){

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

  # Stop if measure is Risk Difference and log scale is TRUE
  if(measure == "RD" && log.scale == TRUE){stop("Risk Difference cannot be computed in log scale", call. = FALSE)}

  # Validate values of inputs and create list needed for next part of function
  if(is.null(data)){

    # Validate inputs
    dat <- validate_vector_inputs(ai = ai, bi = bi, ci = ci, di = di, n1i = n1i, n2i = n2i, has_n = has_n_totals, has_t = has_table_entries)
    if(!missing(slab)) {slab_vals <- slab}

  } else if(is.data.frame(data)){

    # Check: is variable in dataset?
    validate_variables(names = names(data), ai = deparse(substitute(ai)), bi = deparse(substitute(bi)),
                        ci = deparse(substitute(ci)), di = deparse(substitute(di)), n1i = deparse(substitute(n1i)),
                        n2i = deparse(substitute(n2i)), has_n = has_n_totals, has_t = has_table_entries)

    # Validate inputs
    dat <- validate_vector_inputs(ai = data[[deparse(substitute(ai))]], bi = data[[deparse(substitute(bi))]],
                                   ci = data[[deparse(substitute(ci))]], di = data[[deparse(substitute(di))]],
                                   n1i = data[[deparse(substitute(n1i))]], n2i = data[[deparse(substitute(n2i))]],
                                   has_n = has_n_totals, has_t = has_table_entries)
    if(!missing(slab)) {slab_vals <- data[[deparse(substitute(slab))]]}

  } else if(is.matrix(data) || inherits(data, "Matrix")){

    # Check: is variable in matrix?
    validate_variables(names = colnames(data), ai = deparse(substitute(ai)), bi = deparse(substitute(bi)),
                        ci = deparse(substitute(ci)), di = deparse(substitute(di)), n1i = deparse(substitute(n1i)),
                        n2i = deparse(substitute(n2i)), has_n = has_n_totals, has_t = has_table_entries)

    # Validate inputs
    dat <- validate_vector_inputs(ai = data[, deparse(substitute(ai))], bi = data[, deparse(substitute(bi))],
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
        thetai <- mu_1_k - mu_0_k
        vi <- (mu_0_k*(1- mu_0_k) / n1i_vals) + (mu_1_k*(1 - mu_1_k) / n2i_vals)
        ci.lb <- thetai - 1.96*sqrt(vi)
        ci.ub <- thetai + 1.96*sqrt(vi)
      }
      else if(measure == "RR"){
        ratio_i <- mu_1_k / mu_0_k
        vi_log <- (n1i_vals - ai_vals) / (n1i_vals * ai_vals) + (n2i_vals - ci_vals) / (n2i_vals * ci_vals)
        ci_log_i <- 1.96 * sqrt(vi_log)

        if(log.scale){ thetai <- log(ratio_i); vi <- vi_log; ci.lb <- thetai - ci_log_i; ci.ub <- thetai + ci_log_i}
        else if(!log.scale){ thetai <- ratio_i; vi <- vi_log*ratio_i^2; ef <- exp(ci_log_i); ci.lb <- thetai / ef; ci.ub <- thetai * ef}

      }
      else if(measure == "OR"){
        ratio_i <- (ai_vals * (n2i_vals - ci_vals)) / (ci_vals * (n1i_vals - ai_vals))
        vi_log <- (1/ai_vals) + (1/(n1i_vals - ai_vals)) + (1/ci_vals) + (1/(n2i_vals - ci_vals))
        ci_log_i <- 1.96 * sqrt(vi_log)

        if(log.scale){ thetai <- log(ratio_i); vi <- vi_log; ci.lb <- thetai - ci_log_i; ci.ub <- thetai + ci_log_i}
        else if(!log.scale){ thetai <- ratio_i; vi <- vi_log*ratio_i^2; ef <- exp(ci_log_i); ci.lb <- thetai / ef; ci.ub <- thetai * ef}

      }
      else if(measure == "SR"){
        ratio_i <- (1 - mu_1_k) / (1 - mu_0_k)
        vi_log <- (ai_vals / (n1i_vals * (n1i_vals - ai_vals))) + (ci_vals / (n2i_vals * (n2i_vals - ci_vals)))
        ci_log_i <- 1.96 * sqrt(vi_log)

        if(log.scale){ thetai <- log(ratio_i); vi <- vi_log; ci.lb <- thetai - ci_log_i; ci.ub <- thetai + ci_log_i}
        else if(!log.scale){ thetai <- ratio_i; vi <- vi_log*ratio_i^2; ef <- exp(ci_log_i); ci.lb <- thetai / ef; ci.ub <- thetai * ef}

      }

      tibble::tibble(
        proba_k = proba_k,
        mu_0_k = mu_0_k,
        mu_1_k = mu_1_k,
        thetai = thetai,
        sei = sqrt(vi),
        ci.lb = ci.lb,
        ci.ub = ci.ub,
        n1i = n1i_vals,
        n2i = n2i_vals,
        n.k = n1i_vals + n2i_vals
      )
  })

  if(!missing(slab)) {results$study <- slab_vals}
  else {results$study <- paste0("Study ", seq_along(results$thetai))}

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

  if(!is.na(final_result)){

  e_k <- results$n1i / results$n.k # propensity scores

  sigma2_0 <- sum(results$proba_k * results$mu_0_k^2) - EY0^2 + sum((results$n.k * results$mu_0_k * (1 - results$mu_0_k)) / ((1 - e_k) * n_total))

  sigma2_1 <- sum(results$proba_k * results$mu_1_k^2) - EY1^2 + sum((results$n.k * results$mu_1_k * (1 - results$mu_1_k)) / (e_k * n_total))

  switch(measure,
         "RD" = {
           d_psi_0 <- -1
           d_psi_1 <- 1
         },
         "RR" = {
           d_psi_0 <- - EY1 / EY0^2
           d_psi_1 <- 1 / EY0
         },
         "OR" = {
           d_psi_0 <- (- EY1) / ((1 - EY1) * EY0^2)
           d_psi_1 <- (1 - EY0) / ((1 - EY1)^2 * EY0)
         },
         "SR" = {
           d_psi_0 <-  (1 - EY1) / (1 - EY0)^2
           d_psi_1 <- - 1 / (1 - EY0)
         })

  variance <- (sigma2_0 * d_psi_0^2 + sigma2_1 * d_psi_1^2) / n_total

  if(log.scale){
    variance <- variance / final_result^2 # delta method
  }

  } else NA_real_

  if(log.scale){final_result <- log(final_result)} # Change it after computing the variance since we need ratios to compute variance

  result <- tibble::tibble(estimate = final_result, se = sqrt(variance))

  # Random effects model

  if(random.effects && measure != "SR"){
    dat <- metafor::escalc(measure=measure, ai=ai_vals, n1i=n1i_vals, ci=ci_vals, n2i=n2i_vals, data=dat)
    res_random_effects <- metafor::rma(dat$yi, dat$vi, method="DL")

    res_random_effects_beta <- if(is.element(measure,c("OR","RR")) && log.scale == FALSE) exp(as.numeric(res_random_effects$beta)) else as.numeric(res_random_effects$beta)
    res_random_effects_se <- if(is.element(measure,c("OR","RR")) && log.scale == FALSE) res_random_effects_beta*res_random_effects$se else res_random_effects$se
  }
  else if(random.effects && is.element(measure,c("SR"))){
    dat <- metafor::escalc(measure="RR", ai=ci_vals, n1i=n2i_vals, ci=ai_vals, n2i=n1i_vals, data=dat)
    res_random_effects <- metafor::rma(dat$yi, dat$vi, method="DL")

    res_random_effects_beta <- if(log.scale == FALSE) exp(as.numeric(res_random_effects$beta)) else as.numeric(res_random_effects$beta)
    res_random_effects_se <- if(log.scale == FALSE) res_random_effects_beta*res_random_effects$se else res_random_effects$se
  }

  if(plot){

    res <- metafor::rma(yi = result$estimate, sei = result$se, method="FE")

    if(log.scale == FALSE && random.effects){ # is.element(measure,c("OR","RR")) &&
      #res_random_effects_beta <- exp(res_random_effects$beta) # estimate of the measure with random effects model
      ef_random_effects <- exp(1.96 * res_random_effects$se) # error factor
      # upper and lower 95% confidence bounds for the estimate
      res_random_effects_ci.lb <- res_random_effects_beta / ef_random_effects
      res_random_effects_ci.ub <- res_random_effects_beta * ef_random_effects
    }

    measure_name <- if(log.scale) paste("log(",measure,")", sep = "") else measure
    refline <- if(is.element(measure,c("RD")) || log.scale == TRUE) 0 else 1

    par(mar=c(2,2,4,2))
    if(!missing(slab)){
      metafor::forest(x = results$thetai,
                      ci.lb = results$ci.lb,
                      ci.ub = results$ci.ub,
                      header = c("Study",paste(measure_name,"[95% CI]")),
                      top=2,
                      ylim=c(-2, n_study+2),
                      refline = refline,
                      cex = 0.8,
                      slab = slab_vals)
      }
    else {
      metafor::forest(x = results$thetai,
                      ci.lb = results$ci.lb,
                      ci.ub = results$ci.ub,
                      header = c("Study",paste(measure_name,"[95% CI]")),
                      top=2,
                      ylim=c(-2, n_study+2),
                      refline = refline,
                      cex = 0.8)
    }

    metafor::addpoly(res, row = 0, mlab="Causal meta-analysis", cex = 0.8)

    if(random.effects){ # measure!="SR" &&
      if(is.element(measure,c("OR","RR")) && log.scale == FALSE){
        metafor::addpoly(x = res_random_effects_beta,
                         ci.lb = res_random_effects_ci.lb,
                         ci.ub = res_random_effects_ci.ub,
                         row = -1,
                         mlab="Random-effects model",
                         cex = 0.8)
        }
    else {
      metafor::addpoly(res_random_effects, row = -1, mlab="Random-effects model", cex = 0.8)}
    }

  }

  if(random.effects) { # && measure != "SR"
    random_effects_result <- tibble::tibble(estimate = res_random_effects_beta, se = res_random_effects_se)
    return(list(study_results = dplyr::select(results,c(study,thetai,sei)), final_result = result, random_effects_model = random_effects_result))
  }
  else {return(list(study_results = dplyr::select(results,c(study,thetai,sei)), final_result = result))}

}



validate_variables <- function(names,ai,bi,ci,di,n1i,n2i,has_n,has_t){

  if(has_n){ vector_variables <- list(ai = ai, ci = ci, n1i= n1i, n2i = n2i) }
  else if(has_t){ vector_variables <- list(ai = ai, bi = bi, ci = ci, di = di) }

  # Check if all variables are present in dataset/matrix
  not_in_names <- setdiff(vector_variables,names)
  if(length(not_in_names) > 0L){
    stop("Variable(s) not present in 'data': ",paste(not_in_names, collapse = ", "), call. = FALSE)
  }

  invisible(TRUE)

}


validate_vector_inputs <- function(ai,bi,ci,di,n1i,n2i,has_n,has_t){

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

