#' Compare Zero-Inflation Models for Virome Data
#'
#' Performs a statistical comparison of fixed versus variable zero-inflation models
#' for virome data. This function helps researchers determine which approach better
#' fits their data and understand the nature of zero-inflation in their dataset.
#'
#' @param count_matrix A numeric matrix of viral counts with taxa as rows and samples as columns
#' @param presence_threshold Count threshold to consider a taxon as present (default: 1)
#' @param min_prevalence Minimum prevalence (proportion of samples) to include a taxon (default: 0.05)
#' @param variable_zi_method Method for fitting variable zero-inflation: "beta", "power_law", or "mixed" (default: "beta")
#' @param n_bootstrap Number of bootstrap samples for computing confidence intervals (default: 100)
#' @param visualize Whether to create comparison plots (default: TRUE)
#' @param detailed_output Whether to return detailed model fit results (default: FALSE)
#'
#' @return A list containing:
#'   \item{best_model}{Character string indicating the better-fitting model ("fixed" or "variable")}
#'   \item{likelihood_ratio}{Likelihood ratio between the models}
#'   \item{p_value}{P-value from likelihood ratio test (if applicable)}
#'   \item{aic_values}{AIC values for both models}
#'   \item{bic_values}{BIC values for both models}
#'   \item{variable_zi_params}{Estimated parameters for the variable zero-inflation model}
#'   \item{fixed_zi_params}{Estimated parameters for the fixed zero-inflation model}
#'   \item{plots}{List of comparison plots (if visualize=TRUE)}
#'   \item{detailed_results}{Detailed model fit results (if detailed_output=TRUE)}
#'
#' @details
#' This function performs a formal statistical comparison between two approaches to
#' modeling zero-inflation in virome data:
#' 
#' 1. Fixed zero-inflation: All viral taxa have the same structural zero probability
#' 2. Variable zero-inflation: Structural zero probabilities vary across taxa according to
#'    a parametric distribution (e.g., Beta distribution)
#' 
#' The function:
#' - Fits both models to the observed data
#' - Compares them using likelihood ratio tests and information criteria (AIC/BIC)
#' - Provides parameter estimates and uncertainty measures
#' - Creates visualizations to aid interpretation
#'
#' For variable zero-inflation, several distribution models can be used:
#' \itemize{
#'   \item "beta": Beta distribution (flexible shape, bounded between 0 and 1)
#'   \item "power_law": Power law relationship between abundance and zero-inflation
#'   \item "mixed": Mixture of Beta distributions (for multimodal patterns)
#' }
#'
#' @examples
#' \dontrun{
#' # Basic model comparison
#' model_comparison <- compare_zi_models(
#'   count_matrix = my_virome_counts,
#'   presence_threshold = 5
#' )
#' 
#' # Print summary of results
#' cat("Better model:", model_comparison$best_model, "\n")
#' cat("AIC difference:", model_comparison$aic_values$variable - model_comparison$aic_values$fixed, "\n")
#' 
#' # Plot comparison
#' print(model_comparison$plots$observed_vs_predicted)
#' 
#' # Detailed comparison with custom settings
#' detailed_comparison <- compare_zi_models(
#'   count_matrix = my_virome_counts,
#'   presence_threshold = 10,
#'   variable_zi_method = "mixed",
#'   n_bootstrap = 200,
#'   detailed_output = TRUE
#' )
#' }
#'
#' @export
compare_zi_models <- function(count_matrix,
                             presence_threshold = 1,
                             min_prevalence = 0.05,
                             variable_zi_method = "beta",
                             n_bootstrap = 100,
                             visualize = TRUE,
                             detailed_output = FALSE) {
  
  # Check if required packages are available
  req_pkgs <- c("fitdistrplus", "MASS")
  missing_pkgs <- req_pkgs[!sapply(req_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Required packages not available: ", paste(missing_pkgs, collapse = ", "), 
         ". Please install with install.packages(c('", 
         paste(missing_pkgs, collapse = "', '"), "'))")
  }
  
  if (visualize && !requireNamespace("ggplot2", quietly = TRUE)) {
    warning("Package 'ggplot2' is required for visualization. No plots will be generated.")
    visualize <- FALSE
  }
  
  # Convert to presence/absence matrix
  presence_absence <- count_matrix >= presence_threshold
  
  # Calculate prevalence of each taxon
  n_samples <- ncol(count_matrix)
  prevalence <- rowSums(presence_absence) / n_samples
  
  # Calculate zero-inflation rates (1 - prevalence)
  zero_inflation <- 1 - prevalence
  
  # Filter out rare taxa (those with very low prevalence)
  valid_taxa <- prevalence >= min_prevalence & prevalence <= (1 - min_prevalence)
  
  if (sum(valid_taxa) < 10) {
    warning("Very few taxa with sufficient prevalence for model comparison (n=", sum(valid_taxa), 
            "). Results may be unreliable.")
  }
  
  filtered_zi <- zero_inflation[valid_taxa]
  filtered_counts <- count_matrix[valid_taxa, , drop = FALSE]
  
  # Ensure filtered values are bounded away from 0 and 1 for Beta fitting
  filtered_zi <- pmin(pmax(filtered_zi, 0.001), 0.999)
  
  # Calculate abundance metrics for each taxon
  abundance <- rowMeans(filtered_counts)
  
  # 1. Fit fixed zero-inflation model
  fixed_zi_rate <- mean(filtered_zi)
  
  # Calculate log-likelihood for fixed ZI model (Beta with very high concentration)
  fixed_loglik <- sum(dbeta(filtered_zi, 
                           shape1 = 1000 * fixed_zi_rate, 
                           shape2 = 1000 * (1 - fixed_zi_rate), 
                           log = TRUE))
  
  # Number of parameters for fixed model: 1 (just the mean)
  fixed_params <- 1
  
  # 2. Fit variable zero-inflation model based on specified method
  if (variable_zi_method == "beta") {
    # Fit Beta distribution to ZI rates
    beta_fit <- try(fitdistrplus::fitdist(filtered_zi, "beta", method = "mle"), silent = TRUE)
    
    if (inherits(beta_fit, "try-error")) {
      warning("Beta distribution fitting failed")
      variable_zi_params <- list(
        method = "beta",
        fit_failed = TRUE,
        alpha = NA,
        beta = NA
      )
      variable_loglik <- -Inf
      variable_params <- 2  # alpha and beta
    } else {
      alpha <- beta_fit$estimate["shape1"]
      beta <- beta_fit$estimate["shape2"]
      
      variable_zi_params <- list(
        method = "beta",
        alpha = alpha,
        beta = beta,
        mean = alpha / (alpha + beta),
        variance = (alpha * beta) / ((alpha + beta)^2 * (alpha + beta + 1)),
        fit_details = beta_fit
      )
      
      variable_loglik <- beta_fit$loglik
      variable_params <- 2  # alpha and beta
    }
    
  } else if (variable_zi_method == "power_law") {
    # Fit power law model: ZI ~ a * abundance^b
    log_abundance <- log(abundance + 1)  # Add 1 to handle zeros
    log_zi <- log(filtered_zi)
    
    # Linear regression on log scale
    power_model <- lm(log_zi ~ log_abundance)
    
    # Calculate predicted values on original scale
    a <- exp(coef(power_model)[1])
    b <- coef(power_model)[2]
    predicted_zi <- a * (abundance + 1)^b
    
    # Bound predictions to valid range for ZI
    predicted_zi <- pmin(pmax(predicted_zi, 0.001), 0.999)
    
    # Calculate log-likelihood assuming Beta distribution around predicted values
    # We use a high concentration to model values close to predictions
    concentration <- 50  # Higher value = tighter distribution around prediction
    alpha_pred <- concentration * predicted_zi
    beta_pred <- concentration * (1 - predicted_zi)
    
    variable_loglik <- sum(dbeta(filtered_zi, 
                               shape1 = alpha_pred, 
                               shape2 = beta_pred, 
                               log = TRUE))
    
    variable_zi_params <- list(
      method = "power_law",
      a = a,
      b = b,
      r_squared = summary(power_model)$r.squared,
      fit_details = power_model
    )
    
    variable_params <- 2  # a and b parameters
    
  } else if (variable_zi_method == "mixed") {
    # Fit mixture of two Beta distributions (for bimodal data)
    if (!requireNamespace("mixtools", quietly = TRUE)) {
      warning("Package 'mixtools' is required for mixed model. Falling back to beta model.")
      return(compare_zi_models(
        count_matrix, presence_threshold, min_prevalence, 
        variable_zi_method = "beta", n_bootstrap, visualize, detailed_output
      ))
    }
    
    # Try to fit a 2-component beta mixture model
    beta_mix <- try(mixtools::betamix(filtered_zi, k = 2), silent = TRUE)
    
    if (inherits(beta_mix, "try-error")) {
      warning("Beta mixture model fitting failed")
      variable_zi_params <- list(
        method = "mixed",
        fit_failed = TRUE
      )
      variable_loglik <- -Inf
      variable_params <- 5  # 2 alphas, 2 betas, and mixture weight
    } else {
      # Extract parameters
      alpha1 <- beta_mix$alpha[1]
      beta1 <- beta_mix$beta[1]
      alpha2 <- beta_mix$alpha[2]
      beta2 <- beta_mix$beta[2]
      lambda <- beta_mix$lambda
      
      # Calculate log-likelihood
      variable_loglik <- sum(log(
        lambda[1] * dbeta(filtered_zi, shape1 = alpha1, shape2 = beta1) +
          lambda[2] * dbeta(filtered_zi, shape1 = alpha2, shape2 = beta2)
      ))
      
      variable_zi_params <- list(
        method = "mixed",
        alpha1 = alpha1,
        beta1 = beta1,
        alpha2 = alpha2,
        beta2 = beta2,
        lambda = lambda,
        mean1 = alpha1 / (alpha1 + beta1),
        mean2 = alpha2 / (alpha2 + beta2),
        overall_mean = lambda[1] * alpha1 / (alpha1 + beta1) + 
                       lambda[2] * alpha2 / (alpha2 + beta2),
        fit_details = beta_mix
      )
      
      variable_params <- 5  # 2 alphas, 2 betas, and mixture weight
    }
  } else {
    stop("Invalid variable_zi_method. Use 'beta', 'power_law', or 'mixed'.")
  }
  
  # 3. Calculate information criteria
  aic_fixed <- -2 * fixed_loglik + 2 * fixed_params
  bic_fixed <- -2 * fixed_loglik + log(length(filtered_zi)) * fixed_params
  
  aic_variable <- -2 * variable_loglik + 2 * variable_params
  bic_variable <- -2 * variable_loglik + log(length(filtered_zi)) * variable_params
  
  # 4. Likelihood ratio test (only valid for nested models)
  # These models are not strictly nested, but the test can still be informative
  lrt_statistic <- 2 * (variable_loglik - fixed_loglik)
  lrt_df <- variable_params - fixed_params
  lrt_p_value <- pchisq(lrt_statistic, df = lrt_df, lower.tail = FALSE)
  
  # 5. Bootstrap confidence intervals for parameters
  if (n_bootstrap > 0) {
    bootstrap_results <- matrix(NA, nrow = n_bootstrap, ncol = variable_params)
    colnames(bootstrap_results) <- paste0("param", 1:variable_params)
    
    for (i in 1:n_bootstrap) {
      # Bootstrap resample
      bootstrap_indices <- sample(length(filtered_zi), replace = TRUE)
      bootstrap_zi <- filtered_zi[bootstrap_indices]
      
      # Fit model to bootstrap sample
      if (variable_zi_method == "beta") {
        boot_fit <- try(fitdistrplus::fitdist(bootstrap_zi, "beta", method = "mle"), silent = TRUE)
        if (!inherits(boot_fit, "try-error")) {
          bootstrap_results[i, 1] <- boot_fit$estimate["shape1"]  # alpha
          bootstrap_results[i, 2] <- boot_fit$estimate["shape2"]  # beta
        }
      } else if (variable_zi_method == "power_law") {
        bootstrap_abundance <- abundance[bootstrap_indices]
        log_bootstrap_abundance <- log(bootstrap_abundance + 1)
        log_bootstrap_zi <- log(bootstrap_zi)
        
        boot_model <- try(lm(log_bootstrap_zi ~ log_bootstrap_abundance), silent = TRUE)
        if (!inherits(boot_model, "try-error")) {
          bootstrap_results[i, 1] <- exp(coef(boot_model)[1])  # a
          bootstrap_results[i, 2] <- coef(boot_model)[2]       # b
        }
      } else if (variable_zi_method == "mixed") {
        # Mixture models are complex to bootstrap - simplified approach
        # This is a placeholder; a full implementation would be more involved
        bootstrap_results[i, ] <- NA
      }
    }
    
    # Calculate bootstrap confidence intervals (excluding NAs)
    bootstrap_ci <- apply(bootstrap_results, 2, function(x) {
      if (all(is.na(x))) {
        return(c(NA, NA))
      }
      quantile(x, c(0.025, 0.975), na.rm = TRUE)
    })
    
    # Add to parameters
    if (variable_zi_method == "beta") {
      variable_zi_params$alpha_ci <- bootstrap_ci[, 1]
      variable_zi_params$beta_ci <- bootstrap_ci[, 2]
    } else if (variable_zi_method == "power_law") {
      variable_zi_params$a_ci <- bootstrap_ci[, 1]
      variable_zi_params$b_ci <- bootstrap_ci[, 2]
    }
  }
  
  # 6. Create visualizations if requested
  if (visualize) {
    plots <- list()
    
    # Generate histogram of observed ZI rates
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      # Basic histogram
      p_hist <- ggplot2::ggplot(data.frame(zi = filtered_zi), ggplot2::aes(x = zi)) +
        ggplot2::geom_histogram(bins = 30, fill = "steelblue", color = "white", alpha = 0.7) +
        ggplot2::geom_vline(xintercept = fixed_zi_rate, color = "firebrick", 
                            linetype = "dashed", size = 1) +
        ggplot2::labs(
          title = "Distribution of Zero-Inflation Rates",
          subtitle = sprintf("Fixed ZI rate: %.3f", fixed_zi_rate),
          x = "Zero-Inflation Rate",
          y = "Count"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold"),
          plot.subtitle = ggplot2::element_text(face = "italic")
        )
      plots$histogram <- p_hist
      
      # Observed vs. model predictions
      x_seq <- seq(0.01, 0.99, length.out = 100)
      
      # Fixed model density
      fixed_density <- dbeta(x_seq, 
                            shape1 = 1000 * fixed_zi_rate, 
                            shape2 = 1000 * (1 - fixed_zi_rate))
      
      # Variable model density
      if (variable_zi_method == "beta" && !is.na(variable_zi_params$alpha)) {
        variable_density <- dbeta(x_seq, 
                                shape1 = variable_zi_params$alpha, 
                                shape2 = variable_zi_params$beta)
        
        density_df <- data.frame(
          x = rep(x_seq, 2),
          density = c(fixed_density, variable_density),
          model = rep(c("Fixed", "Variable"), each = length(x_seq))
        )
        
        p_density <- ggplot2::ggplot() +
          ggplot2::geom_histogram(
            data = data.frame(zi = filtered_zi),
            ggplot2::aes(x = zi, y = ..density..),
            bins = 20, fill = "gray80", color = "white", alpha = 0.5
          ) +
          ggplot2::geom_line(
            data = density_df,
            ggplot2::aes(x = x, y = density, color = model, linetype = model),
            size = 1.2
          ) +
          ggplot2::scale_color_manual(values = c("Fixed" = "firebrick", "Variable" = "steelblue")) +
          ggplot2::labs(
            title = "Model Comparison",
            subtitle = sprintf("AIC: Fixed = %.1f, Variable = %.1f", aic_fixed, aic_variable),
            x = "Zero-Inflation Rate",
            y = "Density",
            color = "Model",
            linetype = "Model"
          ) +
          ggplot2::theme_minimal() +
          ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold"),
            plot.subtitle = ggplot2::element_text(face = "italic"),
            legend.position = "bottom"
          )
        plots$model_comparison <- p_density
      }
      
      # Relationship with abundance (scatter plot)
      relationship_df <- data.frame(
        abundance = abundance,
        zi_rate = filtered_zi
      )
      
      p_relationship <- ggplot2::ggplot(relationship_df, 
                                       ggplot2::aes(x = abundance, y = zi_rate)) +
        ggplot2::geom_point(alpha = 0.7, color = "steelblue") +
        ggplot2::geom_hline(yintercept = fixed_zi_rate, 
                           color = "firebrick", linetype = "dashed") +
        ggplot2::labs(
          title = "Relationship Between Abundance and Zero-Inflation",
          x = "Mean Abundance",
          y = "Zero-Inflation Rate"
        ) +
        ggplot2::scale_x_log10() +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold")
        )
      
      # Add power law fit if applicable
      if (variable_zi_method == "power_law" && !is.na(variable_zi_params$a)) {
        # Generate prediction curve
        x_curve <- exp(seq(log(min(abundance + 1)), log(max(abundance + 1)), length.out = 100)) - 1
        y_curve <- variable_zi_params$a * (x_curve + 1) ^ variable_zi_params$b
        
        p_relationship <- p_relationship +
          ggplot2::geom_line(
            data = data.frame(x = x_curve, y = y_curve),
            ggplot2::aes(x = x, y = y),
            color = "darkred", size = 1.2
          ) +
          ggplot2::labs(
            subtitle = sprintf("Power law fit: ZI = %.3f × abundance^%.3f, R² = %.2f",
                             variable_zi_params$a, variable_zi_params$b, 
                             variable_zi_params$r_squared)
          )
      }
      
      plots$abundance_relationship <- p_relationship
    }
  } else {
    plots <- NULL
  }
  
  # 7. Determine the best model based on AIC
  aic_diff <- aic_fixed - aic_variable
  bic_diff <- bic_fixed - bic_variable
  
  if (aic_diff > 2) {
    best_model <- "variable"
  } else if (aic_diff < -2) {
    best_model <- "fixed"
  } else {
    best_model <- "inconclusive"
  }
  
  # 8. Prepare results
  fixed_zi_params <- list(
    fixed_rate = fixed_zi_rate,
    n_observations = length(filtered_zi)
  )
  
  # Add bootstrapped CI for fixed rate if available
  if (n_bootstrap > 0) {
    boot_means <- replicate(n_bootstrap, mean(sample(filtered_zi, replace = TRUE)))
    fixed_zi_params$ci <- quantile(boot_means, c(0.025, 0.975))
  }
  
  results <- list(
    best_model = best_model,
    likelihood_ratio = exp(variable_loglik - fixed_loglik),
    lrt_statistic = lrt_statistic,
    lrt_p_value = lrt_p_value,
    aic_values = list(fixed = aic_fixed, variable = aic_variable, difference = aic_diff),
    bic_values = list(fixed = bic_fixed, variable = bic_variable, difference = bic_diff),
    variable_zi_params = variable_zi_params,
    fixed_zi_params = fixed_zi_params,
    plots = plots
  )
  
  if (detailed_output) {
    # Add detailed model results
    results$detailed_results <- list(
      filtered_zi = filtered_zi,
      filtered_counts = filtered_counts,
      abundance = abundance,
      n_valid_taxa = sum(valid_taxa),
      n_total_taxa = length(zero_inflation),
      fixed_loglik = fixed_loglik,
      variable_loglik = variable_loglik
    )
    
    if (n_bootstrap > 0) {
      results$detailed_results$bootstrap_results <- bootstrap_results
    }
  }
  
  return(results)
}