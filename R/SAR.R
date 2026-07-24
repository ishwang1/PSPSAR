SAR <- local({
  prepare_X <- function(X, T_periods, N_units) {
    if (is.data.frame(X)) {
      X <- data.matrix(X)
    }
    if (is.null(dim(X))) {
      X <- matrix(as.numeric(X), nrow = N_units, ncol = 1)
    }
    if (length(dim(X)) == 2L) {
      return(as.matrix(X))
    }
    X
  }

  subset_X <- function(X, keep_units) {
    if (length(dim(X)) == 2L) {
      return(X[keep_units, , drop = FALSE])
    }
    X[, keep_units, , drop = FALSE]
  }

  get_X_cov_t <- function(X, t) {
    if (length(dim(X)) == 2L) {
      return(X)
    }
    X_t <- X[t, , , drop = FALSE]
    dim(X_t) <- dim(X)[2:3]
    X_t
  }

  get_covariate_names <- function(X) {
    if (length(dim(X)) == 2L) {
      x_names <- colnames(X)
      p <- ncol(X)
    } else {
      x_names <- dimnames(X)[[3]]
      p <- dim(X)[3]
    }
    if (is.null(x_names)) {
      x_names <- paste0("x", seq_len(p))
    }
    x_names
  }

  augment_X <- function(X_t) {
    cbind("(Intercept)" = 1, X_t)
  }

  solve_checked <- function(a, b = NULL, step) {
    tryCatch(
      {
        if (is.null(b)) {
          solve(a)
        } else {
          solve(a, b)
        }
      },
      error = function(e) {
        if (!is.matrix(a) || nrow(a) != ncol(a)) {
          stop(
            sprintf("Matrix inversion failed at %s: %s", step, conditionMessage(e)),
            call. = FALSE
          )
        }
        a_ginv <- MASS::ginv(a)
        if (is.null(b)) {
          return(a_ginv)
        }
        result <- a_ginv %*% b
        if (is.matrix(b)) {
          return(result)
        }
        as.numeric(result)
      }
    )
  }

  is_ready_W <- function(W) {
    W <- as.matrix(W)
    tol <- sqrt(.Machine$double.eps)
    all(abs(diag(W)) <= tol) && all(abs(rowSums(W) - 1) <= tol)
  }

  prepare_W_data <- function(Y, X, W, preprocess_W = FALSE) {
    W <- as.matrix(W)
    if (isTRUE(preprocess_W)) {
      return(list(
        Y = Y,
        X = X,
        W = W,
        dropped_units = integer(0),
        W_processed = FALSE
      ))
    }

    if (is_ready_W(W)) {
      return(list(
        Y = Y,
        X = X,
        W = W,
        dropped_units = integer(0),
        W_processed = FALSE
      ))
    }

    diag(W) <- 0
    kept_units <- seq_len(nrow(W))
    dropped_units <- integer(0)

    repeat {
      row_sums <- rowSums(W)
      zero_rows <- which(abs(row_sums) <= sqrt(.Machine$double.eps))
      if (length(zero_rows) == 0L) {
        break
      }
      dropped_units <- c(dropped_units, kept_units[zero_rows])
      keep_local <- setdiff(seq_len(nrow(W)), zero_rows)
      kept_units <- kept_units[keep_local]
      W <- W[keep_local, keep_local, drop = FALSE]
      diag(W) <- 0
    }

    row_sums <- rowSums(W)
    if (any(row_sums <= 0)) {
      stop("W still contains zero rows after preprocessing.", call. = FALSE)
    }

    Y <- Y[, kept_units, drop = FALSE]
    X <- subset_X(X, kept_units)
    W <- W / row_sums
    W[!is.finite(W)] <- 0
    diag(W) <- 0

    list(
      Y = Y,
      X = X,
      W = W,
      dropped_units = as.integer(dropped_units),
      W_processed = TRUE
    )
  }

  resolve_hac_bandwidth <- function(score_matrix, method) {
    score_matrix <- as.matrix(score_matrix)
    max_lag <- max(0L, nrow(score_matrix) - 1L)
    max_bandwidth <- max(1, nrow(score_matrix) - 1L)

    if (identical(method, "Andrews")) {
      candidate <- tryCatch(
        suppressWarnings(as.numeric(sandwich::bwAndrews(
          score_matrix,
          kernel = "Quadratic Spectral",
          prewhite = 0
        ))),
        error = function(e) NA_real_
      )
      if (!is.finite(candidate) || candidate <= 0) {
        candidate <- max(1, nrow(score_matrix)^(1 / 5))
      }
      return(min(as.numeric(candidate), max_bandwidth))
    }

    candidate <- tryCatch(
      suppressWarnings(as.numeric(sandwich::bwNeweyWest(
        score_matrix,
        kernel = "Bartlett",
        prewhite = 0
      ))),
      error = function(e) NA_real_
    )
    if (!is.finite(candidate) || candidate < 0) {
      candidate <- nrow(score_matrix)^(1 / 3)
    }
    min(max_lag, max(0L, floor(candidate)))
  }

  estimate_omega <- function(score_matrix, method, bandwidth) {
    if (identical(method, "Andrews")) {
      return(suppressWarnings(sandwich::kernHAC(
        stats::lm(score_matrix ~ 1),
        bw = bandwidth,
        kernel = "Quadratic Spectral",
        prewhite = 0,
        adjust = FALSE,
        sandwich = FALSE
      )))
    }

    suppressWarnings(sandwich::NeweyWest(
      stats::lm(score_matrix ~ 1),
      lag = bandwidth,
      prewhite = FALSE,
      adjust = FALSE,
      sandwich = FALSE
    ))
  }

  build_time_objects <- function(X, W, T_periods, iv_lag, covariate_names) {
    X_cov <- vector("list", T_periods)
    X_aug <- vector("list", T_periods)
    instruments <- vector("list", T_periods)

    for (t in seq_len(T_periods)) {
      X_cov_t <- as.matrix(get_X_cov_t(X, t))
      colnames(X_cov_t) <- covariate_names
      X_aug_t <- augment_X(X_cov_t)
      instr_t <- X_aug_t

      if (iv_lag > 0L && ncol(X_cov_t) > 0L) {
        WX <- X_cov_t
        for (lag in seq_len(iv_lag)) {
          WX <- W %*% WX
          colnames(WX) <- paste0("W", lag, "_", covariate_names)
          instr_t <- cbind(instr_t, WX)
        }
      }

      X_cov[[t]] <- X_cov_t
      X_aug[[t]] <- X_aug_t
      instruments[[t]] <- instr_t
    }

    list(
      X_cov = X_cov,
      X_aug = X_aug,
      instruments = instruments
    )
  }

  compute_first_stage_diagnostics <- function(time_objects, WY) {
    T_periods <- nrow(WY)
    X_stack <- do.call(rbind, time_objects$X_aug)
    Z_stack <- do.call(rbind, time_objects$instruments)
    y_stack <- as.numeric(t(WY))

    restricted_fit <- stats::lm.fit(x = X_stack, y = y_stack)
    full_fit <- stats::lm.fit(x = Z_stack, y = y_stack)

    rss_restricted <- sum(restricted_fit$residuals^2)
    rss_full <- sum(full_fit$residuals^2)
    df_num <- full_fit$rank - restricted_fit$rank
    df_denom <- nrow(Z_stack) - full_fit$rank

    if (df_num <= 0L || df_denom <= 0L || !is.finite(rss_restricted) || !is.finite(rss_full)) {
      return(list(
        statistic = NA_real_,
        p_value = NA_real_,
        df_num = as.integer(df_num),
        df_denom = as.integer(df_denom),
        partial_r2 = NA_real_,
        method = "Pooled excluded-instrument F test for Wy on X and WX"
      ))
    }

    # Guard against tiny negative differences from floating-point noise.
    ss_excluded <- max(rss_restricted - rss_full, 0)
    mse_full <- rss_full / df_denom
    statistic <- if (mse_full > 0) (ss_excluded / df_num) / mse_full else NA_real_
    p_value <- if (is.finite(statistic)) {
      stats::pf(statistic, df_num, df_denom, lower.tail = FALSE)
    } else {
      NA_real_
    }
    partial_r2 <- if (rss_restricted > 0) ss_excluded / rss_restricted else NA_real_

    list(
      statistic = as.numeric(statistic),
      p_value = as.numeric(p_value),
      df_num = as.integer(df_num),
      df_denom = as.integer(df_denom),
      partial_r2 = as.numeric(partial_r2),
      method = "Pooled excluded-instrument F test for Wy on X and WX"
    )
  }

  compute_gmm_components <- function(time_objects, Y, WY, weight_matrix, step) {
    T_periods <- nrow(Y)
    N_units <- ncol(Y)
    delta_dim <- 1L + ncol(time_objects$X_aug[[1]])
    instrument_dim <- ncol(time_objects$instruments[[1]])

    eta <- numeric(instrument_dim)
    G <- matrix(0, instrument_dim, delta_dim)
    ZZ <- matrix(0, instrument_dim, instrument_dim)

    for (t in seq_len(T_periods)) {
      Z_t <- time_objects$instruments[[t]]
      R_t <- cbind(WY[t, ], time_objects$X_aug[[t]])
      eta <- eta + as.numeric(crossprod(Z_t, Y[t, ]))
      G <- G + crossprod(Z_t, R_t)
      ZZ <- ZZ + crossprod(Z_t)
    }

    denom <- N_units * T_periods
    eta <- eta / denom
    G <- G / denom
    ZZ <- ZZ / denom

    A <- if (identical(weight_matrix, "2SLS")) {
      solve_checked(ZZ, step = sprintf("%s weighting matrix", step))
    } else {
      weight_matrix
    }

    info <- crossprod(G, A %*% G)
    rhs <- crossprod(G, A %*% eta)
    delta <- as.numeric(solve_checked(
      info,
      rhs,
      step = sprintf("%s coefficient matrix", step)
    ))

    list(
      delta = delta,
      eta = eta,
      G = G,
      ZZ = ZZ,
      A = A,
      info = info
    )
  }

  compute_residuals <- function(time_objects, Y, WY, delta) {
    T_periods <- nrow(Y)
    residuals <- matrix(NA_real_, nrow(Y), ncol(Y))

    for (t in seq_len(T_periods)) {
      fitted_t <- delta[1] * WY[t, ] + as.numeric(time_objects$X_aug[[t]] %*% delta[-1])
      residuals[t, ] <- Y[t, ] - fitted_t
    }

    residuals
  }

  compute_score_matrix <- function(time_objects, residuals) {
    T_periods <- nrow(residuals)
    N_units <- ncol(residuals)
    instrument_dim <- ncol(time_objects$instruments[[1]])
    score <- matrix(NA_real_, T_periods, instrument_dim)

    for (t in seq_len(T_periods)) {
      score[t, ] <- as.numeric(crossprod(time_objects$instruments[[t]], residuals[t, ]) / sqrt(N_units))
    }

    score
  }

  finalize_output <- function(delta, vcov_delta, ci_level, beta_names) {
    z_alpha <- stats::qnorm((1 + ci_level) / 2)
    se <- sqrt(pmax(diag(vcov_delta), 0))
    lower <- delta - z_alpha * se
    upper <- delta + z_alpha * se
    statistic <- delta / se
    p <- 2 * stats::pnorm(abs(statistic), lower.tail = FALSE)

    beta_output <- data.frame(
      term = beta_names,
      estimate = delta[-1],
      se = se[-1],
      lower = lower[-1],
      upper = upper[-1],
      p = p[-1],
      stringsAsFactors = FALSE
    )
    # Keep covariates first and place the intercept last in the public table.
    beta_output <- beta_output[c(
      which(beta_output$term != "(Intercept)"),
      which(beta_output$term == "(Intercept)")
    ), , drop = FALSE]
    rownames(beta_output) <- NULL

    list(
      rho = data.frame(
        term = "rho",
        estimate = delta[1],
        se = se[1],
        lower = lower[1],
        upper = upper[1],
        p = p[1],
        stringsAsFactors = FALSE
      ),
      beta = beta_output
    )
  }

  function(
    Y, X, W,
    iv_lag = 1,
    estimator = c("2SLS", "GMM"),
    ci_level = 0.95,
    hac_method = c("Newey-West", "Andrews"),
    preprocess_W = FALSE
  ) {
    estimator <- match.arg(estimator)
    hac_method <- match.arg(hac_method)

    Y <- as.matrix(Y)
    ci_level <- as.numeric(ci_level[[1]])

    T_periods <- nrow(Y)
    N_units <- ncol(Y)
    X <- prepare_X(X, T_periods, N_units)

    prepared <- prepare_W_data(Y, X, W, preprocess_W = preprocess_W)
    Y <- prepared$Y
    X <- prepared$X
    W <- prepared$W
    dropped_units <- prepared$dropped_units
    W_processed <- prepared$W_processed

    N_units <- ncol(Y)
    covariate_names <- get_covariate_names(X)
    beta_names <- c("(Intercept)", covariate_names)
    time_objects <- build_time_objects(X, W, T_periods, iv_lag, covariate_names)
    WY <- Y %*% t(W)

    fit_2sls <- compute_gmm_components(
      time_objects = time_objects,
      Y = Y,
      WY = WY,
      weight_matrix = "2SLS",
      step = "2SLS"
    )
    first_stage <- compute_first_stage_diagnostics(
      time_objects = time_objects,
      WY = WY
    )

    final_fit <- fit_2sls
    final_residuals <- compute_residuals(time_objects, Y, WY, fit_2sls$delta)

    if (identical(estimator, "GMM")) {
      score_init <- compute_score_matrix(time_objects, final_residuals)
      hac_bandwidth <- resolve_hac_bandwidth(score_init, hac_method)
      omega_init <- estimate_omega(score_init, hac_method, hac_bandwidth)
      final_fit <- compute_gmm_components(
        time_objects = time_objects,
        Y = Y,
        WY = WY,
        weight_matrix = solve_checked(omega_init, step = "GMM feasible weighting matrix"),
        step = "GMM"
      )
      final_residuals <- compute_residuals(time_objects, Y, WY, final_fit$delta)
      score_final <- compute_score_matrix(time_objects, final_residuals)
    } else {
      score_final <- compute_score_matrix(time_objects, final_residuals)
      hac_bandwidth <- resolve_hac_bandwidth(score_final, hac_method)
    }

    if (!exists("hac_bandwidth", inherits = FALSE)) {
      hac_bandwidth <- resolve_hac_bandwidth(score_final, hac_method)
    }

    omega_final <- estimate_omega(score_final, hac_method, hac_bandwidth)
    info_inv <- solve_checked(
      final_fit$info,
      step = sprintf("%s sandwich factor", estimator)
    )
    V_delta <- info_inv %*% crossprod(
      final_fit$G,
      final_fit$A %*% omega_final %*% final_fit$A %*% final_fit$G
    ) %*% t(info_inv)
    V_delta <- V_delta / (N_units * T_periods)

    output <- finalize_output(
      delta = final_fit$delta,
      vcov_delta = V_delta,
      ci_level = ci_level,
      beta_names = beta_names
    )

    result <- list(
      beta = output$beta,
      rho = output$rho,
      residuals = final_residuals
    )
    if (isTRUE(W_processed)) {
      result$W <- W
    }
    result$diagnostics <- list(
      estimator = estimator,
      iv_lag = as.integer(iv_lag),
      first_stage = first_stage,
      hac_method = hac_method,
      hac_bandwidth = hac_bandwidth,
      W_processed = W_processed,
      dropped_units = dropped_units,
      N_used = N_units
    )
    result
  }
})
