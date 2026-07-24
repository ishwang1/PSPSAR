PSPSAR <- local({
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

  robust_scale_z <- function(z) {
    z_sd <- stats::sd(z)
    z_iqr <- stats::IQR(z) / 1.34
    scale <- min(z_sd, z_iqr, na.rm = TRUE)
    if (!is.finite(scale) || scale <= 0) {
      scale <- max(z_sd, z_iqr, na.rm = TRUE)
    }
    if (!is.finite(scale) || scale <= 0) {
      scale <- 1
    }
    scale
  }

  kernel_weight <- function(u, kernel) {
    if (identical(kernel, "gaussian")) {
      return(stats::dnorm(u))
    }
    if (identical(kernel, "epanechnikov")) {
      return(ifelse(abs(u) <= 1, 0.75 * (1 - u^2), 0))
    }
    if (identical(kernel, "triangular")) {
      return(ifelse(abs(u) <= 1, 1 - abs(u), 0))
    }
    if (identical(kernel, "uniform")) {
      return(ifelse(abs(u) <= 1, 0.5, 0))
    }
    ifelse(abs(u) <= 1, 15 / 16 * (1 - u^2)^2, 0)
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

  build_time_objects <- function(X, W, T_periods, iv_lag, covariate_names) {
    X_cov_list <- vector("list", T_periods)
    X_aug_list <- vector("list", T_periods)
    base_instr_list <- vector("list", T_periods)

    for (t in seq_len(T_periods)) {
      X_cov_t <- as.matrix(get_X_cov_t(X, t))
      colnames(X_cov_t) <- covariate_names
      X_aug_t <- augment_X(X_cov_t)
      base_instr_t <- X_aug_t

      if (iv_lag > 0L && ncol(X_cov_t) > 0L) {
        WX <- X_cov_t
        for (lag in seq_len(iv_lag)) {
          WX <- W %*% WX
          colnames(WX) <- paste0("W", lag, "_", covariate_names)
          base_instr_t <- cbind(base_instr_t, WX)
        }
      }

      X_cov_list[[t]] <- X_cov_t
      X_aug_list[[t]] <- X_aug_t
      base_instr_list[[t]] <- base_instr_t
    }

    list(
      X_cov = X_cov_list,
      X_aug = X_aug_list,
      base_instr = base_instr_list
    )
  }

  build_Q_H <- function(base_instr_t, WY_t, u, order) {
    basis <- u^(0:order)
    Q_t <- do.call(cbind, lapply(basis, function(power) base_instr_t * power))
    H_t <- do.call(cbind, lapply(basis, function(power) matrix(WY_t * power, ncol = 1)))

    q_names <- colnames(base_instr_t)
    colnames(Q_t) <- unlist(lapply(seq_along(basis), function(idx) {
      suffix <- if (idx == 1L) "" else paste0("_u", idx - 1L)
      paste0(q_names, suffix)
    }), use.names = FALSE)
    colnames(H_t) <- paste0("theta", seq_len(order + 1L))

    list(Q = Q_t, H = H_t)
  }

  compute_local_moments <- function(time_objects, Y, WY, z, z0, bandwidth, kernel, order) {
    T_periods <- nrow(Y)
    N_units <- ncol(Y)
    k <- ncol(time_objects$X_aug[[1]])
    q0 <- ncol(time_objects$base_instr[[1]])
    q <- q0 * (order + 1L)
    m <- order + 1L

    eta <- numeric(q)
    Xi <- matrix(0, q, k)
    Gamma <- matrix(0, q, m)
    QQ <- matrix(0, q, q)
    components <- vector("list", T_periods)
    denom <- N_units * T_periods * bandwidth

    for (t in seq_len(T_periods)) {
      u <- (z[t] - z0) / bandwidth
      k_t <- kernel_weight(u, kernel)
      qh <- build_Q_H(time_objects$base_instr[[t]], WY[t, ], u, order)
      Q_t <- qh$Q
      H_t <- qh$H

      eta <- eta + k_t * as.numeric(crossprod(Q_t, Y[t, ]))
      Xi <- Xi + k_t * crossprod(Q_t, time_objects$X_aug[[t]])
      Gamma <- Gamma + k_t * crossprod(Q_t, H_t)
      QQ <- QQ + k_t * crossprod(Q_t)

      components[[t]] <- list(
        weight = k_t,
        Q = Q_t,
        H = H_t,
        u = u
      )
    }

    list(
      eta = eta / denom,
      Xi = Xi / denom,
      Gamma = Gamma / denom,
      QQ = QQ / denom,
      components = components
    )
  }

  build_local_operator <- function(local_moments, A, step) {
    info_matrix <- crossprod(local_moments$Gamma, A %*% local_moments$Gamma)
    info_inv <- solve_checked(info_matrix, step = step)
    P_mat <- info_inv %*% crossprod(local_moments$Gamma, A)
    gamma_y <- as.numeric(P_mat[1, , drop = FALSE] %*% local_moments$eta)
    gamma_X <- as.numeric(P_mat[1, , drop = FALSE] %*% local_moments$Xi)

    list(
      A = A,
      P = P_mat,
      info_inv = info_inv,
      gamma_y = gamma_y,
      gamma_X = gamma_X
    )
  }

  compute_theta_from_operator <- function(local_moments, operator, beta) {
    as.numeric(operator$P %*% (local_moments$eta - local_moments$Xi %*% beta))
  }

  compute_beta_profile <- function(time_objects, Y, WY, gamma_y_list, gamma_X_list, z_index, step) {
    T_periods <- nrow(Y)
    k <- ncol(time_objects$X_aug[[1]])
    G_beta <- matrix(0, k, k)
    rhs <- numeric(k)
    denom <- ncol(Y) * T_periods

    for (t in seq_len(T_periods)) {
      X_t <- time_objects$X_aug[[t]]
      gamma_y_t <- gamma_y_list[[z_index[t]]]
      gamma_X_t <- gamma_X_list[[z_index[t]]]
      G_beta <- G_beta + crossprod(X_t, X_t - WY[t, ] %o% gamma_X_t)
      rhs <- rhs + as.numeric(crossprod(X_t, Y[t, ] - WY[t, ] * gamma_y_t))
    }

    beta <- solve_checked(
      G_beta / denom,
      rhs / denom,
      step = step
    )

    list(
      beta = as.numeric(beta),
      G_beta = G_beta / denom,
      rhs = rhs / denom
    )
  }

  fit_stage <- function(grid, beta_for_scores, time_objects, Y, WY, z, bandwidth, kernel,
                        A_mode, omega_list = NULL, stage_label) {
    local_moments <- vector("list", length(grid))
    operators <- vector("list", length(grid))
    theta_list <- vector("list", length(grid))
    rho <- numeric(length(grid))
    score_list <- vector("list", length(grid))
    omega_out <- vector("list", length(grid))

    for (g in seq_along(grid)) {
      local_moments[[g]] <- compute_local_moments(
        time_objects = time_objects,
        Y = Y,
        WY = WY,
        z = z,
        z0 = grid[g],
        bandwidth = bandwidth,
        kernel = kernel,
        order = 1L
      )

      if (identical(A_mode, "2SLS")) {
        A_g <- solve_checked(
          local_moments[[g]]$QQ,
          step = sprintf("%s 2SLS weighting matrix at z=%s", stage_label, format(grid[g], digits = 6))
        )
      } else {
        A_g <- solve_checked(
          omega_list[[g]],
          step = sprintf("%s GMM weighting matrix at z=%s", stage_label, format(grid[g], digits = 6))
        )
        omega_out[[g]] <- omega_list[[g]]
      }

      operators[[g]] <- build_local_operator(
        local_moments = local_moments[[g]],
        A = A_g,
        step = sprintf("%s coefficient matrix at z=%s", stage_label, format(grid[g], digits = 6))
      )
      theta_list[[g]] <- compute_theta_from_operator(local_moments[[g]], operators[[g]], beta_for_scores)
      rho[g] <- theta_list[[g]][1]

      score_list[[g]] <- compute_local_score(
        local_moments = local_moments[[g]],
        time_objects = time_objects,
        Y = Y,
        z = z,
        beta = beta_for_scores,
        theta = theta_list[[g]],
        bandwidth = bandwidth
      )
    }

    list(
      local_moments = local_moments,
      operators = operators,
      theta = theta_list,
      rho = rho,
      score = score_list,
      omega = omega_out
    )
  }

  compute_local_score <- function(local_moments, time_objects, Y, z, beta, theta, bandwidth) {
    T_periods <- nrow(Y)
    N_units <- ncol(Y)
    q <- nrow(local_moments$Gamma)
    Xi_score <- matrix(NA_real_, T_periods, q)

    for (t in seq_len(T_periods)) {
      H_t <- local_moments$components[[t]]$H
      Q_t <- local_moments$components[[t]]$Q
      k_t <- local_moments$components[[t]]$weight
      e_t <- Y[t, ] - as.numeric(time_objects$X_aug[[t]] %*% beta) - as.numeric(H_t %*% theta)
      Xi_score[t, ] <- as.numeric(k_t * crossprod(Q_t, e_t) / sqrt(N_units * bandwidth))
    }

    Xi_score
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

  resolve_pointwise_hac_bandwidths <- function(score_list, method) {
    bandwidths <- numeric(length(score_list))
    for (g in seq_along(score_list)) {
      bandwidths[g] <- resolve_hac_bandwidth(score_list[[g]], method)
    }
    bandwidths
  }

  build_feasible_omega_list <- function(grid, beta, stage_2sls, time_objects, Y, WY, z,
                                        bandwidth, kernel, hac_method, hac_bandwidths) {
    omega_list <- vector("list", length(grid))

    for (g in seq_along(grid)) {
      score_g <- stage_2sls$score[[g]]
      if (is.null(score_g)) {
        score_g <- compute_local_score(
          local_moments = stage_2sls$local_moments[[g]],
          time_objects = time_objects,
          Y = Y,
          z = z,
          beta = beta,
          theta = stage_2sls$theta[[g]],
          bandwidth = bandwidth
        )
      }
      omega_list[[g]] <- estimate_omega(score_g, hac_method, hac_bandwidths[g])
    }

    omega_list
  }

  evaluate_raw_estimator <- function(grid, beta, time_objects, Y, WY, z, bandwidth, kernel,
                                     estimator, hac_method, hac_bandwidths) {
    stage_2sls <- fit_stage(
      grid = grid,
      beta_for_scores = beta,
      time_objects = time_objects,
      Y = Y,
      WY = WY,
      z = z,
      bandwidth = bandwidth,
      kernel = kernel,
      A_mode = "2SLS",
      stage_label = "profile 2SLS"
    )

    if (identical(estimator, "2SLS")) {
      return(stage_2sls)
    }

    omega_list <- build_feasible_omega_list(
      grid = grid,
      beta = beta,
      stage_2sls = stage_2sls,
      time_objects = time_objects,
      Y = Y,
      WY = WY,
      z = z,
      bandwidth = bandwidth,
      kernel = kernel,
      hac_method = hac_method,
      hac_bandwidths = hac_bandwidths
    )

    fit_stage(
      grid = grid,
      beta_for_scores = beta,
      time_objects = time_objects,
      Y = Y,
      WY = WY,
      z = z,
      bandwidth = bandwidth,
      kernel = kernel,
      A_mode = "GMM",
      omega_list = omega_list,
      stage_label = "profile GMM"
    )
  }

  estimate_beta_score_matrix <- function(time_objects, Y, WY, rho_at_t, beta_raw) {
    T_periods <- nrow(Y)
    k <- length(beta_raw)
    phi <- matrix(NA_real_, T_periods, k)

    for (t in seq_len(T_periods)) {
      residual_t <- Y[t, ] - as.numeric(time_objects$X_aug[[t]] %*% beta_raw) - rho_at_t[t] * WY[t, ]
      phi[t, ] <- as.numeric(crossprod(time_objects$X_aug[[t]], residual_t) / sqrt(ncol(Y)))
    }

    phi
  }

  compute_raw_residuals <- function(time_objects, Y, WY, rho_at_t, beta_raw) {
    T_periods <- nrow(Y)
    residuals <- matrix(NA_real_, nrow(Y), ncol(Y))

    for (t in seq_len(T_periods)) {
      fitted_t <- as.numeric(time_objects$X_aug[[t]] %*% beta_raw) + rho_at_t[t] * WY[t, ]
      residuals[t, ] <- Y[t, ] - fitted_t
    }

    residuals
  }

  compute_beta_bias <- function(time_objects, Y, WY, z_index, gamma_upsilon_sample, G_beta) {
    T_periods <- nrow(Y)
    k <- ncol(time_objects$X_aug[[1]])
    rhs <- numeric(k)
    denom <- ncol(Y) * T_periods

    for (t in seq_len(T_periods)) {
      rhs <- rhs + as.numeric(crossprod(
        time_objects$X_aug[[t]],
        WY[t, ] * gamma_upsilon_sample[[z_index[t]]]
      ))
    }

    as.numeric(-solve_checked(
      G_beta,
      rhs / denom,
      step = "beta bias correction"
    ))
  }

  compute_local_quadratic_fit <- function(grid, beta, time_objects, Y, WY, z, bandwidth, kernel) {
    fits <- vector("list", length(grid))

    for (g in seq_along(grid)) {
      local_moments <- compute_local_moments(
        time_objects = time_objects,
        Y = Y,
        WY = WY,
        z = z,
        z0 = grid[g],
        bandwidth = bandwidth,
        kernel = kernel,
        order = 2L
      )
      A_g <- solve_checked(
        local_moments$QQ,
        step = sprintf("local quadratic 2SLS weighting matrix at z=%s", format(grid[g], digits = 6))
      )
      operator <- build_local_operator(
        local_moments = local_moments,
        A = A_g,
        step = sprintf("local quadratic coefficient matrix at z=%s", format(grid[g], digits = 6))
      )
      theta <- compute_theta_from_operator(local_moments, operator, beta)
      fits[[g]] <- list(local_moments = local_moments, operator = operator, theta = theta)
    }

    fits
  }

  compute_gamma_upsilon <- function(grid, raw_stage, quad_stage, WY, z, bandwidth, pilot_bandwidth, kernel) {
    gamma_upsilon <- vector("list", length(grid))

    for (g in seq_along(grid)) {
      z0 <- grid[g]
      local_moments <- raw_stage$local_moments[[g]]
      theta_quad <- quad_stage[[g]]$theta
      theta3_hat <- theta_quad[3]
      upsilon <- numeric(length(local_moments$eta))
      denom <- ncol(WY) * length(z) * bandwidth

      for (t in seq_along(z)) {
        u_main <- (z[t] - z0) / bandwidth
        k_t <- kernel_weight(u_main, kernel)
        u_pilot <- (z[t] - z0) / pilot_bandwidth
        r_t <- theta3_hat * (u_pilot^2) * WY[t, ]
        upsilon <- upsilon + k_t * as.numeric(crossprod(local_moments$components[[t]]$Q, r_t))
      }

      gamma_upsilon[[g]] <- as.numeric(raw_stage$operators[[g]]$P[1, , drop = FALSE] %*% (upsilon / denom))
    }

    gamma_upsilon
  }

  finalize_beta_output <- function(beta_raw, beta_bias, beta_se, ci_level, beta_names, inference) {
    z_alpha <- stats::qnorm((1 + ci_level) / 2)
    beta_center <- if (identical(inference, "debiased")) {
      beta_raw - beta_bias
    } else {
      beta_raw
    }
    # Use the inference center so debiased p-values match the reported interval.
    beta_statistic <- beta_center / beta_se
    beta_p <- 2 * stats::pnorm(abs(beta_statistic), lower.tail = FALSE)

    # The public point estimate remains raw; debiased inference only shifts the CI center.
    beta_output <- data.frame(
      term = beta_names,
      estimate = beta_raw,
      debiased = beta_center,
      bias = beta_bias,
      se = beta_se,
      lower = beta_center - z_alpha * beta_se,
      upper = beta_center + z_alpha * beta_se,
      p = beta_p,
      stringsAsFactors = FALSE
    )

    # Keep covariates first and place the intercept last in the public table.
    beta_output <- beta_output[c(
      which(beta_output$term != "(Intercept)"),
      which(beta_output$term == "(Intercept)")
    ), , drop = FALSE]
    rownames(beta_output) <- NULL
    beta_output
  }

  finalize_rho_output <- function(grid, rho_raw, rho_bias, rho_se, ci_level, inference) {
    z_alpha <- stats::qnorm((1 + ci_level) / 2)
    rho_center <- if (identical(inference, "debiased")) {
      rho_raw - rho_bias
    } else {
      rho_raw
    }
    rho_statistic <- rho_center / rho_se
    rho_p <- 2 * stats::pnorm(abs(rho_statistic), lower.tail = FALSE)

    data.frame(
      z = grid,
      estimate = rho_raw,
      debiased = rho_center,
      bias = rho_bias,
      se = rho_se,
      lower = rho_center - z_alpha * rho_se,
      upper = rho_center + z_alpha * rho_se,
      p = rho_p,
      stringsAsFactors = FALSE
    )
  }

  function(
    Y, X, W, z,
    grid = sort(unique(z)),
    iv_lag = 1,
    estimator = c("2SLS", "GMM"),
    inference = c("debiased", "undersmooth"),
    bandwidth = NULL,
    pilot_bandwidth = NULL,
    kernel = c("gaussian", "epanechnikov", "triangular", "uniform", "quartic"),
    ci_level = 0.95,
    hac_method = c("Newey-West", "Andrews"),
    preprocess_W = FALSE
  ) {
    estimator <- match.arg(estimator)
    inference <- match.arg(inference)
    kernel <- match.arg(kernel)
    hac_method <- match.arg(hac_method)

    Y <- as.matrix(Y)
    z <- as.numeric(z)
    grid <- as.numeric(grid)

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

    z_scale <- robust_scale_z(z)
    default_bandwidth <- 1.06 * z_scale * (N_units * T_periods)^(-1 / 5)
    if (!is.finite(default_bandwidth) || default_bandwidth <= 0) {
      default_bandwidth <- 1
    }
    bandwidth_source <- if (is.null(bandwidth)) "default_mse" else "user_supplied"
    if (identical(inference, "undersmooth") && is.null(bandwidth)) {
      bandwidth_source <- "default_mse_for_undersmooth"
    }
    if (is.null(bandwidth)) {
      bandwidth <- default_bandwidth
    } else {
      bandwidth <- as.numeric(bandwidth[[1]])
    }

    if (identical(inference, "debiased")) {
      if (is.null(pilot_bandwidth)) {
        pilot_bandwidth <- 1.06 * z_scale * (N_units * T_periods)^(-1 / 7)
      } else {
        pilot_bandwidth <- as.numeric(pilot_bandwidth[[1]])
      }
    } else {
      pilot_bandwidth <- NA_real_
    }

    covariate_names <- get_covariate_names(X)
    beta_names <- c("(Intercept)", covariate_names)
    time_objects <- build_time_objects(X, W, T_periods, iv_lag, covariate_names)
    WY <- Y %*% t(W)
    sample_grid <- sort(unique(z))
    z_index <- match(z, sample_grid)

    preliminary_operators <- fit_stage(
      grid = sample_grid,
      beta_for_scores = rep(0, length(beta_names)),
      time_objects = time_objects,
      Y = Y,
      WY = WY,
      z = z,
      bandwidth = bandwidth,
      kernel = kernel,
      A_mode = "2SLS",
      stage_label = "preliminary profile 2SLS"
    )
    gamma_y_init <- lapply(preliminary_operators$operators, `[[`, "gamma_y")
    gamma_X_init <- lapply(preliminary_operators$operators, `[[`, "gamma_X")
    beta_init_fit <- compute_beta_profile(
      time_objects = time_objects,
      Y = Y,
      WY = WY,
      gamma_y_list = gamma_y_init,
      gamma_X_list = gamma_X_init,
      z_index = z_index,
      step = "preliminary profiled beta"
    )
    beta_init <- beta_init_fit$beta

    preliminary_2sls_sample <- evaluate_raw_estimator(
      grid = sample_grid,
      beta = beta_init,
      time_objects = time_objects,
      Y = Y,
      WY = WY,
      z = z,
      bandwidth = bandwidth,
      kernel = kernel,
      estimator = "2SLS",
      hac_method = hac_method,
      hac_bandwidths = numeric(0)
    )
    rho_sample_hac_bandwidths <- if (identical(estimator, "GMM")) {
      resolve_pointwise_hac_bandwidths(preliminary_2sls_sample$score, hac_method)
    } else {
      numeric(0)
    }

    preliminary_stage_sample <- if (identical(estimator, "GMM")) {
      evaluate_raw_estimator(
        grid = sample_grid,
        beta = beta_init,
        time_objects = time_objects,
        Y = Y,
        WY = WY,
        z = z,
        bandwidth = bandwidth,
        kernel = kernel,
        estimator = "GMM",
        hac_method = hac_method,
        hac_bandwidths = rho_sample_hac_bandwidths
      )
    } else {
      preliminary_2sls_sample
    }
    gamma_y_sample <- lapply(preliminary_stage_sample$operators, `[[`, "gamma_y")
    gamma_X_sample <- lapply(preliminary_stage_sample$operators, `[[`, "gamma_X")
    beta_raw_fit <- compute_beta_profile(
      time_objects = time_objects,
      Y = Y,
      WY = WY,
      gamma_y_list = gamma_y_sample,
      gamma_X_list = gamma_X_sample,
      z_index = z_index,
      step = sprintf("profile %s beta", estimator)
    )
    beta_raw <- beta_raw_fit$beta

    raw_stage_sample <- evaluate_raw_estimator(
      grid = sample_grid,
      beta = beta_raw,
      time_objects = time_objects,
      Y = Y,
      WY = WY,
      z = z,
      bandwidth = bandwidth,
      kernel = kernel,
      estimator = estimator,
      hac_method = hac_method,
      hac_bandwidths = rho_sample_hac_bandwidths
    )
    raw_stage_grid_2sls <- evaluate_raw_estimator(
      grid = grid,
      beta = beta_raw,
      time_objects = time_objects,
      Y = Y,
      WY = WY,
      z = z,
      bandwidth = bandwidth,
      kernel = kernel,
      estimator = "2SLS",
      hac_method = hac_method,
      hac_bandwidths = numeric(0)
    )
    rho_grid_hac_bandwidths <- resolve_pointwise_hac_bandwidths(raw_stage_grid_2sls$score, hac_method)
    raw_stage_output <- if (identical(estimator, "GMM")) {
      evaluate_raw_estimator(
        grid = grid,
        beta = beta_raw,
        time_objects = time_objects,
        Y = Y,
        WY = WY,
        z = z,
        bandwidth = bandwidth,
        kernel = kernel,
        estimator = "GMM",
        hac_method = hac_method,
        hac_bandwidths = rho_grid_hac_bandwidths
      )
    } else {
      raw_stage_grid_2sls
    }

    rho_raw_output <- raw_stage_output$rho
    rho_raw_sample <- raw_stage_sample$rho
    rho_at_t <- rho_raw_sample[z_index]

    phi_matrix <- estimate_beta_score_matrix(time_objects, Y, WY, rho_at_t, beta_raw)
    beta_hac_bandwidth <- resolve_hac_bandwidth(phi_matrix, hac_method)
    omega_phi <- estimate_omega(phi_matrix, hac_method, beta_hac_bandwidth)
    G_beta <- beta_raw_fit$G_beta
    V_beta_const <- solve_checked(
      G_beta,
      step = "beta sandwich left factor"
    ) %*% omega_phi %*% t(solve_checked(
      G_beta,
      step = "beta sandwich right factor"
    ))
    beta_se <- sqrt(pmax(diag(V_beta_const), 0) / (N_units * T_periods))

    rho_se <- numeric(length(grid))
    for (g in seq_along(grid)) {
      local_fit <- raw_stage_output
      omega_rho <- estimate_omega(local_fit$score[[g]], hac_method, rho_grid_hac_bandwidths[g])
      operator_g <- local_fit$operators[[g]]
      local_moments_g <- local_fit$local_moments[[g]]
      meat_g <- crossprod(
        local_moments_g$Gamma,
        operator_g$A %*% omega_rho %*% operator_g$A %*% local_moments_g$Gamma
      )
      V_theta_const <- operator_g$info_inv %*% meat_g %*% operator_g$info_inv
      rho_se[g] <- sqrt(max(V_theta_const[1, 1], 0) / (N_units * T_periods * bandwidth))
    }

    if (identical(inference, "debiased")) {
      quadratic_sample <- compute_local_quadratic_fit(
        grid = sample_grid,
        beta = beta_raw,
        time_objects = time_objects,
        Y = Y,
        WY = WY,
        z = z,
        bandwidth = pilot_bandwidth,
        kernel = kernel
      )
      quadratic_output <- compute_local_quadratic_fit(
        grid = grid,
        beta = beta_raw,
        time_objects = time_objects,
        Y = Y,
        WY = WY,
        z = z,
        bandwidth = pilot_bandwidth,
        kernel = kernel
      )
      gamma_upsilon_sample <- compute_gamma_upsilon(
        grid = sample_grid,
        raw_stage = raw_stage_sample,
        quad_stage = quadratic_sample,
        WY = WY,
        z = z,
        bandwidth = bandwidth,
        pilot_bandwidth = pilot_bandwidth,
        kernel = kernel
      )
      gamma_upsilon_output <- compute_gamma_upsilon(
        grid = grid,
        raw_stage = raw_stage_output,
        quad_stage = quadratic_output,
        WY = WY,
        z = z,
        bandwidth = bandwidth,
        pilot_bandwidth = pilot_bandwidth,
        kernel = kernel
      )
      beta_bias <- compute_beta_bias(
        time_objects = time_objects,
        Y = Y,
        WY = WY,
        z_index = z_index,
        gamma_upsilon_sample = gamma_upsilon_sample,
        G_beta = G_beta
      )
      rho_bias <- numeric(length(grid))
      for (g in seq_along(grid)) {
        rho_bias[g] <- gamma_upsilon_output[[g]] - sum(raw_stage_output$operators[[g]]$gamma_X * beta_bias)
      }
    } else {
      beta_bias <- rep(0, length(beta_raw))
      rho_bias <- rep(0, length(grid))
    }

    residuals_raw <- compute_raw_residuals(
      time_objects = time_objects,
      Y = Y,
      WY = WY,
      rho_at_t = rho_raw_sample[z_index],
      beta_raw = beta_raw
    )

    rho_sample_hac_df <- if (identical(estimator, "GMM")) {
      data.frame(
        z = sample_grid,
        bandwidth = rho_sample_hac_bandwidths,
        stringsAsFactors = FALSE
      )
    } else {
      # Keep the public shape stable even though sample-z GMM weighting is unused under 2SLS.
      data.frame(
        z = numeric(0),
        bandwidth = numeric(0),
        stringsAsFactors = FALSE
      )
    }
    rho_grid_hac_df <- data.frame(
      z = grid,
      bandwidth = rho_grid_hac_bandwidths,
      stringsAsFactors = FALSE
    )

    result <- list(
      beta = finalize_beta_output(
        beta_raw = beta_raw,
        beta_bias = beta_bias,
        beta_se = beta_se,
        ci_level = ci_level,
        beta_names = beta_names,
        inference = inference
      ),
      rho = finalize_rho_output(
        grid = grid,
        rho_raw = rho_raw_output,
        rho_bias = rho_bias,
        rho_se = rho_se,
        ci_level = ci_level,
        inference = inference
      ),
      residuals = residuals_raw,
      hac_bandwidths = list(
        beta = data.frame(
          bandwidth = beta_hac_bandwidth,
          stringsAsFactors = FALSE
        ),
        rho_sample = rho_sample_hac_df,
        rho_grid = rho_grid_hac_df
      )
    )
    if (isTRUE(W_processed)) {
      result$W <- W
    }
    result$diagnostics <- list(
      estimator = estimator,
      inference = inference,
      kernel = kernel,
      main_bandwidth = bandwidth,
      pilot_bandwidth = pilot_bandwidth,
      bandwidth_source = bandwidth_source,
      hac_method = hac_method,
      iv_lag = as.integer(iv_lag),
      W_processed = W_processed,
      dropped_units = dropped_units,
      N_used = N_units
    )
    result
  }
})
