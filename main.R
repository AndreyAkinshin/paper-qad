source("utils.R")
source("mad-factors.R")

knitr::opts_chunk$set(echo = FALSE, fig.align = "center", fig.pos = "ht!", fig.height = 3.5)
options(knitr.kable.NA = "-")
theme_set(theme_bw())
options(scipen = 999)

# Settings ---------------------------------------------------------------------

efficiency_ns <- c(
  2:100,
  c(109, 110) + rep(seq(0, 90, by = 10), each = 2),
  c(249, 250) + rep(seq(0, 250, by = 50), each = 2),
  600, 700, 800, 900, 1000, 1500, seq(2000, 10000, by = 1000),
  25000, 50000, 100000
)
settings <- list(
  constants = list(
    rebuild = FALSE,
    filename = "data-consistency-constants.csv",
    repetitions = 25 * 1000 * 1000,
    ns = c(
      2:100,
      c(109, 110) + rep(seq(0, 90, by = 10), each = 2),
      c(249, 250) + rep(seq(0, 250, by = 50), each = 2),
      seq(600, 1000, by = 100),
      seq(1500, 5000, by = 500),
      seq(6000, 10000, by = 1000)
    )
  ),
  scale_efficiency = list(
    rebuild = FALSE,
    filename = "data-scale-efficiency.csv",
    repetitions = 10 * 1000 * 1000,
    ns = efficiency_ns
  ),
  location_efficiency = list(
    rebuild = FALSE,
    filename = "data-location-efficiency.csv",
    repetitions = 10 * 1000 * 1000,
    ns = efficiency_ns
  )
)

# Functions --------------------------------------------------------------------

format_int_latex <- function(n) {
  s <- format(n, scientific = FALSE)
  res <- ""
  while (nchar(s) > 0) {
    if (nchar(res) > 0) {
      res <- paste0("\\,", res)
    }
    res <- paste0(substr(s, nchar(s) - 2, nchar(s)), res)
    s <- substr(s, 1, nchar(s) - 3)
  }
  res
}

erfinv <- function(x) sqrt(qchisq(abs(x), 1) / 2) * sign(x)

mad <- function(x) median(abs(x - median(x)))
qad <- function(x, p) as.numeric(quantile(abs(x - median(x)), p))

d_trimodal <- list(
  title = "A trimodal distribution",
  d = function(x) dbeta(x, 1.5, 1.5) / 4 + dbeta(x - 4, 1.5, 1.5) / 2 + dbeta(x - 8, 1.5, 1.5) / 4,
  q = function(p) {
    ifelse(p <= 0.25,
      qbeta(p * 4, 1.5, 1.5),
      ifelse(p <= 0.75,
        4 + qbeta((p - 0.25) * 2, 1.5, 1.5),
        8 + qbeta((p - 0.75) * 4, 1.5, 1.5)
      )
    )
  },
  p = function(q) {
    ifelse(q < 4,
      pbeta(q, 1.5, 1.5) / 4,
      ifelse(q < 8,
        0.25 + pbeta(q - 4, 1.5, 1.5) / 2,
        0.75 + pbeta(q - 8, 1.5, 1.5) / 4
      )
    )
  },
  r = function(n) {
    r1 <- rbeta(n, 1.5, 1.5)
    r2 <- rbeta(n, 1.5, 1.5) + 4
    r3 <- rbeta(n, 1.5, 1.5) + 8
    m <- rbind(r1, r2, r2, r3)
    ind <- sample(1:nrow(m), n, T)
    sapply(1:n, function(i) m[ind[i], i])
  }
)

reload_constants <- function() {
  constants <<- if (file.exists(settings$constants$filename)) read.csv(settings$constants$filename) else data.frame()
  constants
}
invisible(reload_constants())

Pm <- 0.5
Ps <- pnorm(1) - pnorm(-1)
Po <- 0.861678977787423

get_asymptotic_constant <- function(p) {
  if (p == "sqad") {
    return(1)
  }
  if (p == "oqad") {
    return(1 / (sqrt(2) * erfinv(Po)))
  }
  1 / (sqrt(2) * erfinv(p))
}

get_bias_coefficients <- function(estimator) { # estimator is "sqad" or "oqad"
  df <- constants
  df <- df[df$n >= 100 & df$n <= 1000, ]
  df$bias <- df[, estimator] / get_asymptotic_constant(estimator) - 1
  fit <- lm(bias ~ 0 + I(n^(-1)) + I(n^(-2)), data = df)
  fit$coefficients
}

get_constant <- function(estimator, n, always_predict = FALSE) {
  if (n %in% constants$n && !always_predict) {
    return(constants[constants$n == n, estimator])
  }
  coef <- get_bias_coefficients(estimator)
  as.numeric(1 + coef[1] / n + coef[2] / n^2) * get_asymptotic_constant(estimator)
}

sqad_constant <- function(n, always_predict = FALSE) get_constant("sqad", n, always_predict)
oqad_constant <- function(n, always_predict = FALSE) get_constant("oqad", n, always_predict)

sqad <- function(x, constant = NULL) {
  n <- length(x)
  if (is.null(constant)) {
    constant <- sqad_constant(n)
  }
  constant * qad(x, Ps)
}

oqad <- function(x, constant = NULL) {
  n <- length(x)
  if (is.null(constant)) {
    constant <- oqad_constant(n)
  }
  constant * qad(x, Po)
}

c4 <- function(n) ifelse(n < 300, sqrt(2 / (n - 1)) * gamma(n / 2) / gamma((n - 1) / 2), 1)
sd_unbiased <- function(x) sd(x) / c4(length(x))

sthdme <- function(x) quantile.thd(x, 0.5, Ps)
othdme <- function(x) quantile.thd(x, 0.5, Po)

qad_asymptotic_efficiency <- function(ps) {
  unname(sapply(ps, function(p) {
    if (p == "mad") {
      return(qad_asymptotic_efficiency(Pm))
    }
    if (p == "sqad") {
      return(qad_asymptotic_efficiency(Ps))
    }
    if (p == "oqad") {
      return(qad_asymptotic_efficiency(Po))
    }
    2 * erfinv(p)^2 / (pi * p * (1 - p) * exp(2 * erfinv(p)^2))
  }))
}

get_po <- function(eps = 1e-8) {
  l <- 0.1
  r <- 0.9
  while (abs(r - l) > eps) {
    m <- (l + r) / 2
    if (qad_asympt_eff(m) <= qad_asympt_eff(m + eps)) {
      l <- m
    } else {
      r <- m
    }
  }
  (l + r) / 2
}

# Simulations ------------------------------------------------------------------

simulation_constants <- function(rebuild = NULL, filename = NULL, repetitions = NULL, ns = NULL) {
  apply_settings(settings$constants)

  process <- function(n) {
    sqad_constant <- 1 / mean(future_replicate(repetitions, sqad(rnorm(n), 1)))
    oqad_constant <- 1 / mean(future_replicate(repetitions, oqad(rnorm(n), 1)))
    data.frame(n = n, sqad = round(sqad_constant, 6), oqad = round(oqad_constant, 6))
  }

  multi_estimate(rebuild, filename, ns, process)
  reload_constants()
}

simulation_scale_efficiency <- function(rebuild = NULL, filename = NULL, repetitions = NULL, ns = NULL) {
  apply_settings(settings$scale_efficiency)

  estimate <- function(x) {
    c(
      n = length(x),
      sd = sd_unbiased(x),
      mad = mad.sm(x),
      sqad = sqad(x),
      oqad = oqad(x)
    )
  }

  process <- function(n) {
    df_n <- data.frame(t(future_replicate(repetitions, estimate(rnorm(n)))))

    df <- data.frame(
      n = n,
      bias.sd = mean(df_n$sd),
      bias.mad = mean(df_n$mad),
      bias.sqad = mean(df_n$sqad),
      bias.oqad = mean(df_n$oqad),
      svar.sd = n * var(df_n$sd) / mean(df_n$sd)^2,
      svar.mad = n * var(df_n$mad) / mean(df_n$mad)^2,
      svar.sqad = n * var(df_n$sqad) / mean(df_n$sqad)^2,
      svar.oqad = n * var(df_n$oqad) / mean(df_n$oqad)^2
    )
    df$eff.mad <- df$svar.sd / df$svar.mad
    df$eff.sqad <- df$svar.sd / df$svar.sqad
    df$eff.oqad <- df$svar.sd / df$svar.oqad

    round(df, 6)
  }

  multi_estimate(rebuild, filename, ns, process)
}

simulation_location_efficiency <- function(rebuild = NULL, filename = NULL, repetitions = NULL, ns = NULL) {
  apply_settings(settings$location_efficiency)

  estimate <- function(x) {
    c(
      mean = mean(x),
      median = median(x),
      sthdme = sthdme(x),
      othdme = othdme(x)
    )
  }

  process <- function(n) {
    df0 <- data.frame(t(future_replicate(repetitions, estimate(rnorm(n)))))

    df <- data.frame(
      n = n,
      median = (var(df0$mean) / mean(df0$mean)^2) / (var(df0$median) / mean(df0$median)^2),
      sthdme = (var(df0$mean) / mean(df0$mean)^2) / (var(df0$sthdme) / mean(df0$sthdme)^2),
      othdme = (var(df0$mean) / mean(df0$mean)^2) / (var(df0$othdme) / mean(df0$othdme)^2)
    )

    round(df, 6)
  }

  multi_estimate(rebuild, filename, ns, process)
}

data_constant_compare <- function(estimator) {
  df <- simulation_constants()
  df$predicted <- sapply(df$n, function(n) get_constant(estimator, n, n >= 100))
  df$diff <- abs(df[, estimator] - df$predicted)
  df
}

# Inlines ----------------------------------------------------------------------

inline_constant_diff <- function(estimator) {
  df <- data_constant_compare(estimator)
  round(max(df$diff), 6)
}

inline_pareto_3mad_coverage <- function() {
  M <- qpareto(0.5, 1)
  MAD <- uniroot.all(function(mad) ppareto(M + mad, 1) - ppareto(M - mad, 1) - 0.5,
    c(0, 10),
    tol = 1e-15
  )
  c <- 1 / qnorm(0.75)
  k <- 3 * c
  round(ppareto(M + k * MAD, 1) - ppareto(M - k * MAD, 1), 3) * 100
}

# Tables -----------------------------------------------------------------------

table_constants <- function() {
  df <- reload_constants()

  size <- 50
  caption <- "Finite-sample consistency constants for the $\\SQAD$ and the $\\OQAD$."

  if (nrow(df[df$n == 1, ]) == 0) {
    df <- rbind(data.frame(n = 1, sqad = NA, oqad = NA), df)
  }

  slice <- function(index) df[((index - 1) * size + 1):(index * size), ]

  slice_count <- ceiling(nrow(df) / size)
  df2 <- slice(1)
  for (i in 2:slice_count) {
    df2 <- cbind(df2, slice(i))
  }

  header <- rep(c("n", "$\\SQAD$", "$\\OQAD$"), slice_count)
  knitr::kable(df2, caption = caption, col.names = header, escape = FALSE, digits = 4) %>%
    kable_styling(latex_options = "hold_position")
}

sliced_table_efficiency <- function(df, size, caption, header) {
  fake_row <- function(n) {
    df_fake <- data.frame(matrix(ncol = ncol(df), nrow = 0))
    names(df_fake) <- names(df)
    df_fake[1, 1] <- n
    df_fake
  }
  if (nrow(df[df$n == 1, ]) == 0) {
    df <- rbind(fake_row(1), df)
  }

  slice <- function(index) df[((index - 1) * size + 1):(index * size), ]

  slice_count <- ceiling(nrow(df) / size)
  df2 <- slice(1)
  for (i in 2:slice_count) {
    df2 <- cbind(df2, slice(i))
  }

  header <- rep(header, slice_count)
  knitr::kable(df2, caption = caption, col.names = header, escape = FALSE, digits = 4) %>%
    kable_styling(latex_options = "hold_position") %>%
    row_spec(0, font_size = 6)
}

table_scale_efficiency <- function() {
  df <- simulation_scale_efficiency() %>% extract_coumns("eff")
  sliced_table_efficiency(
    df, 50,
    "Finite-sample Gaussian efficiency of MAD, SQAD, OQAD.",
    c("n", "MAD", "SQAD", "OQAD")
  )
}

table_location_efficiency <- function() {
  df <- simulation_location_efficiency()
  sliced_table_efficiency(
    df, 50,
    "Finite-sample Gaussian efficiency of SM, STHDME, OTHDME.",
    TeX(c("n", "SM", "STHDME", "OTHDME"))
  )
}

table_mad_coverage <- function() {
  c <- 1 / qnorm(0.75)
  k <- round(c(1, c, 2 * c, 3 * c), 2)
  cover.normal <- function(k) {
    M <- qnorm(0.5)
    MAD <- uniroot.all(function(mad) pnorm(M + mad) - pnorm(M - mad) - 0.5,
      c(0, 10),
      tol = 1e-15
    )
    round(pnorm(M + k * MAD) - pnorm(M - k * MAD), 3)
  }
  cover.pareto <- function(k) {
    M <- qpareto(0.5, 1)
    MAD <- uniroot.all(function(mad) ppareto(M + mad, 1) - ppareto(M - mad, 1) - 0.5,
      c(0, 10),
      tol = 1e-15
    )
    round(ppareto(M + k * MAD, 1) - ppareto(M - k * MAD, 1), 3)
  }
  df <- data.frame(
    k = round(c(1, c, 2 * c, 3 * c), 2),
    normal = cover.normal(k),
    pareto = cover.pareto(k)
  )
  caption <- "Coverage of the Normal(0, 1) and Pareto(1, 1) distributions by various intervals."
  header <- c("", "Normal", "Pareto(1, 1)")
  kbl(df, caption = caption, col.names = header, escape = FALSE, align = "c") %>%
    add_header_above(c(
      "k",
      "$\\\\mathbb{P}(M - k\\\\cdot \\\\mathrm{MAD}) \\\\leq X \\\\leq \\\\mathbb{P}(M + k\\\\cdot \\\\mathrm{MAD})$" = 2
    ), escape = FALSE) %>%
    column_spec(2, width = "4cm") %>%
    column_spec(3, width = "4cm") %>%
    kable_styling(latex_options = "hold_position")
}

# Figures (MAD) ----------------------------------------------------------------

figure_mad_efficiency <- function(maxn = 100) {
  df <- simulation_scale_efficiency() %>% extract_coumns("eff")
  df <- df[df$n > 2 & df$n <= maxn, ]
  df$parity <- factor(ifelse(df$n %% 2 == 0, "Even", "Odd"), levels = c("Even", "Odd"))

  ggplot(df, aes(n, mad, shape = parity)) +
    geom_hline(yintercept = 0.3675) +
    geom_line(alpha = 0.2, col = cbp$red) +
    geom_point(size = 1, col = cbp$red) +
    labs(
      title = "Finite-sample Gaussian efficiency of MAD",
      x = "Sample size (n)",
      y = "Gaussian efficiency",
      col = "Parity of n",
      shape = "Parity of n"
    ) +
    scale_color_manual(values = cbp$values) +
    scale_shape_manual(values = c(3, 4)) +
    scale_y_continuous(breaks = c(0.3675, seq(0.4, 0.55, by = 0.05)))
}

figure_mad_multimodal <- function() {
  x <- seq(0, 9, by = 0.01)
  y <- d_trimodal$d(x)
  df.density <- data.frame(x, y)
  df.labels <- data.frame(
    x = c(0.5, 4.5, 8.5),
    y = c(0.1, 0.1, 0.1),
    label = c("25%", "50%", "25%"),
    size = c(5, 7, 5)
  )
  p1 <- ggplot(df.density, aes(x, y)) +
    geom_line() +
    geom_text(data = df.labels, size = df.labels$size, aes(x, y, label = label)) +
    scale_x_continuous(breaks = seq(0, 9, by = 1)) +
    ggtitle("(a) A trimodal distribution") +
    labs(y = "Density") +
    theme(legend.position = "none")

  set.seed(1729)
  mads <- replicate(1000, mad(d_trimodal$r(100)))
  p2 <- ggplot(data.frame(x = mads), aes(x)) +
    geom_density(bw = "SJ") +
    scale_x_continuous(breaks = seq(0, 6, by = 1), limits = c(0, 5)) +
    ggtitle("(b) Distribution of MAD estimations from a trimodal distribution") +
    labs(y = "Density")

  ggarrange(p1, p2, ncol = 1)
}

figure_mad_discrete <- function() {
  lambda <- 0.6
  x <- 0:6
  y <- dpois(0:6, lambda)
  df <- data.frame(
    x = x,
    y = y,
    yend = rep(0, length(x)),
    ytext = y + 0.05,
    label = signif(y, 2),
    lambda = lambda
  )
  ggplot(df, aes(x = x, y = y, xend = x, yend = yend)) +
    geom_point() +
    geom_segment() +
    scale_x_continuous(breaks = seq(min(df$x), max(df$x), by = 1), expand = c(0, 0.4)) +
    scale_y_continuous(limits = c(0, max(df$y) + 0.2)) +
    ggtitle(TeX(paste0("Poisson distribution for $\\lambda$=", lambda))) +
    labs(y = "Probability") +
    geom_text(aes(x, ytext, label = label))
}

figure_mad_heavy <- function() {
  M <- qpareto(0.5, 1)
  MAD <- uniroot.all(function(mad) ppareto(M + mad, 1) - ppareto(M - mad, 1) - 0.5,
    c(0, 10),
    tol = 1e-15
  )

  build.df <- function() {
    x <- seq(1, 10, by = 0.01)
    data.frame(
      x = x,
      y = dpareto(x, 1)
    )
  }
  build.df.seg <- function() {
    x <- c(
      M - 4.45 * MAD, M - 2.97 * MAD, M - 1.48 * MAD, M - MAD, M,
      M + MAD, M + 1.48 * MAD, M + 2.97 * MAD, M + 4.45 * MAD
    )
    ytext <- dpareto(x, 1) + 0.01
    df.seg <- data.frame(
      x = x,
      y = dpareto(x, 1),
      ytext = ytext,
      yend = rep(0, length(x)),
      angle = c(rep(-45, 3), rep(0, 2), rep(45, 4)),
      hjust = c(rep(1.1, 3), rep(0, 2), rep(-0.1, 4)),
      vjust = c(rep(0.3, 3), rep(-0.5, 6)),
      linetype = c(rep("dashed", 4), "solid", rep("dashed", 4)),
      label = c(
        "M - 4.45 MAD", "M - 2.97 MAD", "M - 1.48 MAD", "M - MAD", "M",
        "M + MAD", "M + 1.48 MAD", "M + 2.97 MAD", "M + 4.45 MAD"
      )
    )
  }
  df <- build.df()
  df.seg <- build.df.seg()
  color <- cbp$navy
  braket <- function(x1, x2, depth) {
    df <- data.frame(
      x = c(x1, x1, x2, x2),
      y = c(0, depth, depth, 0) * -0.03
    )
    geom_path(data = df, aes(x, y), col = color, linetype = "dashed")
  }
  ggplot(df, aes(x, y)) +
    geom_hline(yintercept = 0, col = cbp$grey) +
    geom_area(fill = cbp$grey, alpha = 0.3) +
    geom_line() +
    geom_segment(aes(x, y, xend = x, yend = yend),
      df.seg,
      col = color, linetype = df.seg$linetype
    ) +
    geom_point(data = df.seg, col = color) +
    geom_text(aes(x, ytext, label = label, hjust = hjust, vjust = vjust, angle = angle),
      df.seg,
      size = 2.5, col = color
    ) +
    scale_x_continuous(limits = c(min(df.seg$x) - 0.7, max(df$x)), breaks = seq(-3, 10, by = 1)) +
    scale_y_continuous(limits = c(NA, 1), breaks = seq(0, 1, by = 0.1)) +
    ggtitle("Pareto(1, 1) distribution") +
    labs(y = "Density") +
    braket(M - MAD, M + MAD, 1) +
    braket(M - 1.48 * MAD, M + 1.48 * MAD, 2) +
    braket(M - 2.97 * MAD, M + 2.97 * MAD, 3) +
    braket(M - 4.45 * MAD, M + 4.45 * MAD, 4)
}

# Figures (QAD) ----------------------------------------------------------------

figure_qad_functions <- function() {
  build_df <- function(pD, qD, caption) {
    ps <- seq(0, 0.99, by = 0.001)
    qad <- function(p) {
      tail(
        uniroot.all(function(qad) pD(qD(0.5) + qad) - pD(qD(0.5) - qad) - p,
          c(0, qD(0.995) - qD(0.005)),
          tol = 1e-7
        ), 1
      )
    }
    qads <- sapply(ps, function(p) qad(p))
    df <- data.frame(ps, qads, caption)
  }
  df <- rbind(
    build_df(pnorm, qnorm, "(a) Normal(0, 1)"),
    build_df(punif, qunif, "(b) Uniform(0, 1)"),
    build_df(pexp, qexp, "(c) Exponential(1)"),
    build_df(d_trimodal$p, d_trimodal$q, "(d) Trimodal"),
    build_df(function(q) ppois(q, 0.6), function(p) qpois(p, 0.6), "(e) Pois(0.6)"),
    build_df(function(q) ppareto(q, 1), function(p) qpareto(p, 1), "(f) Pareto(1, 1)")
  )
  df$caption <- factor(df$caption)
  df_mad <- df[df$ps == 0.5, ]
  df_trimodal <- data.frame(
    x = 0.5,
    y = df[df$caption == "(d) Trimodal" & df$ps == 0.49, ]$qads,
    xend = 0.5,
    yend = df[df$caption == "(d) Trimodal" & df$ps == 0.50, ]$qads,
    caption = "(d) Trimodal"
  )
  ggplot(df, aes(ps, qads)) +
    geom_line() +
    geom_point(data = df_mad, col = cbp$navy) +
    geom_text(data = df_mad, col = cbp$navy, label = "MAD", vjust = -1, hjust = 0.7) +
    geom_segment(data = df_trimodal, col = cbp$navy, mapping = aes(x, y, xend = xend, yend = yend), size = 1) +
    facet_wrap(vars(caption), scales = "free_y", ncol = 2) +
    labs(title = "QAD function for various distributions", x = "p", y = "QAD")
}

figure_qad_efficiency <- function() {
  step <- 0.001
  p <- seq(step, 1 - step, by = step)
  e <- qad_asymptotic_efficiency(p)
  p <- c(0, p, 1)
  e <- c(0, e, 0)
  e0 <- qad_asymptotic_efficiency(Po)
  ggplot(data.frame(p, e), aes(p, e)) +
    geom_line() +
    geom_point(x = Po, y = e0, col = cbp$navy, shape = 8) +
    labs(y = "Gaussian efficiency", title = "Asymptotic Gaussian efficiency of QAD(X, p)")
}

# Figures (Simulations) --------------------------------------------------------

figure_scale_constants1 <- function(estimator, color) {
  df <- simulation_constants()
  df <- df[df$n <= 100, ]
  df$value <- df[, estimator]

  letter <- substr(estimator, 1, 1)
  y_breaks <- round(c(get_asymptotic_constant(estimator), seq(round(min(df$value) + 0.1, 1), round(max(df$value), 1), by = 0.1)), 4)
  y_breaks_labels <- c(sprintf("%0.4f", y_breaks[1]), sprintf("%0.3f", tail(y_breaks, -1)))

  ggplot(df, aes(n, value)) +
    geom_hline(yintercept = get_asymptotic_constant(estimator), linetype = "dashed", col = color) +
    geom_point(col = color) +
    scale_color_manual(values = c(cbp$blue, cbp$green), labels = c("SQAD", "OQAD")) +
    scale_y_continuous(breaks = y_breaks, labels = y_breaks_labels) +
    labs(
      title = TeX(paste0("(a) ", toupper(estimator), " consistency constants ($n\\leq 100$)")),
      y = TeX(paste0("$K_{", letter, ",n}$")),
      x = "Sample size (n)",
      col = "Estimator"
    )
}

figure_scale_constants2 <- function(estimator, color) {
  df <- data_constant_compare(estimator)
  df$value <- df[, estimator]
  df <- df[df$n >= 100, ]

  letter <- substr(estimator, 1, 1)

  y_breaks <- c(get_asymptotic_constant(estimator), seq(round(min(df$value) + 0.001, 3), round(max(df$value), 3), by = 0.001))
  y_breaks_labels <- c(sprintf("%0.4f", y_breaks[1]), sprintf("%0.3f", tail(y_breaks, -1)))

  ggplot() +
    geom_hline(yintercept = get_asymptotic_constant(estimator), linetype = "dashed", col = color) +
    geom_line(data = df, aes(n, predicted), col = cbp$grey) +
    geom_point(data = df, aes(n, value), col = color) +
    scale_y_continuous(breaks = y_breaks, labels = y_breaks_labels) +
    labs(
      title = TeX(paste0("(b) ", toupper(estimator), " consistency constants: actual and predicted ($n\\geq 100$)")),
      x = "Sample size (n)",
      y = TeX(paste0("$K_{", letter, ",n}$"))
    )
}

figure_scale_constants <- function(estimator, color) {
  p1 <- figure_scale_constants1(estimator, color)
  p2 <- figure_scale_constants2(estimator, color)
  ggarrange(p1, p2, ncol = 1)
}
figure_sqad_constants <- function() figure_scale_constants("sqad", cbp$blue)
figure_oqad_constants <- function() figure_scale_constants("oqad", cbp$green)

figure_scale_efficiencyX <- function(maxn, prefix) {
  df <- simulation_scale_efficiency() %>% extract_coumns("eff")
  df <- df[df$n <= maxn, ]
  df2 <- df %>% gather("type", "value", -n)
  df2$type <- factor(df2$type, levels = c("mad", "sqad", "oqad"))

  y_breaks <- round(sort(c(
    qad_asymptotic_efficiency("mad"), qad_asymptotic_efficiency("sqad"), qad_asymptotic_efficiency("oqad"),
    seq(round(min(df2$value), 1), round(max(df2$value), 1), by = 0.1)
  )), 4)

  ggplot(df2, aes(n, value, col = type, shape = type)) +
    geom_point(size = 1) +
    geom_hline(yintercept = qad_asymptotic_efficiency("mad"), col = cbp$red, alpha = 0.5, linetype = "dashed") +
    geom_hline(yintercept = qad_asymptotic_efficiency("sqad"), col = cbp$blue, alpha = 0.5, linetype = "dashed") +
    geom_hline(yintercept = qad_asymptotic_efficiency("oqad"), col = cbp$green, alpha = 0.5, linetype = "dashed") +
    labs(
      title = TeX(paste0(prefix, " Gaussian efficiency of MAD, SQAD, OQAD (n \\leq ", maxn, ")")),
      x = "Sample size (n)",
      y = "Relative efficiency",
      col = "Estimator",
      shape = "Estimator"
    ) +
    scale_color_manual(values = c(cbp$red, cbp$blue, cbp$green), labels = c("MAD", "SQAD", "OQAD")) +
    scale_shape_manual(values = c(16, 17, 18), labels = c("MAD", "SQAD", "OQAD")) +
    scale_y_continuous(breaks = y_breaks)
}

figure_scale_efficiency1 <- function() figure_scale_efficiencyX(100, "(a)")

figure_scale_efficiency2 <- function() figure_scale_efficiencyX(1000, "(b)")

figure_scale_efficiency <- function() {
  p1 <- figure_scale_efficiency1()
  p2 <- figure_scale_efficiency2()
  ggarrange(p1, p2, ncol = 1)
}

figure_location_efficiency1 <- function() {
  df <- simulation_location_efficiency() %>% gather("estimator", "value", -n)
  df <- df[df$n <= 100, ]
  df$estimator <- factor(df$estimator, levels = c("median", "sthdme", "othdme"))
  df$parity <- factor(ifelse(df$n %% 2 == 0, "Even", "Odd"), levels = c("Even", "Odd"))

  y_breaks <- c(round(2 / pi, 3), seq(0.6, max(df$value), by = 0.1))
  ggplot(df, aes(n, value, col = estimator, shape = estimator)) +
    geom_line(aes(group = interaction(estimator, parity)), alpha = 0.3) +
    geom_point() +
    scale_color_manual(values = c(cbp$red, cbp$blue, cbp$green), labels = c("SM", "STHDME", "OTHDME")) +
    scale_shape_manual(values = c(16, 17, 18), labels = c("SM", "STHDME", "OTHDME")) +
    scale_y_continuous(breaks = y_breaks, labels = TeX(c("$2/\\pi$", tail(y_breaks, -1)))) +
    geom_hline(yintercept = 2 / pi, col = cbp$red, alpha = 0.5, linetype = "dashed") +
    labs(
      title = TeX("(a) Gaussian efficiency of SM, STHDME, OTHDME ($n \\leq 100$)"),
      col = "Estimator", shape = "Estimator",
      y = "Gaussian efficiency",
      x = "Sample size (n)"
    )
}

figure_location_efficiency2 <- function(maxn = 10000) {
  df <- simulation_location_efficiency() %>% gather("estimator", "value", -n)
  df$estimator <- factor(df$estimator, levels = c("median", "sthdme", "othdme"))
  df <- df[df$n > 100 & df$n <= 10000, ]

  y_breaks <- c(round(2 / pi, 3), seq(0.64, max(df$value), by = 0.01))
  ggplot(df, aes(n, value, col = estimator, shape = estimator)) +
    geom_point() +
    scale_y_continuous(breaks = y_breaks, labels = TeX(c("$2/\\pi$", tail(y_breaks, -1)))) +
    scale_color_manual(values = c(cbp$red, cbp$blue, cbp$green), labels = c("SM", "STHDME", "OTHDME")) +
    scale_shape_manual(values = c(16, 17, 18), labels = c("SM", "STHDME", "OTHDME")) +
    geom_hline(yintercept = 2 / pi, col = cbp$red, alpha = 0.5, linetype = "dashed") +
    labs(
      title = TeX(paste0("(b) Gaussian efficiency of SM, STHDME, OTHDME ($100 < n \\leq ", max(df$n), "$)")),
      col = "Estimator", shape = "Estimator",
      y = "Gaussian efficiency",
      x = "Sample size (n)"
    )
}

figure_location_efficiency <- function() {
  p1 <- figure_location_efficiency1()
  p2 <- figure_location_efficiency2()
  ggarrange(p1, p2, ncol = 1)
}

simulation_constants()
simulation_scale_efficiency()
simulation_location_efficiency()
