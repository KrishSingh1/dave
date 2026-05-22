#!/usr/bin/env Rscript

# Plot mortality equations used in mortality_guess().
# Base-R only; no package dependencies.

# ----------------------------
# User-tunable settings
# ----------------------------

# Equation defaults from code/documentation.
f_survive <- 0.001
q_shape <- 2

# Longevity values (years) to compare.
longevity_values <- c(100, 200, 400, 800)

# Growth-efficiency mortality settings.
k_mort <- 0.3
greff_min_values <- c(0.01, 0.02, 0.05, 0.1)

# Output settings.
out_dir <- "docs/science/img/mortality"
out_file_age_light <- file.path(out_dir, "age_mortality_curve-light.svg")
out_file_age_dark <- file.path(out_dir, "age_mortality_curve-dark.svg")
out_file_greff_smooth_light <- file.path(out_dir, "greff_mortality_smooth-light.svg")
out_file_greff_smooth_dark <- file.path(out_dir, "greff_mortality_smooth-dark.svg")
out_file_greff_step_light <- file.path(out_dir, "greff_mortality_step-light.svg")
out_file_greff_step_dark <- file.path(out_dir, "greff_mortality_step-dark.svg")
out_file_greff_by_min_light <- file.path(out_dir, "greff_mortality_vs_greff_by_greffmin-light.svg")
out_file_greff_by_min_dark <- file.path(out_dir, "greff_mortality_vs_greff_by_greffmin-dark.svg")
out_file_crowd_light <- file.path(out_dir, "crowd_mortality_surface-light.svg")
out_file_crowd_dark <- file.path(out_dir, "crowd_mortality_surface-dark.svg")
bg_col <- "transparent"

# Domain for self-thinning mortality plot.
crowd_ac_min <- 0.05
crowd_ac_max <- 4.0
crowd_d_min <- 0.05
crowd_d_max <- 4.0
crowd_points <- 240

# ----------------------------
# Functions
# ----------------------------

mortality_age <- function(age, longevity, f = f_survive, q = q_shape) {
    -log(f) * ((q + 1) / longevity) * (age / longevity)^q
}

mortality_greff_smooth <- function(greff_ratio, k = k_mort) {
    k / (1 + greff_ratio^5)
}

mortality_greff_step <- function(greff_ratio, k = k_mort) {
    ifelse(greff_ratio < 1, k, 0)
}

mortality_crowd <- function(crown_area, density) {
    crowding <- crown_area * density
    ifelse(crowding > 1, (crowding - 1) / crowding, 0)
}

line_cols <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02")
line_cols <- rep(line_cols, length.out = length(longevity_values))

if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

plot_age_variant <- function(out_file, ink_col, grid_col) {
    svg(out_file, width = 10, height = 6, bg = bg_col)
    par(mar = c(5.5, 5.5, 5.0, 2.5),
        xaxs = "i",
        yaxs = "i",
        bg = bg_col,
        fg = ink_col,
        col.axis = ink_col,
        col.lab = ink_col,
        col.main = ink_col)

    # Use a common relative-age domain so curves are comparable.
    rel_age_vals <- seq(0, 1, length.out = 400)

    curves <- sapply(longevity_values, function(longevity) {
        ages <- rel_age_vals * longevity
        mortality_age(ages, longevity)
    })

    if (is.null(dim(curves))) {
        curves <- matrix(curves, ncol = 1)
    }

    ylim_max <- max(curves, na.rm = TRUE) * 1.05

    plot(rel_age_vals,
         curves[, 1],
         type = "l",
         lwd = 3,
         col = line_cols[1],
         xlab = "Relative age (age / longevity)",
         ylab = "Fraction of cohort killed per year",
         main = "Age-related mortality",
         ylim = c(0, ylim_max))

    if (ncol(curves) > 1) {
        for (i in 2:ncol(curves)) {
            lines(rel_age_vals, curves[, i], lwd = 3, col = line_cols[i])
        }
    }

    abline(h = 0, col = grid_col, lty = 2)
    box(col = ink_col)

    legend("topleft",
           legend = paste0("longevity = ", longevity_values, " yr"),
           col = line_cols[seq_along(longevity_values)],
           lwd = 3,
           bty = "n",
           cex = 1.0,
           text.col = ink_col)

    caption <- paste(sprintf("F=%.3g", f_survive),
                     sprintf("Q=%g", q_shape),
                     "m_age = -ln(F) * (Q+1)/L * (a/L)^Q",
                     sep = "  ")
    mtext(caption, side = 1, line = 4.0, cex = 0.85)

    dev.off()
    message("Wrote: ", out_file)
}

plot_greff_smooth_variant <- function(out_file, ink_col, grid_col) {
    svg(out_file, width = 10, height = 6, bg = bg_col)
    par(mar = c(5.5, 5.5, 5.0, 2.5),
        xaxs = "i",
        yaxs = "i",
        bg = bg_col,
        fg = ink_col,
        col.axis = ink_col,
        col.lab = ink_col,
        col.main = ink_col)

    greff_ratio <- seq(0, 3, length.out = 500)
    mort_vals <- mortality_greff_smooth(greff_ratio)

    plot(greff_ratio,
         mort_vals,
         type = "l",
         lwd = 3,
         col = "#1b9e77",
         xlab = "Relative growth efficiency (greff_5 / greff_min)",
         ylab = "Fraction of cohort killed per year",
         main = "Growth-efficiency mortality (smooth mode)",
         ylim = c(0, max(mort_vals) * 1.05))

    abline(v = 1, col = grid_col, lty = 3)
    abline(h = 0, col = grid_col, lty = 2)
    box(col = ink_col)

    caption <- paste(sprintf("K=%.3g", k_mort),
                     "m_greff = K / (1 + (greff_5 / greff_min)^5)",
                     sep = "  ")
    mtext(caption, side = 1, line = 4.0, cex = 0.85)

    dev.off()
    message("Wrote: ", out_file)
}

plot_greff_step_variant <- function(out_file, ink_col, grid_col) {
    svg(out_file, width = 10, height = 6, bg = bg_col)
    par(mar = c(5.5, 5.5, 5.0, 2.5),
        xaxs = "i",
        yaxs = "i",
        bg = bg_col,
        fg = ink_col,
        col.axis = ink_col,
        col.lab = ink_col,
        col.main = ink_col)

    greff_ratio <- seq(0, 3, length.out = 500)
    mort_vals <- mortality_greff_step(greff_ratio)

    plot(greff_ratio,
         mort_vals,
         type = "l",
         lwd = 3,
         col = "#d95f02",
         xlab = "Relative growth efficiency (greff_5 / greff_min)",
         ylab = "Fraction of cohort killed per year",
         main = "Growth-efficiency mortality (non-smooth mode)",
         ylim = c(0, max(k_mort, mort_vals) * 1.05))

    abline(v = 1, col = grid_col, lty = 3)
    abline(h = 0, col = grid_col, lty = 2)
    box(col = ink_col)

    caption <- paste(sprintf("K=%.3g", k_mort),
                     "m_greff = K if greff_5/greff_min < 1, else 0",
                     sep = "  ")
    mtext(caption, side = 1, line = 4.0, cex = 0.85)

    dev.off()
    message("Wrote: ", out_file)
}

plot_greff_by_min_variant <- function(out_file, ink_col, grid_col) {
    svg(out_file, width = 10, height = 6, bg = bg_col)
    par(mar = c(5.5, 5.5, 5.0, 2.5),
        xaxs = "i",
        yaxs = "i",
        bg = bg_col,
        fg = ink_col,
        col.axis = ink_col,
        col.lab = ink_col,
        col.main = ink_col)

    greff_vals <- seq(0, max(greff_min_values) * 3, length.out = 600)
    curves <- sapply(greff_min_values, function(greff_min_val) {
        mortality_greff_smooth(greff_vals / greff_min_val)
    })

    if (is.null(dim(curves))) {
        curves <- matrix(curves, ncol = 1)
    }

    plot(greff_vals,
         curves[, 1],
         type = "l",
         lwd = 3,
         col = line_cols[1],
         xlab = "Growth efficiency, greff",
         ylab = expression(m[greff]),
         main = expression("Growth-efficiency mortality vs greff (varying " * greff[min] * ")"),
         ylim = c(0, k_mort * 1.05))

    if (ncol(curves) > 1) {
        for (i in 2:ncol(curves)) {
            lines(greff_vals, curves[, i], lwd = 3, col = line_cols[i])
        }
    }

    abline(h = 0, col = grid_col, lty = 2)
    box(col = ink_col)

    legend("topright",
           legend = paste0("greff_min = ", format(greff_min_values, trim = TRUE)),
           col = line_cols[seq_along(greff_min_values)],
           lwd = 3,
           bty = "n",
           cex = 0.95,
           text.col = ink_col)

    caption <- paste(sprintf("K=%.3g", k_mort),
                     "Smooth mode: m_greff = K / (1 + (greff / greff_min)^5)",
                     sep = "  ")
    mtext(caption, side = 1, line = 4.0, cex = 0.85)

    dev.off()
    message("Wrote: ", out_file)
}

plot_crowd_variant <- function(out_file, ink_col, grid_col, palette_name) {
    svg(out_file, width = 10, height = 7, bg = bg_col)
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)

    # Layout with main panel and color bar panel.
    layout(matrix(c(1, 2), nrow = 1), widths = c(4.8, 0.8))

    par(mar = c(4, 5.2, 4.0, 0.8),
        xaxs = "i",
        yaxs = "i",
        bg = bg_col,
        fg = ink_col,
        col.axis = ink_col,
        col.lab = ink_col,
        col.main = ink_col)

    crown_area_vals <- seq(crowd_ac_min, crowd_ac_max, length.out = crowd_points)
    density_vals <- seq(crowd_d_min, crowd_d_max, length.out = crowd_points)

    z_vals <- outer(crown_area_vals, density_vals, mortality_crowd)
     palette_cols <- hcl.colors(80, palette_name)

     image(crown_area_vals,
           density_vals,
           z_vals,
           col = palette_cols,
           xlab = expression("Crown area, " * A[c] * " (" * m^2 * ")"),
           ylab = expression("Density, " * D * " (individuals per " * m^2 * ")"),
           main = expression("Self-thinning mortality surface (" * m[crowd] * ")"),
           zlim = c(0, 1),
           useRaster = TRUE)

    contour(crown_area_vals,
            density_vals,
            z_vals,
            add = TRUE,
            levels = seq(0.1, 0.9, by = 0.1),
            drawlabels = FALSE,
            col = adjustcolor(ink_col, alpha.f = 0.45),
            lwd = 0.8)

    # Threshold where crowding starts: A_c * D = 1.
    d_threshold <- 1 / crown_area_vals
    valid <- d_threshold >= crowd_d_min & d_threshold <= crowd_d_max
    lines(crown_area_vals[valid], d_threshold[valid], col = ink_col, lwd = 2, lty = 2)

    # One-row horizontal legend above (outside) the main plot area.
    legend("topleft",
           inset = c(0, -0.06),
           legend = as.expression(bquote(A[c] * D == 1 ~ "threshold")),
           lty = 2,
           lwd = 2,
           col = ink_col,
           bty = "n",
           text.col = ink_col,
           cex = 0.95,
           horiz = TRUE,
           xpd = NA)

    # Color bar panel.
    par(mar = c(4, 0.2, 4.0, 4.6),
        xaxs = "i",
        yaxs = "i",
        bg = bg_col,
        fg = ink_col,
        col.axis = ink_col,
        col.lab = ink_col,
        col.main = ink_col)

    zbar <- matrix(seq(0, 1, length.out = 200), nrow = 1)
    image(x = c(0, 1),
          y = seq(0, 1, length.out = ncol(zbar)),
          z = zbar,
          col = palette_cols,
          axes = FALSE,
          xlab = "",
          ylab = "",
          useRaster = TRUE)
    axis(4, at = seq(0, 1, by = 0.2), labels = sprintf("%.1f", seq(0, 1, by = 0.2)), las = 1)
    mtext(expression(m[crowd]), side = 4, line = 1.9)
    box(col = ink_col)

    dev.off()
    message("Wrote: ", out_file)
}

# Light-theme variants.
plot_age_variant(
    out_file = out_file_age_light,
    ink_col = "#2e2e2e",
    grid_col = "#9a9a9a"
)

plot_greff_smooth_variant(
    out_file = out_file_greff_smooth_light,
    ink_col = "#2e2e2e",
    grid_col = "#9a9a9a"
)

plot_greff_step_variant(
    out_file = out_file_greff_step_light,
    ink_col = "#2e2e2e",
    grid_col = "#9a9a9a"
)

plot_greff_by_min_variant(
    out_file = out_file_greff_by_min_light,
    ink_col = "#2e2e2e",
    grid_col = "#9a9a9a"
)

plot_crowd_variant(
    out_file = out_file_crowd_light,
    ink_col = "#2e2e2e",
    grid_col = "#9a9a9a",
    palette_name = "YlOrRd"
)

# Dark-theme variants.
plot_age_variant(
    out_file = out_file_age_dark,
    ink_col = "#e8e8e8",
    grid_col = "#bdbdbd"
)

plot_greff_smooth_variant(
    out_file = out_file_greff_smooth_dark,
    ink_col = "#e8e8e8",
    grid_col = "#bdbdbd"
)

plot_greff_step_variant(
    out_file = out_file_greff_step_dark,
    ink_col = "#e8e8e8",
    grid_col = "#bdbdbd"
)

plot_greff_by_min_variant(
    out_file = out_file_greff_by_min_dark,
    ink_col = "#e8e8e8",
    grid_col = "#bdbdbd"
)

plot_crowd_variant(
    out_file = out_file_crowd_dark,
    ink_col = "#e8e8e8",
    grid_col = "#bdbdbd",
    palette_name = "Inferno"
)
