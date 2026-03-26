#!/usr/bin/env Rscript

# Plot expected saplings vs F = anetps_ff / anetps_ff_max for multiple alphar
# values. Base-R only; no package dependencies.
#
# This script defaults to "no SME, propagule pool present", which isolates
# alphar effects while staying close to common parameter values in pft_dave.ins.

# ----------------------------
# User-tunable settings
# ----------------------------

# Alphar values to compare (values seen in benchmarks/ozflux/pft_dave.ins).
alphar_values <- c(2, 7, 11, 50)

# Other recruitment parameters (representative values from pft_dave.ins groups).
est_max <- 0.2
kest_repr <- 200
kest_pres <- 1.0
kest_bg <- 0.1

# Set cmass_repr so kest_repr * cmass_repr ~ 1 (as intended by model comments).
cmass_repr <- 0.005

# Patch area in m2.
patcharea <- 1.0

# Choose recruitment branch.
ifsme <- FALSE
propagule_pool_present <- TRUE

# Optional woody-PFT normalization from establishment_guess:
# est *= (3 / nwoodypfts_estab)
apply_woody_normalization <- FALSE
nwoodypfts_estab <- 3

# F range: avoid F=0 exactly (division by zero in mu).
f_min <- 0.02
f_max <- 1.0
n_points <- 500

# anetps_ff_max context (for optional top axis label conversion).
anetps_ff_max <- 4.0  # kgC m^-2 yr^-1 (typical upper range context)

# Output.
out_dir <- "docs/science/img/establishment"
out_file_light <- file.path(out_dir, "saplings_vs_F_alphar-light.svg")
out_file_dark <- file.path(out_dir, "saplings_vs_F_alphar-dark.svg")

# Transparent background and theme-specific ink colors.
bg_col <- "transparent"

# ----------------------------
# Functions
# ----------------------------

mu_fun <- function(f, alphar) {
    exp(alphar * (1 - 1 / f))
}

expected_saplings <- function(f, alphar) {
    c_val <- mu_fun(f, alphar) * est_max * patcharea

    if (ifsme) {
        est <- c_val * (kest_repr * cmass_repr + kest_bg)
    } else if (propagule_pool_present) {
        est <- c_val * (kest_pres + kest_bg)
    } else {
        est <- c_val * kest_bg
    }

    if (apply_woody_normalization) {
        est <- est * (3 / nwoodypfts_estab)
    }

    est
}

# ----------------------------
# Compute curves
# ----------------------------

f_vals <- seq(f_min, f_max, length.out = n_points)
est_curves <- sapply(alphar_values, function(a) expected_saplings(f_vals, a))

# Ensure matrix shape even for one alphar value.
if (is.null(dim(est_curves))) {
    est_curves <- matrix(est_curves, ncol = 1)
}

line_cols <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02")
line_cols <- rep(line_cols, length.out = length(alphar_values))

# ----------------------------
# Plot
# ----------------------------

if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

plot_variant <- function(out_file, ink_col, grid_col) {
    svg(out_file, width = 10, height = 6, bg = bg_col)
    par(mar = c(5.5, 5.5, 5.0, 2.5),
        xaxs = "i",
        yaxs = "i",
        bg = bg_col,
        fg = ink_col,
        col.axis = ink_col,
        col.lab = ink_col,
        col.main = ink_col)

    ylim_max <- max(est_curves, na.rm = TRUE) * 1.05
    plot(f_vals,
         est_curves[, 1],
         type = "l",
         lwd = 3,
         col = line_cols[1],
         xlab = "F = anetps_ff / anetps_ff_max",
         ylab = "Expected # saplings (/year)",
         main = "Effect of alphar on recruitment rate",
         ylim = c(0, ylim_max))

    if (ncol(est_curves) > 1) {
        for (i in 2:ncol(est_curves)) {
            lines(f_vals, est_curves[, i], lwd = 3, col = line_cols[i])
        }
    }

    abline(h = 0, col = grid_col, lty = 2)
    box(col = ink_col)

    legend("topleft",
           legend = paste0("alphar = ", alphar_values),
           col = line_cols[seq_along(alphar_values)],
           lwd = 3,
           bty = "n",
           cex = 1.0,
           text.col = ink_col)

    caption <- paste(sprintf("est_max=%.3g", est_max),
                     sprintf("kest_bg=%.3g", kest_bg),
                     sprintf("kest_pres=%.3g", kest_pres),
                     sprintf("kest_repr=%.3g", kest_repr),
                     sprintf("cmass_repr=%.3g", cmass_repr),
                     sprintf("patcharea=%.3g", patcharea),
                     sprintf("ifsme=%s", ifsme),
                     sprintf("propagule_pool_present=%s",
                             propagule_pool_present),
                     sep = "  ")
    mtext(caption, side = 1, line = 4.0, cex = 0.8)

    dev.off()
    message("Wrote: ", out_file)
}

# Light-theme variant.
plot_variant(
    out_file = out_file_light,
    ink_col = "#2e2e2e",
    grid_col = "#9a9a9a"
)

# Dark-theme variant.
plot_variant(
    out_file = out_file_dark,
    ink_col = "#e8e8e8",
    grid_col = "#bdbdbd"
)
