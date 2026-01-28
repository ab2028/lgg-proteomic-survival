# 2. Parallel Analysis
# Note: this can be skipped as the number of factors is directly entered into step 3

library(psych)

# ----- Run parallel analysis -----
pa <- fa.parallel(
  protein_data_final,
  fm = "minres",
  fa = "fa",
  n.iter = 50,
  show.legend = TRUE,
  main = "Parallel Analysis (Full Plot)"
)

# Extract eigenvalues
observed <- pa$fa.values              # Observed eigenvalues
simulated <- pa$fa.sim               # Simulated eigenvalues (mean of random data)

# ----- Zoomed plot (first 40 components) -----

zoom_range <- 1:40

plot(zoom_range, observed[zoom_range],
     type = "b", pch = 19, col = "black", lwd = 2,
     xlab = "Component Number",
     ylab = "Eigenvalue",
     main = "Parallel Analysis Scree Plot (First 40 Components)")

lines(zoom_range, simulated[zoom_range],
      type = "b", pch = 17, col = "blue", lwd = 2, lty = 2)

# Vertical cutoff at factor = 18
abline(v = 18, col = "red", lwd = 3, lty = 2)

# Add text label
text(
  x = 8.5,                                 # move far left
  y = max(observed[zoom_range]) * 0.85,   # vertical placement unchanged
  labels = "18 factors",
  col = "red",
  cex = 1.1,
  pos = 4
)

legend("topright",
       legend = c("Observed eigenvalues", "Simulated eigenvalues"),
       col = c("black", "blue"), lwd = 2, lty = c(1, 2), pch = c(19, 17))


# The result is 18, manually entered in step 3.