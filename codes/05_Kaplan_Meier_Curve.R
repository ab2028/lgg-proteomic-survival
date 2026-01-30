# 5. Kaplan Meier Curve

library(survival)
library(survminer)
library(ggplot2)

# make sure directory exists
if (!dir.exists("figures")) dir.create("figures")

km <- survfit(Surv(OS.time, OS) ~ 1, data = df_model)

pdf(
  file = "figures/km_curve.pdf",
  width = 6,
  height = 4.5
)

p <- ggsurvplot(
  km,
  data = df_model,
  xlab = "Survival Time (days)",
  ylab = "Overall survival probability",
  title = "Kaplan-Meier Estimate of Overall Survival",
  conf.int = TRUE,
  conf.int.style = "ribbon",
  conf.int.alpha = 0.25,
  palette = "#3d7ca8",
  conf.int.fill = "#2596be",
  size = 1.2,
  censor = TRUE,
  censor.shape = "|",
  censor.size = 2.5,
  ggtheme = theme_classic(),
  legend = "none"
)

# Modify the ggplot
p$plot <- p$plot +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5
    ),
    axis.title.x = element_text(
      margin = margin(t = 10)
    ),
    axis.title.y = element_text(
      margin = margin(r = 10)
    )
  )

# save
print(p$plot)

dev.off()

# print to plots tab
print(p$plot)