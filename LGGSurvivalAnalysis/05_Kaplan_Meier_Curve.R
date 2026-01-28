# 5. Kaplan Meier Curve

library(survival)
library(survminer)

km <- survfit(Surv(OS.time, OS) ~ 1, data = df_model)

#pdf("km_curve.pdf", width = 6, height = 4.5)

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

p$plot <- p$plot +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      size = 18,      
      face = "bold",
      hjust = 0.5
    )
  )

p$plot <- p$plot +
  ggplot2::theme(
    axis.title.x = ggplot2::element_text(
      margin = ggplot2::margin(t = 10)  # space above x-axis title
    ),
    axis.title.y = ggplot2::element_text(
      margin = ggplot2::margin(r = 10)  # space to the right of y-axis title
    )
  )



# dev.off()

print(p)
