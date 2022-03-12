library(readr)
library(dplyr)
library(ggplot2)

# Environment preparation

if (!dir.exists("output")) {
  dir.create("output")
} else {
  unlink("output", recursive = TRUE)
  dir.create("output")
}

quadratic_formula <- function(a, b, c, mod) {
  up_term <- -b + (sqrt(b^2 - 4*a*c) * mod)
  down_term <- 2 * a
  return(up_term / down_term)
}

# Nominal ranks preparation

medium_ranks <- c("MS", "VW", "KC")
pgr_ranks <- c("TDZ", "BA/BAP", "Kn")
explant_ranks <- c("PLB", "Shoot")

# Thesis data preparation

thesis <- read_csv("data/thesis.csv", show_col_types = FALSE)
thesis_survival_rate <- data.frame(
  survival_rate = c(100, 84.3, 95, 80, 90, 84.3, 89, 100, 84, 90, 90, 97, 75, 75, 100, 80, 80, 90, 100, 95),
  pgr = ""
)

thesis_recoded <- thesis

thesis_recoded$pgr <- factor(recode(thesis_recoded$pgr, tdz = "TDZ", bap = "BA/BAP", kn = "Kn"))
thesis_recoded$medium <- factor(recode(thesis_recoded$medium, ms = "MS", vw = "VW", kc = "KC"))
thesis_recoded$explant <- factor(recode(thesis_recoded$explant, plb = "PLB", shoot = "Shoot"))

thesis_filtered <- thesis %>%
  filter(pgr_concentration < 6) %>%
  filter(yield > -10)

thesis_filtered$pgr <- recode(thesis_filtered$pgr, tdz = "TDZ", bap = "BA/BAP", kn = "Kn")
thesis_filtered$medium <- recode(thesis_filtered$medium, ms = "MS", vw = "VW", kc = "KC")
thesis_filtered$explant <- recode(thesis_filtered$explant, plb = "PLB", shoot = "Shoot")

thesis_regression_model <- thesis_filtered$yield ~ poly(thesis_filtered$pgr_concentration, 2, raw = TRUE)

# Setup output sink

output_file <- file("output/output.txt", open = "wt")
sink(output_file)
sink(output_file, type = "message")

# PGR Regression plotting method
pgr_regression <- function(data, explant_name, filename) {
  ggplot(data, mapping = aes(x = pgr_concentration, y = yield, color = medium)) +
    labs(x = "PGR Concentration (mg/L)", y = paste(explant_name, "produced per explant"), shape = "PGR", color = "Medium") +
    geom_smooth(inherit.aes = FALSE,
                method = "lm",
                formula = y ~ poly(x, 2, raw = TRUE),
                se = FALSE,
                fill = NA,
                alpha = 0.5,
                aes(x = pgr_concentration, y = yield, color = medium)) +
    geom_point(alpha = 0.5, size = 1.5) +
    facet_grid(pgr ~ .) +
    theme(text = element_text(size = 11, family = "Serif"),
          axis.title.y = element_text(vjust = 1),
          axis.title.x = element_text(vjust = -1),
          axis.text = element_text(size = rel(0.7)))
  ggsave(paste0("output/", filename), width = 6.5, height = 6.5, units = "in", dpi = "print", type = "cairo-png")
}

# PGR Box plotting method
pgr_boxplot <- function(data, explant_name, filename) {
  ggplot(data, mapping = aes(x = medium, y = yield)) +
    labs(x = "Medium", y = paste(explant_name, "produced per explant"), shape = "PGR", color = "Medium") +
    geom_boxplot(show.legend = FALSE, aes(color = medium)) +
    facet_grid(. ~ pgr) +
    theme(text = element_text(size = 11, family = "Serif"),
          axis.text = element_text(size = rel(0.8)))
  ggsave(paste0("output/", filename), width = 6.5, height = 6.5, units = "in", dpi = "print", type = "cairo-png")
}

# Outlier
ggplot(thesis, mapping = aes(x = pgr_concentration, y = yield)) +
  labs(x = "PGR Concentration (mg/L)", y = "Regenerants produced per explant") +
  geom_point(alpha = 0.5, size = 1.5) +
  theme(text = element_text(size = 10, family = "Serif"))
ggsave("output/outlier.png", width = 6.5, height = 3.5, units = "in", dpi = "print", type = "cairo-png")

# PGR Regression
pgr_regression(thesis_filtered, "Regenerants", "pgr_regression_combined.png")
pgr_regression(thesis_filtered %>% filter(explant == "PLB"), "PLBs", "pgr_regression_plb.png")
pgr_regression(thesis_filtered %>% filter(explant == "Shoot"), "Shoots", "pgr_regression_shoot.png")

# PGR Box plot
pgr_boxplot(thesis_filtered, "Regenerants", "pgr_boxplot_combined.png")
pgr_boxplot(thesis_filtered %>% filter(explant == "PLB"), "PLBs", "pgr_boxplot_plb.png")
pgr_boxplot(thesis_filtered %>% filter(explant == "Shoot"), "Shoots", "pgr_boxplot_shoot.png")

# Regression detail report
regression_detail_report <- function(data, medium_filter, pgr_filter, explant_filter) {
  medium_filter_text <- "*"
  pgr_filter_text <- "*"
  explant_filter_text <- "*"
  if (medium_filter != FALSE) {
    data <- data %>% filter(medium == medium_filter)
    medium_filter_text <- medium_filter
  }
  if (pgr_filter != FALSE) {
    data <- data %>% filter(pgr == pgr_filter)
    pgr_filter_text <- pgr_filter
  }
  if (explant_filter != FALSE) {
    data <- data %>% filter(explant == explant_filter)
    explant_filter_text <- explant_filter
  }

  message(paste0(medium_filter_text, "-", pgr_filter_text, "-", explant_filter_text))
  print(summary(data))
  try({
    data_regression = lm(yield ~ poly(pgr_concentration, 2, raw = TRUE), data = data)
    print(summary(data_regression))
    a <- coef(data_regression)[3]
    b <- coef(data_regression)[2]
    c <- coef(data_regression)[1]
    cat(paste0("Model: y = ", a, "x^2 + (", b, "x) + (", c, ")\n"))
    discriminant <- b^2 - 4 * a * c
    if (discriminant == 0) {
      cat(paste0("Roots: ", quadratic_formula(a, b, c, 1), "\n"))
    }
    else if (discriminant > 0) {
      cat(paste0("Roots: ", quadratic_formula(a, b, c, 1), ", ", quadratic_formula(a, b, c, -1), "\n"))
    }
  })
}

process_report <- function(
  medium_filters = vector(length = 1),
  pgr_filters = vector(length = 1),
  explant_filters = vector(length = 1)
) {
  for (medium_filter in medium_filters) {
    for (pgr_filter in pgr_filters) {
      for (explant_filter in explant_filters) {
        regression_detail_report(thesis_filtered, medium_filter, pgr_filter, explant_filter)
        cat("\n")
      }
    }
  }
}

process_report()
process_report(pgr_filters = pgr_ranks)
process_report(pgr_filters = pgr_ranks, medium_filters = medium_ranks)
process_report(pgr_filters = pgr_ranks, explant_filters = explant_ranks)
process_report(medium_rank)
process_report(pgr_filters = pgr_ranks, medium_filters = medium_ranks, explant_filters = explant_ranks)

cat("\nSignificance Model: \n")
significance_model <- lm(yield ~ poly(pgr_concentration, 2, raw = TRUE) + pgr + medium + explant, data = thesis_recoded)
print(summary(significance_model))

cost_efficiency <- data.frame(
  name=c("MS/TDZ", "MS/6-BAP", "MS/KN", "VW/TDZ", "VW/6-BAP", "KC/TDZ", "KC/6-BAP"),
  value=c(105.526177, 168.9452874, 147.8822365, 33.22136364, 32.75779221, 45.32821916, 157.25)
)

ggplot(cost_efficiency, aes(x=name, y=value)) +
  labs(x = "Medium/PGR", y = "Cost-efficiency (lower is better)") +
  geom_bar(stat = "identity", fill=rgb(0.1, 0.4, 0.5, 0.7)) +
  theme(text = element_text(size = 11, family = "Serif"),
        axis.title.y = element_text(vjust = 1),
        axis.title.x = element_text(vjust = -1),
        axis.text = element_text(size = rel(0.7)))

ggsave("output/cost_efficiency.png", width = 6.5, height = 4.0, units = "in", dpi = "print", type = "cairo-png") 

sink()
sink(type = "message")
close(output_file)