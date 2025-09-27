#!/usr/bin/env Rscript

### TRANSCRIPT PROCESS SUMMARIZATION ####
# Date: May 29th, 2025
# Author: MA Yaqin
# Description: regression plot

#======== Libraries 
library(ggplot2)
library(dplyr)
library(scales)

#======== Libraries 
# Directory(es)
RESDIR <- paste0(getwd(),"/result/")

#======== Data import
summary <- read.delim("processed/kallisto_summary.tsv") %>%
  mutate(log_n_processed = log10(n_processed))

#======== Linear model fitting
# Fit linear model
fit <- lm(runtime_seconds ~ log_n_processed, data = summary)
slope <- coef(fit)[2]
intercept <- coef(fit)[1]
r2 <- summary(fit)$r.squared

# Make plot
p1 <- ggplot(summary, aes(x = log_n_processed, y = runtime_seconds)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", col = "red", se = TRUE) +
  labs(title = "Kallisto Runtime vs Total Reads",
       x = "Number of Reads (log10 scale)",
       y = "Runtime (seconds)") +
  annotate("text",
           x = max(summary$log_n_processed) * 0.9,
           y = max(summary$runtime_seconds) * 0.9,
           label = paste0("Slope = ", round(slope, 3),
                          "\nR² = ", round(r2, 3)),
           hjust = 0, size = 4) +
  theme_minimal() 
p1

# Save plot
ggsave(paste0(RESDIR,"kallisto_runtime_vs_reads.png"), p1, width=6, height=4)
