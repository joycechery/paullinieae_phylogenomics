###############################################################################
# Plotting the number of transitions to vascular variants through time 
# Author: Joyce Onyenedum with help from ChatGPT
# Description: extracts transition rates for vascular categories, and visualizes them as a line plot
###############################################################################

library(phytools)
library(corHMM)
library(grDevices)
library(geiger)
library(ggplot2)
library(dplyr)

setwd("~/Desktop")


simmap.trees_vc <- readRDS("simmap.trees_vc.RDS")

# === TRANSITIONS THROUGH TIME ===
root_age <- max(nodeHeights(simmap.trees_vc[[1]]))
bin_width <- 2
time_bins <- seq(0, root_age + bin_width, by = bin_width)

get_transitions_from_typical <- function(tree) {
  transition_times <- c()
  for (i in seq_len(nrow(tree$edge))) {
    segments <- tree$maps[[i]]
    if (length(segments) < 2) next
    state_seq <- names(segments)
    for (j in seq_len(length(state_seq) - 1)) {
      if (state_seq[j] == "1" && state_seq[j + 1] != "1") {
        start_time <- nodeHeights(tree)[i, 1]
        trans_time <- sum(segments[1:j])
        transition_times <- c(transition_times, start_time + trans_time)
      }
    }
  }
  return(transition_times)
}

# Extract transitions across all simmaps
all_transitions <- lapply(simmap.trees_vc, get_transitions_from_typical)

# Bin transitions
binned_counts <- lapply(all_transitions, function(times) {
  if (length(times) == 0) return(rep(0, length(time_bins) - 1))
  transition_ages <- root_age - times
  hist(transition_ages, breaks = time_bins, plot = FALSE)$counts
})

binned_matrix <- do.call(rbind, binned_counts)
avg_per_bin <- colMeans(binned_matrix)
se_per_bin <- apply(binned_matrix, 2, sd) / sqrt(nrow(binned_matrix))

# === PLOT TRANSITION HISTOGRAM ===
plot_data <- data.frame(
  Time = time_bins[-length(time_bins)],
  AvgTransitions = avg_per_bin,
  SE = se_per_bin
) %>%
  mutate(Upper = AvgTransitions + 2 * SE,
         Lower = pmax(0, AvgTransitions - 2 * SE))

ggplot(plot_data, aes(x = Time, y = AvgTransitions)) +
  annotate("rect", xmin = 0, xmax = 17, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.4) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "#33a02c", alpha = 0.5) +
  geom_line(color = "#1b7837", size = 0.5) +
  geom_point(color = "#00441b", size = 1, alpha = 0.8) +
  scale_x_continuous(limits = c(0, 65), breaks = seq(0, 60, by = 10)) +
  labs(
    x = "Time before present (Ma)",
    y = "Mean transitions per 1 My",
    title = "Transitions from Typical to Vascular Variants Through Time"
    # subtitle = paste("Mean ± 2 SE across", n_sim, "SIMMAP trees")
  ) +
  theme_bw(base_size = 20) +
  theme(
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 20, face = "bold"),
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 20)
  )