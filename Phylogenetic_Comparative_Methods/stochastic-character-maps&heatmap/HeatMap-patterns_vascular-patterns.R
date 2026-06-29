###############################################################################
# Phylogenetic Trait Evolution: Transition Rate Visualization
# Author: Joyce Onyenedum with help from ChatGPT
# Description: extracts transition rates for vascular patterns, and visualizes them as a heatmap.
###############################################################################

# ---- Load required packages ----
library(phytools)
library(corHMM)
library(grDevices)
library(geiger)
library(reshape2)
library(ggplot2)

# ---- Set working directory ----
setwd("Desktop/stochastic-character-maps&heatmap-FINAL/")

# ---- Format for corHMM ----
# Define your intended biological order of states
state_levels <- c(
  "Typical", 
  "Divided", 
  "Compound", 
  "Phloem-wedges", 
  "Fissured", 
  "Lobed", 
  "Ectopic-cambia", 
  "Divided+Ectopic-cambia", 
  "Compound+Ectopic-cambia", 
  "Compound+Fissured+Ectopic-cambia", 
  "Compound+Phloem-wedges", 
  "Compound+Phloem-wedges+Ectopic-cambia", 
  "Phloem-wedges+Ectopic-cambia", 
  "Lobed+Phloem-wedges"
)

# ---- Summarize SIMMAP results ----
simmap_trees <- readRDS("simmap.trees_vp.RDS")
summary_stm <- describe.simmap(simmap_trees)
print(summary_stm)

###############################################################################
# Transition Matrix Parsing and Heatmap Visualization
###############################################################################

# ---- Raw text output copied from summary ----
raw_text <- "changes are of the following types:
      1,10 1,11  1,12  1,13  1,14   1,2   1,3   1,4   1,5   1,6   1,7   1,8   1,9  10,1 10,11 10,12 10,13 10,14  10,2  10,3  10,4  10,5  10,6 10,7  10,8 10,9  11,1
x->y 0.348  0.5 0.326 0.719 2.459 0.896 2.738 7.716 0.881 5.324 5.241 0.899 1.279 0.224 0.017 0.011 0.014  0.06 0.016 0.039 0.142 0.019 0.088 0.07 0.015 0.02 0.343
     11,10 11,12 11,13 11,14  11,2  11,3  11,4  11,5  11,6  11,7  11,8  11,9  12,1 12,10 12,11 12,13 12,14 12,2  12,3  12,4  12,5  12,6 12,7  12,8 12,9  13,1 13,10
x->y 0.023 0.015 0.016 0.052 0.018 0.053 0.129 0.026 0.076 0.069 0.034 0.016 0.242 0.027 0.029 0.023 0.057 0.02 0.039 0.141 0.014 0.075 0.06 0.023 0.02 0.284 0.027
     13,11 13,12 13,14  13,2  13,3  13,4  13,5  13,6  13,7  13,8  13,9  14,1 14,10 14,11 14,12 14,13  14,2  14,3  14,4  14,5  14,6  14,7  14,8  14,9   2,1  2,10
x->y 0.027 0.021 0.068 0.021 0.048 0.272 0.025 0.069 0.077 0.022 0.022 0.271 0.039 0.024 0.031  0.03 0.018 0.058 0.365 0.228 0.445 0.083 0.032 0.017 0.272 0.021
      2,11  2,12  2,13  2,14   2,3   2,4   2,5   2,6  2,7   2,8   2,9   3,1  3,10  3,11  3,12  3,13  3,14   3,2   3,4   3,5   3,6   3,7   3,8  3,9   4,1  4,10  4,11
x->y 0.049 0.027 0.023 0.068 0.052 0.162 0.022 0.079 0.07 0.828 0.023 9.276 0.967 0.827 1.085 0.526 0.447 0.208 4.795 0.216 1.284 0.274 0.188 1.13 1.383 0.031 0.047
      4,12  4,13  4,14   4,2   4,3   4,5   4,6   4,7   4,8   4,9   5,1  5,10 5,11  5,12 5,13  5,14   5,2   5,3   5,4   5,6   5,7   5,8   5,9   6,1  6,10  6,11 6,12
x->y 0.028 0.189 0.268 0.041 0.157 0.035 0.198 0.976 0.039 0.021 0.268 0.027 0.02 0.023 0.02 0.272 0.021 0.043 0.132 0.135 0.069 0.015 0.019 3.393 0.075 0.077 0.06
      6,13  6,14   6,2   6,3   6,4  6,5   6,7   6,8   6,9   7,1  7,10  7,11  7,12  7,13  7,14   7,2   7,3   7,4   7,5   7,6   7,8   7,9   8,1  8,10  8,11  8,12  8,13
x->y 0.072 1.672 0.079 2.434 1.243 0.13 0.191 0.081 0.081 1.576 0.037 0.024 0.036 0.022 0.085 0.033 0.058 0.186 0.027 0.114 0.035 0.027 0.247 0.024 0.032 0.019 0.021
      8,14   8,2  8,3   8,4   8,5   8,6   8,7   8,9   9,1  9,10  9,11  9,12 9,13  9,14   9,2   9,3   9,4   9,5  9,6   9,7   9,8
x->y 0.069 0.858 0.04 0.158 0.026 0.084 0.073 0.019 0.251 0.023 0.011 0.012 0.02 0.061 0.023 0.047 0.141 0.029 0.07 0.062 0.024"

# ---- Clean and parse the raw text ----
lines <- trimws(unlist(strsplit(raw_text, "\n")))
lines <- lines[lines != "" & !grepl("^changes are", lines, ignore.case = TRUE)]

transition_lines <- lines[!grepl("^x->y", lines)]
rate_lines <- lines[grepl("^x->y", lines)]

transitions <- unlist(strsplit(paste(transition_lines, collapse = " "), "\\s+"))
transitions <- transitions[transitions != ""]

rates <- unlist(strsplit(paste(rate_lines, collapse = " "), "\\s+"))
rates <- as.numeric(rates[!grepl("x->y", rates)])

# ---- Build transition data frame ----
transition_df <- do.call(rbind, strsplit(transitions, ",")) |> as.data.frame()
colnames(transition_df) <- c("From", "To")
transition_df$Rate <- rates
transition_df$From <- as.integer(transition_df$From)
transition_df$To <- as.integer(transition_df$To)

# ---- Construct full 14x14 matrix ----
states <- 1:14
rate_matrix <- matrix(0, nrow = length(states), ncol = length(states),
                      dimnames = list(states, states))

for (i in seq_len(nrow(transition_df))) {
  f <- as.character(transition_df$From[i])
  t <- as.character(transition_df$To[i])
  rate_matrix[f, t] <- transition_df$Rate[i]
}

###############################################################################
# ---- HEATMAP VISUALIZATION ----
###############################################################################


rownames(rate_matrix) <- state_levels
colnames(rate_matrix) <- state_levels

# ---- Melt matrix for plotting ----
df <- melt(rate_matrix, varnames = c("From", "To"), value.name = "Rate")
df$fill_color <- with(df, ifelse(From == To, NA, ifelse(Rate == 0, "lightgrey", Rate)))

# ---- Plot heatmap ----
# ---- Melt matrix for plotting ----
df <- melt(rate_matrix, varnames = c("From", "To"), value.name = "Rate")

# Make fill_color purely numeric; leave diagonals as NA
df$fill_color <- ifelse(df$From == df$To, NA, df$Rate)

# ---- Plot heatmap ----
ggplot(df, aes(x = To, y = From, fill = fill_color)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(is.na(fill_color), "", round(Rate, 2))), size = 5) +
  scale_fill_gradientn(
    colours = c("lightgrey", "lightyellow", "tomato"),
    values = scales::rescale(c(0, 0.6, max(df$Rate, na.rm = TRUE))),
    limits = c(0, max(df$Rate, na.rm = TRUE)),
    name = "Rate",
    na.value = "darkgrey"
  ) +
  coord_fixed() +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev(state_levels)) +
  labs(x = "To State", y = "From State") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 15, angle = 45, hjust = 0, vjust = 0.02,
                               margin = margin(b = -5)),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(margin = margin(t = -5)),
    axis.title.y = element_text(margin = margin(r = 5)),
    panel.grid = element_blank(),
    plot.title = element_blank()
    )

###############################################################################
# ---- GENERATE & EXPORT TABLE OF RESULTS ----
###############################################################################

# 1. Convert the matrix into a structured data frame
heatmap_table <- as.data.frame(rate_matrix)

# 2. Print the raw table cleanly to the R Console
cat("\n--- TRANSITION RATE MATRIX TABLE ---\n")
print(round(heatmap_table, 3))

# 3.Export  table to a CSV 
write.csv(heatmap_table, "heatmap_vascular_patterns_rates.csv", row.names = TRUE)

