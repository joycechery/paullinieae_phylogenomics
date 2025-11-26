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
setwd("/Users/joyce.onyenedum/Library/CloudStorage/Box-Box/Onyenedum_Lab/Joyce_Onyenedum/Projects/Paullinieae_Phylogeny/*MANUSCRIPT/10.DRAFTS/Draft_V9/Revision_v1/V1/ASR_revision")

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
     1,10  1,11  1,12  1,13  1,14  1,2   1,3   1,4   1,5   1,6   1,7   1,8   1,9  10,1 10,11 10,12 10,13 10,14  10,2  10,3  10,4  10,5 10,6  10,7  10,8  10,9  11,1
x->y 0.32 0.495 0.308 0.725 2.398 0.88 2.677 7.725 0.894 5.343 5.256 0.922 1.288 0.212 0.014 0.023 0.013 0.053 0.019 0.045 0.131 0.024 0.07 0.068 0.018 0.012 0.396
     11,10 11,12 11,13 11,14  11,2  11,3  11,4  11,5  11,6 11,7  11,8  11,9  12,1 12,10 12,11 12,13 12,14  12,2  12,3  12,4  12,5  12,6 12,7  12,8  12,9  13,1 13,10
x->y 0.015 0.017  0.02 0.059 0.019 0.039 0.125 0.025 0.065 0.06 0.027 0.027 0.231 0.013 0.014 0.018 0.084 0.019 0.044 0.129 0.023 0.072 0.07 0.018 0.018 0.297 0.026
     13,11 13,12 13,14  13,2  13,3  13,4  13,5  13,6  13,7  13,8  13,9  14,1 14,10 14,11 14,12 14,13  14,2  14,3  14,4  14,5  14,6  14,7  14,8  14,9   2,1 2,10  2,11
x->y 0.021 0.027 0.053 0.022 0.037 0.287 0.026 0.069 0.063 0.015 0.016 0.275 0.027 0.035 0.035 0.022 0.026 0.057 0.346 0.249 0.424 0.069 0.024 0.035 0.266 0.03 0.033
      2,12  2,13  2,14   2,3   2,4   2,5   2,6   2,7   2,8   2,9   3,1  3,10 3,11  3,12  3,13  3,14   3,2   3,4 3,5   3,6   3,7   3,8   3,9   4,1  4,10  4,11  4,12
x->y 0.029 0.018 0.072 0.052 0.139 0.029 0.087 0.078 0.795 0.025 9.297 0.984 0.85 1.063 0.499 0.428 0.204 4.829 0.2 1.301 0.258 0.198 1.105 1.432 0.043 0.037 0.038
      4,13  4,14  4,2   4,3   4,5   4,6   4,7   4,8   4,9   5,1  5,10  5,11  5,12  5,13  5,14   5,2  5,3   5,4   5,6   5,7   5,8   5,9   6,1  6,10  6,11 6,12  6,13
x->y 0.182 0.283 0.05 0.134 0.039 0.205 0.977 0.056 0.025 0.284 0.027 0.024 0.014 0.016 0.275 0.026 0.04 0.159 0.139 0.083 0.023 0.022 3.382 0.081 0.091 0.07 0.076
      6,14   6,2   6,3   6,4   6,5   6,7   6,8   6,9   7,1  7,10  7,11  7,12  7,13 7,14   7,2   7,3   7,4   7,5   7,6   7,8  7,9  8,1  8,10  8,11  8,12  8,13 8,14   8,2
x->y 1.716 0.086 2.412 1.239 0.147 0.208 0.065 0.071 1.528 0.034 0.034 0.019 0.033 0.07 0.036 0.051 0.198 0.029 0.094 0.023 0.03 0.26 0.028 0.024 0.029 0.023 0.06 0.886
       8,3   8,4   8,5   8,6   8,7  8,9   9,1  9,10  9,11  9,12  9,13  9,14   9,2   9,3   9,4   9,5   9,6   9,7   9,8
x->y 0.041 0.155 0.026 0.081 0.065 0.03 0.251 0.024 0.012 0.012 0.014 0.079 0.026 0.046 0.151 0.021 0.064 0.067 0.023
"

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
  geom_text(aes(label = ifelse(is.na(fill_color), "", round(Rate, 2))), size = 3) +
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
    axis.text.x = element_text(size = 10, angle = 45, hjust = 0, vjust = 0.02,
                               margin = margin(b = -5)),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(margin = margin(t = -5)),
    axis.title.y = element_text(margin = margin(r = 5)),
    panel.grid = element_blank(),
    plot.title = element_blank()
  )

