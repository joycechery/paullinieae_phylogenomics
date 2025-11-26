###############################################################################
# Phylogenetic Trait Evolution: Transition Rate Visualization
# Author: Joyce Onyenedum with help from ChatGPT
# Description: extracts transition rates for vascular categories, and visualizes them as a heatmap.
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

# Define your intended biological order of states
# Factorize the trait data
state_levels <- c("Typical" ,"Procambial" ,"Cambial" ,"Ectopic-cambia" ,"Procambial+Cambial" ,"Procambial+Cambial+Ectopic-cambia" ,"Procambial+Ectopic-cambia" ,"Cambial+Ectopic-cambia")

#Read in and summarize simmap RDS
simmap_trees <- readRDS("simmap.trees_vc.RDS")
summary_stm <- describe.simmap(simmap_trees)
print(summary_stm)

###############################################################################
# Transition Matrix Parsing and Heatmap Visualization
###############################################################################

# ---- Raw text output copied from summary ----
raw_text <- "changes are of the following types:
      1,2    1,3  1,4   1,5 1,6   1,7   1,8    2,1  2,3   2,4 2,5   2,6   2,7   2,8    3,1   3,2   3,4   3,5 3,6
x->y 2.85 17.255 5.62 0.015   0 0.892 0.058 11.415 6.33 0.007   0 0.745 5.487 0.017 16.906 6.626 1.545 0.012   0
       3,7   3,8   4,1 4,2   4,3 4,5 4,6 4,7 4,8   5,1 5,2 5,3 5,4 5,6   5,7 5,8 6,1   6,2 6,3 6,4 6,5   6,7 6,8
x->y 0.004 0.936 2.692   0 0.289   0   0   0   0 0.001   0   0   0   0 0.211   0   0 0.008   0   0   0 0.551   0
       7,1   7,2 7,3 7,4   7,5   7,6 7,8 8,1 8,2   8,3 8,4 8,5 8,6 8,7
x->y 0.055 0.218   0   0 1.181 1.802   0   0   0 0.014   0   0   0   0"


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
states <- 1:8
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
df$fill_color <- with(df, ifelse(From == To, NA, Rate))

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
