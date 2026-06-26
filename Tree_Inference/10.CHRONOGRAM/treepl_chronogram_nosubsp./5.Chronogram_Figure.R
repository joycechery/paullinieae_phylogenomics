###############################################################################
# Generate publication figure of chronogram (Updated for SumTrees)
# By Joyce G. Onyenedum
###############################################################################

library(ape)
library(treeio)
library(ggtree)
library(ggplot2)
library(deeptime)

## TREE WITH NODE VALUES
setwd("~/Desktop/")
# Load your new SumTrees file
tree <- read.beast("final_summarized_chronogram.nex")

# Extract min and max values from the SumTrees 'age_hpd95' field safely
tree@data$age_hpd95_min <- sapply(tree@data$age_hpd95, function(x) {
  if (is.null(x) || length(x) == 0 || any(is.na(x))) return(NA)
  return(as.numeric(x[1]))
})

tree@data$age_hpd95_max <- sapply(tree@data$age_hpd95, function(x) {
  if (is.null(x) || length(x) == 0 || any(is.na(x))) return(NA)
  return(as.numeric(x[length(x)])) # Grabs the second value if present
})

# Plot the tree
p <- ggtree(tree) +
  # Updated to use the SumTrees HPD variable name
  geom_range("age_hpd95", color = "steelblue", size = 1.7, alpha = 0.8) +
  geom_tiplab(size = 1.5) +
  geom_text2(
    aes(label = ifelse(
      !isTip & !is.na(age_mean) & age_mean > 0,
      paste0(
        sprintf("%.2f", age_mean), " [",
        sprintf("%.2f", age_hpd95_min), "-",
        sprintf("%.2f", age_hpd95_max), "]"
      ),
      NA
    )),
    hjust = -0.15,
    size = .65,
    na.rm = TRUE
  ) +
  # Adjust xlim if your root age (~68.8 Mya) is cut off by the geological timeline
  coord_geo(
    xlim = c(-75, 5), 
    ylim = c(0, 230), # Adjusted ylim to fit screme
    dat = list("epochs"),
    abbrv = FALSE,
    size = list(4, 3),
    neg = TRUE
  ) +
  scale_x_continuous(
    breaks = seq(-70, 0, 10),
    labels = abs(seq(-70, 0, 10))
  ) +
  theme_tree2()

# Reverse the time scale so 0 is present-day on the right
revts(p)
