###############################################################################
# Generate publication figure of chronogram 
# By Joyce G. Onyenedum
###############################################################################

library(ape)
library(treeio)
library(ggtree)
library(ggplot2)
library(deeptime)

## TREE WITH NODE VALUES
setwd("~/Desktop/10.CHRONOGRAM_Final/")
# Load your new SumTrees file
tree <- read.beast("CI-combined_chronograms.nex")

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
  geom_tiplab(size = 1.8) +
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
    size = 1,
    na.rm = TRUE
  ) +
# Adjust xlim 
coord_geo(
  xlim = c(-75, 10), 
  ylim = c(0, 227), # Adjusted ylim to fit screme
  dat = list("epochs"),
  fill = list(c("grey90")),
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

#exported 9 x 13
