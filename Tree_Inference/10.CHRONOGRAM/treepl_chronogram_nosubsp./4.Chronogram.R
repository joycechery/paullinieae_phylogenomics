## TREE WITH NODE VALUES
setwd("/Users/joyce.onyenedum/Library/CloudStorage/Box-Box/Onyenedum_Lab/Joyce_Onyenedum/Projects/Paullinieae_Phylogeny/*MANUSCRIPT/5.CHRONOGRAM.pl/tree.pl/pau333s_351g/treepl_chronogram_nosubsp./chronograms/")
tree <- read.beast("CI-combined_chronograms.tre")

# Extract min and max values from height_0.95_HPD
tree@data$height_0.95_HPD_min <- sapply(tree@data$height_0.95_HPD, function(x) x[1])
tree@data$height_0.95_HPD_max <- sapply(tree@data$height_0.95_HPD, function(x) x[2])

# Plot the tree
p <- ggtree(tree) +
  geom_range("height_0.95_HPD", color = "steelblue", size = 1.7, alpha = 0.8) +
  geom_tiplab(size = 1.5) +
  geom_text2(
    aes(label = ifelse(
      !isTip,
      paste0(
        sprintf("%.2f", height), " [",
        sprintf("%.2f", height_0.95_HPD_min), "-",
        sprintf("%.2f", height_0.95_HPD_max), "]"
      ),
      NA
    )),
    hjust = -0.15,
    size = .65,
    #  nudge_y = 2,
    na.rm = TRUE
  ) +
  coord_geo(
    xlim = c(-90, 12),
    ylim = c(0, 225),
    dat = list("epochs", "periods"),
    abbrv = FALSE,
    size = list(3, 3),
    neg = TRUE
  ) +
  scale_x_continuous(
    breaks = seq(-80, 0, 10),
    labels = abs(seq(-80, 0, 10))
  ) +
  theme_tree2()

revts(p)
Results
