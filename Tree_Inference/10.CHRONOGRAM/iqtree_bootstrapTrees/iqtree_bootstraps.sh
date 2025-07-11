#!/bin/bash

# Run IQ-TREE with specified parameters
./iqtree2 \
  -s pau_333s_351g.fasta \
  -m GTR \
  -bb 1000 \
  -spp pau_333s_351g_partitions \
  -t pau_333s_351g_partitions.contree \
  -wbtl pau_333s_351g_partitions.BS_BL.tre

