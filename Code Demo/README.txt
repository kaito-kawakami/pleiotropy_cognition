#---------------------------------------------------------------------------------------------------------
# Code demo for: "Multi-Scale Characterization of Pleiotropy and Unique Genetic Influences on Cognitive Abilities"
# This script walksthrough how to estimate an undirected network model from LDSC-derived genetic correlations
# using our implementation of the ggmModSelect algorithm.
#
# Written by: Kaito Kawakami and Jacob Knypsel
#---------------------------------------------------------------------------------------------------------

We performed all analyses using R version 4.3.0. The following R libraries are required to run this code demo,
we provide installation code in the script:

GNA (version 0.0.1)
GenomicSEM (version 0.0.5)
Psychonetrics (version 0.13)
glasso (version 1.11)
Matrix (version 1.7.0)
gdata (version 3.0.1)

#----------------------------------------------------
# Installation Guide
#----------------------------------------------------

Step 1: Download the R programming language (https://www.r-project.org/)

Step 2: Download R studio (https://posit.co/downloads/)

Step 3: Download the Code folder. Open the cognitive_network_demo.R script in R studio, and install dependencies (about 2-3 minutes on a normal desktop computer), and set the folder as the working directory.

Step 4: Run the analysis script to estimate undirected genomic network (about 1-2 minutes on a normal desktop computer)
