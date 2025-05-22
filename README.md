# Code for "Multi-Scale Characterization of Pleiotropy and Unique Genetic Influences on Cognitive Abilities"
Code accompanying the paper: "Multi-Scale Characterization of Pleiotropy and Unique Genetic Influences on Cognitive Abilities"

We adapted the ggmModSelect algorithm from the qgraph R package to estimate an undirected genomic network model. A working code demo (including the real LDSC covariance matrix used in our paper) in available in the Code Demo folder. To run it, please download the folder, and set it as your working directory in R (instructions provided in README file within it). This is the only novel algorithm used in the paper, the analytic pipeline otherwise entirely makes use of existing package/functions written by other people.

Primary Analyses.R contains the main analytic pipeline used for this paper

Code MiXeR, LAVA, and S-LDSC analyses are available in their respective folders. 
