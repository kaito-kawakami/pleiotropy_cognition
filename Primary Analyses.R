#---------------------------------------------------------------------------------------------------------
# Code for: "Genomic network modelling reveals unique genetic influences on specific cognitive abilities"
# Written by: Kaito Kawakami, Jacob Knypsel, and Johan Källberg Zvrskovec
#---------------------------------------------------------------------------------------------------------
setwd("/cephfs/volumes/hpc_data_usr/k20051211/8e486e0b-f00f-4534-8dcc-d66a07c10fd5/scratch_tmp/network_final/")

### Step 1: reverse code beta column for GWAS sumstats 
reverse_code() <- function(files) {
  # Prepare a list to hold the data.frames
  out_list <- vector("list", length(files))
  names(out_list) <- files
  
  for (i in seq_along(files)) {
    f <- files[i]
    
    # 1) Read in
    con_in <- gzfile(f, "rt")
    df     <- read.table(con_in, header = TRUE, sep = "\t")
    close(con_in)
    
    # 2) Flip the sign of BETA
    df$BETA <- -df$BETA
    
    # 3) Write it back out (gzipped, overwriting the original)
    con_out <- gzfile(f, "wt")
    write.table(df,
                con_out,
                sep = "\t",
                row.names = FALSE,
                quote     = FALSE)
    close(con_out)
    
    # Store in list
    out_list[[i]] <- df
  }
  
  return(out_list)
}

files <- c("20157.tmtb.1.fastGWA.gz", "20023.rt.1.fastGWA.gz", "399.memory.1.fastGWA.gz")
corrected_dfs <- reverse_code(files)

### Load in libraries 

# load igraph
library(igraph, lib.loc="/cephfs/volumes/hpc_data_usr/k20051211/8e486e0b-f00f-4534-8dcc-d66a07c10fd5/scratch_tmp/network_final/")

# load jpeg (dependency)
dyn.load("/users/k20051211/lib/libjpeg.so.9.5.0")
dyn.load("/users/k20051211/R/x86_64-pc-linux-gnu-library/4.3/jpeg/libs/jpeg.so")
library(jpeg)

library(Matrix); library(dplyr); library(glasso)
library(psychonetrics); library(writexl); library(gdata)

library(GNA)
library(GenomicSEM)

##########################################################################################
# Step 2: Preprocessing
# Munge sumstats
# Create genetic covariance matrix using ldsc() from GenomicSEM
##########################################################################################

## Munging

#create vector of the summary statistics files
files <- c("399.pairs.1.fastGWA.gz",
           "20023.rt.1.fastGWA.gz",
           "20018.prospective.1.fastGWA.gz",
           "20159.symbol.1.fastGWA.gz",
           "20240.numeric.1.fastGWA.gz",
           "20016.vnr.1.fastGWA.gz",
           "20157.tmtb.1.fastGWA.gz")

#define the reference file being used to allign alleles across summary stats
#here we are using hapmap3
hm3<-"eur_w_ld_chr/w_hm3.snplist"

#name the traits
trait.names<-c("memory", "RT", "prospective", "symbol", "numeric",
               "VNR", "TMTB")

#list the sample sizes. Eising has SNP-specific N, while De La Fuente has one N for all SNPs.
N=c(NA,NA,NA,NA,NA,NA,NA)

#definte the imputation quality filter
info.filter=0.9

#define the MAF filter
maf.filter=0.01

#run munge
munge(files=files,hm3=hm3,trait.names=trait.names,N=N,info.filter=info.filter,maf.filter=maf.filter,
      parallel = TRUE)

## LDSC

#vector of munged summary statisitcs
traits<-c("memory.sumstats.gz", "RT.sumstats.gz", "prospective.sumstats.gz",
          "symbol.sumstats.gz", "numeric.sumstats.gz",
          "TMTB.sumstats.gz", "VNR.sumstats.gz")

# Set sample and population prevalence to NA for all traits (because continuous)
trait_length<-length(traits)
sample.prev<-rep(NA, trait_length)
population.prev<-rep(NA, trait_length)

#the folder of LD scores
ld<-"eur_w_ld_chr/"

#the folder of LD weights [typically the same as folder of LD scores]
wld<-"eur_w_ld_chr/"

#name the traits
trait.names<- c("memory", "RT", "prospective", "symbol", "numeric", "TMTB", "VNR")

#run LDSC
LDSCoutput_final<-ldsc(traits=traits,sample.prev=sample.prev,population.prev=population.prev,
                       ld=ld,wld=wld,trait.names=trait.names, stand = TRUE)

# saveRDS(LDSCoutput_final, file="LDSCoutput_final.RDS")
LDSCoutput_final <- readRDS("LDSCoutput_final.RDS")

corrs <- LDSCoutput_final$S_Stand
rownames(corrs) <- colnames(corrs)

##########################################################################################
# Step 2: Genetic correlation analyses
##########################################################################################

# Correlation heatmap
library(ggcorrplot)
library(ggplot2)
library(grid)

p <- ggcorrplot(corrs,
                type = "lower",
                insig = "blank",
                lab = TRUE,
                digits = 3
) +
  guides(fill = guide_colorbar(barheight = unit(10, "cm"),   # adjust this value as needed
                               barwidth = unit(1, "cm"))) +
  theme(legend.position = "right")


#-------------------------------
# average rg correlations
#-------------------------------

# 1. Overall average correlation (exclude diagonal)
# Use only the upper triangle (or lower triangle) since the matrix is symmetric.
upper_corr <- corrs[upper.tri(corrs)]
overall_avg <- mean(upper_corr)
print(paste("Overall average correlation:", round(overall_avg, 3)))
# 0.479

# 2. Average correlation for Reaction Time (RT)
# Exclude the self-correlation
RT_corr <- corrs["RT", ]
RT_corr_no_self <- RT_corr[names(RT_corr) != "RT"]
avg_RT <- mean(RT_corr_no_self)
print(paste("Average correlation for RT:", round(avg_RT, 3)))
# 0.220

# 3. Average correlation for TMTB
TMTB_corr <- corrs["TMTB", ]
TMTB_corr_no_self <- TMTB_corr[names(TMTB_corr) != "TMTB"]
avg_TMTB <- mean(TMTB_corr_no_self)
print(paste("Average correlation for TMTB:", round(avg_TMTB, 3)))
# 0.639

# check Z-score h2
S <- LDSCoutput_final$S
k<-nrow(LDSCoutput_final$S)
SE<-matrix(0, k, k)
SE[lower.tri(SE,diag=TRUE)] <-sqrt(diag(LDSCoutput_final$V))
S / SE

Stand <- LDSCoutput_final$S_Stand
rownames(Stand) <- colnames(Stand)

#-------------------------------
# common factor model
#-------------------------------
common <- commonfactor(LDSCoutput_final)
fit <- as.data.frame(common$modelfit)
results <- as.data.frame(common$results)

##########################################################################################
# LAVA
##########################################################################################
setwd("/cephfs/volumes/hpc_data_usr/k20051211/8e486e0b-f00f-4534-8dcc-d66a07c10fd5/scratch_tmp/network_final/lava/results")

#----------------------------
# 1. Chord plot
#----------------------------

# Load required packages
library(dplyr)
library(circlize)

# List all files with .bivar extension in the working directory
files <- list.files(pattern = "\\.bivar$")

# Read each file, parse trait names from filename, and combine
data_list <- lapply(files, function(f) {
  df <- read.table(f, header = TRUE)  
  
  # Extract trait1, trait2 from filename "phenotype1.phenotype2.bivar"
  parts <- strsplit(f, "\\.")[[1]]
  df$trait1 <- parts[1]
  df$trait2 <- parts[2]
  
  return(df)
})

# Combine into one data frame
data_all <- do.call(rbind, data_list)
write.csv(data_all, "lava_all.csv", row.names=FALSE)

# Number of total tests
total_tests <- nrow(data_all)
bonf_threshold <- 0.05 / total_tests
cat("Bonferroni threshold:", bonf_threshold, "\n")

# Filter for significant results
sig_data <- data_all %>%
  filter(p < bonf_threshold)

write.csv(sig_data, "sig_data.csv")

# Count how many significant tests per trait pair
sig_counts <- sig_data %>%
  count(trait1, trait2)

# Get the unique traits from both trait1 and trait2
all_traits <- unique(c(sig_counts$trait1, sig_counts$trait2))

# Initialize a square matrix of zeros
matrix_data <- matrix(
  0, 
  nrow = length(all_traits), 
  ncol = length(all_traits),
  dimnames = list(all_traits, all_traits)
)

# Fill the matrix with counts of significant tests
for(i in seq_len(nrow(sig_counts))) {
  t1 <- sig_counts$trait1[i]
  t2 <- sig_counts$trait2[i]
  matrix_data[t1, t2] <- sig_counts$n[i]
  matrix_data[t2, t1] <- sig_counts$n[i]  
}

# 5A) Specify the order of traits around the circle, if desired.
trait_order <- sort(all_traits)

# 5B) Define colors for each trait.
n_traits <- length(all_traits)
my_colors <- colorRampPalette(c("#03045e", "#0077b6", "#00b4d8", "#ffba08", "#e85d04"))(n_traits)

# Create a named vector: names = trait, value = color
trait_colors <- setNames(my_colors, trait_order)

# Clear any existing circos plots and reset device margins
circos.clear()
par(mai = rep(0, 4))  # set all margins to zero

# Adjust circos parameters with a smaller track height
circos.par(
  start.degree = 90,
  gap.after = rep(2, length(trait_order)),
  cell.padding = c(0, 0, 0, 0),
  track.margin = c(0, 0),
  track.height = 0.05,             
  canvas.xlim = c(-1.2, 1.2),     
  canvas.ylim = c(-1.2, 1.2)
)

# Draw the chord diagram
chordDiagram(
  x = matrix_data,
  order = trait_order,
  grid.col = trait_colors,
  transparency = 0.3,
  directional = 0,
  annotationTrack = "grid"
)

# Add bigger labels to each sector
circos.trackPlotRegion(
  track.index = get.current.track.index(),
  panel.fun = function(x, y) {
    sector.name <- get.cell.meta.data("sector.index")
    circos.text(
      x = CELL_META$xcenter,
      y = CELL_META$ylim[1] + 1.1,
      labels = sector.name,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 1.1   
    )
  },
  bg.border = NA
)

# Reset circos so it doesn't affect future plots
circos.clear()

### volcano plot for numeric and RT
### https://bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html

# import data
library(data.table)
setwd("/scratch_tmp/users/k20051211/network_final/lava/results")
lava_numeric_rt <- fread("RT.numeric.bivar") 

# plot
library(ggplot2)
library(ggrepel)

# Prepare the data: mark significant points
lava_numeric_rt$signif <- ifelse(lava_numeric_rt$p < 7.05e-06, "Significant", "Not Significant")

# Create a surreal volcano plot inspired by Salvador Dalí
ggplot(lava_numeric_rt, aes(x = rho, y = -log10(p), color = signif)) +
  geom_point(alpha = 0.8, size = 3) +
  geom_text_repel(data = subset(lava_numeric_rt, p < 7.05e-06),
                  aes(label = locus),
                  size = 4, fontface = "bold", color = "black",
                  box.padding = 0.35, point.padding = 0.5, segment.color = "grey50") +
  scale_color_manual(values = c("Significant" = "#FF69B4", "Not Significant" = "#1E90FF")) +
  labs(
    title = "Local Genetic Correlations (RT-Numeric)",
    subtitle = "",
    x = "Genetic Correlation (ρ)",
    y = "-log10(p-value)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5, color = "black"),
    plot.caption = element_text(face = "italic", hjust = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "none"
  )


#----------------------------
# 2. Network plot from LAVA
#----------------------------
# Load required libraries
library(dplyr)
library(qgraph)

# 1. Identify the locus with the most connections
locus_counts <- sig_data %>% 
  count(locus) %>% 
  arrange(desc(n))
max_locus <- locus_counts$locus[1]
cat("Locus with most connections:", max_locus, "\n")

write.csv(locus_counts, "locus_counts.csv")

# 2. Subset the data for that locus
subset_data <- sig_data %>% 
  filter(locus == max_locus)

# 3. Create a list of unique original node labels 
original_nodes <- unique(c(as.character(subset_data$phen1), as.character(subset_data$phen2)))

# 4. Create numeric labels for nodes (as characters)
numeric_nodes <- as.character(seq_along(original_nodes))

# 5. Create a mapping for the legend 
node_mapping <- paste(numeric_nodes, original_nodes, sep=": ")

# 6. Initialize an empty symmetric matrix for edge weights (omega)
omega <- matrix(0, nrow = length(original_nodes), ncol = length(original_nodes),
                dimnames = list(original_nodes, original_nodes))

# 7. Populate the matrix with rho values from the subset data
for(i in 1:nrow(subset_data)) {
  p1 <- as.character(subset_data$phen1[i])
  p2 <- as.character(subset_data$phen2[i])
  weight <- subset_data$rho[i]
  omega[p1, p2] <- weight
  omega[p2, p1] <- weight  # Ensure the matrix is symmetric
}

# 8. Convert the omega matrix to use numeric node names for plotting.
omega_numeric <- omega[original_nodes, original_nodes]
rownames(omega_numeric) <- numeric_nodes
colnames(omega_numeric) <- numeric_nodes

labs <- c('memory', 'numeric', 'symbol', 'TMTB', 'VNR')

# 9. Create the network plot using qgraph with numeric node names
qgraph(omega_numeric,
       layout = "spring",
       nodeNames = labs,
       legend.cex = 0.6,
       layoutScale = c(0.9, 0.9),
       fade = FALSE,
       esize = 15,
       color = c("#440154FF"),
       border.color = "white",
       border.width = 2,
       label.color = "white",
       vsize = 7,
       curve = 0.1,
       curveAll = TRUE,
       title.cex = 1.3,
       title = paste("Network Plot for Locus", max_locus))

#----------------------------
# 6. Bar Chart: Significant Connections per Trait
#----------------------------
library(ggplot2)
library(dplyr)
library(tidyr)

# Compute connection counts per trait by pivoting trait1 and trait2 into one column
trait_connection_counts <- sig_data %>%
  pivot_longer(cols = c(trait1, trait2), names_to = "role", values_to = "trait") %>%
  count(trait, name = "connections") %>%
  arrange(desc(connections))

# Create the bar chart using ggplot2
ggplot(trait_connection_counts, aes(x = reorder(trait, -connections), y = connections)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = connections), vjust = -0.5, size = 4) +
  labs(title = "Number of Significant Local Genetic Correlations per Trait",
       x = "Trait",
       y = "Number of Connections") +
  theme_minimal(base_size = 20)

##########################################################################################
# MiXeR
##########################################################################################
library(readxl)

### univariate mixer
mixer_univariate <- read_excel('univariate_mixer.xlsx')

# Load required libraries
library(tidyr)
library(dplyr)
library(ggplot2)

# Pivot the data from wide to long format 
mixer_long <- mixer_univariate %>%
  pivot_longer(
    cols = c(`nc@p9 (mean)`, `sig2_beta (mean)`, `h2 (mean)`,
             `nc@p9 (std)`, `sig2_beta (std)`, `h2 (std)`),
    names_to = c("measure", "stat"),
    names_pattern = "(.*) \\((.*)\\)"
  ) %>%
  pivot_wider(
    names_from = stat,
    values_from = value
  )

mixer_long <- mixer_long %>%
  mutate(measure = recode(measure,
                          "nc@p9" = "Polygenicity",
                          "sig2_beta" = "Discoverability",
                          "h2" = "Heritability"))

# Set the factor level order to control the panel order: nc@p9, then Discoverability, then Heritability
mixer_long$measure <- factor(mixer_long$measure, levels = c("Polygenicity", "Discoverability", "Heritability"))

# plot
ggplot(mixer_long, aes(x = fname, y = mean, fill = fname)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = mean - std, ymax = mean + std), width = 0.2) +
  facet_wrap(~ measure, scales = "free_y", nrow = 1) +
  labs(
    title = "Bar Charts for Polygenicity, Discoverability, and Heritability",
    x = "Trait",
    y = "Value"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )


#import Excel file into R
mixer <- read_excel('bivariate_mixer.xlsx')
mixer <- mixer %>% filter(best_vs_min_AIC > 0 & best_vs_max_AIC > 0)

# LDSC correlations
ldsc_cor <- LDSCoutput_final$S_Stand
rownames(ldsc_cor) <- colnames(LDSCoutput_final$S_Stand)

# LDSC SE
k<-nrow(LDSCoutput_final$S_Stand)
ldsc_SE<-matrix(0, k, k)
ldsc_SE[lower.tri(ldsc_SE,diag=TRUE)] <-sqrt(diag(LDSCoutput_final$V_Stand))
rownames(ldsc_SE) <- colnames(ldsc_SE) <- colnames(LDSCoutput_final$S_Stand)
ldsc_SE <- (ldsc_SE + t(ldsc_SE)) / 2

#### scatter
# Load required libraries
library(tidyr)
library(dplyr)
library(ggplot2)

# 1. Convert the LDSC matrices wide to long
ldsc_cor_df <- as.data.frame(as.table(ldsc_cor))
colnames(ldsc_cor_df) <- c("trait1", "trait2", "ldsc_cor")

ldsc_se_df <- as.data.frame(as.table(ldsc_SE))
colnames(ldsc_se_df) <- c("trait1", "trait2", "ldsc_se")

# 2. Merge the LDSC correlation and standard error data frames 
ldsc_df <- inner_join(ldsc_cor_df, ldsc_se_df, by = c("trait1", "trait2"))

ldsc_df <- ldsc_df %>%
  mutate(trait1 = tolower(trait1),
         trait2 = tolower(trait2))

mixer_dice <- mixer %>%
  select(trait1, trait2, `dice (mean)`, `dice (std)`, `fraction_concordant_within_shared (mean)`) %>%
  mutate(trait1 = tolower(trait1),
         trait2 = tolower(trait2))

# 4. Merge the dice values from the mixer data with the LDSC data
merged_data <- inner_join(mixer_dice, ldsc_df, by = c("trait1", "trait2"))

# Check the merged data
print(merged_data)

# 5. Create the scatter plot with error bars using ggplot2
ggplot(merged_data, aes(x = `dice (mean)`, y = ldsc_cor)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ldsc_cor - ldsc_se, ymax = ldsc_cor + ldsc_se), width = 0.01) +
  geom_errorbarh(aes(xmin = `dice (mean)` - `dice (std)`, xmax = `dice (mean)` + `dice (std)`), height = 0.01) +
  theme_minimal(base_size = 15) +
  labs(x = "Dice Coefficient",
       y = "LDSC Correlation",
       title = "MiXeR Dice Coefficient vs LDSC Correlation")+
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dotted") +
  xlim(c(0,1))+
  ylim(c(0,1))


mean(merged_data$ldsc_cor)
mean(merged_data$`dice (mean)`)
mean(merged_data$`fraction_concordant_within_shared (mean)`)


##########################################################################################
# Network analysis
##########################################################################################

# GNA, signifiance based model
test <- traitNET(covstruc = LDSCoutput_final, graph_layout="circle")
fit1 <- test$model_results$sparse$modelfit # poor model fit
omega_test <- test$model_results$sparse$parameters

# Import GNA functions to do EBICglasso
.sandwichSE<-function(Model_Results,V,toler){
  
  ### get necessary model matrices / information
  omega <- getmatrix(Model_Results, "omega") #weight matrix
  delta <- getmatrix(Model_Results, "delta") #scaling matrix
  sigma <- getmatrix(Model_Results, "sigma") #mod-implied cov matrix
  
  # extra matrices
  mat <- lapply(Model_Results@extramatrices, as.matrix)
  L <- mat$L
  D <- mat$D #duplication matrix
  Dstar <- mat$Dstar
  In <- mat$In #k x k identity matrix
  A <- mat$A
  
  ### wls weight matrix from stage 2
  estimation <- Model_Results@estimator
  
  if(estimation == "ML"){
    sigma_inv <- solve(sigma)
    S2.W <- 0.5 * t(D) %*% kronecker(sigma_inv, sigma_inv) %*% D
  }
  
  if(estimation == "DWLS"){
    S2.W <- Model_Results@sample@WLS.W[[1]]
  }
  
  ### Jacobian matrix of sigma in respect to theta (all model parameters) i.e "Delta" in genomic sem
  IminOinv <- solve(diag(ncol(omega)) - omega) #inverse of (identify matrix - omega)
  delta_IminOinv <- delta %*% IminOinv
  
  # model derivatives sigma in respect to omega
  d_sigma_omega <- L %*% (delta_IminOinv %x% delta_IminOinv) %*% Dstar
  
  # model derivatives sigma in respect to delta
  d_sigma_delta <- L %*% ((delta_IminOinv %x% In) + (In %x% delta_IminOinv)) %*% A
  
  # full jacobian
  S2.delt <- cbind(d_sigma_omega, d_sigma_delta)
  
  ### sandwich correction
  #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specified
  bread <- solve(t(S2.delt) %*% S2.W %*% S2.delt, tol = toler)
  
  #create the "lettuce" part of the sandwich
  lettuce <- S2.W %*% S2.delt
  
  #ohm-hat-theta-tilde is the corrected sampling covariance matrix of the model parameters
  Ohtt <- bread %*% t(lettuce) %*% V %*% lettuce %*% bread
  
  #the lettuce plus inner "meat" (V) of the sandwich adjusts the naive covariance matrix by using the correct sampling covariance matrix of the observed covariance matrix in the computation
  SE <- as.vector(sqrt(diag(Ohtt)))
  
  return(SE)
}
.runGGM <- function(covstruc,fix_omega="full",saturated=NULL,estimation="ML",toler=NULL) {
  
  # read in V (sampling covariance) and S (covariance) matrices
  V_LD <- as.matrix(covstruc[[1]]) 
  S_LD <- as.matrix(covstruc[[2]])
  
  # If fixed omega matrix, recode - edges fixed to zero = 0; edge free = 1 (other integers encode equality constraints)
  if (is.matrix(fix_omega)) {
    fix_omega[fix_omega != 0] <- 1
    diag(fix_omega) <- 0
  }
  
  ### run the GGM
  if(estimation == "ML"){
    model <- varcov(type = "ggm", covs = S_LD, omega = fix_omega, nobs = 200, covtype = "ML", estimator = "ML", optimizer = "nlminb")
  }
  if(estimation == "DWLS"){
    W <- diag(diag(solve(V_LD)),nrow(V_LD)) #diagonal weight matrix
    model <- varcov(type = "ggm", covs = S_LD, omega = fix_omega, nobs = 200, covtype = "ML", estimator = "DWLS", WLS.W = W, optimizer = "nlminb")
  }
  
  Model_Results <- runmodel(model)
  
  ### get necessary model matrices / information
  omega <- getmatrix(Model_Results, "omega") #weight matrix
  delta <- getmatrix(Model_Results, "delta") #scaling matrix
  sigma <- getmatrix(Model_Results, "sigma") #mod-implied cov matrix
  
  #get sandwich corrected SEs from internal function
  SE<-.sandwichSE(Model_Results,V_LD,toler)
  
  ### extract results
  params <- data.frame(Model_Results@parameters)
  params$se <- SE
  params$se[params$par == 0] <- NA
  params$z <- params$est / params$se
  params$p <- round(2 * pnorm(abs(params$est / params$se), lower.tail = FALSE),5)
  params <- params[,c("var1","op","var2","est","se","z","p","par","matrix","row","col")]
  colnames(params) <- c("trait1","op","trait2","est","se","z","p","free","matrix","row","col")
  
  #calculate model fit if there are pruned edges (otherwise fully saturated and not relevant)
  if(is.matrix(fix_omega)){
    
    #calculate model chi-square:
    
    #eigen values of the V matrix  
    Eig<-as.matrix(eigen(V_LD)$values)
    Eig2<-diag(ncol(V_LD))
    diag(Eig2)<-Eig
    
    #Pull P1 (the eigen vectors of V_eta)
    P1<-eigen(V_LD)$vectors
    
    #residual matrix: difference between model implied matrix and observed matrix
    resid<-S_LD-sigma
    
    #eta: the vector of unique elements of the residual matrix
    eta<-as.vector(lowerTriangle(resid,diag=TRUE))
    
    #matrix algebra weighting the vector of residuals by the precision of those residuals (i.e., P1 and Eig)
    model_chi<-t(eta)%*%P1%*%solve(Eig2)%*%t(P1)%*%eta
    
    #degrees of freedom of model (how many unique edges are fixed to 0)
    df<-Model_Results@fitmeasures$df
    
    #calculate p-value for model chi-square
    model_chi_p<-pchisq(model_chi,df,lower.tail=FALSE)
    
    #calculate SRMR
    
    #unique values of genetic correlation  matrix
    obs <-  cov2cor(S_LD)[!lower.tri(cov2cor(S_LD))]
    
    #unique values of model implied correlation matrix
    imp <-  cov2cor(sigma)[!lower.tri(cov2cor(sigma))]
    
    #square root of the average squared residual 
    SRMR<-sqrt(mean((imp - obs)^2))
    
    #calculate CFI:
    
    #the difference between the model implied matrix of the independence model [only models variances]
    #and the observed matrix is the observed matrix with diagonals set to 0
    resid_CFI<-S_LD
    diag(resid_CFI)<-0
    
    #eta: the vector of unique elements of the residual matrix
    eta_CFI<-as.vector(lowerTriangle(resid_CFI,diag=TRUE))
    
    #matrix algebra weighting the vector of residuals by the precision of those residuals (i.e., P1 and Eig)
    CFI_chi<-t(eta_CFI)%*%P1%*%solve(Eig2)%*%t(P1)%*%eta_CFI
    
    ##df of independence Model
    k<-ncol(S_LD)
    dfCFI <- (((k * (k + 1))/2) - k)
    
    #calculate CFI
    CFI<-as.numeric(((CFI_chi-dfCFI)-(model_chi-df))/(CFI_chi-dfCFI))
    
    #calculate AIC
    AIC<-(model_chi + 2*Model_Results@fitmeasures$df)
    
    ###calculate eBIC - TESTING, NEED TO CHECK
    npar = Model_Results@fitmeasures$npar
    nvar = Model_Results@fitmeasures$nvar
    # effective N of weight matrix. based on formula for estiming sampling variance of partial correlation (see van Aert & Goos, 2023).. mean of edge-wise Ns 
    neff <- mean(((((1-(saturated$parameters$est[saturated$parameters$matrix=="omega"]^2))^2)/(saturated$parameters$se[saturated$parameters$matrix=="omega"]^2)) + nvar + 1), na.rm = TRUE)
    
    # calc eBICs (based on Foygel & Drton, 2010)
    BIC <- model_chi + (npar*log(neff))
    eBIC.25 <- model_chi + (npar*log(neff)) + (4*npar*0.25*log(nvar))
    eBIC.50 <- model_chi + (npar*log(neff)) + (4*npar*0.50*log(nvar))
    eBIC.75 <- model_chi + (npar*log(neff)) + (4*npar*0.75*log(nvar))
    eBIC1 <- model_chi + (npar*log(neff)) + (4*npar*1*log(nvar))
    
    #combine model fit indices
    modelfit<-cbind(model_chi,df,model_chi_p, SRMR,CFI,AIC,BIC,eBIC.25,eBIC.50,eBIC.75,eBIC1)
    colnames(modelfit)=c("model_chisquare","df","modelchi_pvalue", "SRMR", "CFI","AIC","BIC","eBIC.25","eBIC.50","eBIC.75","eBIC1")
    
  }else{
    modelfit <- NULL
  }
  
  traitnames <- list(colnames(S_LD),colnames(S_LD))
  dimnames(omega) <- traitnames
  dimnames(delta) <- traitnames
  dimnames(sigma) <- traitnames
  
  return(list(parameters=params,modelfit=modelfit,omega=omega,delta=delta,sigma=sigma))
}
wi2net <- function(x) {
  x <- -cov2cor(x)
  diag(x) <- 0
  x <- forceSymmetric(x)
  return(x)
}
.ebic_unregularized <- function(LDSCoutput, gamma = 0.5, lambda.min.ratio = 0.01, nlambda = 100, estimation = "ML", toler = NULL, penalize.diagonal = FALSE) {
  
  # Estimate saturated model
  model_out <- .runGGM(LDSCoutput, fix_omega = "full", saturated = NULL, estimation = estimation, toler = toler)
  model_results <- list(saturated = model_out)
  
  # Load LDSC matrices
  S <- LDSCoutput$S  # Observed covariance matrix
  V <- LDSCoutput$V  # Sampling covariance matrix
  
  # Compute lambda sequence (code taken from huge package):
  lambda.max = max(max(S - diag(nrow(S))), -min(S - diag(nrow(S))))
  lambda.min = lambda.min.ratio*lambda.max
  lambda = exp(seq(log(lambda.min), log(lambda.max), length = nlambda))
  
  # Run glasso path:
  glas_path <- glassopath(S, lambda, trace = 0,
                          penalize.diagonal=penalize.diagonal)
  
  # Evaluate EBIC for each rho
  ebic_list <- numeric(length(lambda))
  
  for (i in seq_along(lambda)) {
    K <- glas_path$wi[,,i]
    omega <- as.matrix(wi2net(K))
    diag(omega) <- 0
    fix_omega <- omega != 0
    
    # Fit unregularized network model using the saturated model as a starting point
    model_test <- .runGGM(LDSCoutput, fix_omega = fix_omega, saturated = model_results$saturated, estimation = estimation, toler = toler)
    
    # Store the model and its EBIC
    ebic_list[i] <- model_test$modelfit[9] }
  
  # Identify the best model
  best_index <- which.min(ebic_list)
  best_model <- list(
    lambda = lambda[best_index],
    lambda_all = lambda,
    precision_matrix = glas_path$wi[,,best_index],
    partial_correlation_matrix = wi2net(glas_path$wi[,,best_index]),
    ebic_all = ebic_list,
    ebic = ebic_list[best_index]         # Store all EBIC values for review
  )
  
  # Return the best model with all tracked models and EBICs
  return(best_model)
}
mod_select <- function(LDSCoutput, max_iterations = 100, estimation = "ML", toler = NULL, verbose = TRUE) {
  
  # Estimate saturated model
  model_out <- .runGGM(LDSCoutput, fix_omega = "full", saturated = NULL, estimation = estimation, toler = toler)
  model_saturated <- list(saturated = model_out)
  
  message("Running unregularized EBICglasso.")
  
  # Run unregularized ebicglasso
  best_model <- .ebic_unregularized(LDSCoutput)
  
  # Initialize omega and fix_omega
  omega <- as.matrix(best_model$partial_correlation_matrix)
  omega_free <- omega != 0
  diag(omega_free) <- FALSE
  
  # Initialize variables for stepwise search
  best_ebic <- best_model$ebic
  current_model <- best_model$full_model
  iteration <- 0
  improved <- TRUE
  
  while (improved && iteration < max_iterations) {
    iteration <- iteration + 1
    improved <- FALSE
    
    # List all edges to test (upper triangle)
    all_edges <- which(upper.tri(omega_free), arr.ind = TRUE)
    if (verbose) message(sprintf("Iteration %d: Testing %d edges.", iteration, nrow(all_edges)))
    
    # Serial processing of edge testing
    results <- lapply(seq_len(nrow(all_edges)), function(edge) {
      edge_indices <- all_edges[edge, ]
      modified_omega_free <- omega_free
      modified_omega_free[edge_indices[1], edge_indices[2]] <- !modified_omega_free[edge_indices[1], edge_indices[2]]
      modified_omega_free[edge_indices[2], edge_indices[1]] <- modified_omega_free[edge_indices[1], edge_indices[2]]
      
      # Compute the proposed model with modified edges
      model_test <- tryCatch(
        .runGGM(LDSCoutput, fix_omega = modified_omega_free, saturated = model_saturated$saturated, 
                estimation = estimation, toler = toler),
        error = function(e) NULL
      )
      
      # Safely extract EBIC
      if (!is.null(model_test) && !is.null(model_test$modelfit)) {
        proposed_ebic <- model_test$modelfit[9]  # Adjust index for EBIC column
        list(edge1 = edge_indices[1], edge2 = edge_indices[2], ebic = proposed_ebic)
      } else {
        list(edge1 = edge_indices[1], edge2 = edge_indices[2], ebic = Inf)
      }
    })
    
    # Flatten the list of results and convert to a data frame
    results_df <- do.call(rbind, lapply(results, as.data.frame))
    colnames(results_df) <- c("edge1", "edge2", "ebic")
    best_edge_idx <- which.min(results_df$ebic)
    best_edge <- as.numeric(results_df[best_edge_idx, c("edge1", "edge2")])
    best_edge_ebic <- results_df$ebic[best_edge_idx]
    
    # Update if EBIC improved
    if (best_edge_ebic < best_ebic) {
      best_ebic <- best_edge_ebic
      omega_free[best_edge[1], best_edge[2]] <- !omega_free[best_edge[1], best_edge[2]]
      omega_free[best_edge[2], best_edge[1]] <- omega_free[best_edge[1], best_edge[2]]
      current_model <- .runGGM(LDSCoutput, fix_omega = omega_free, saturated = model_saturated$saturated, 
                               estimation = estimation, toler = toler)
      improved <- TRUE
      
      if (verbose) {
        action <- if (omega_free[best_edge[1], best_edge[2]]) "Added" else "Removed"
        message(sprintf("Iteration %d: %s edge %d-%d. EBIC: %.4f", 
                        iteration, action, best_edge[1], best_edge[2], best_ebic))
      }
    } else {
      if (verbose) message("No further improvements found. Stopping.")
    }
  }
  
  return(list(model = current_model, ebic = best_ebic))
}

# estimate network
network_model <- mod_select(LDSCoutput_final, max_iterations = 100, verbose = TRUE)

params <- network_model$model$parameters
fit <- network_model$model$modelfit

omega <- network_model$model$omega 
omega <- ifelse(omega > 0, TRUE, FALSE)

test <- traitNET(covstruc = LDSCoutput_final, graph_layout="spring", fix_omega = omega, reestimate = FALSE, recursive = FALSE, prune = FALSE)
test$model_results
test$model_results$base$modelfit


omega <- as.matrix(network_model$model$omega)
# saveRDS(omega, file="omega.RDS")
omega <- readRDS("omega.RDS")

colnames(omega) <- c(1:7)

# plot
library(qgraph)

# Trait names with numbers
nodeNames <- c("memory", "RT", "prospective", "symbol", "numeric", "TMTB", "VNR")

stepwise_graph <- qgraph(omega, layout = "spring", nodeNames = nodeNames, legend.cex = 0.7,
                         layoutScale = c(0.9,0.9), fade = FALSE, esize = 18,
                         color = c("#440154FF"), border.color = "white", border.width = 2, label.color = "white", 
                         vsize = 9, curve = 0.1, curveAll = TRUE)

layout_old <- stepwise_graph$layout

# Network centrality statisics
library(ggplot2)
omega <- readRDS("omega.RDS")
colnames(omega) <- rownames(omega)
centrality <- centralityTable(omega, standardized = FALSE)

library(igraph) # eigenvector centrality
g <- graph_from_adjacency_matrix(omega, mode = "undirected", weighted = TRUE, diag = FALSE)
ev_cent <- eigen_centrality(g, weights = E(g)$weight)$vector
print(ev_cent)

# Create a data frame for eigenvector centrality
ev_df <- data.frame(
  graph = "graph 1",
  type = NA,
  node = names(ev_cent),
  measure = "Eigenvector",
  value = as.numeric(ev_cent),
  stringsAsFactors = FALSE
)

all_centrality <- bind_rows(centrality, ev_df)

# Faceted bar chart (one panel per centrality measure)
library(grid)  

centrality_plot <- ggplot(all_centrality, aes(x = node, y = value, fill = node)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ measure, scales = "free_y") +
  labs(title = "Centrality Measures by Node", y = "Centrality Value") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

grid.newpage()
grid.text("b", x = unit(0, "npc"), y = unit(1, "npc"), 
          just = c("left", "top"), gp = gpar(fontsize = 16))
print(centrality_plot, newpage = FALSE)

#------------------------------------------------
# Network with AM and ASC correction
#-----------------------------------------------

## Munge
munge(files=c("GCST90267294_buildGRCh37.tsv","parter_choice_index.txt.gz"),
      hm3=hm3,
      trait.names=c("ascertainment", "AM"),
      N=c(765283,NA),
      info.filter=0.9,maf.filter=0.01, parallel = TRUE)

## LDSC
traits<-c("memory.sumstats.gz", "RT.sumstats.gz", "prospective.sumstats.gz",
          "symbol.sumstats.gz", "numeric.sumstats.gz",
          "TMTB.sumstats.gz", "VNR.sumstats.gz",
          "AM.sumstats.gz", "ascertainment.sumstats.gz")

trait.names<-c("memory", "RT", "prospective", "symbol", "numeric",
               "VNR", "TMTB", "AM", "ASC")

trait_length<-length(traits)
sample.prev<-rep(NA, trait_length)
population.prev<-rep(NA, trait_length)

LDSCoutput_corrected<-ldsc(traits=traits,sample.prev=sample.prev,population.prev=population.prev,
                    ld=ld,wld=wld,trait.names=trait.names, stand = TRUE)

saveRDS(LDSCoutput_corrected, file="LDSCoutput_corrected.RDS")
LDSCoutput_corrected <- readRDS("LDSCoutput_corrected.RDS")
rownames(LDSCoutput_corrected$S_Stand) <- colnames(LDSCoutput_corrected$S_Stand)

# estimate network
ea_model <- mod_select(LDSCoutput_corrected, max_iterations = 100, verbose = TRUE)
params <- ea_model$model$parameters
fit <- ea_model$model$modelfit
omega <- as.matrix(ea_model$model$omega)

params_df <- as.data.frame(params)
fit_df <- as.data.frame(fit)

# Create an empty layout matrix for 10 nodes
layout_new <- matrix(NA, nrow = 9, ncol = 2)
layout_new[1:7, ] <- layout_old

layout_new[8, ]  <- c(min(layout_old[,1]) + 0.5, min(layout_old[,2]) + 1.1)
layout_new[9, ]  <- c(max(layout_old[,1]) - 0.4, min(layout_old[,2]) + 1.1)

# New node names 
colnames(omega) <- c(1:9)
nodeNames <- c("memory", "RT", "prospective", "symbol", "numeric", "VNR", "TMTB", "Assortative Mating", "Ascertainment")

# Define node colors: first 7 nodes original color, next 3 nodes highlighted 
node_colors <- c(rep("#440154FF", 7), rep("#2a9d8f", 2))

# Plot 
stepwise_graph_corrected <- qgraph(omega, layout = layout_new, nodeNames = nodeNames,
                                   legend.cex = 0.6, layoutScale = c(0.9, 0.9),
                                   fade = FALSE, esize = 15, color = node_colors,
                                   border.color = "white", border.width = 2,
                                   label.color = "white", vsize = 7, curve = 0.1,
                                   curveAll = TRUE)

#------------------------------------------------
# Exploratory network with EA and Income
#-----------------------------------------------

LDSCoutput_ea <- readRDS("LDSCoutput_ea.RDS")

# estimate network
ea_model <- mod_select(LDSCoutput_ea, max_iterations = 100, verbose = TRUE)
params <- ea_model$model$parameters
fit <- ea_model$model$modelfit
omega <- as.matrix(ea_model$model$omega)

params_df <- as.data.frame(params)
fit_df <- as.data.frame(fit)  

# Create an empty layout matrix for 10 nodes
layout_new <- matrix(NA, nrow = 10, ncol = 2)
layout_new[1:7, ] <- layout_old

layout_new[8, ]  <- c(min(layout_old[,1]) + 0.5, min(layout_old[,2]) + 1.1)
layout_new[9, ]  <- c(max(layout_old[,1]) - 0.4, min(layout_old[,2]) + 1.1)
layout_new[10, ]  <- c(max(layout_old[,1]) - 0.6, min(layout_old[,2]) + 1.4)

# New node names 
colnames(omega) <- c(1:10)
nodeNames <- c("memory", "RT", "prospective", "symbol", "numeric", "TMTB", "VNR", "Assortative Mating", "Ascertainment", "EA")

# Define node colors
node_colors <- c(rep("#440154FF", 7), rep("#2a9d8f", 2), rep("orange", 2))

# Plot 
stepwise_graph_corrected <- qgraph(omega, layout = layout_new, nodeNames = nodeNames,
                                   legend.cex = 0.6, layoutScale = c(0.9, 0.9),
                                   fade = FALSE, esize = 15, color = node_colors,
                                   border.color = "white", border.width = 2,
                                   label.color = "white", vsize = 7, curve = 0.1,
                                   curveAll = TRUE)


#------------------------------------------------
# Network GWAS
#-----------------------------------------------

# # the GWAS files provided in the same order that they were listed for LDSC
files <- c("399.pairs.1.fastGWA.gz",
           "20023.rt.1.fastGWA.gz",
           "20018.prospective.1.fastGWA.gz",
           "20159.symbol.1.fastGWA.gz",
           "20240.numeric.1.fastGWA.gz",
           "20157.tmtb.1.fastGWA.gz",
           "20016.vnr.1.fastGWA.gz")

# #
# # #the reference file that matches the genetic ancestry of the included traits
ref <- "reference.1000G.maf.0.005.txt"

# # #the name of the traits
trait.names <- c("memory", "RT", "prospective", "symbol", "numeric", "TMTB", "VNR")

# # #list the sample size of each trait that has TRUE for OLS or linprob argument
# # #should reflect sum effective N for case control traits; total N for continuous traits
# # #argument only necessary when this data is not included as a GWAS column
N <- c(455695, 453043, 150165, 111471, 104620, 98075, 146808)
#
# # #whether the standard error of the GWAS estimate is on a logistic scale
# # #(only relevant for binary outcomes)
se.logit <- c(F,F,F,F,F,F,F)
#
# # #whether the phenotype was a continuous outcome analyzed using a linear estimator
OLS <- c(T,T,T,T,T,T,T)
#
# # #whether the GWAS data only includes Z-statistics -or- is a binary outcome analyzed using a linear estimator
linprob <- c(F,F,F,F,F,F,F)

SNPDATA_MET <- sumstats(files=files, ref=ref, trait.names=trait.names, se.logit=se.logit, OLS=OLS, linprob=linprob, N=N,
                     parallel = TRUE)

# #optional code to save SNP output for future use
# #saveRDS(SNPDATA_MET, file="SNPDATA_MET.RDS")
# 
SNPDATA_MET <- readRDS("SNPDATA_MET.RDS")

# ### Split SNP data into chunks for common factor GWAS
setwd("/cephfs/volumes/hpc_data_usr/k20051211/8e486e0b-f00f-4534-8dcc-d66a07c10fd5/scratch_tmp/network_final/net_gwas_rerun")

library(future)
library(future.apply)

# Number of chunks
num_chunks <- 20
n_total <- nrow(SNPDATA_MET)
chunk_size <- ceiling(n_total / num_chunks)

# Define chunk indices
chunk_indices <- rep(1:num_chunks, each = chunk_size)[1:n_total]

# Split data into chunks
data_chunks <- split(SNPDATA_MET, chunk_indices)

# Set up parallel backend
plan(multisession) # Use multiple sessions for parallel processing
options(future.globals.maxSize = 2 * 1024^3) # Increase to 2 GiB

# Save each chunk in parallel
future_lapply(1:num_chunks, function(i) {
  chunk_data <- data_chunks[[i]]
  saveRDS(chunk_data, file = paste0("SNPDATA_MET_chunk_", i, ".RDS"))
})

# ### combine networkGWAS chunks
setwd("/cephfs/volumes/hpc_data_usr/k20051211/8e486e0b-f00f-4534-8dcc-d66a07c10fd5/scratch_tmp/network_final/net_gwas_rerun")

# List all files matching the pattern "GWASNET_MET_chunk_*.RDS"
input_dir <- "/cephfs/volumes/hpc_data_usr/k20051211/8e486e0b-f00f-4534-8dcc-d66a07c10fd5/scratch_tmp/network_final/net_gwas_rerun"
rds_files <- list.files(input_dir, pattern = "GWASNET_MET_chunk_\\d+\\.RDS", full.names = TRUE)
combined_data <- list()

# Loop through each file, read the RDS, and store it in the list
for (file in rds_files) {
  message("Reading file: ", file)  # Display progress
  data <- readRDS(file)           # Read the RDS file
  combined_data[[file]] <- data   # Add it to the list
}

# Combine all data frames into one (assuming all have the same structure)
final_data <- do.call(rbind, combined_data)

#save the combined data to a new RDS or CSV file
combined_GWASNET_data <- readRDS("combined_GWASNET_data.RDS")

### create net gwas supplementary tables
library(dplyr)

# 1) define your traits
traits <- c("memory", "numeric", "prospective", "RT", "symbol", "TMTB", "VNR")

# 2) build a named list of filtered data.frames
trait_dfs <- lapply(traits, function(trait) {
  combined_GWASNET_data %>%
    # pick core cols + trait-specific cols
    select(
      SNP, CHR, BP, MAF, A1, A2,
      est    = !!sym(paste0("SNPpcor_est_", trait)),
      se     = !!sym(paste0("SNPpcor_se_",  trait)),
      Zstat  = !!sym(paste0("Zstat_",        trait)),
      p      = !!sym(paste0("p_",            trait))
    ) %>%
    # filter on p-value
    filter(p < 1e-5)
})

# 3) give the list elements names so you can pull them out easily
names(trait_dfs) <- traits

for(trait in traits) {
  df <- combined_GWASNET_data %>%
    select(
      SNP, CHR, BP, MAF, A1, A2,
      est    = !!sym(paste0("SNPpcor_est_", trait)),
      se     = !!sym(paste0("SNPpcor_se_",  trait)),
      Zstat  = !!sym(paste0("Zstat_",        trait)),
      p      = !!sym(paste0("p_",            trait))
    ) %>%
    filter(p < 1e-5)
  assign(paste0("df_", trait), df)
}

for(trait in names(trait_dfs)) {
  write.csv(
    trait_dfs[[trait]],
    file      = paste0("df_", trait, ".csv"),
    row.names = FALSE
  )
}


### Save network gwas sumstats
output_dir <- "/cephfs/volumes/hpc_data_usr/k20051211/8e486e0b-f00f-4534-8dcc-d66a07c10fd5/scratch_tmp/network_final/net_gwas_rerun"

final_data_df <- as.data.frame(final_data)

# List of traits to process
traits <- c("memory", "numeric", "prospective", "RT", "symbol", "TMTB", "VNR")

# Columns to always include
common_cols <- c("SNP", "CHR", "BP", "MAF", "A1", "A2")

for (trait in traits) {
  trait_cols <- c(paste0("SNPpcor_est_", trait),
                  paste0("SNPpcor_se_", trait),
                  paste0("Zstat_", trait),
                  paste0("p_", trait))

  required_cols <- c(common_cols, trait_cols)

  # Check if all required columns exist
  missing_cols <- setdiff(required_cols, colnames(final_data_df))
  if (length(missing_cols) > 0) {
    warning("Missing columns for trait '", trait, "': ",
            paste(missing_cols, collapse = ", "))
    next
  }

  # Subset the data.frame
  export_data <- final_data_df[, required_cols, drop = FALSE]

  # Rename trait-specific columns to a consistent format
  new_colnames <- c(
    common_cols,         # unchanged
    "SNPpcor_est",       # from SNPpcor_est_<trait>
    "SNPpcor_se",        # from SNPpcor_se_<trait>
    "Zstat",             # from Zstat_<trait>
    "p"                  # from p_<trait>
  )
  colnames(export_data) <- new_colnames

  # Write gzipped tab-separated file via a gzfile connection
  output_file <- file.path(output_dir, paste0("GWAS_", trait, ".csv.gz"))
  gz_con <- gzfile(output_file, "wt")  # write text mode
  write.table(export_data, file = gz_con, sep = "\t", row.names = FALSE, quote = FALSE)
  close(gz_con)

  message("Exported: ", output_file)
}

### Investigate significant network gwas hits
net_gwas_p <- final_data_df %>% select(SNP, p_memory, p_numeric, p_prospective, p_RT,
                                       p_symbol, p_TMTB, p_VNR)

# Apply FDR correction to each p-value column
net_gwas_p <- net_gwas_p %>%
  mutate(fdr_memory = p.adjust(p_memory, method = "fdr"),
         fdr_numeric = p.adjust(p_numeric, method = "fdr"),
         fdr_prospective = p.adjust(p_prospective, method = "fdr"),
         fdr_RT = p.adjust(p_RT, method = "fdr"),
         fdr_symbol = p.adjust(p_symbol, method = "fdr"),
         fdr_TMTB = p.adjust(p_TMTB, method = "fdr"),
         fdr_VNR = p.adjust(p_VNR, method = "fdr"))

# Determine significance (e.g., FDR < 0.05)
net_gwas_p <- net_gwas_p %>%
  mutate(sig_memory = fdr_memory < 0.05,
         sig_numeric = fdr_numeric < 0.05,
         sig_prospective = fdr_prospective < 0.05,
         sig_RT = fdr_RT < 0.05,
         sig_symbol = fdr_symbol < 0.05,
         sig_TMTB = fdr_TMTB < 0.05,
         sig_VNR = fdr_VNR < 0.05)

# View a summary of significant results
summary_significant <- net_gwas_p %>%
  summarise(across(starts_with("sig_"), sum))

## examine suggestive threshold 
suggestive_memory <- combined_GWASNET_data %>% filter(p_memory < 1e-05)
suggestive_prospective<- combined_GWASNET_data %>% filter(p_prospective < 1e-05)
suggestive_numeric <- combined_GWASNET_data %>% filter(p_numeric < 1e-05)
suggestive_symbol<- combined_GWASNET_data %>% filter(p_symbol < 1e-05)
suggestive_TMTB<- combined_GWASNET_data %>% filter(p_TMTB < 1e-05)
suggestive_VNR<- combined_GWASNET_data %>% filter(p_VNR < 1e-05)

## clumping
# LD clump

# memory
library(ieugwasr)
net_memory_vars <-suggestive_memory %>%  select(rsid = SNP, pval = p_memory)
clumpnet_memory_sig <- ld_clump(net_memory_vars,
                            clump_kb = 500,       # Clump window size reduced from 10,000 kb to 250 kb
                            clump_r2 = 0.1,       # Increase LD threshold from 0.001 to 0.1
                            opengwas_jwt = 'eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJrMjAwNTEyMTFAa2NsLmFjLnVrIiwiaWF0IjoxNzQ1OTYxNDI2LCJleHAiOjE3NDcxNzEwMjZ9.WY3RbvpEa7p-5nCzHQW7QXqWrlcN8zVcKr7uS0J7ob6_Spe1Ft0GIHf-OH6bwJnS-XbRIVQf3H-7Eq0MEmDFFUHQfouGwzaSceBTiMkN59Xzy7KwHIVbTlq4g4jQFhgb0-sdSzR_eBe8Xv9yHyCJ8PCFfqlfGgwuYG39bdwMZ7weQcBy4O1VCdurs-CkOzWd_T0G-gqadr42uxMNI1jBqOuKUWQmkoEMA0aakzwIsli103Zuy6tst1dpG_grgE18jLB06Ym91JFLuTBuTzUl0kwVeeN0JB_x5lWxdt1qMm8Vl80iiqjtmd8rqOIJG8AVH9t3kCDIKFwJmkAUsuGKCQ')
# prospective
library(ieugwasr)
net_prospective_vars <-suggestive_prospective %>%  select(rsid = SNP, pval = p_prospective)
clumpnet_prospective_sig <- ld_clump(net_prospective_vars,
                                clump_kb = 500,       # Clump window size reduced from 10,000 kb to 250 kb
                                clump_r2 = 0.1,       # Increase LD threshold from 0.001 to 0.1
                                opengwas_jwt = 'eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJrMjAwNTEyMTFAa2NsLmFjLnVrIiwiaWF0IjoxNzQ1OTYxNDI2LCJleHAiOjE3NDcxNzEwMjZ9.WY3RbvpEa7p-5nCzHQW7QXqWrlcN8zVcKr7uS0J7ob6_Spe1Ft0GIHf-OH6bwJnS-XbRIVQf3H-7Eq0MEmDFFUHQfouGwzaSceBTiMkN59Xzy7KwHIVbTlq4g4jQFhgb0-sdSzR_eBe8Xv9yHyCJ8PCFfqlfGgwuYG39bdwMZ7weQcBy4O1VCdurs-CkOzWd_T0G-gqadr42uxMNI1jBqOuKUWQmkoEMA0aakzwIsli103Zuy6tst1dpG_grgE18jLB06Ym91JFLuTBuTzUl0kwVeeN0JB_x5lWxdt1qMm8Vl80iiqjtmd8rqOIJG8AVH9t3kCDIKFwJmkAUsuGKCQ')
# numeric
library(ieugwasr)
net_numeric_vars <-suggestive_numeric %>%  select(rsid = SNP, pval = p_numeric)
clumpnet_numeric_sig <- ld_clump(net_numeric_vars,
                                     clump_kb = 500,       # Clump window size reduced from 10,000 kb to 250 kb
                                     clump_r2 = 0.1,       # Increase LD threshold from 0.001 to 0.1
                                 opengwas_jwt = 'eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJrMjAwNTEyMTFAa2NsLmFjLnVrIiwiaWF0IjoxNzQ1OTYxNDI2LCJleHAiOjE3NDcxNzEwMjZ9.WY3RbvpEa7p-5nCzHQW7QXqWrlcN8zVcKr7uS0J7ob6_Spe1Ft0GIHf-OH6bwJnS-XbRIVQf3H-7Eq0MEmDFFUHQfouGwzaSceBTiMkN59Xzy7KwHIVbTlq4g4jQFhgb0-sdSzR_eBe8Xv9yHyCJ8PCFfqlfGgwuYG39bdwMZ7weQcBy4O1VCdurs-CkOzWd_T0G-gqadr42uxMNI1jBqOuKUWQmkoEMA0aakzwIsli103Zuy6tst1dpG_grgE18jLB06Ym91JFLuTBuTzUl0kwVeeN0JB_x5lWxdt1qMm8Vl80iiqjtmd8rqOIJG8AVH9t3kCDIKFwJmkAUsuGKCQ')
# symbol
# library(ieugwasr)
# net_symbol_vars <-suggestive_symbol %>%  select(rsid = SNP, pval = p_symbol)
# clumpnet_symbol_sig <- ld_clump(net_symbol_vars,
#                                      clump_kb = 500,       # Clump window size reduced from 10,000 kb to 250 kb
#                                      clump_r2 = 0.1,       # Increase LD threshold from 0.001 to 0.1
#                                      opengwas_jwt = 'eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJrMjAwNTEyMTFAa2NsLmFjLnVrIiwiaWF0IjoxNzQ0NzE2NTAxLCJleHAiOjE3NDU5MjYxMDF9.n-A87uePkYW0J3ckBADjS5S24ACgWHKB1hMZOTR91MVFKNx5QTXh2DriF_24IVT_vs4_rhJ6slCb_tjx-HjnS85D7uVKBz7Zy5l1-W81rLyTlxBjBFe8gAum1gsWtUKdrHQhnGjrZV3uKqhqOTVBiqXwANWEuA21QZhWFG4yyeEg4Zc_uKfcX_1up13cw9Vmwp7_TbU88RC9XyErSpLwfhWSrD99LUTkf9p6XGmPKnVuT9j_r8mu3QO0Qowwng9bivkgo_S_ZEMjJ5uz5Z_aqiwORSplbaWu9J1rXoaZ7eLmPSe9jyYvSqNeREd3XOBv6buOQBN3P3gxHFSMSPfHaw')

# TMTB
# library(ieugwasr)
# net_TMTB_vars <-suggestive_TMTB %>%  select(rsid = SNP, pval = p_TMTB)
# clumpnet_TMTB_sig <- ld_clump(net_TMTB_vars,
#                                      clump_kb = 500,       # Clump window size reduced from 10,000 kb to 250 kb
#                                      clump_r2 = 0.1,       # Increase LD threshold from 0.001 to 0.1
#                                      opengwas_jwt = 'eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJrMjAwNTEyMTFAa2NsLmFjLnVrIiwiaWF0IjoxNzQ0NzE2NTAxLCJleHAiOjE3NDU5MjYxMDF9.n-A87uePkYW0J3ckBADjS5S24ACgWHKB1hMZOTR91MVFKNx5QTXh2DriF_24IVT_vs4_rhJ6slCb_tjx-HjnS85D7uVKBz7Zy5l1-W81rLyTlxBjBFe8gAum1gsWtUKdrHQhnGjrZV3uKqhqOTVBiqXwANWEuA21QZhWFG4yyeEg4Zc_uKfcX_1up13cw9Vmwp7_TbU88RC9XyErSpLwfhWSrD99LUTkf9p6XGmPKnVuT9j_r8mu3QO0Qowwng9bivkgo_S_ZEMjJ5uz5Z_aqiwORSplbaWu9J1rXoaZ7eLmPSe9jyYvSqNeREd3XOBv6buOQBN3P3gxHFSMSPfHaw')

# VNR
library(ieugwasr)
net_VNR_vars <-suggestive_VNR %>%  select(rsid = SNP, pval = p_VNR)
clumpnet_VNR_sig <- ld_clump(net_VNR_vars,
                              clump_kb = 500,       # Clump window size reduced from 10,000 kb to 250 kb
                              clump_r2 = 0.1,       # Increase LD threshold from 0.001 to 0.1
                             opengwas_jwt = 'eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJrMjAwNTEyMTFAa2NsLmFjLnVrIiwiaWF0IjoxNzQ1OTYxNDI2LCJleHAiOjE3NDcxNzEwMjZ9.WY3RbvpEa7p-5nCzHQW7QXqWrlcN8zVcKr7uS0J7ob6_Spe1Ft0GIHf-OH6bwJnS-XbRIVQf3H-7Eq0MEmDFFUHQfouGwzaSceBTiMkN59Xzy7KwHIVbTlq4g4jQFhgb0-sdSzR_eBe8Xv9yHyCJ8PCFfqlfGgwuYG39bdwMZ7weQcBy4O1VCdurs-CkOzWd_T0G-gqadr42uxMNI1jBqOuKUWQmkoEMA0aakzwIsli103Zuy6tst1dpG_grgE18jLB06Ym91JFLuTBuTzUl0kwVeeN0JB_x5lWxdt1qMm8Vl80iiqjtmd8rqOIJG8AVH9t3kCDIKFwJmkAUsuGKCQ')

# sig rt
rt_sig <- combined_GWASNET_data %>% filter(p_RT < 5e-08)
net_rt_vars <-rt_sig %>%  select(rsid = SNP, pval = p_memory)
clumpnet_rt_sig <- ld_clump(net_rt_vars,
                                clump_kb = 500,       # Clump window size reduced from 10,000 kb to 250 kb
                                clump_r2 = 0.1,       # Increase LD threshold from 0.001 to 0.1
                                opengwas_jwt = 'eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJrMjAwNTEyMTFAa2NsLmFjLnVrIiwiaWF0IjoxNzQ1OTYxNDI2LCJleHAiOjE3NDcxNzEwMjZ9.WY3RbvpEa7p-5nCzHQW7QXqWrlcN8zVcKr7uS0J7ob6_Spe1Ft0GIHf-OH6bwJnS-XbRIVQf3H-7Eq0MEmDFFUHQfouGwzaSceBTiMkN59Xzy7KwHIVbTlq4g4jQFhgb0-sdSzR_eBe8Xv9yHyCJ8PCFfqlfGgwuYG39bdwMZ7weQcBy4O1VCdurs-CkOzWd_T0G-gqadr42uxMNI1jBqOuKUWQmkoEMA0aakzwIsli103Zuy6tst1dpG_grgE18jLB06Ym91JFLuTBuTzUl0kwVeeN0JB_x5lWxdt1qMm8Vl80iiqjtmd8rqOIJG8AVH9t3kCDIKFwJmkAUsuGKCQ')


#######################################################
# Correlation between net and regular sumstats
#####################################################

# ### Munge sumstats
setwd("/cephfs/volumes/hpc_data_usr/k20051211/8e486e0b-f00f-4534-8dcc-d66a07c10fd5/scratch_tmp/network_final/net_gwas_rerun")

#create vector of the summary statistics files
files <- c("GWAS_memory.csv.gz",
           "GWAS_RT.csv.gz",
           "GWAS_prospective.csv.gz",
           "GWAS_symbol.csv.gz",
           "GWAS_numeric.csv.gz",
           "GWAS_TMTB.csv.gz",
           "GWAS_VNR.csv.gz")

# #define the reference file being used to allign alleles across summary stats
# #here we are using hapmap3
hm3<-"/cephfs/volumes/hpc_data_usr/k20051211/8e486e0b-f00f-4534-8dcc-d66a07c10fd5/scratch_tmp/network_final/eur_w_ld_chr/w_hm3.snplist"

#the name of the traits
trait.names <- c("network_memory", "network_RT", "network_prospective", "network_symbol",
                 "network_numeric", "network_TMTB", "network_VNR")

#list the sample size of each trait that has TRUE for OLS or linprob argument
#should reflect sum effective N for case control traits; total N for continuous traits
#argument only necessary when this data is not included as a GWAS column
N <- c(455695, 453043, 150165, 111471, 104620, 98075, 146808)

#definte the imputation quality filter
info.filter=0.9

#define the MAF filter
maf.filter=0.01

#run munge
munge(files=files,hm3=hm3,trait.names=trait.names,N=N,info.filter=info.filter,maf.filter=maf.filter,
      parallel = TRUE)

### Calculate h2 network
devtools::install_github("mglev1n/ldscr", force = TRUE)
library(ldscr)

h2_memory_net <- ldsc_h2(munged_sumstats = "network_memory.sumstats.gz", ancestry = "EUR")
h2_rt_net <- ldsc_h2(munged_sumstats = "network_RT.sumstats.gz", ancestry = "EUR")
h2_prospective_net <- ldsc_h2(munged_sumstats = "network_prospective.sumstats.gz", ancestry = "EUR")
h2_symbol_net <- ldsc_h2(munged_sumstats = "network_symbol.sumstats.gz", ancestry = "EUR")
h2_numeric_net <- ldsc_h2(munged_sumstats = "network_numeric.sumstats.gz", ancestry = "EUR")
h2_tmtb_net <- ldsc_h2(munged_sumstats = "network_TMTB.sumstats.gz", ancestry = "EUR")
h2_vnr_net <- ldsc_h2(munged_sumstats = "network_VNR.sumstats.gz", ancestry = "EUR")

### Calculate h2 regular
setwd("/cephfs/volumes/hpc_data_usr/k20051211/8e486e0b-f00f-4534-8dcc-d66a07c10fd5/scratch_tmp/network_final/")

h2_memory <- ldsc_h2(munged_sumstats = "memory.sumstats.gz", ancestry = "EUR")
h2_rt <- ldsc_h2(munged_sumstats = "RT.sumstats.gz", ancestry = "EUR")
h2_prospective <- ldsc_h2(munged_sumstats = "prospective.sumstats.gz", ancestry = "EUR")
h2_symbol <- ldsc_h2(munged_sumstats = "symbol.sumstats.gz", ancestry = "EUR")
h2_numeric <- ldsc_h2(munged_sumstats = "numeric.sumstats.gz", ancestry = "EUR")
h2_tmtb <- ldsc_h2(munged_sumstats = "TMTB.sumstats.gz", ancestry = "EUR")
h2_vnr <- ldsc_h2(munged_sumstats = "VNR.sumstats.gz", ancestry = "EUR")

# barplot
# Load required libraries
library(tidyr)
library(dplyr)
library(ggplot2)

# # Define trait names in the order you want them to appear
traits <- c("memory", "RT", "prospective", "symbol", "numeric", "TMTB", "VNR")

# Create a data frame with the h2 estimates and standard errors
df_h2 <- data.frame(
  trait = rep(traits, 2),
  type = rep(c("network", "standard"), each = length(traits)),
  h2 = c(
    h2_memory_net$h2_observed,
    h2_rt_net$h2_observed,
    h2_prospective_net$h2_observed,
    h2_symbol_net$h2_observed,
    h2_numeric_net$h2_observed,
    h2_tmtb_net$h2_observed,
    h2_vnr_net$h2_observed,
    h2_memory$h2_observed,
    h2_rt$h2_observed,
    h2_prospective$h2_observed,
    h2_symbol$h2_observed,
    h2_numeric$h2_observed,
    h2_tmtb$h2_observed,
    h2_vnr$h2_observed
  ),
  h2_se = c(
    h2_memory_net$h2_observed_se,
    h2_rt_net$h2_observed_se,
    h2_prospective_net$h2_observed_se,
    h2_symbol_net$h2_observed_se,
    h2_numeric_net$h2_observed_se,
    h2_tmtb_net$h2_observed_se,
    h2_vnr_net$h2_observed_se,
    h2_memory$h2_observed_se,
    h2_rt$h2_observed_se,
    h2_prospective$h2_observed_se,
    h2_symbol$h2_observed_se,
    h2_numeric$h2_observed_se,
    h2_tmtb$h2_observed_se,
    h2_vnr$h2_observed_se
  )
)

saveRDS(df_h2, "df_h2.rds")
df_h2 <- readRDS("df_h2.rds")

# Create the barplot with error bars
ggplot(df_h2, aes(x = trait, y = h2, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +  # Barplot for h2 values
  geom_errorbar(aes(ymin = h2 - h2_se, ymax = h2 + h2_se),
                position = position_dodge(width = 0.9), width = 0.25) +  # Error bars for SE
  labs(x = "Trait", y = "Heritability (h2)", fill = "Version") +
  theme_minimal(base_size = 12) +
  ggtitle("Heritability Estimates (h2) with Standard Errors\n(Network vs. Standard)")


### calculate rg between regular and non-network
rg_memory <- ldsc_rg(
  munged_sumstats = list(
    "memory" = "memory.sumstats.gz",
    "memory_net" = "net_gwas_rerun/network_memory.sumstats.gz"
  ),
  ancestry = "EUR"
)

rg_RT <- ldsc_rg(
  munged_sumstats = list(
    "RT" = "RT.sumstats.gz",
    "RT_net" = "net_gwas_rerun/network_RT.sumstats.gz"
  ),
  ancestry = "EUR"
)

rg_prospective <- ldsc_rg(
  munged_sumstats = list(
    "prospective" = "prospective.sumstats.gz",
    "prospective_net" = "net_gwas_rerun/network_prospective.sumstats.gz"
  ),
  ancestry = "EUR"
)

rg_symbol <- ldsc_rg(
  munged_sumstats = list(
    "symbol" = "symbol.sumstats.gz",
    "symbol_net" = "net_gwas_rerun/network_symbol.sumstats.gz"
  ),
  ancestry = "EUR"
)

rg_numeric <- ldsc_rg(
  munged_sumstats = list(
    "numeric" = "numeric.sumstats.gz",
    "numeric_net" = "net_gwas_rerun/network_numeric.sumstats.gz"
  ),
  ancestry = "EUR"
)

rg_TMTB <- ldsc_rg(
  munged_sumstats = list(
    "TMTB" = "TMTB.sumstats.gz",
    "TMTB_net" = "net_gwas_rerun/network_TMTB.sumstats.gz"
  ),
  ancestry = "EUR"
)

rg_VNR <- ldsc_rg(
  munged_sumstats = list(
    "VNR" = "VNR.sumstats.gz",
    "VNR_net" = "net_gwas_rerun/network_VNR.sumstats.gz"
  ),
  ancestry = "EUR"
)

# forest plot of rg
# Load required package
library(ggplot2)

# Create a data frame with your genetic correlation estimates and standard errors.
# For each trait, we compute the lower and upper bounds of the 95% CI.
df <- data.frame(
  trait = c("memory", "RT", "prospective", "symbol", "numeric", "TMTB", "VNR"),
  rg = c(
    rg_memory$rg$rg,
    rg_RT$rg$rg,
    rg_prospective$rg$rg,
    rg_symbol$rg$rg,
    rg_numeric$rg$rg,
    rg_TMTB$rg$rg,
    rg_VNR$rg$rg
  ),
  se = c(
    rg_memory$rg$rg_se,
    rg_RT$rg$rg_se,
    rg_prospective$rg$rg_se,
    rg_symbol$rg$rg_se,
    rg_numeric$rg$rg_se,
    rg_TMTB$rg$rg_se,
    rg_VNR$rg$rg_se
  )
)

# Compute the 95% confidence intervals
df$lower <- df$rg - 1.96 * df$se
df$upper <- df$rg + 1.96 * df$se

# For a nice ordering on the plot, make 'trait' a factor.
df$trait <- factor(df$trait, levels = df$trait)

# Create the forest plot.
# We use coord_flip() to swap axes so that trait names appear on the y-axis.
ggplot(df, aes(x = rg, y = trait)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Genetic Correlation (rg)",
       y = "Trait",
       title = "Forest Plot of Genetic Correlations (Regular vs. Network GWAS)") +
  theme_minimal(base_size = 12)

#-------------------------------------------------------------------------
# S-LDSC analysis
#-------------------------------------------------------------------------
results_dir <- "/cephfs/volumes/hpc_data_usr/k20051211/8e486e0b-f00f-4534-8dcc-d66a07c10fd5/scratch_tmp/s_ldsc/ldsc/results_updated/"

# 2) list all .results files
files <- list.files(results_dir, pattern="\\.results$", full.names=TRUE)

# 3) for each file: parse its name, read it in, grab the L2_1 row, add metadata
res_list <- lapply(files, function(fpath) {
  fname     <- basename(fpath)
  
  # determine 'network' vs 'standard'
  type      <- if (substr(fname, 1, 8) == "network_") "network" else "standard"
  
  # strip leading "network_" and trailing ".results"
  base      <- sub("^network_", "", fname)
  base      <- sub("\\.results$",   "", base)
  
  # first underscore‐piece is trait; the rest joined are cell_type
  parts     <- strsplit(base, "_", fixed=TRUE)[[1]]
  trait     <- parts[1]
  cell_type <- paste(parts[-1], collapse = "_")
  
  # read in with read.delim()
  df        <- read.delim(fpath, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  
  # subset to the Category == "L2_1" row
  row       <- df[df$Category == "L2_1", ]
  
  # add our metadata columns
  row$trait     <- trait
  row$type      <- type
  row$cell_type <- cell_type
  
  row
})

# 4) bind into one big data.frame
combined_df <- do.call(rbind, res_list)

# 5) reorder so metadata come first
meta_cols   <- c("trait","type","cell_type")
other_cols  <- setdiff(names(combined_df), meta_cols)
combined_df <- combined_df[ , c(meta_cols, other_cols) ]

### create forest plot 
library(ggplot2)

# 1) compute FDR and a “star” indicator
combined_df$FDR <- p.adjust(combined_df$Enrichment_p, method = "fdr")
combined_df$star <- ifelse(combined_df$FDR < 0.05, "*", "")

# 2) order factor levels so the smallest enrichment appears at the bottom
combined_df$cell_type <- with(combined_df, 
                              reorder(cell_type, Enrichment)
)

library(dplyr)
library(ggplot2)
library(forcats)

combined_df %>%
  filter(trait == "symbol") %>%
  
  mutate(cell_type = fct_reorder(cell_type, Enrichment)) %>%
  
  ggplot(aes(x = Enrichment, y = cell_type, color = type)) +
  
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey70") +
  
 
  geom_errorbarh(aes(
    xmin = Enrichment - Enrichment_std_error,
    xmax = Enrichment + Enrichment_std_error
  ),
  height = 0.2,           
  size   = 0.5,           
  color  = "grey50") +
  

  geom_point(size = 3) +
  
 
  geom_text(aes(label = star),
            hjust = -0.3,     
            vjust = 0.5,      
            size  = 5,
            show.legend = FALSE) +
  
 
  scale_color_brewer(
    palette = "Set1",
    name    = "Enrichment type",
    labels  = c("Network", "Standard")
  ) +
  

  labs(
    x        = "Enrichment (± SE)",
    y        = NULL,                            
    title    = "S-LDSC Enrichment by Neuronal Cell Type\n(Symbol)",
    subtitle = "Dashed line = no enrichment; “*” = FDR < 0.05"
  ) +
  
 
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),      
    panel.grid.minor   = element_blank(),
    axis.text.y        = element_text(size = 10),
    legend.position    = "bottom",
    plot.title         = element_text(face = "bold", size = 14),
    plot.subtitle      = element_text(size = 10)
  ) +
  

  coord_cartesian(
    xlim = c(
      0,
      max(combined_df$Enrichment + combined_df$Enrichment_std_error, na.rm = TRUE) * 1.1
    )
  )
