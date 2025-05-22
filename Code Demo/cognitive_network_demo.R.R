#---------------------------------------------------------------------------------------------------------
# Code demo for: "Multi-Scale Characterization of Pleiotropy and Unique Genetic Influences on Cognitive Abilities"
# This script walksthrough how to estimate an undirected network model from LDSC-derived genetic correlations
# using our implementation of the ggmModSelect algorithm.
#
# Written by: Kaito Kawakami and Jacob Knypsel
#---------------------------------------------------------------------------------------------------------

# Set working directory
setwd("C:/Users/kaito/Desktop/NHB forms/Code Demo")

# Install packages 
install.packages("devtools") # install devtools if needed
library(devtools)

install_github("GenomicNetworkAnalysis/GNA")
install_github("GenomicSEM/GenomicSEM")
install.packages('psychonetrics')
install.packages('glasso')
install.packages('Matrix')
install.packages('gdata')

# load libraries
library(GNA)
library(GenomicSEM)
library(psychonetrics)
library(glasso)
library(Matrix) # for forceSymmetric() function
library(gdata) # for lowerTriangle() function

##########################################################################################
# Step 1. Load in LDSC output 
# Using the 7 GWAS summary statistics of cognitive abilities, we ran them through the
# munge() function from the GenomicSEM package to format the summary statistics,
# and then ran them through the ldsc() function to get the genetic correlation matrix.
##########################################################################################

# Load the LDSC output
LDSCoutput_final <- readRDS("LDSCoutput_final.RDS")

corrs <- LDSCoutput_final$S_Stand
rownames(corrs) <- colnames(corrs)


##########################################################################################
# Step 2. Define functions for genomic network analysis
# First, we define some functions from the GNA package to run the network analysis.
##########################################################################################

# Import GNA function to get robust standard errors
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

# Import GNA function to estimate Gaussian Graphical Model (network model) from LDSC object
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

# Import function to convert precision matrix to partial correlation matrix from the qgraph R package
wi2net <- function(x) {
  x <- -cov2cor(x)
  diag(x) <- 0
  x <- forceSymmetric(x)
  return(x)
}

### Next is our implementation of the ggmodselect() algorithm developed by Sasha Epskamp from the
### qgraph R package.
### See: https://rdrr.io/cran/qgraph/man/ggmModSelect.html

## Step 1: .ebic_unregularized: initialize GGM via EBIC‐guided glasso path 
##
## 1) Fit fully saturated GGM to get baseline model. 
## 2) Compute 100 logarthimically spaced rho values to explore using glasso (graphical lasso).
## 3) Run glassopath on candidate rho values to obtain precision matrix estimates.
## 4) For each glasso-estimated precision matrix, re‐fit unregularized GGM by creating a constraining
##    matrix (fix_omega) that sets all edges estimated as zero to zero, and freely estimates the others.    
##    This is done using the runGGM function.
## 5) Compute EBIC_γ = χ² + |E|·log(N_eff) + 4γ|E|·log(P); choose λ minimizing EBIC (taken from the GGM package)
## 6) Return best λ, precision & partial‐corr matrices, plus EBIC path.
## 7) This weight matrix serves as a starting point for the next step!
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

## Step 2: mod_select: stepwise search for best GGM
##
## 1) Using the weight matrix obtained by ebic_unregularized, we run a stepwise search for the best GGM.
##    by iteratively adding and removing each possible edge until no further improvement is found in EBIC.
##
##    Note: I added an edge counter but it currently does not work! Displays the same number of edges irrespective of how many are added/removed
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

##########################################################################################
# Step 3. Estimate network model
##########################################################################################

# Estimate network model using ggmModSelect
network_model <- mod_select(LDSCoutput_final, max_iterations = 100, verbose = TRUE)

# Obtain parameters and model fit
params <- network_model$model$parameters
fit <- network_model$model$modelfit

# Check to see if GNA estimates the same parameters/fit using our weight matrix
omega <- network_model$model$omega 
omega <- ifelse(omega == 0, FALSE, TRUE) # This is the weight matrix

test <- traitNET(covstruc = LDSCoutput_final, graph_layout="spring", 
                 fix_omega = omega, # this is the weight matrix from mod_select
                 reestimate = FALSE, 
                 recursive = FALSE, 
                 prune = FALSE)

# Same results and fit
test$model_results
test$model_results$base$modelfit

# compare to p-value thresholding
p_val_threshold <- traitNET(covstruc = LDSCoutput_final, graph_layout="spring")
p_val_threshold$model_results$sparse$modelfit # poorer fit in our case

# Plot
library(qgraph)

omega <- as.matrix(network_model$model$omega)
colnames(omega) <- c(1:7)

# Trait names with numbers
nodeNames <- c("memory", "RT", "prospective", "symbol", "numeric", "TMTB", "VNR")

# This is Figure 3A in the paper
stepwise_graph <- qgraph(omega, layout = "spring", nodeNames = nodeNames, legend.cex = 0.7,
                         layoutScale = c(0.9,0.9), fade = FALSE, esize = 18,
                         color = c("#440154FF"), border.color = "white", border.width = 2, label.color = "white", 
                         vsize = 9, curve = 0.1, curveAll = TRUE)

