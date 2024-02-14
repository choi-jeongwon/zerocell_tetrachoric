# -------------------------------------------------------------------------------------------------------------------------------#
# R code for tetrachoric correlation estimation and efficient computations                                                       #
# for the paper 'The Impact of Different Zero-Frequency Cell Correction Strategies on the Tetrachoric Correlation Estimation'    #
# Last updated: January 17, 2024                                                                                                 #
# -------------------------------------------------------------------------------------------------------------------------------#

# This R script demonstrates codes for computing and estimating tetrachoric correlations, as outlined in the paper.
# The code for estimation is a modified version of the polychor function from the polycor package (Fox, 2022).

# References
# Fox, J. (2022). polycor: Polychoric and polyserial correlations. (R package version 0.8-1). Retrieved from https://CRAN.R-project.org/package=polycor


#------------------#
# Section 0: Setup
#------------------#

#---  Load libraries---#
if(!require(tidyverse)){install.packages('tidyverse')}
if(!require(mvtnorm)){install.packages('mvtnorm')}
if(!require(MASS)){install.packages('MASS')}
if(!require(numDeriv)){install.packages('numDeriv')}

#--- Set seed for reproducibility ---#
set.seed(2023)

#--- example table for demonstration purpose ---#
x <- matrix(c(40,5,3,2),ncol=2,byrow=FALSE) 


# ===========================================================================
# Part I: Bivariate Simulation
# ===========================================================================

#-----------------------------#
# Section 1: Data Generation
#-----------------------------#

#--- Data generation based on sample size, correlation, and thresholds ---#

gen_data <- function(sample=100, correlation=0.3, thresh=c(1.5, 1.5),replications=1000){
  
  cormat <- matrix(c(1, correlation, correlation, 1), nrow=2, ncol=2)
  p00_p <- pmvnorm(lower = c(-Inf, -Inf), upper = thresh, corr = cormat)
  p10_p <- pmvnorm(lower=c(thresh[1],-Inf), upper = c(Inf,thresh[2]), corr = cormat)
  p01_p <- pmvnorm(lower=c(-Inf, thresh[2]), upper = c(thresh[1],Inf), corr = cormat)
  p11_p <- pmvnorm(lower = thresh, upper = c(Inf, Inf), corr = cormat)
  
  p00 <- p00_p/(p00_p+p10_p+p01_p+p11_p)
  p10 <- p10_p/(p00_p+p10_p+p01_p+p11_p)
  p01 <- p01_p/(p00_p+p10_p+p01_p+p11_p)
  p11 <- p11_p/(p00_p+p10_p+p01_p+p11_p)
  
  ctb <- matrix(c(p00, p10, p01, p11), 2, 2, byrow=TRUE)
  
  ctv_v <- c(ctb)
  
  samp <- rmultinom(replications, sample, ctv_v)
  results <- t(samp)
  
  return(results)
}

#---------------------------------------------------#
# Section 2: Estimation of Tetrachoric Correlations
#---------------------------------------------------#

#--- Define the function f to be minimized ---#

f <- function(pars) {
  
  if (length(pars) == 1) {
    rho <- pars
    row.cuts <- rc
    col.cuts <- cc
  } else {
    rho <- pars[1]
    row.cuts <- pars[2:r]
    col.cuts <- pars[(r + 1):(r + c - 1)]
    
    if (any(diff(row.cuts) < 0) || any(diff(col.cuts) < 0)) {
      return(Inf)
    }
    
    if (abs(rho)>=1) {
       return(NA) #Assign NA if correlation goes over [-1,1]
    }
  }   
  P <- newbinBvn(rho, row.cuts, col.cuts) 
  N<-sum(tab)

  -sum(tab[tab>0] * log(P[tab>0]/tab[tab>0]*N))/N # loss function
}


#--- calculate thresholds (rc, cc) ---#

# To use unadjusted thresholds in the second stage, tab should be uncorrected (raw) table.
tab <- x
n <- sum(tab)
r <- nrow(tab)
c <- ncol(tab)
rc <- qnorm(cumsum(rowSums(tab))/n)[-r]
cc <- qnorm(cumsum(colSums(tab))/n)[-c]


#--- Internal function (original function name used in polycor package: binBvn) ---#

newbinBvn <- function (rho, row.cuts, col.cuts, bins = 4) 
{
  row.cuts <- if (missing(row.cuts)) 
    c(-Inf, 1:(bins - 1)/bins, Inf)
  else c(-Inf, row.cuts, Inf)
  col.cuts <- if (missing(col.cuts)) 
    c(-Inf, 1:(bins - 1)/bins, Inf)
  else c(-Inf, col.cuts, Inf)
  r <- length(row.cuts) - 1
  c <- length(col.cuts) - 1
  P <- matrix(0, r, c)
  R <- matrix(c(1, rho, rho, 1), 2, 2)
  for (i in 1:r) {
    for (j in 1:c) {
      P[i, j] <- mvtnorm::pmvnorm(lower = c(row.cuts[i], col.cuts[j]), upper = c(row.cuts[i + 1], col.cuts[j + 1]), corr = R)
    }
  }
  P
}


#--- optimize f to get tetrachoric correlation estimates ---#
# Set the starting value with the rho estimated with one-dimensional optimization

rho_2s <- optimise(f, interval = c(-1, 1))$minimum

# Get correlation estimate with maximum likelihood

rho <- optim(par=c(rho_2s, rc, cc), fn=f, method="Nelder-Mead", hessian = FALSE)$par[1]

#--- estimate standard error using numDeriv package ---#

se <- sqrt(solve(hessian(f, x=c(rho,rc,cc),  method="Richardson", method.args=list(eps=0.04, d=(1-abs(rho))/2, zero.tol=0.9)))[1,1])/sqrt(n)


#--- print the result ---#

result <- list('estimate' = rho, 'row threshold' = rc, 'column threshold' =cc, 'standard error' = se)
print(result)


#------------------------------------------------------------------#
# Section 3: Demonstration of Strategies for Efficient Computations
#------------------------------------------------------------------#
# This section demonstrates ways to make the simulations more efficient when there are many redundant tables need to be analyzed.
  # It consists of 1) eliminating the unnecessary simulation conditions, and 2) using prototypes for similar tables. 


#--- Simplifying the conditions by reducing the number of thresholds ---#

# make thresholds sets with the following rule:the first threshold>0 & abs(first threshold)>=abs(second threshold), rather than using every possible threshold combination.

threshset <- function(thresh){
  th.list <- list()
  counter <- 1 
  for (i in 1:length(thresh)){
    for (j in 1:length(thresh)){
      this_threshs <- c(thresh[i], thresh[j])
      if((abs(this_threshs[1])>=abs(this_threshs[2]))& this_threshs[1]>0){
        th.list[[counter]] <- this_threshs
        counter <- counter +1 
      }
    }
  }
  return(th.list)
}

  #  For example, threshset(c(1.5, 1, 0.8, -0.8, -1, -1.5)) produces 12 possible threshold sets.


#--- Using prototypes for similar tables ---#

# 1) Generate 8 different types of tables for a given table.
do.perm <- function(arr){
  # Arguments: arr is a vector of 4 elements
  # Return: mat (matrix of 8 tables that yield similar or identical rho, t1, and t2)
  a <- arr[1]; b <- arr[2]; c <- arr[3]; d <- arr[4]
  
  mat <- matrix(
    c(a,b,c,d,
      a,c,b,d,
      d,c,b,a,
      d,b,c,a,
      c,a,d,b,
      c,d,a,b,
      b,d,a,c,
      b,a,d,c),
    nrow=8, ncol=4,byrow=T)
  return(mat)
}


# 2) make a prototype for a given table
make.proto <- function(this_row){
  # Arguments: this row: This function takes a row from raw data (this_row)
  # Returns: It returns a vector c(proto,var.type)
  #  proto: prototype table
  #  var.type: variant type of given table (row index in do.perm(proto))
  
  this_row <- as.numeric(this_row[1:4])
  # make a list of all possible permutations
  mat<- do.perm(this_row)
  
  # determine which is the prototypical form
  idx <- mat[,1]>=mat[,2] & mat[,2]>=mat[,3] & mat[,1]>=mat[,4]
  proto <- as.numeric(mat[match(T,idx),])

  # mark vartype for the given table (determining vartype based on do.perm(proto))
  mat2 <- do.perm(proto)
  var.idx <- mat2[,1]==this_row[1] & mat2[,2]==this_row[2] & mat2[,3]==this_row[3] & mat2[,4]==this_row[4] 
  var.type <- match(T, var.idx)
  
  return(c(proto,var.type))
}


# 3) update correlation and thresholds based on the vartype
change_by_var <- function(df){
  # Arguments: df: a dataframe with estimation results only for prototype tables
  # Returns: updated dataframe with all estimates
  
  # a list consisting of index vectors (for each vartype value. from 1 to 8)
  list_idx <- lapply(X=1:8,FUN=function(x){ return(x==df$vartype) })
  
  rho_idx <- which(startsWith(colnames(df),"rho_"))
  tau1_idx <- which(startsWith(colnames(df),"tau1_"))
  tau2_idx <- which(startsWith(colnames(df),"tau2_"))
  
# change estimates based on vartype values
  # vartype==2, rho=rho; tau1=tau2; tau2=tau1
  temp <- df[list_idx[[2]], tau2_idx]
  df[list_idx[[2]], tau2_idx] <- df[list_idx[[2]], tau1_idx]
  df[list_idx[[2]], tau1_idx] <- temp
  
  # vartype==3, rho=rho; tau1=-tau1; tau2=-tau2
  df[list_idx[[3]], tau1_idx] <- -df[list_idx[[3]], tau1_idx]
  df[list_idx[[3]], tau2_idx] <- -df[list_idx[[3]], tau2_idx]
  
  # vartype==4, rho=rho; tau1=-tau2; tau2=-tau1
  temp <- df[list_idx[[4]], tau2_idx]
  df[list_idx[[4]], tau2_idx] <- -df[list_idx[[4]], tau1_idx]
  df[list_idx[[4]], tau1_idx] <- -temp
  
  # vartype==5, rho=-rho; tau1=tau2; tau2=-tau1
  temp <- df[list_idx[[5]], tau2_idx]
  df[list_idx[[5]], rho_idx] <- -df[list_idx[[5]], rho_idx]
  df[list_idx[[5]], tau2_idx] <- -df[list_idx[[5]], tau1_idx]
  df[list_idx[[5]], tau1_idx] <- temp
  
  # vartype==6, rho=-rho; tau1=-tau1; tau2=tau2
  df[list_idx[[6]], rho_idx] <- -df[list_idx[[6]], rho_idx]
  df[list_idx[[6]], tau1_idx] <- -df[list_idx[[6]], tau1_idx]
  
  # vartype==7, rho=-rho; tau1=-tau2; tau2=tau1
  temp <- df[list_idx[[7]], tau2_idx]
  df[list_idx[[7]], rho_idx] <- -df[list_idx[[7]], rho_idx]
  df[list_idx[[7]], tau2_idx] <- df[list_idx[[7]], tau1_idx]
  df[list_idx[[7]], tau1_idx] <- -temp
  
  # vartype==8, rho=-rho; tau1=tau1; tau2=-tau2
  df[list_idx[[8]], rho_idx] <- -df[list_idx[[8]], rho_idx]
  df[list_idx[[8]], tau2_idx] <- -df[list_idx[[8]], tau2_idx]
  
  return(df)
}


# ===========================================================================
# Part II: Multivariate Simulation Data Generation
# ===========================================================================

multi_gen_data <- function(n_samples, means, covmat, thresholds, returntable=FALSE){
  # Step 1. Generate data from multivariate (6) variate normal distribution
  # Step 2. Discretize the variables with thresholds
  # Step 3. Make contingency tables for each pair of variables. (total 15 tables)
  
  # if returntable=TRUE, return the contingency table. If F, return raw data
  
  samples <- mvrnorm(n_samples, means, covmat)
  colnames(samples) <- c('V1','V2','V3','V4','V5','V6')
  samples_th <- matrix(0, nrow = nrow(samples), ncol = ncol(samples))
  colnames(samples_th) <- c('V1','V2','V3','V4','V5','V6')
  
  # discretize with thresholds
  for (i in 1:ncol(samples)){
    samples_th[, i] <- as.numeric(samples[,i] > thresholds[i])
  }
  
  if (returntable){
    comb_pairs <- combn(1:ncol(samples),m=2)
    df_tb <- tibble()
    for (i in 1:ncol(comb_pairs)){
      this_rawdata <- cbind(samples_th[,comb_pairs[1,i]],samples_th[,comb_pairs[2,i]])
      this_table <- table(factor(this_rawdata[,1],levels=c(0,1)),
                          factor(this_rawdata[,2],levels=c(0,1)),
                          dnn=c('1)','2)'))
      df_tb <- rbind(df_tb,tibble(V1=this_table[1,1],V2=this_table[2,1],V3=this_table[1,2],V4=this_table[2,2]))
    }

    df_tb <- df_tb %>% add_column(pairs = paste(comb_pairs[1,],sep="&",comb_pairs[2,]))
    
    return(df_tb)
  }
  return(samples_th)
}


#--- example n_samples, means, covmat, thresholds for cor 0.4 & mixed threshold condition ---#
  # n_samples <- 50
  # means <- rep(0,6)
  # covmat <- matrix(0.4,nrow=6,ncol=6)
  # diag(covmat) <- 1
  # thresholds <- c(1.5, 1.0, 0.8, -1.5, -1.0, -0.8)


# For estimation of multivariate simulation, we first search for previously computed rho and standard errors in the bivariate simulation. 
# If the table was not present in the bivariate result, we conducted estimation using the codes used for bivariate analysis, with adjusted thresholds.