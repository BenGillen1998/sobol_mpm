library(sensitivity)


## Define 2 samples

N <- 5000 # number of Monte Carlo Samples for Quasi Monte Carlo method

# Define a function that will generate a sample matrix for all our parameters
sample_params <- function(N = 5000) {
  out <- data.frame(Fec = runif(N, min = 0, max = 2),   ##### sample for fecundity
                    P_AJ = runif(N, min = 0, max = 1),         ##### sample for prob juvenile grow up
                    P_AA = runif(N, min = 0, max = 1))          ##### sample for prob adults survive
  
  return(out)
}

# Take 2 samples
sample1 <- sample_params()
sample2 <- sample_params()


## Define function for output values of interest (if needed)
eigen_val <- function(sample_matrix) {
  # Initialize a vector to store results
  output <- numeric(nrow(sample_matrix))

  # Iterate through the rows of the sample matrix
  for (i in 1:nrow(sample_matrix)) {
    
    # add a 0
    params <- c(0, as.matrix(sample_matrix[i, ]))
    
    # convert vector of parameters into a matrix
    mat <- matrix(params, byrow = TRUE, nrow = 2)
  
    # find eigenvalues of matrix
    out <- eigen(mat)
    
    # store the dominant eigenvalue
   output[i] <- out$values[1]
  }
  
  # return the vector of dominant eigenvalues
  return(output)
}


## 
sobol_indices <- sobol(model = eigen_val,
                       X1 = sample1,
                       X2 = sample2,
                       order = 3)

sobol_indices
