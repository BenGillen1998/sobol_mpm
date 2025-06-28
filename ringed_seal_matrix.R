library(sensitivity)


## Define 2 samples

N <- 5000 # number of Monte Carlo Samples for Quasi Monte Carlo method

# Define a function that will generate a sample matrix for all our parameters
sample_params <- function(N = 5000) {
  
  out <- matrix(NA, nrow = N, ncol = 12)
  colnames(out) <- c("P0", "P1","P2","P3","P4","P5","P6","P7","m4","m5","m6","m7")
  
  # Vector of mean transition probabilities (P0, P1, ..., P7 in Jody's paper)
  P_means <- c(0.65, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92) 
  
  # Vector of mean Fecundity parameters (m4, m5, m6, m7 in Jody's paper)
  # from a normal year
  m_means <- c(0.098, 0.144, 0.195, 0.406)
  
  all_means <- c(P_means, m_means)
  
  for(i in 1:12) {
    out[ , i] <- rnorm(N, mean = all_means[i], sd = 0.1)
    
    # Truncate the distribution at 1 for the probabilities
    if (i <= 8) {
      out[ , i] <- pmin(out[ , i], 1)
    }
    
    # set negative values to 0
    out[ , i] <- pmax(out[ , i], 0)
  }
  
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
    
    # set up a matrix of zeros
    mat <- matrix(0, nrow = 8, ncol = 8)
    
    # add sampled values to the appropriate spot in the matrix
    for (j in 1:7) {
      mat[j+1,j] <- sample_matrix[i, j]
    }
    
    mat[8,8] <- sample_matrix[i, 8]
    
    
    mat[1, 5] <- sample_matrix[i, 9]*sample_matrix[i, 5]
    mat[1, 6] <- sample_matrix[i, 10]*sample_matrix[i, 6]
    mat[1, 7] <- sample_matrix[i, 11]*sample_matrix[i, 7]
    mat[1, 8] <- sample_matrix[i, 12]*sample_matrix[i, 8]
    
    # find eigenvalues of matrix
    out <- eigen(mat)
    
    # store the dominant eigenvalue
    output[i] <- out$values[1]
  }
  
  # return the vector of dominant eigenvalues
  return(output)
}






## compute the first order  and total order sobol indices
sobol_indices_tot <- soboljansen(model = eigen_val,
                          X1 = sample1,
                          X2 = sample2)

print(sobol_indices_tot)


## Make a histogram of computed dominant eigenvalues
vect_eigen <- eigen_val(sample1)
vect_eigen <- c(vect_eigen, eigen_val(sample2))

hist(Re(vect_eigen), main = "Histogram of Dominant Eigenvalue", xlab = "lambda")

hist(Re(vect_eigen), breaks = c(0.8, 1, 1.2), main = "Histogram of Dominant Eigenvalue", xlab = "lambda")


