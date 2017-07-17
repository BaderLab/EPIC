# This code implements R code for Bayesian correlation computations for sequence count
# data, as described in the paper "Bayesian Correlation Analysis for Sequence Count Data"
# To use, simply 'source' the file, and then run whichever functions you like on the data.
# The first function allows you to specify your own condition-specific priors, while the
# subsequent functions implement the three priors described in the paper.
# Questions or comments: rpara26@gmail.com ; theodore.j.perkins@gmail.com


# The following is an R function to compute the Bayesian correlation coefficients 
# for a given data matrix of counts, X.
# The rows represent entities and the columns represent conditions
# INPUTs: alpha0 and beta0 are row vectors giving condition-specific priors
# OUTPUTs: Bayesian correlation square matrix

#######################
# Defining FUNCTION to compute the Bayesian correlation coefficients
Bayes_Corr <- function(alpha0, beta0, X){
        nrowsX <- nrow(X)
        k <- ncol(X)
	cs <- colSums(X)
        alphas <- matrix(rep(alpha0,nrowsX), nrow=nrowsX, byrow=TRUE) + X
        betas  <- matrix(rep(beta0,nrowsX), nrow=nrowsX, byrow=TRUE) + matrix(rep(cs,nrowsX), nrow=nrowsX, byrow=TRUE) - X
        alphasPLUSbetas <- alphas + betas
        
        # First BIG product term for covariance formula
        Psi <- alphas/alphasPLUSbetas - matrix(rep(rowSums(alphas/alphasPLUSbetas)/k, k), ncol=k, byrow=FALSE) 
        
        # Covariance matrix
        cov_mtrx <- Psi %*% t(Psi) / k
        
        # Variances (this is a column vector of length = nrowsX)
        var_vec <- as.matrix( ( rowSums( (alphas*betas)/( (alphasPLUSbetas^2)*(alphasPLUSbetas+1) ) ) + rowSums(Psi^2) )/k )
        
        Bcorrvals <- cov_mtrx / sqrt( var_vec %*% t(var_vec) )
	diag(Bcorrvals) <- 1
        return(Bcorrvals)
}
# End of FUNCTION
#######################

#######################
# Computing the Bayesian correlations assuming first (uniform) prior
Bayes_Corr_Prior1 <- function(X){
	d <- dim(X)
	alpha0 <- rep(1,d[2])
	beta0 <- rep(1,d[2])
	Bcorrvals <- Bayes_Corr(alpha0,beta0,X)
	return(Bcorrvals)
}
# End of FUNCTION
#######################

#######################
# Computing the Bayesian correlations assuming second (Dirichlet-marginalized) prior
Bayes_Corr_Prior2 <- function(X){
	d <- dim(X)
	alpha0 <- rep(1/d[1],d[2])
	beta0 <- rep(1-1/d[1],d[2])
	Bcorrvals <- Bayes_Corr(alpha0,beta0,X)
	return(Bcorrvals)
}
# End of FUNCTION
#######################

#######################
# Computing the Bayesian correlations assuming third (zero count-motivated) prior
Bayes_Corr_Prior3 <- function(X){
	d <- dim(X)
	cs <- colSums(X)
	alpha0 <- (cs+1)/(max(cs)+1)
	beta0 <- rep(1,d[2])
	Bcorrvals <- Bayes_Corr(alpha0,beta0,X)
	return(Bcorrvals)
}
# End of FUNCTION
#######################

