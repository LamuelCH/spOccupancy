PGOcc_LC <- function(occ.formula, det.formula, data, inits, priors, n.samples,
                  n.omp.threads = 1, verbose = TRUE, n.report = 100,
                  n.burn = round(.10 * n.samples), n.thin = 1, n.chains = 1,
                  k.fold, k.fold.threads = 1, k.fold.seed, k.fold.only = FALSE,
                  folds) {
  
  ptm <- proc.time()
  
  # Make it look nice
  if (verbose) {
    cat("----------------------------------------\n")
    cat("\tPreparing to run the model\n")
    cat("----------------------------------------\n")
  }
  
  # Functions ---------------------------------------------------------------
  logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}
  rigamma <- function(n, a, b){
    1/rgamma(n = n, shape = a, rate = b)
  }
  
  # Check for unused arguments ------------------------------------------
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
    if(! i %in% formal.args)
      warning("'",i, "' is not an argument")
  }
  
  # Check function arguments ------------------------------------------------
  # Call
  cl <- match.call()
  
  # Data ------------------------------------------------------------------
  if (missing(data)) {
    stop("error: data must be specified")
  }
  if (!is.list(data)) {
    stop("error: data must be a list")
  }
  names(data) <- tolower(names(data))
  if (missing(occ.formula)) {
    stop("error: occ.formula must be specified")
  }
  if (missing(det.formula)) {
    stop("error: det.formula must be specified")
  }
  if (!'y' %in% names(data)) {
    stop("error: detection-nondetection data y must be specified in data")
  }
  y <- data$y
  if (length(dim(y)) != 2) {
    stop("error: detection-nondetection data y must be a two-dimensional array with dimensions corresponding to sites and replicates.")
  }
  J <- nrow(y)
  n.rep <- ncol(y)
  if (!'occ.covs' %in% names(data)) {
    if (occ.formula == ~ 1) {
      if (verbose) {
        message("occupancy covariates (occ.covs) not specified in data.\nAssuming intercept only occupancy model.\n")
      }
      data$occ.covs <- matrix(1, J, 1)
    } else {
      stop("error: occ.covs must be specified in data for an occupancy model with covariates")
    }
  }
  if (!is.matrix(data$occ.covs) & !is.data.frame(data$occ.covs)) {
    stop("error: occ.covs must be a matrix or data frame")
  }
  if (is.matrix(data$occ.covs)) {
    data$occ.covs <- data.frame(data$occ.covs)
  }
  if (!'det.covs' %in% names(data)) {
    if (det.formula == ~ 1) {
      if (verbose) {
        message("detection covariates (det.covs) not specified in data.\nAssuming interept only detection model.\n")
      }
      data$det.covs <- list(int = matrix(1, J, n.rep))
    } else {
      stop("error: det.covs must be specified in data for a detection model with covariates")
    }
  }
  if (!is.list(data$det.covs)) {
    stop("error: det.covs must be a list of matrices, data frames, and/or vectors")
  }
  if (missing(n.samples)) {
    stop("error: n.samples must be specified")
  }
  if (n.burn > n.samples) {
    stop("n.burn must be less than n.samples")
  }
  if (n.thin > n.samples) {
    stop("n.thin must be less than n.samples")
  }
  # Check if n.burn and n.thin are divisible by n.samples
  if (((n.samples - n.burn) / n.thin) %% 1 != 0) {
    stop("the number of posterior samples to save ((n.samples - n.burn) / n.thin) is not a whole number.")
  }
  if (!missing(k.fold)) {
    if (!is.numeric(k.fold) | k.fold < 2) {
      stop("k.fold must be an integer greater than 1")
    }
  }
  
  # Formula -----------------------------------------------------------------
  # Occupancy
  X <- model.matrix(occ.formula, data$occ.covs)
  p <- ncol(X)
  # Detection
  if (length(data$det.covs) > 1) {
    for (i in 1:length(data$det.covs)) {
      if (is.matrix(data$det.covs[[i]])) {
        data$det.covs[[i]] <- data.frame(data$det.covs[[i]])
      }
    }
    X.p <- model.matrix(det.formula, data.frame(lapply(data$det.covs, function(a) unlist(c(a)))))
  } else {
    X.p <- model.matrix(det.formula, data$det.covs)
  }
  p.det <- ncol(X.p)
  
  # Priors ------------------------------------------------------------------
  if (missing(priors)) {
    priors <- list()
  }
  names(priors) <- tolower(names(priors))
  # beta
  if ("beta.normal" %in% names(priors)) {
    if (!is.list(priors$beta.normal) | length(priors$beta.normal) != 2) {
      stop("error: beta.normal must be a list of length 2")
    }
    mu.beta <- priors$beta.normal[[1]]
    sigma.beta <- priors$beta.normal[[2]]
    if (length(mu.beta) != p & length(mu.beta) != 1) {
      if (p == 1) {
        stop(paste("error: beta.normal[[1]] must be a vector of length ", p, sep = ""))
      } else {
        stop(paste("error: beta.normal[[1]] must be a vector of length ", p, " or 1", sep = ""))
      }
    }
    if (length(sigma.beta) != p & length(sigma.beta) != 1) {
      if (p == 1) {
        stop(paste("error: beta.normal[[2]] must be a vector of length ", p, sep = ""))
      } else {
        stop(paste("error: beta.normal[[2]] must be a vector of length ", p, " or 1", sep = ""))
      }
    }
    if (length(mu.beta) != p) {
      mu.beta <- rep(mu.beta, p)
    }
    if (length(sigma.beta) != p) {
      sigma.beta <- rep(sigma.beta, p)
    }
  } else {
    if (verbose) {
      message("No prior specified for beta.normal.\nSetting prior mean to 0 and prior variance to 2.72\n")
    }
    mu.beta <- rep(0, p)
    sigma.beta <- rep(2.72, p)
  }
  # alpha
  if ("alpha.normal" %in% names(priors)) {
    if (!is.list(priors$alpha.normal) | length(priors$alpha.normal) != 2) {
      stop("error: alpha.normal must be a list of length 2")
    }
    mu.alpha <- priors$alpha.normal[[1]]
    sigma.alpha <- priors$alpha.normal[[2]]
    if (length(mu.alpha) != p.det & length(mu.alpha) != 1) {
      if (p.det == 1) {
        stop(paste("error: alpha.normal[[1]] must be a vector of length ", p.det, sep = ""))
      } else {
        stop(paste("error: alpha.normal[[1]] must be a vector of length ", p.det, " or 1", sep = ""))
      }
    }
    if (length(sigma.alpha) != p.det & length(sigma.alpha) != 1) {
      if (p.det == 1) {
        stop(paste("error: alpha.normal[[2]] must be a vector of length ", p.det, sep = ""))
      } else {
        stop(paste("error: alpha.normal[[2]] must be a vector of length ", p.det, " or 1", sep = ""))
      }
    }
    if (length(mu.alpha) != p.det) {
      mu.alpha <- rep(mu.alpha, p.det)
    }
    if (length(sigma.alpha) != p.det) {
      sigma.alpha <- rep(sigma.alpha, p.det)
    }
  } else {
    if (verbose) {
      message("No prior specified for alpha.normal.\nSetting prior mean to 0 and prior variance to 2.72\n")
    }
    mu.alpha <- rep(0, p.det)
    sigma.alpha <- rep(2.72, p.det)
  }
  
  # Starting values -------------------------------------------------------
  if (missing(inits)) {
    inits <- list()
  }
  names(inits) <- tolower(names(inits))
  # z
  if ("z" %in% names(inits)) {
    z.inits <- inits$z
    if (length(z.inits) != J) {
      stop(paste("z must be a vector of length ", J, " with initial values for the latent occupancy state", sep = ""))
    }
  } else {
    z.inits <- apply(y, 1, max, na.rm = TRUE)
    z.inits[z.inits == '-Inf'] <- 0
    if (verbose) {
      message("z is not specified in initial values.\nSetting initial values based on observed data\n")
    }
  }
  # beta
  if ("beta" %in% names(inits)) {
    beta.inits <- inits[["beta"]]
    if (length(beta.inits) != p) {
      stop(paste("error: initial values for beta must be of length ", p, sep = ""))
    }
  } else {
    beta.inits <- rnorm(p, mu.beta, 0.5)
    if (verbose) {
      message("beta is not specified in initial values.\nSetting initial values to random values from the prior\n")
    }
  }
  # alpha
  if ("alpha" %in% names(inits)) {
    alpha.inits <- inits[["alpha"]]
    if (length(alpha.inits) != p.det) {
      stop(paste("error: initial values for alpha must be of length ", p.det, sep = ""))
    }
  } else {
    alpha.inits <- rnorm(p.det, mu.alpha, 0.5)
    if (verbose) {
      message("alpha is not specified in initial values.\nSetting initial values to random values from the prior\n")
    }
  }
  
  # Prep for model fitting ------------------------------------------------
  # Reformat detection covariates
  X.p.reordered <- array(NA, dim = c(J, n.rep, p.det))
  for (j in 1:J) {
    X.p.reordered[j, , ] <- X.p[(j - 1) * n.rep + 1:n.rep, ]
  }
  
  # Get detection covariates in list form
  X.p.list <- list()
  for (i in 1:p.det) {
    X.p.list[[i]] <- X.p.reordered[, , i]
  }
  
  # y is ordered by site, survey
  y.long <- c(y)
  # Occupancy covariates
  X.long <- X
  
  # Detection covariates
  # Covariates are site by survey
  X.p.long <- do.call(rbind, X.p.list)
  
  # Used to subset the data for the BVS portion of the model.
  beta.zero.indx <- which(colSums(abs(X)) == J)
  alpha.zero.indx <- which(colSums(abs(X.p)) == J * n.rep)
  
  # Number of pseudoreplicates for fitting the model.
  n.obs <- sum(!is.na(y))
  
  # Number of occupancy parameters
  n.occ.params <- p
  # Number of detection parameters
  n.det.params <- p.det
  
  # Run model ---------------------------------------------------------------
  # Grab a slice of the arrays for the first chain
  z.inits.i <- z.inits
  beta.inits.i <- beta.inits
  alpha.inits.i <- alpha.inits
  
  # Create a list of initial values for each chain
  inits.list <- list()
  for (i in 1:n.chains) {
    inits.list[[i]] <- list(z = z.inits.i,
                            beta = beta.inits.i,
                            alpha = alpha.inits.i)
  }
  
  # If not running cross-validation
  if (missing(k.fold)) {
    # Fit model
    out <- PGOccModel(y, X, X.p, p, p.det, J, n.rep,
                      n.obs, n.occ.params, n.det.params,
                      mu.beta, sigma.beta, mu.alpha, sigma.alpha,
                      n.samples, n.omp.threads, verbose, n.report,
                      n.burn, n.thin, n.chains, inits.list)
    
    # Return results --------------------------------------------------------
    # Convert lists to matrices
    out$beta.samples <- do.call(rbind, out$beta.samples)
    out$alpha.samples <- do.call(rbind, out$alpha.samples)
    out$z.samples <- do.call(rbind, out$z.samples)
    out$psi.samples <- do.call(rbind, out$psi.samples)
    # Get rid of initial value lists
    out$z.inits <- NULL
    out$beta.inits <- NULL
    out$alpha.inits <- NULL
    
    # Create MCMC list
    out$beta.samples <- mcmc(out$beta.samples)
    out$alpha.samples <- mcmc(out$alpha.samples)
    out$z.samples <- mcmc(out$z.samples)
    out$psi.samples <- mcmc(out$psi.samples)
    # Give parameter names
    # Occupancy covariates
    if (is.null(colnames(X))) {
      occ.cov.names <- paste('occ.cov', 1:p, sep = '')
    } else {
      occ.cov.names <- colnames(X)
    }
    # Detection covariates
    if (is.null(colnames(X.p))) {
      det.cov.names <- paste('det.cov', 1:p.det, sep = '')
    } else {
      det.cov.names <- colnames(X.p)
    }
    colnames(out$beta.samples) <- occ.cov.names
    colnames(out$alpha.samples) <- det.cov.names
    colnames(out$z.samples) <- paste('z[', 1:J, ']', sep = '')
    colnames(out$psi.samples) <- paste('psi[', 1:J, ']', sep = '')
    
    # Add other things to the list
    out$n.samples <- n.samples
    out$n.burn <- n.burn
    out$n.thin <- n.thin
    out$n.chains <- n.chains
    out$p <- p
    out$p.det <- p.det
    out$J <- J
    out$n.rep <- n.rep
    out$occ.formula <- occ.formula
    out$det.formula <- det.formula
    out$call <- cl
    out$X <- X
    out$X.p <- X.p
    out$y <- y
    class(out) <- "PGOcc"
  } else {
    # Run k-fold cross-validation -----------------------------------------
    if (verbose) {
      message(paste("Performing ", k.fold, "-fold cross-validation with ", k.fold.threads,
                    " thread(s).\n", sep = ''))
    }
    if (k.fold.only) {
      if (verbose) {
        message("k.fold.only = TRUE. Skipping full model fit.\n")
      }
    }
    # Set seed for all threads
    if (!missing(k.fold.seed)) {
      set.seed(k.fold.seed)
    }
    
    ## MODIFICATION START
    # Custom folds
    if (missing(folds)) {
      fold.indx <- sample(rep(1:k.fold, length.out = J))
    } else {
      fold.indx <- folds
      k.fold <- max(folds)
    }
    ## MODIFICATION END
    
    # Object to hold of deviance values
    dev.vals <- rep(NA, k.fold)
    
    ## MODIFICATION START
    # Additional evaluation metrics
    auc.vals <- rep(NA, k.fold)
    rmse.vals <- rep(NA, k.fold)
    tjur.vals <- rep(NA, k.fold)
    ## MODIFICATION END
    
    # Set up threads for parallel programming
    # Show progress
    pb <- txtProgressBar(min = 0, max = k.fold, style = 3)
    progress <- function(n) {
      setTxtProgressBar(pb, n)
    }
    opts <- list(progress = progress)
    # Set up the cluster
    cl.i <- makeCluster(k.fold.threads, outfile = "")
    registerDoParallel(cl.i)
    # Run the model in parallel
    out.list <- foreach(i = 1:k.fold, .packages = c('spOccupancy', 'coda'),
                        .options.snow = opts) %dopar% {
                          # Create data objects for held out fold
                          # Held out sites
                          sites.hold <- which(fold.indx == i)
                          # Observed sites
                          sites.fit <- which(fold.indx != i)
                          # Number of held out sites
                          J.hold <- length(sites.hold)
                          # Number of observed sites
                          J.fit <- length(sites.fit)
                          # Data
                          y.fit <- y[sites.fit, , drop = FALSE]
                          y.hold <- y[sites.hold, , drop = FALSE]
                          # Occupancy covariates
                          X.fit <- X[sites.fit, , drop = FALSE]
                          X.hold <- X[sites.hold, , drop = FALSE]
                          # Detection covariates
                          X.p.fit <- X.p[do.call(c, lapply(sites.fit, function(a) (a-1)*n.rep + 1:n.rep)), , drop = FALSE]
                          X.p.hold <- X.p[do.call(c, lapply(sites.hold, function(a) (a-1)*n.rep + 1:n.rep)), , drop = FALSE]
                          # Initial values
                          # Grab a slice of the arrays for the first chain
                          z.inits.i <- z.inits[sites.fit]
                          beta.inits.i <- beta.inits
                          alpha.inits.i <- alpha.inits
                          # Need to make sure initial values for z are correct.
                          z.inits.i[z.inits.i == -Inf] <- 0
                          # Create a list of initial values for each chain
                          inits.list.i <- list()
                          for (j in 1:n.chains) {
                            inits.list.i[[j]] <- list(z = z.inits.i,
                                                      beta = beta.inits.i,
                                                      alpha = alpha.inits.i)
                          }
                          # Reformat detection covariates
                          X.p.reordered <- array(NA, dim = c(J.fit, n.rep, p.det))
                          for (j in 1:J.fit) {
                            X.p.reordered[j, , ] <- X.p.fit[(j - 1) * n.rep + 1:n.rep, ]
                          }
                          # Get detection covariates in list form
                          X.p.list <- list()
                          for (j in 1:p.det) {
                            X.p.list[[j]] <- X.p.reordered[, , j]
                          }
                          # Fit the model
                          out.i <- PGOccModel(y.fit, X.fit, X.p.fit, p, p.det, J.fit, n.rep,
                                              n.obs, n.occ.params, n.det.params,
                                              mu.beta, sigma.beta, mu.alpha, sigma.alpha,
                                              n.samples, n.omp.threads, verbose = FALSE, n.report = n.report,
                                              n.burn, n.thin, n.chains, inits.list.i)
                          # Convert lists to matrices
                          out.i$beta.samples <- do.call(rbind, out.i$beta.samples)
                          out.i$alpha.samples <- do.call(rbind, out.i$alpha.samples)
                          n.post.samples <- nrow(out.i$beta.samples)
                          
                          # Predict occurrence at held out locations
                          psi.vals <- matrix(NA, J.hold, n.post.samples)
                          for (j in 1:n.post.samples) {
                            psi.vals[, j] <- logit.inv(X.hold %*% out.i$beta.samples[j, ])
                          }
                          psi.means <- apply(psi.vals, 1, mean)
                          z.obs.hold <- apply(y.hold, 1, max, na.rm = TRUE)
                          # Replace -Inf with 0
                          z.obs.hold[z.obs.hold == -Inf] <- 0
                          dev.i <- -2 * sum(dbinom(z.obs.hold, 1, psi.means, log = TRUE), na.rm = TRUE)
                          
                          ## MODIFICATION START
                          # Calculate AUC, RMSE, and Tjur's R2
                          pred <- prediction(psi.means, z.obs.hold)
                          auc.i <- performance(pred, "auc")@y.values[[1]]
                          rmse.i <- sqrt(mean((z.obs.hold - psi.means)^2))
                          tjur.i <- mean(psi.means[which(z.obs.hold == 1)]) - mean(psi.means[which(z.obs.hold == 0)])
                          ## MODIFICATION END
                          
                          list(dev = dev.i,
                               auc = auc.i,
                               rmse = rmse.i,
                               tjur = tjur.i,
                               psi.means = psi.means,
                               z.obs.hold = z.obs.hold)
                        } # i
    # Stop the cluster
    stopCluster(cl.i)
    
    # Separate out the results
    dev.vals <- sapply(out.list, function(a) a$dev)
    auc.vals <- sapply(out.list, function(a) a$auc)
    rmse.vals <- sapply(out.list, function(a) a$rmse)
    tjur.vals <- sapply(out.list, function(a) a$tjur)
    
    if (k.fold.only) {
      out <- list(k.fold.deviance = dev.vals,
                  k.fold.auc = auc.vals,
                  k.fold.rmse = rmse.vals,
                  k.fold.tjur = tjur.vals,
                  call = cl)
      class(out) <- "PGOcc"
    } else {
      # Fit model to the full data set
      if (verbose) {
        message("Finished cross-validation. Fitting model to full data set.\n")
      }
      # Fit model
      out <- PGOccModel(y, X, X.p, p, p.det, J, n.rep,
                        n.obs, n.occ.params, n.det.params,
                        mu.beta, sigma.beta, mu.alpha, sigma.alpha,
                        n.samples, n.omp.threads, verbose, n.report,
                        n.burn, n.thin, n.chains, inits.list)
      # Convert lists to matrices
      out$beta.samples <- do.call(rbind, out$beta.samples)
      out$alpha.samples <- do.call(rbind, out$alpha.samples)
      out$z.samples <- do.call(rbind, out$z.samples)
      out$psi.samples <- do.call(rbind, out$psi.samples)
      # Get rid of initial value lists
      out$z.inits <- NULL
      out$beta.inits <- NULL
      out$alpha.inits <- NULL
      
      # Create MCMC list
      out$beta.samples <- mcmc(out$beta.samples)
      out$alpha.samples <- mcmc(out$alpha.samples)
      out$z.samples <- mcmc(out$z.samples)
      out$psi.samples <- mcmc(out$psi.samples)
      # Give parameter names
      # Occupancy covariates
      if (is.null(colnames(X))) {
        occ.cov.names <- paste('occ.cov', 1:p, sep = '')
      } else {
        occ.cov.names <- colnames(X)
      }
      # Detection covariates
      if (is.null(colnames(X.p))) {
        det.cov.names <- paste('det.cov', 1:p.det, sep = '')
      } else {
        det.cov.names <- colnames(X.p)
      }
      colnames(out$beta.samples) <- occ.cov.names
      colnames(out$alpha.samples) <- det.cov.names
      colnames(out$z.samples) <- paste('z[', 1:J, ']', sep = '')
      colnames(out$psi.samples) <- paste('psi[', 1:J, ']', sep = '')
      
      # Add other things to the list
      out$n.samples <- n.samples
      out$n.burn <- n.burn
      out$n.thin <- n.thin
      out$n.chains <- n.chains
      out$p <- p
      out$p.det <- p.det
      out$J <- J
      out$n.rep <- n.rep
      out$occ.formula <- occ.formula
      out$det.formula <- det.formula
      out$call <- cl
      out$X <- X
      out$X.p <- X.p
      out$y <- y
      out$k.fold.deviance <- dev.vals
      out$k.fold.auc <- auc.vals
      out$k.fold.rmse <- rmse.vals
      out$k.fold.tjur <- tjur.vals
      class(out) <- "PGOcc"
    }
  }
  out$run.time <- proc.time() - ptm
  out
}