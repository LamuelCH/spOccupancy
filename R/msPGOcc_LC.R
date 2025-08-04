msPGOcc <- function(occ.formula, det.formula, data, inits, priors, n.samples,
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
  if (length(dim(y)) != 3) {
    stop("error: detection-nondetection data y must be a three-dimensional array with dimensions corresponding to species, sites, and replicates.")
  }
  N <- dim(y)[1]
  J <- dim(y)[2]
  n.rep <- dim(y)[3]
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
  # beta.comm
  if ("beta.comm.normal" %in% names(priors)) {
    if (!is.list(priors$beta.comm.normal) | length(priors$beta.comm.normal) != 2) {
      stop("error: beta.comm.normal must be a list of length 2")
    }
    mu.beta.comm <- priors$beta.comm.normal[[1]]
    sigma.beta.comm <- priors$beta.comm.normal[[2]]
    if (length(mu.beta.comm) != p & length(mu.beta.comm) != 1) {
      if (p == 1) {
        stop(paste("error: beta.comm.normal[[1]] must be a vector of length ", p, sep = ""))
      } else {
        stop(paste("error: beta.comm.normal[[1]] must be a vector of length ", p, " or 1", sep = ""))
      }
    }
    if (length(sigma.beta.comm) != p & length(sigma.beta.comm) != 1) {
      if (p == 1) {
        stop(paste("error: beta.comm.normal[[2]] must be a vector of length ", p, sep = ""))
      } else {
        stop(paste("error: beta.comm.normal[[2]] must be a vector of length ", p, " or 1", sep = ""))
      }
    }
    if (length(mu.beta.comm) != p) {
      mu.beta.comm <- rep(mu.beta.comm, p)
    }
    if (length(sigma.beta.comm) != p) {
      sigma.beta.comm <- rep(sigma.beta.comm, p)
    }
  } else {
    if (verbose) {
      message("No prior specified for beta.comm.normal.\nSetting prior mean to 0 and prior variance to 2.72\n")
    }
    mu.beta.comm <- rep(0, p)
    sigma.beta.comm <- rep(2.72, p)
  }
  # alpha.comm
  if ("alpha.comm.normal" %in% names(priors)) {
    if (!is.list(priors$alpha.comm.normal) | length(priors$alpha.comm.normal) != 2) {
      stop("error: alpha.comm.normal must be a list of length 2")
    }
    mu.alpha.comm <- priors$alpha.comm.normal[[1]]
    sigma.alpha.comm <- priors$alpha.comm.normal[[2]]
    if (length(mu.alpha.comm) != p.det & length(mu.alpha.comm) != 1) {
      if (p.det == 1) {
        stop(paste("error: alpha.comm.normal[[1]] must be a vector of length ", p.det, sep = ""))
      } else {
        stop(paste("error: alpha.comm.normal[[1]] must be a vector of length ", p.det, " or 1", sep = ""))
      }
    }
    if (length(sigma.alpha.comm) != p.det & length(sigma.alpha.comm) != 1) {
      if (p.det == 1) {
        stop(paste("error: alpha.comm.normal[[2]] must be a vector of length ", p.det, sep = ""))
      } else {
        stop(paste("error: alpha.comm.normal[[2]] must be a vector of length ", p.det, " or 1", sep = ""))
      }
    }
    if (length(mu.alpha.comm) != p.det) {
      mu.alpha.comm <- rep(mu.alpha.comm, p.det)
    }
    if (length(sigma.alpha.comm) != p.det) {
      sigma.alpha.comm <- rep(sigma.alpha.comm, p.det)
    }
  } else {
    if (verbose) {
      message("No prior specified for alpha.comm.normal.\nSetting prior mean to 0 and prior variance to 2.72\n")
    }
    mu.alpha.comm <- rep(0, p.det)
    sigma.alpha.comm <- rep(2.72, p.det)
  }
  # tau.sq.beta
  if ("tau.sq.beta.ig" %in% names(priors)) {
    if (!is.list(priors$tau.sq.beta.ig) | length(priors$tau.sq.beta.ig) != 2) {
      stop("error: tau.sq.beta.ig must be a list of length 2")
    }
    tau.sq.beta.a <- priors$tau.sq.beta.ig[[1]]
    tau.sq.beta.b <- priors$tau.sq.beta.ig[[2]]
    if (length(tau.sq.beta.a) != p & length(tau.sq.beta.a) != 1) {
      if (p == 1) {
        stop(paste("error: tau.sq.beta.ig[[1]] must be a vector of length ", p, sep = ""))
      } else {
        stop(paste("error: tau.sq.beta.ig[[1]] must be a vector of length ", p, " or 1", sep = ""))
      }
    }
    if (length(tau.sq.beta.b) != p & length(tau.sq.beta.b) != 1) {
      if (p == 1) {
        stop(paste("error: tau.sq.beta.ig[[2]] must be a vector of length ", p, sep = ""))
      } else {
        stop(paste("error: tau.sq.beta.ig[[2]] must be a vector of length ", p, " or 1", sep = ""))
      }
    }
    if (length(tau.sq.beta.a) != p) {
      tau.sq.beta.a <- rep(tau.sq.beta.a, p)
    }
    if (length(tau.sq.beta.b) != p) {
      tau.sq.beta.b <- rep(tau.sq.beta.b, p)
    }
  } else {
    if (verbose) {
      message("No prior specified for tau.sq.beta.ig.\nSetting prior shape to 0.1 and prior scale to 0.1\n")
    }
    tau.sq.beta.a <- rep(0.1, p)
    tau.sq.beta.b <- rep(0.1, p)
  }
  # tau.sq.alpha
  if ("tau.sq.alpha.ig" %in% names(priors)) {
    if (!is.list(priors$tau.sq.alpha.ig) | length(priors$tau.sq.alpha.ig) != 2) {
      stop("error: tau.sq.alpha.ig must be a list of length 2")
    }
    tau.sq.alpha.a <- priors$tau.sq.alpha.ig[[1]]
    tau.sq.alpha.b <- priors$tau.sq.alpha.ig[[2]]
    if (length(tau.sq.alpha.a) != p.det & length(tau.sq.alpha.a) != 1) {
      if (p.det == 1) {
        stop(paste("error: tau.sq.alpha.ig[[1]] must be a vector of length ", p.det, sep = ""))
      } else {
        stop(paste("error: tau.sq.alpha.ig[[1]] must be a vector of length ", p.det, " or 1", sep = ""))
      }
    }
    if (length(tau.sq.alpha.b) != p.det & length(tau.sq.alpha.b) != 1) {
      if (p.det == 1) {
        stop(paste("error: tau.sq.alpha.ig[[2]] must be a vector of length ", p.det, sep = ""))
      } else {
        stop(paste("error: tau.sq.alpha.ig[[2]] must be a vector of length ", p.det, " or 1", sep = ""))
      }
    }
    if (length(tau.sq.alpha.a) != p.det) {
      tau.sq.alpha.a <- rep(tau.sq.alpha.a, p.det)
    }
    if (length(tau.sq.alpha.b) != p.det) {
      tau.sq.alpha.b <- rep(tau.sq.alpha.b, p.det)
    }
  } else {
    if (verbose) {
      message("No prior specified for tau.sq.alpha.ig.\nSetting prior shape to 0.1 and prior scale to 0.1\n")
    }
    tau.sq.alpha.a <- rep(0.1, p.det)
    tau.sq.alpha.b <- rep(0.1, p.det)
  }
  
  # Starting values -------------------------------------------------------
  if (missing(inits)) {
    inits <- list()
  }
  names(inits) <- tolower(names(inits))
  # z
  if ("z" %in% names(inits)) {
    z.inits <- inits$z
    if (length(dim(z.inits)) != 2 | dim(z.inits)[1] != N | dim(z.inits)[2] != J) {
      stop(paste("z must be a ", N, " x ", J, " matrix of initial values", sep = ""))
    }
  } else {
    z.inits <- apply(y, c(1, 2), max, na.rm = TRUE)
    z.inits[z.inits == '-Inf'] <- 0
    if (verbose) {
      message("z is not specified in initial values.\nSetting initial values based on observed data\n")
    }
  }
  # beta.comm
  if ("beta.comm" %in% names(inits)) {
    beta.comm.inits <- inits[["beta.comm"]]
    if (length(beta.comm.inits) != p & length(beta.comm.inits) != 1) {
      if (p == 1) {
        stop(paste("error: initial values for beta.comm must be of length ", p, sep = ""))
      } else {
        stop(paste("error: initial values for beta.comm must be of length ", p, " or 1", sep = ""))
      }
    }
    if (length(beta.comm.inits) != p) {
      beta.comm.inits <- rep(beta.comm.inits, p)
    }
  } else {
    beta.comm.inits <- rnorm(p, mu.beta.comm, 0.5)
    if (verbose) {
      message("beta.comm is not specified in initial values.\nSetting initial values to random values from the prior\n")
    }
  }
  # alpha.comm
  if ("alpha.comm" %in% names(inits)) {
    alpha.comm.inits <- inits[["alpha.comm"]]
    if (length(alpha.comm.inits) != p.det & length(alpha.comm.inits) != 1) {
      if (p.det == 1) {
        stop(paste("error: initial values for alpha.comm must be of length ", p.det, sep = ""))
      } else {
        stop(paste("error: initial values for alpha.comm must be of length ", p.det, " or 1", sep = ""))
      }
    }
    if (length(alpha.comm.inits) != p.det) {
      alpha.comm.inits <- rep(alpha.comm.inits, p.det)
    }
  } else {
    alpha.comm.inits <- rnorm(p.det, mu.alpha.comm, 0.5)
    if (verbose) {
      message("alpha.comm is not specified in initial values.\nSetting initial values to random values from the prior\n")
    }
  }
  # tau.sq.beta
  if ("tau.sq.beta" %in% names(inits)) {
    tau.sq.beta.inits <- inits[["tau.sq.beta"]]
    if (length(tau.sq.beta.inits) != p & length(tau.sq.beta.inits) != 1) {
      if (p == 1) {
        stop(paste("error: initial values for tau.sq.beta must be of length ", p, sep = ""))
      } else {
        stop(paste("error: initial values for tau.sq.beta must be of length ", p, " or 1", sep = ""))
      }
    }
    if (length(tau.sq.beta.inits) != p) {
      tau.sq.beta.inits <- rep(tau.sq.beta.inits, p)
    }
  } else {
    tau.sq.beta.inits <- runif(p, 0.05, 2)
    if (verbose) {
      message("tau.sq.beta is not specified in initial values.\nSetting initial values to random values between 0.05 and 2\n")
    }
  }
  # tau.sq.alpha
  if ("tau.sq.alpha" %in% names(inits)) {
    tau.sq.alpha.inits <- inits[["tau.sq.alpha"]]
    if (length(tau.sq.alpha.inits) != p.det & length(tau.sq.alpha.inits) != 1) {
      if (p.det == 1) {
        stop(paste("error: initial values for tau.sq.alpha must be of length ", p.det, sep = ""))
      } else {
        stop(paste("error: initial values for tau.sq.alpha must be of length ", p.det, " or 1", sep = ""))
      }
    }
    if (length(tau.sq.alpha.inits) != p.det) {
      tau.sq.alpha.inits <- rep(tau.sq.alpha.inits, p.det)
    }
  } else {
    tau.sq.alpha.inits <- runif(p.det, 0.05, 2)
    if (verbose) {
      message("tau.sq.alpha is not specified in initial values.\nSetting initial values to random values between 0.05 and 2\n")
    }
  }
  # beta
  if ("beta" %in% names(inits)) {
    beta.inits <- inits[["beta"]]
    if (is.matrix(beta.inits)) {
      if (nrow(beta.inits) != N | ncol(beta.inits) != p) {
        stop(paste("error: initial values for beta must be a matrix with dimensions ", N, "x", p, " or a single vector of length ", p, sep = ""))
      }
    } else {
      if (length(beta.inits) != p) {
        stop(paste("error: initial values for beta must be a matrix with dimensions ", N, "x", p, " or a single vector of length ", p, sep = ""))
      }
      beta.inits <- matrix(beta.inits, N, p, byrow = TRUE)
    }
  } else {
    beta.inits <- matrix(rnorm(N * p, beta.comm.inits, sqrt(tau.sq.beta.inits)), N, p, byrow = TRUE)
    if (verbose) {
      message("beta is not specified in initial values.\nSetting initial values to random values from the community-level normal distribution\n")
    }
  }
  # alpha
  if ("alpha" %in% names(inits)) {
    alpha.inits <- inits[["alpha"]]
    if (is.matrix(alpha.inits)) {
      if (nrow(alpha.inits) != N | ncol(alpha.inits) != p.det) {
        stop(paste("error: initial values for alpha must be a matrix with dimensions ", N, "x", p.det, " or a single vector of length ", p.det, sep = ""))
      }
    } else {
      if (length(alpha.inits) != p.det) {
        stop(paste("error: initial values for alpha must be a matrix with dimensions ", N, "x", p.det, " or a single vector of length ", p.det, sep = ""))
      }
      alpha.inits <- matrix(alpha.inits, N, p.det, byrow = TRUE)
    }
  } else {
    alpha.inits <- matrix(rnorm(N * p.det, alpha.comm.inits, sqrt(tau.sq.alpha.inits)), N, p.det, byrow = TRUE)
    if (verbose) {
      message("alpha is not specified in initial values.\nSetting initial values to random values from the community-level normal distribution\n")
    }
  }
  
  # Get species names -----------------------------------------------------
  if (is.null(dimnames(y))) {
    sp.names <- paste('sp', 1:N, sep = '')
  } else {
    sp.names <- dimnames(y)[[1]]
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
  
  # y is ordered by species, site, survey
  # We want it ordered by species, site, survey
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
  beta.comm.inits.i <- beta.comm.inits
  alpha.comm.inits.i <- alpha.comm.inits
  tau.sq.beta.inits.i <- tau.sq.beta.inits
  tau.sq.alpha.inits.i <- tau.sq.alpha.inits
  
  # Create a list of initial values for each chain
  inits.list <- list()
  for (i in 1:n.chains) {
    inits.list[[i]] <- list(z = z.inits.i,
                            beta = beta.inits.i,
                            alpha = alpha.inits.i,
                            beta.comm = beta.comm.inits.i,
                            alpha.comm = alpha.comm.inits.i,
                            tau.sq.beta = tau.sq.beta.inits.i,
                            tau.sq.alpha = tau.sq.alpha.inits.i)
  }
  
  # If not running cross-validation
  if (missing(k.fold)) {
    # Fit model
    out <- msPGOccModel(y, X, X.p, p, p.det, N, J, n.rep,
                        n.obs, n.occ.params, n.det.params,
                        mu.beta.comm, sigma.beta.comm, mu.alpha.comm, sigma.alpha.comm,
                        tau.sq.beta.a, tau.sq.beta.b, tau.sq.alpha.a, tau.sq.alpha.b,
                        n.samples, n.omp.threads, verbose, n.report,
                        n.burn, n.thin, n.chains, inits.list)
    
    # Return results --------------------------------------------------------
    # Convert lists to matrices
    out$beta.samples <- do.call(rbind, out$beta.samples)
    out$alpha.samples <- do.call(rbind, out$alpha.samples)
    out$beta.comm.samples <- do.call(rbind, out$beta.comm.samples)
    out$alpha.comm.samples <- do.call(rbind, out$alpha.comm.samples)
    out$tau.sq.beta.samples <- do.call(rbind, out$tau.sq.beta.samples)
    out$tau.sq.alpha.samples <- do.call(rbind, out$tau.sq.alpha.samples)
    out$z.samples <- do.call(abind::abind, c(out$z.samples, along = 0))
    # Get rid of initial value lists
    out$z.inits <- NULL
    out$beta.inits <- NULL
    out$alpha.inits <- NULL
    
    # Create MCMC list
    out$beta.samples <- mcmc(out$beta.samples)
    out$alpha.samples <- mcmc(out$alpha.samples)
    out$beta.comm.samples <- mcmc(out$beta.comm.samples)
    out$alpha.comm.samples <- mcmc(out$alpha.comm.samples)
    out$tau.sq.beta.samples <- mcmc(out$tau.sq.beta.samples)
    out$tau.sq.alpha.samples <- mcmc(out$tau.sq.alpha.samples)
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
    colnames(out$beta.comm.samples) <- occ.cov.names
    colnames(out$tau.sq.beta.samples) <- occ.cov.names
    colnames(out$alpha.comm.samples) <- det.cov.names
    colnames(out$tau.sq.alpha.samples) <- det.cov.names
    # Species-specific effects
    beta.names <- paste(rep(sp.names, each = p), occ.cov.names, sep = '-')
    alpha.names <- paste(rep(sp.names, each = p.det), det.cov.names, sep = '-')
    colnames(out$beta.samples) <- beta.names
    colnames(out$alpha.samples) <- alpha.names
    
    # Add other things to the list
    out$n.samples <- n.samples
    out$n.burn <- n.burn
    out$n.thin <- n.thin
    out$n.chains <- n.chains
    out$p <- p
    out$p.det <- p.det
    out$N <- N
    out$J <- J
    out$n.rep <- n.rep
    out$occ.formula <- occ.formula
    out$det.formula <- det.formula
    out$call <- cl
    out$X <- X
    out$X.p <- X.p
    out$y <- y
    class(out) <- "msPGOcc"
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
    dev.vals <- matrix(NA, k.fold, N)
    colnames(dev.vals) <- sp.names
    
    ## MODIFICATION START
    # Additional evaluation metrics
    auc.vals <- matrix(NA, k.fold, N)
    rmsev.vals <- matrix(NA, k.fold, N)
    tjur.vals <- matrix(NA, k.fold, N)
    colnames(auc.vals) <- sp.names
    colnames(rmsev.vals) <- sp.names
    colnames(tjur.vals) <- sp.names
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
                          y.fit <- y[, sites.fit, , drop = FALSE]
                          y.hold <- y[, sites.hold, , drop = FALSE]
                          # Occupancy covariates
                          X.fit <- X[sites.fit, , drop = FALSE]
                          X.hold <- X[sites.hold, , drop = FALSE]
                          # Detection covariates
                          X.p.fit <- X.p[do.call(c, lapply(sites.fit, function(a) (a-1)*n.rep + 1:n.rep)), , drop = FALSE]
                          X.p.hold <- X.p[do.call(c, lapply(sites.hold, function(a) (a-1)*n.rep + 1:n.rep)), , drop = FALSE]
                          # Initial values
                          # Grab a slice of the arrays for the first chain
                          z.inits.i <- z.inits[, sites.fit]
                          beta.inits.i <- beta.inits
                          alpha.inits.i <- alpha.inits
                          beta.comm.inits.i <- beta.comm.inits
                          alpha.comm.inits.i <- alpha.comm.inits
                          tau.sq.beta.inits.i <- tau.sq.beta.inits
                          tau.sq.alpha.inits.i <- tau.sq.alpha.inits
                          # Need to make sure initial values for z are correct.
                          z.inits.i[z.inits.i == -Inf] <- 0
                          # Create a list of initial values for each chain
                          inits.list.i <- list()
                          for (j in 1:n.chains) {
                            inits.list.i[[j]] <- list(z = z.inits.i,
                                                      beta = beta.inits.i,
                                                      alpha = alpha.inits.i,
                                                      beta.comm = beta.comm.inits.i,
                                                      alpha.comm = alpha.comm.inits.i,
                                                      tau.sq.beta = tau.sq.beta.inits.i,
                                                      tau.sq.alpha = tau.sq.alpha.inits.i)
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
                          out.i <- msPGOccModel(y.fit, X.fit, X.p.fit, p, p.det, N, J.fit, n.rep,
                                                n.obs, n.occ.params, n.det.params,
                                                mu.beta.comm, sigma.beta.comm, mu.alpha.comm, sigma.alpha.comm,
                                                tau.sq.beta.a, tau.sq.beta.b, tau.sq.alpha.a, tau.sq.alpha.b,
                                                n.samples, n.omp.threads, verbose = FALSE, n.report = n.report,
                                                n.burn, n.thin, n.chains, inits.list.i)
                          # Convert lists to matrices
                          out.i$beta.samples <- do.call(rbind, out.i$beta.samples)
                          out.i$alpha.samples <- do.call(rbind, out.i$alpha.samples)
                          n.post.samples <- nrow(out.i$beta.samples)
                          
                          # Predict occurrence at held out locations
                          psi.vals <- array(NA, dim = c(N, J.hold, n.post.samples))
                          for (j in 1:n.post.samples) {
                            psi.vals[, , j] <- logit.inv(X.hold %*% t(matrix(out.i$beta.samples[j, ], N, p, byrow = TRUE)))
                          }
                          psi.means <- apply(psi.vals, c(1, 2), mean)
                          z.obs.hold <- apply(y.hold, c(1, 2), max, na.rm = TRUE)
                          # Replace -Inf with 0
                          z.obs.hold[z.obs.hold == -Inf] <- 0
                          dev.i <- -2 * sum(dbinom(z.obs.hold, 1, psi.means, log = TRUE), na.rm = TRUE)
                          
                          ## MODIFICATION START
                          # Calculate AUC, RMSE, and Tjur's R2
                          auc.i <- numeric(N)
                          rmse.i <- numeric(N)
                          tjur.i <- numeric(N)
                          for (k in 1:N) {
                            pred <- prediction(psi.means[k, ], z.obs.hold[k, ])
                            auc.i[k] <- performance(pred, "auc")@y.values[[1]]
                            rmse.i[k] <- sqrt(mean((z.obs.hold[k, ] - psi.means[k, ])^2))
                            tjur.i[k] <- mean(psi.means[k, which(z.obs.hold[k, ] == 1)]) - mean(psi.means[k, which(z.obs.hold[k, ] == 0)])
                          }
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
      class(out) <- "msPGOcc"
    } else {
      # Fit model to the full data set
      if (verbose) {
        message("Finished cross-validation. Fitting model to full data set.\n")
      }
      # Fit model
      out <- msPGOccModel(y, X, X.p, p, p.det, N, J, n.rep,
                          n.obs, n.occ.params, n.det.params,
                          mu.beta.comm, sigma.beta.comm, mu.alpha.comm, sigma.alpha.comm,
                          tau.sq.beta.a, tau.sq.beta.b, tau.sq.alpha.a, tau.sq.alpha.b,
                          n.samples, n.omp.threads, verbose, n.report,
                          n.burn, n.thin, n.chains, inits.list)
      # Convert lists to matrices
      out$beta.samples <- do.call(rbind, out$beta.samples)
      out$alpha.samples <- do.call(rbind, out$alpha.samples)
      out$beta.comm.samples <- do.call(rbind, out$beta.comm.samples)
      out$alpha.comm.samples <- do.call(rbind, out$alpha.comm.samples)
      out$tau.sq.beta.samples <- do.call(rbind, out$tau.sq.beta.samples)
      out$tau.sq.alpha.samples <- do.call(rbind, out$tau.sq.alpha.samples)
      out$z.samples <- do.call(abind::abind, c(out$z.samples, along = 0))
      # Get rid of initial value lists
      out$z.inits <- NULL
      out$beta.inits <- NULL
      out$alpha.inits <- NULL
      
      # Create MCMC list
      out$beta.samples <- mcmc(out$beta.samples)
      out$alpha.samples <- mcmc(out$alpha.samples)
      out$beta.comm.samples <- mcmc(out$beta.comm.samples)
      out$alpha.comm.samples <- mcmc(out$alpha.comm.samples)
      out$tau.sq.beta.samples <- mcmc(out$tau.sq.beta.samples)
      out$tau.sq.alpha.samples <- mcmc(out$tau.sq.alpha.samples)
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
      colnames(out$beta.comm.samples) <- occ.cov.names
      colnames(out$tau.sq.beta.samples) <- occ.cov.names
      colnames(out$alpha.comm.samples) <- det.cov.names
      colnames(out$tau.sq.alpha.samples) <- det.cov.names
      # Species-specific effects
      beta.names <- paste(rep(sp.names, each = p), occ.cov.names, sep = '-')
      alpha.names <- paste(rep(sp.names, each = p.det), det.cov.names, sep = '-')
      colnames(out$beta.samples) <- beta.names
      colnames(out$alpha.samples) <- alpha.names
      
      # Add other things to the list
      out$n.samples <- n.samples
      out$n.burn <- n.burn
      out$n.thin <- n.thin
      out$n.chains <- n.chains
      out$p <- p
      out$p.det <- p.det
      out$N <- N
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
      class(out) <- "msPGOcc"
    }
  }
  out$run.time <- proc.time() - ptm
  out
}