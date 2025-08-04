spPGOcc <- function(occ.formula, det.formula, data, inits, priors, tuning,
                    cov.model = 'exponential', NNGP = TRUE, n.neighbors = 15,
                    search.type = "cb", n.batch, batch.length,
                    accept.rate = 0.43, n.omp.threads = 1, verbose = TRUE,
                    n.report = 100, n.burn = round(.10 * n.batch * batch.length),
                    n.thin = 1, n.chains = 1, k.fold, k.fold.threads = 1,
                    k.fold.seed, k.fold.only = FALSE, folds) {
  
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
  if (!'coords' %in% names(data)) {
    stop("error: coords must be specified in data for a spatial occupancy model")
  }
  coords <- as.matrix(data$coords)
  if (missing(n.batch)) {
    stop("error: n.batch must be specified")
  }
  if (missing(batch.length)) {
    stop("error: batch.length must be specified")
  }
  n.samples <- n.batch * batch.length
  if (n.burn > n.samples) {
    stop("n.burn must be less than n.samples")
  }
  if (n.thin > n.samples) {
    stop("n.thin must be less than n.samples")
  }
  if (n.thin %% 1 != 0) {
    stop("n.thin must be a whole number")
  }
  if (!missing(k.fold)) {
    if (!is.numeric(k.fold) | k.fold < 2) {
      stop("k.fold must be an integer greater than 1")
    }
  }
  if (NNGP) {
    if (!(search.type %in% c('brute', 'cb'))) {
      stop("search.type must be either 'brute' or 'cb'")
    }
  }
  cov.model.names <- c("exponential", "spherical", "matern", "gaussian")
  if(! cov.model %in% cov.model.names){
    stop("error: specified cov.model '", cov.model, "' is not a valid option; choose from ",
         paste(cov.model.names, collapse=", ", sep="") ,".")
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
  # sigma.sq
  if ("sigma.sq.ig" %in% names(priors)) {
    if (!is.list(priors$sigma.sq.ig) | length(priors$sigma.sq.ig) != 2) {
      stop("error: sigma.sq.ig must be a list of length 2")
    }
    sigma.sq.a <- priors$sigma.sq.ig[[1]]
    sigma.sq.b <- priors$sigma.sq.ig[[2]]
    if (length(sigma.sq.a) != 1 | length(sigma.sq.b) != 1) {
      stop("error: sigma.sq.ig must be a list of length 2 with each element a single value")
    }
  } else {
    if (verbose) {
      message("No prior specified for sigma.sq.ig.\nSetting prior shape to 2 and prior scale to 1\n")
    }
    sigma.sq.a <- 2
    sigma.sq.b <- 1
  }
  # phi
  if ("phi.unif" %in% names(priors)) {
    if (!is.list(priors$phi.unif) | length(priors$phi.unif) != 2) {
      stop("error: phi.unif must be a list of length 2")
    }
    phi.a <- priors$phi.unif[[1]]
    phi.b <- priors$phi.unif[[2]]
    if (length(phi.a) != 1 | length(phi.b) != 1) {
      stop("error: phi.unif must be a list of length 2 with each element a single value")
    }
  } else {
    if (verbose) {
      message("No prior specified for phi.unif.\nUsing uniform(3/max(dist.non), 3/min(dist.non)) prior for phi\n")
    }
    coords.D <- dist(coords)
    phi.a <- 3 / max(coords.D)
    phi.b <- 3 / min(coords.D[coords.D > 0])
  }
  # nu
  if (cov.model == "matern") {
    if ("nu.unif" %in% names(priors)) {
      if (!is.list(priors$nu.unif) | length(priors$nu.unif) != 2) {
        stop("error: nu.unif must be a list of length 2")
      }
      nu.a <- priors$nu.unif[[1]]
      nu.b <- priors$nu.unif[[2]]
      if (length(nu.a) != 1 | length(nu.b) != 1) {
        stop("error: nu.unif must be a list of length 2 with each element a single value")
      }
    } else {
      if (verbose) {
        message("No prior specified for nu.unif.\nUsing uniform(0.1, 2) prior for nu\n")
      }
      nu.a <- 0.1
      nu.b <- 2
    }
  } else {
    nu.a <- 0
    nu.b <- 0
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
  # sigma.sq
  if ("sigma.sq" %in% names(inits)) {
    sigma.sq.inits <- inits[["sigma.sq"]]
    if (length(sigma.sq.inits) != 1) {
      stop("error: initial values for sigma.sq must be of length 1")
    }
  } else {
    sigma.sq.inits <- runif(1, 0.05, 2)
    if (verbose) {
      message("sigma.sq is not specified in initial values.\nSetting initial value to random value between 0.05 and 2\n")
    }
  }
  # phi
  if ("phi" %in% names(inits)) {
    phi.inits <- inits[["phi"]]
    if (length(phi.inits) != 1) {
      stop("error: initial values for phi must be of length 1")
    }
  } else {
    phi.inits <- runif(1, phi.a, phi.b)
    if (verbose) {
      message("phi is not specified in initial values.\nSetting initial value to random value from the prior\n")
    }
  }
  # nu
  if (cov.model == "matern") {
    if ("nu" %in% names(inits)) {
      nu.inits <- inits[["nu"]]
      if (length(nu.inits) != 1) {
        stop("error: initial values for nu must be of length 1")
      }
    } else {
      nu.inits <- runif(1, nu.a, nu.b)
      if (verbose) {
        message("nu is not specified in initial values.\nSetting initial value to random value from the prior\n")
      }
    }
  } else {
    nu.inits <- 0
  }
  # w
  if ("w" %in% names(inits)) {
    w.inits <- inits[["w"]]
    if (length(w.inits) != J) {
      stop(paste("error: initial values for w must be of length ", J, sep = ""))
    }
  } else {
    w.inits <- rep(0, J)
    if (verbose) {
      message("w is not specified in initial values.\nSetting initial values to 0\n")
    }
  }
  
  # Tuning ------------------------------------------------------------------
  if (missing(tuning)) {
    phi.tuning <- 1
    if (cov.model == 'matern') {
      nu.tuning <- 1
    } else {
      nu.tuning <- 0
    }
  } else {
    names(tuning) <- tolower(names(tuning))
    if (!'phi' %in% names(tuning)) {
      stop("error: phi must be specified in tuning value list")
    }
    phi.tuning <- tuning$phi
    if (cov.model == 'matern') {
      if (!'nu' %in% names(tuning)) {
        stop("error: nu must be specified in tuning value list for a matern cov model")
      }
      nu.tuning <- tuning$nu
    } else {
      nu.tuning <- 0
    }
  }
  
  # NNGP --------------------------------------------------------------------
  if (NNGP) {
    # Get neighbors
    if (search.type == 'brute') {
      ord <- order(coords[,1])
      neighbor.indx <- brute_force_neighbors(coords, n.neighbors, ord)
    } else {
      neighbor.indx <- cb_neighbors(coords, n.neighbors)
      ord <- neighbor.indx$ord
      neighbor.indx <- neighbor.indx$neighbor.indx
    }
    # Save objects for later
    ord.orig <- order(ord)
    u.search.type <- 2
    if (search.type == 'brute') {
      u.search.type <- 1
    }
    neighbor.indx.orig <- neighbor.indx
    # Reorder everything to be consistent with NNGP
    y <- y[ord, ]
    X <- X[ord, , drop = FALSE]
    coords <- coords[ord, ]
    z.inits <- z.inits[ord]
    w.inits <- w.inits[ord]
    # Reorder detection covariates
    X.p.full <- X.p
    X.p <- X.p[do.call(c, lapply(ord, function(a) (a-1)*n.rep + 1:n.rep)), , drop = FALSE]
  } else {
    ord <- 1:J
    ord.orig <- 1:J
    neighbor.indx <- NULL
    neighbor.indx.orig <- NULL
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
  w.inits.i <- w.inits
  sigma.sq.inits.i <- sigma.sq.inits
  phi.inits.i <- phi.inits
  nu.inits.i <- nu.inits
  
  # Create a list of initial values for each chain
  inits.list <- list()
  for (i in 1:n.chains) {
    inits.list[[i]] <- list(z = z.inits.i,
                            beta = beta.inits.i,
                            alpha = alpha.inits.i,
                            sigma.sq = sigma.sq.inits.i,
                            phi = phi.inits.i,
                            nu = nu.inits.i,
                            w = w.inits.i)
  }
  
  # If not running cross-validation
  if (missing(k.fold)) {
    # Fit model
    out <- spPGOccModel(y, X, X.p, coords, p, p.det, J, n.rep,
                        n.obs, n.occ.params, n.det.params,
                        mu.beta, sigma.beta, mu.alpha, sigma.alpha,
                        phi.a, phi.b, sigma.sq.a, sigma.sq.b, nu.a, nu.b,
                        phi.tuning, nu.tuning,
                        n.batch, batch.length, accept.rate,
                        n.omp.threads, verbose, n.report,
                        n.burn, n.thin, n.chains, inits.list,
                        cov.model, NNGP, neighbor.indx, ord)
    
    # Return results --------------------------------------------------------
    # Convert lists to matrices
    out$beta.samples <- do.call(rbind, out$beta.samples)
    out$alpha.samples <- do.call(rbind, out$alpha.samples)
    out$z.samples <- do.call(rbind, out$z.samples)
    out$psi.samples <- do.call(rbind, out$psi.samples)
    out$w.samples <- do.call(rbind, out$w.samples)
    out$theta.samples <- do.call(rbind, out$theta.samples)
    # Get rid of initial value lists
    out$z.inits <- NULL
    out$beta.inits <- NULL
    out$alpha.inits <- NULL
    out$w.inits <- NULL
    
    # Create MCMC list
    out$beta.samples <- mcmc(out$beta.samples)
    out$alpha.samples <- mcmc(out$alpha.samples)
    out$z.samples <- mcmc(out$z.samples[, ord.orig])
    out$psi.samples <- mcmc(out$psi.samples[, ord.orig])
    out$theta.samples <- mcmc(out$theta.samples)
    out$w.samples <- mcmc(out$w.samples[, ord.orig])
    
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
    if (cov.model == 'matern') {
      colnames(out$theta.samples) <- c('sigma.sq', 'phi', 'nu')
    } else {
      colnames(out$theta.samples) <- c('sigma.sq', 'phi')
    }
    w.names <- paste('w[', 1:J, ']', sep = '')
    colnames(out$w.samples) <- w.names
    out$y <- y[ord.orig, ]
    out$X <- X[ord.orig, , drop = FALSE]
    if (NNGP) {
      out$X.p <- X.p.full
    } else {
      out$X.p <- X.p
    }
    out$coords <- coords[ord.orig, ]
    
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
    out$cov.model <- cov.model
    out$NNGP <- NNGP
    if (NNGP) {
      out$n.neighbors <- n.neighbors
      out$ord <- ord
    }
    class(out) <- "spPGOcc"
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
    
    # Custom folds
    if (missing(folds)) {
      fold.indx <- sample(rep(1:k.fold, length.out = J))
    } else {
      fold.indx <- folds
      k.fold <- max(folds)
    }
    
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
                          # Coordinates
                          coords.fit <- coords[sites.fit, ]
                          coords.hold <- coords[sites.hold, ]
                          # Initial values
                          # Grab a slice of the arrays for the first chain
                          z.inits.i <- z.inits[sites.fit]
                          beta.inits.i <- beta.inits
                          alpha.inits.i <- alpha.inits
                          w.inits.i <- w.inits[sites.fit]
                          sigma.sq.inits.i <- sigma.sq.inits
                          phi.inits.i <- phi.inits
                          nu.inits.i <- nu.inits
                          # Need to make sure initial values for z are correct.
                          z.inits.i[z.inits.i == -Inf] <- 0
                          # Create a list of initial values for each chain
                          inits.list.i <- list()
                          for (j in 1:n.chains) {
                            inits.list.i[[j]] <- list(z = z.inits.i,
                                                      beta = beta.inits.i,
                                                      alpha = alpha.inits.i,
                                                      w = w.inits.i,
                                                      sigma.sq = sigma.sq.inits.i,
                                                      phi = phi.inits.i,
                                                      nu = nu.inits.i)
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
                          # NNGP
                          if (NNGP) {
                            # Get neighbors
                            if (search.type == 'brute') {
                              ord.i <- order(coords.fit[,1])
                              neighbor.indx.i <- brute_force_neighbors(coords.fit, n.neighbors, ord.i)
                            } else {
                              neighbor.indx.i <- cb_neighbors(coords.fit, n.neighbors)
                              ord.i <- neighbor.indx.i$ord
                              neighbor.indx.i <- neighbor.indx.i$neighbor.indx
                            }
                            # Reorder everything to be consistent with NNGP
                            y.fit <- y.fit[ord.i, ]
                            X.fit <- X.fit[ord.i, , drop = FALSE]
                            coords.fit <- coords.fit[ord.i, ]
                            z.inits.i <- z.inits.i[ord.i]
                            w.inits.i <- w.inits.i[ord.i]
                            # Reorder detection covariates
                            X.p.fit.full <- X.p.fit
                            X.p.fit <- X.p.fit[do.call(c, lapply(ord.i, function(a) (a-1)*n.rep + 1:n.rep)), , drop = FALSE]
                          } else {
                            ord.i <- 1:J.fit
                            neighbor.indx.i <- NULL
                          }
                          # Fit the model
                          out.i <- spPGOccModel(y.fit, X.fit, X.p.fit, coords.fit, p, p.det, J.fit, n.rep,
                                                n.obs, n.occ.params, n.det.params,
                                                mu.beta, sigma.beta, mu.alpha, sigma.alpha,
                                                phi.a, phi.b, sigma.sq.a, sigma.sq.b, nu.a, nu.b,
                                                phi.tuning, nu.tuning,
                                                n.batch, batch.length, accept.rate,
                                                n.omp.threads, verbose = FALSE, n.report = n.report,
                                                n.burn, n.thin, n.chains, inits.list.i,
                                                cov.model, NNGP, neighbor.indx.i, ord.i)
                          # Convert lists to matrices
                          out.i$beta.samples <- do.call(rbind, out.i$beta.samples)
                          out.i$alpha.samples <- do.call(rbind, out.i$alpha.samples)
                          out.i$theta.samples <- do.call(rbind, out.i$theta.samples)
                          n.post.samples <- nrow(out.i$beta.samples)
                          
                          # Predict occurrence at held out locations
                          if (NNGP) {
                            coords.all <- rbind(coords.fit, coords.hold)
                            ord.pred <- order(coords.all[, 1])
                            y.pred.indx <- which(ord.pred > J.fit)
                            coords.all <- coords.all[ord.pred, ]
                            if (search.type == 'brute') {
                              neighbor.indx.pred <- brute_force_neighbors(coords.all, n.neighbors, 1:nrow(coords.all))
                            } else {
                              neighbor.indx.pred <- cb_neighbors(coords.all, n.neighbors)
                              neighbor.indx.pred <- neighbor.indx.pred$neighbor.indx
                            }
                            w.pred <- predict(out.i, X.hold, coords.hold, n.omp.threads, verbose = FALSE)$w.samples
                            psi.means <- apply(logit.inv(X.hold %*% t(out.i$beta.samples) + w.pred), 1, mean)
                          } else {
                            w.pred <- predict(out.i, X.hold, coords.hold, n.omp.threads, verbose = FALSE)$w.samples
                            psi.means <- apply(logit.inv(X.hold %*% t(out.i$beta.samples) + w.pred), 1, mean)
                          }
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
      class(out) <- "spPGOcc"
    } else {
      # Fit model to the full data set
      if (verbose) {
        message("Finished cross-validation. Fitting model to full data set.\n")
      }
      # Fit model
      out <- spPGOccModel(y, X, X.p, coords, p, p.det, J, n.rep,
                          n.obs, n.occ.params, n.det.params,
                          mu.beta, sigma.beta, mu.alpha, sigma.alpha,
                          phi.a, phi.b, sigma.sq.a, sigma.sq.b, nu.a, nu.b,
                          phi.tuning, nu.tuning,
                          n.batch, batch.length, accept.rate,
                          n.omp.threads, verbose, n.report,
                          n.burn, n.thin, n.chains, inits.list,
                          cov.model, NNGP, neighbor.indx, ord)
      # Convert lists to matrices
      out$beta.samples <- do.call(rbind, out$beta.samples)
      out$alpha.samples <- do.call(rbind, out$alpha.samples)
      out$z.samples <- do.call(rbind, out$z.samples)
      out$psi.samples <- do.call(rbind, out$psi.samples)
      out$w.samples <- do.call(rbind, out$w.samples)
      out$theta.samples <- do.call(rbind, out$theta.samples)
      # Get rid of initial value lists
      out$z.inits <- NULL
      out$beta.inits <- NULL
      out$alpha.inits <- NULL
      out$w.inits <- NULL
      
      # Create MCMC list
      out$beta.samples <- mcmc(out$beta.samples)
      out$alpha.samples <- mcmc(out$alpha.samples)
      out$z.samples <- mcmc(out$z.samples[, ord.orig])
      out$psi.samples <- mcmc(out$psi.samples[, ord.orig])
      out$theta.samples <- mcmc(out$theta.samples)
      out$w.samples <- mcmc(out$w.samples[, ord.orig])
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
      if (cov.model == 'matern') {
        colnames(out$theta.samples) <- c('sigma.sq', 'phi', 'nu')
      } else {
        colnames(out$theta.samples) <- c('sigma.sq', 'phi')
      }
      w.names <- paste('w[', 1:J, ']', sep = '')
      colnames(out$w.samples) <- w.names
      out$y <- y[ord.orig, ]
      out$X <- X[ord.orig, , drop = FALSE]
      if (NNGP) {
        out$X.p <- X.p.full
      } else {
        out$X.p <- X.p
      }
      out$coords <- coords[ord.orig, ]
      
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
      out$cov.model <- cov.model
      out$NNGP <- NNGP
      if (NNGP) {
        out$n.neighbors <- n.neighbors
        out$ord <- ord
      }
      out$k.fold.deviance <- dev.vals
      out$k.fold.auc <- auc.vals
      out$k.fold.rmse <- rmse.vals
      out$k.fold.tjur <- tjur.vals
      class(out) <- "spPGOcc"
    }
  }
  out$run.time <- proc.time() - ptm
  out
}