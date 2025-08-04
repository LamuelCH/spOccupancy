#' @export
sfMsPGOcc_LC <- function(occ.formula, det.formula, data, inits, priors, tuning,
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
  # lambda
  if ("lambda" %in% names(inits)) {
    lambda.inits <- inits[["lambda"]]
    if (is.matrix(lambda.inits)) {
      if (nrow(lambda.inits) != N | ncol(lambda.inits) != 1) {
        stop(paste("error: initial values for lambda must be a matrix with dimensions ", N, "x", 1, " or a single vector of length ", N, sep = ""))
      }
    } else {
      if (length(lambda.inits) != N) {
        stop(paste("error: initial values for lambda must be a matrix with dimensions ", N, "x", 1, " or a single vector of length ", N, sep = ""))
      }
      lambda.inits <- matrix(lambda.inits, N, 1)
    }
  } else {
    lambda.inits <- matrix(rnorm(N, 0, 1), N, 1)
    if (verbose) {
      message("lambda is not specified in initial values.\nSetting initial values to random values from a standard normal distribution\n")
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
    y <- y[, ord, ]
    X <- X[ord, , drop = FALSE]
    coords <- coords[ord, ]
    z.inits <- z.inits[, ord]
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
  w.inits.i <- w.inits
  phi.inits.i <- phi.inits
  lambda.inits.i <- lambda.inits
  nu.inits.i <- nu.inits
  
  # Create a list of initial values for each chain
  inits.list <- list()
  for (i in 1:n.chains) {
    inits.list[[i]] <- list(z = z.inits.i,
                            beta = beta.inits.i,
                            alpha = alpha.inits.i,
                            beta.comm = beta.comm.inits.i,
                            alpha.comm = alpha.comm.inits.i,
                            tau.sq.beta = tau.sq.beta.inits.i,
                            tau.sq.alpha = tau.sq.alpha.inits.i,
                            phi = phi.inits.i,
                            lambda = lambda.inits.i,
                            nu = nu.inits.i,
                            w = w.inits.i)
  }
  
  # If not running cross-validation
  if (missing(k.fold)) {
    # Fit model
    out <- sfMsPGOccModel(y, X, X.p, coords, p, p.det, N, J, n.rep,
                          n.obs, n.occ.params, n.det.params,
                          mu.beta.comm, sigma.beta.comm, mu.alpha.comm, sigma.alpha.comm,
                          tau.sq.beta.a, tau.sq.beta.b, tau.sq.alpha.a, tau.sq.alpha.b,
                          phi.a, phi.b, nu.a, nu.b,
                          phi.tuning, nu.tuning,
                          n.batch, batch.length, accept.rate,
                          n.omp.threads, verbose, n.report,
                          n.burn, n.thin, n.chains, inits.list,
                          cov.model, NNGP, neighbor.indx, ord)
    
    # Return results --------------------------------------------------------
    # Convert lists to matrices
    out$beta.samples <- do.call(rbind, out$beta.samples)
    out$alpha.samples <- do.call(rbind, out$alpha.samples)
    out$beta.comm.samples <- do.call(rbind, out$beta.comm.samples)
    out$alpha.comm.samples <- do.call(rbind, out$alpha.comm.samples)
    out$tau.sq.beta.samples <- do.call(rbind, out$tau.sq.beta.samples)
    out$tau.sq.alpha.samples <- do.call(rbind, out$tau.sq.alpha.samples)
    out$z.samples <- do.call(abind::abind, c(out$z.samples, along = 0))
    out$psi.samples <- do.call(abind::abind, c(out$psi.samples, along = 0))
    out$w.samples <- do.call(abind::abind, c(out$w.samples, along = 0))
    out$lambda.samples <- do.call(rbind, out$lambda.samples)
    out$theta.samples <- do.call(rbind, out$theta.samples)
    # Get rid of initial value lists
    out$z.inits <- NULL
    out$beta.inits <- NULL
    out$alpha.inits <- NULL
    out$w.inits <- NULL
    out$lambda.inits <- NULL
    
    # Create MCMC list
    out$beta.samples <- mcmc(out$beta.samples)
    out$alpha.samples <- mcmc(out$alpha.samples)
    out$beta.comm.samples <- mcmc(out$beta.comm.samples)
    out$alpha.comm.samples <- mcmc(out$alpha.comm.samples)
    out$tau.sq.beta.samples <- mcmc(out$tau.sq.beta.samples)
    out$tau.sq.alpha.samples <- mcmc(out$tau.sq.alpha.samples)
    out$z.samples <- mcmc(out$z.samples[, , ord.orig])
    out$psi.samples <- mcmc(out$psi.samples[, , ord.orig])
    out$w.samples <- mcmc(out$w.samples[, ord.orig])
    out$lambda.samples <- mcmc(out$lambda.samples)
    out$theta.samples <- mcmc(out$theta.samples)
    
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
    lambda.names <- paste(rep(sp.names, each = 1), 'lambda', sep = '-')
    colnames(out$lambda.samples) <- lambda.names
    if (cov.model == 'matern') {
      colnames(out$theta.samples) <- c('phi', 'nu')
    } else {
      colnames(out$theta.samples) <- 'phi'
    }
    out$y <- y[, ord.orig, ]
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
    out$N <- N
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
    class(out) <- "sfMsPGOcc"
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
                          # Coordinates
                          coords.fit <- coords[sites.fit, ]
                          coords.hold <- coords[sites.hold, ]
                          # Initial values
                          # Grab a slice of the arrays for the first chain
                          z.inits.i <- z.inits[, sites.fit]
                          beta.inits.i <- beta.inits
                          alpha.inits.i <- alpha.inits
                          beta.comm.inits.i <- beta.comm.inits
                          alpha.comm.inits.i <- alpha.comm.inits
                          tau.sq.beta.inits.i <- tau.sq.beta.inits
                          tau.sq.alpha.inits.i <- tau.sq.alpha.inits
                          w.inits.i <- w.inits[sites.fit]
                          phi.inits.i <- phi.inits
                          lambda.inits.i <- lambda.inits
                          nu.inits.i <- nu.inits
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
                                                      tau.sq.alpha = tau.sq.alpha.inits.i,
                                                      phi = phi.inits.i,
                                                      lambda = lambda.inits.i,
                                                      nu = nu.inits.i,
                                                      w = w.inits.i)
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
                            y.fit <- y.fit[, ord.i, ]
                            X.fit <- X.fit[ord.i, , drop = FALSE]
                            coords.fit <- coords.fit[ord.i, ]
                            z.inits.i <- z.inits.i[, ord.i]
                            w.inits.i <- w.inits.i[ord.i]
                            # Reorder detection covariates
                            X.p.fit.full <- X.p.fit
                            X.p.fit <- X.p.fit[do.call(c, lapply(ord.i, function(a) (a-1)*n.rep + 1:n.rep)), , drop = FALSE]
                          } else {
                            ord.i <- 1:J.fit
                            neighbor.indx.i <- NULL
                          }
                          # Fit the model
                          out.i <- sfMsPGOccModel(y.fit, X.fit, X.p.fit, coords.fit, p, p.det, N, J.fit, n.rep,
                                                  n.obs, n.occ.params, n.det.params,
                                                  mu.beta.comm, sigma.beta.comm, mu.alpha.comm, sigma.alpha.comm,
                                                  tau.sq.beta.a, tau.sq.beta.b, tau.sq.alpha.a, tau.sq.alpha.b,
                                                  phi.a, phi.b, nu.a, nu.b,
                                                  phi.tuning, nu.tuning,
                                                  n.batch, batch.length, accept.rate,
                                                  n.omp.threads, verbose = FALSE, n.report = n.report,
                                                  n.burn, n.thin, n.chains, inits.list.i,
                                                  cov.model, NNGP, neighbor.indx.i, ord.i)
                          # Convert lists to matrices
                          out.i$beta.samples <- do.call(rbind, out.i$beta.samples)
                          out.i$alpha.samples <- do.call(rbind, out.i$alpha.samples)
                          out.i$lambda.samples <- do.call(rbind, out.i$lambda.samples)
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
                            psi.means <- apply(logit.inv(X.hold %*% t(matrix(out.i$beta.samples, N, p, byrow = TRUE)) +
                                                           out.i$lambda.samples %*% t(w.pred)), 1, mean)
                          } else {
                            w.pred <- predict(out.i, X.hold, coords.hold, n.omp.threads, verbose = FALSE)$w.samples
                            psi.means <- apply(logit.inv(X.hold %*% t(matrix(out.i$beta.samples, N, p, byrow = TRUE)) +
                                                           out.i$lambda.samples %*% t(w.pred)), 1, mean)
                          }
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
    
    ## MODIFICATION START
    # Removed summary statistics calculation
    if (k.fold.only) {
      out <- list(k.fold.deviance = dev.vals,
                  k.fold.auc = auc.vals,
                  k.fold.rmse = rmse.vals,
                  k.fold.tjur = tjur.vals,
                  call = cl)
      class(out) <- "sfMsPGOcc"
    } else {
      # Fit model to the full data set
      if (verbose) {
        message("Finished cross-validation. Fitting model to full data set.\n")
      }
      # Fit model
      out <- sfMsPGOccModel(y, X, X.p, coords, p, p.det, N, J, n.rep,
                            n.obs, n.occ.params, n.det.params,
                            mu.beta.comm, sigma.beta.comm, mu.alpha.comm, sigma.alpha.comm,
                            tau.sq.beta.a, tau.sq.beta.b, tau.sq.alpha.a, tau.sq.alpha.b,
                            phi.a, phi.b, nu.a, nu.b,
                            phi.tuning, nu.tuning,
                            n.batch, batch.length, accept.rate,
                            n.omp.threads, verbose, n.report,
                            n.burn, n.thin, n.chains, inits.list,
                            cov.model, NNGP, neighbor.indx, ord)
      # Convert lists to matrices
      out$beta.samples <- do.call(rbind, out$beta.samples)
      out$alpha.samples <- do.call(rbind, out$alpha.samples)
      out$beta.comm.samples <- do.call(rbind, out$beta.comm.samples)
      out$alpha.comm.samples <- do.call(rbind, out$alpha.comm.samples)
      out$tau.sq.beta.samples <- do.call(rbind, out$tau.sq.beta.samples)
      out$tau.sq.alpha.samples <- do.call(rbind, out$tau.sq.alpha.samples)
      out$z.samples <- do.call(abind::abind, c(out$z.samples, along = 0))
      out$psi.samples <- do.call(abind::abind, c(out$psi.samples, along = 0))
      out$w.samples <- do.call(abind::abind, c(out$w.samples, along = 0))
      out$lambda.samples <- do.call(rbind, out$lambda.samples)
      out$theta.samples <- do.call(rbind, out$theta.samples)
      # Get rid of initial value lists
      out$z.inits <- NULL
      out$beta.inits <- NULL
      out$alpha.inits <- NULL
      out$w.inits <- NULL
      out$lambda.inits <- NULL
      
      # Create MCMC list
      out$beta.samples <- mcmc(out$beta.samples)
      out$alpha.samples <- mcmc(out$alpha.samples)
      out$beta.comm.samples <- mcmc(out$beta.comm.samples)
      out$alpha.comm.samples <- mcmc(out$alpha.comm.samples)
      out$tau.sq.beta.samples <- mcmc(out$tau.sq.beta.samples)
      out$tau.sq.alpha.samples <- mcmc(out$tau.sq.alpha.samples)
      out$z.samples <- mcmc(out$z.samples[, , ord.orig])
      out$psi.samples <- mcmc(out$psi.samples[, , ord.orig])
      out$w.samples <- mcmc(out$w.samples[, ord.orig])
      out$lambda.samples <- mcmc(out$lambda.samples)
      out$theta.samples <- mcmc(out$theta.samples)
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
      lambda.names <- paste(rep(sp.names, each = 1), 'lambda', sep = '-')
      colnames(out$lambda.samples) <- lambda.names
      if (cov.model == 'matern') {
        colnames(out$theta.samples) <- c('phi', 'nu')
      } else {
        colnames(out$theta.samples) <- 'phi'
      }
      out$y <- y[, ord.orig, ]
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
      out$N <- N
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
      class(out) <- "sfMsPGOcc"
    }
    ## MODIFICATION END
  }
  out$run.time <- proc.time() - ptm
  out
}