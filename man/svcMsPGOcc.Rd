\name{svcMsPGOcc}
\alias{svcMsPGOcc}
\title{Function for Fitting Multi-Species Spatially-Varying Coefficient Occupancy Models}

\description{
  The function \code{svcMsPGOcc} fits multi-species spatially-varying coefficient occupancy models with species correlations (i.e., a spatially-explicit joint species distribution model with imperfect detection). We use Polya-Gamma latent variables and a spatial factor modeling approach. Models are implemented using a Nearest Neighbor Gaussian Process. 
}

\usage{
svcMsPGOcc(occ.formula, det.formula, data, inits, priors, tuning, 
           svc.cols = 1, cov.model = 'exponential', NNGP = TRUE, 
           n.neighbors = 15, search.type = 'cb', std.by.sp = FALSE, 
           n.factors, n.batch, batch.length, 
           accept.rate = 0.43, n.omp.threads = 1, 
           verbose = TRUE, n.report = 100, 
           n.burn = round(.10 * n.batch * batch.length), n.thin = 1, 
           n.chains = 1, ...)
}

\arguments{

  \item{occ.formula}{a symbolic description of the model to be fit
    for the occurrence portion of the model using R's model syntax. Random
    intercepts are allowed using lme4 syntax (Bates et al. 2015). Only
    right-hand side of formula is specified. See example below.}
  
  \item{det.formula}{a symbolic description of the model to be fit
    for the detection portion of the model using R's model syntax. Only
    right-hand side of formula is specified. See example below. Random 
    intercepts are allowed using lme4 syntax (Bates et al. 2015).}

  \item{data}{a list containing data necessary for model fitting.
    Valid tags are \code{y}, \code{occ.covs}, \code{det.covs}, 
    \code{coords}, and \code{range.ind}. 
    \code{y} is a three-dimensional array with first dimension equal to the 
    number of species, second dimension equal to the number of sites, and 
    third dimension equal to the maximum number of replicates at a given site. 
    \code{occ.covs} is a matrix or data frame containing the variables 
    used in the occurrence portion of the model, with \eqn{J}{J} rows for 
    each column (variable). \code{det.covs} is a list of variables included 
    in the detection portion of the model. Each list element is a different 
    detection covariate, which can be site-level or observational-level. 
    Site-level covariates are specified as a vector of length \eqn{J}{J} while 
    observation-level covariates are specified as a matrix or data frame 
    with the number of rows equal to \eqn{J}{J} and number of columns equal 
    to the maximum number of replicates at a given site. \code{coords} is a 
    \eqn{J \times 2}{J x 2} matrix of the observation coordinates. Note that 
    \code{spOccupancy} assumes coordinates are specified in a projected coordinate system.
    \code{range.ind} is a matrix with rows corresponding to species and columns
    corresponding to sites, with each element taking value 1 if that site is 
    within the range of the corresponding species and 0 if it is outside of the 
    range. This matrix is not required, but it can be helpful to restrict the 
    modeled area for each individual species to be within the realistic range
    of locations for that species when estimating the model parameters. This is
    applicable when auxiliary data sources are available on the realistic range
    of the species.}

  \item{inits}{a list with each tag corresponding to a parameter name.
    Valid tags are \code{alpha.comm}, \code{beta.comm}, \code{beta}, 
    \code{alpha}, \code{tau.sq.beta}, \code{tau.sq.alpha}, 
    \code{sigma.sq.psi}, \code{sigma.sq.p}, \code{z}, 
    \code{phi}, \code{lambda}, and \code{nu}. \code{nu} is only specified if 
    \code{cov.model = "matern"}, and \code{sigma.sq.psi} and \code{sigma.sq.p} are
    only specified if random effects are included in \code{occ.formula} or 
    \code{det.formula}, respectively. The value portion of each tag is 
    the parameter's initial value. See \code{priors} description for definition
    of each parameter name. Additionally, the tag \code{fix} can be set to \code{TRUE} 
    to fix the starting values across all chains. If \code{fix} is not specified
    (the default), starting values are varied randomly across chains.}

  \item{priors}{a list with each tag corresponding to a parameter name. 
  Valid tags are \code{beta.comm.normal}, \code{alpha.comm.normal}, 
  \code{tau.sq.beta.ig}, \code{tau.sq.alpha.ig}, \code{sigma.sq.psi}, 
  \code{sigma.sq.p}, \code{phi.unif}, 
  and \code{nu.unif}. Community-level occurrence (\code{beta.comm}) and detection 
  (\code{alpha.comm}) regression coefficients are assumed to follow a 
  normal distribution. The hyperparameters of the normal distribution
  are passed as a list of length two with the first and second elements 
  corresponding to the mean and variance of the normal distribution, 
  which are each specified as vectors of length equal to the number of 
  coefficients to be estimated or of length one if priors are the same for 
  all coefficients. If not specified, prior means are set 
  to 0 and prior variances set to 2.72. By default, community-level variance parameters 
  for occupancy (\code{tau.sq.beta}) and detection (\code{tau.sq.alpha}) are 
  assumed to follow an inverse Gamma distribution. The hyperparameters of 
  the inverse gamma distribution are passed as a list of length two with 
  the first and second elements corresponding to the shape and scale parameters,
  which are each specified as vectors of length equal to the number of 
  coefficients to be estimated or a single value if priors are the same for all 
  parameters. If not specified, prior shape and scale 
  parameters are set to 0.1. 
  The spatial factor model fits \code{n.factors} independent
  spatial processes for each spatially-varying coefficient specified in \code{svc.cols}. 
  The spatial decay \code{phi} and smoothness \code{nu} parameters 
  for each latent factor are assumed to follow Uniform distributions. 
  The hyperparameters of the Uniform are passed as a list with two elements, 
  with both elements being vectors of length \code{n.factors * length(svc.cols)} 
  corresponding to the lower and 
  upper support, respectively, or as a single value if the same value is assigned
  for all factor/SVC combinations. The priors for the factor loadings 
  matrix \code{lambda} for each SVC are fixed
  following the standard spatial factor model to ensure parameter 
  identifiability (Christensen and Amemlya 2002). The
  upper triangular elements of the \code{N x n.factors} matrix are fixed at 0 and the 
  diagonal elements are fixed at 1 for each SVC. The lower triangular elements are assigned a 
  standard normal prior (i.e., mean 0 and variance 1).
  \code{sigma.sq.psi} and \code{sigma.sq.p} are the random 
  effect variances for any occurrence or
  detection random effects, respectively, and are assumed to follow an inverse
  Gamma distribution. The hyperparameters of the inverse-Gamma distribution
  are passed as a list of length two with first and second elements corresponding
  to the shape and scale parameters, respectively, which are each specified as
  vectors of length equal to the number of random intercepts or of length one 
  if priors are the same for all random effect variances.}

  \item{tuning}{a list with each tag corresponding to a parameter
    name. Valid tags are \code{phi} and \code{nu}. The value portion of each
    tag defines the initial variance of the adaptive sampler. We assume the 
    initial variance of the adaptive sampler is the same for each species, 
    although the adaptive sampler will adjust the tuning variances separately
    for each species. See Roberts and Rosenthal (2009) for details.}

  \item{svc.cols}{a vector indicating the variables whose effects will be 
    estimated as spatially-varying coefficients. \code{svc.cols} can be an 
    integer vector with values indicating the order of covariates specified 
    in the model formula (with 1 being the intercept if specified), or it can 
    be specified as a character vector with names corresponding to variable 
    names in \code{occ.covs} (for the intercept, use \code{'(Intercept)'}). \code{svc.cols}
    default argument of 1 results in a spatial occupancy model analogous to 
    \code{sfMsPGOcc} (assuming an intercept is included in the model).} 

  \item{cov.model}{a quoted keyword that specifies the covariance
    function used to model the spatial dependence structure among the
    observations.  Supported covariance model key words are:
    \code{"exponential"}, \code{"matern"}, \code{"spherical"}, and
    \code{"gaussian"}.}
  
  \item{NNGP}{if \code{TRUE}, model is fit with an NNGP. If \code{FALSE}, 
    a full Gaussian process is used. See Datta et al. (2016) and 
    Finley et al. (2019) for more information. Only 
    \code{NNGP = TRUE} is currently supported.}
  
  \item{n.neighbors}{number of neighbors used in the NNGP. Only used if 
  \code{NNGP = TRUE}. Datta et al. (2016) showed that 15 neighbors is usually 
  sufficient, but that as few as 5 neighbors can be adequate for certain data
  sets, which can lead to even greater decreases in run time. We recommend
  starting with 15 neighbors (the default) and if additional gains in computation
  time are desired, subsequently compare the results with a smaller number
  of neighbors using WAIC.}
  
  \item{search.type}{a quoted keyword that specifies the type of nearest
    neighbor search algorithm. Supported method key words are: \code{"cb"} and
    \code{"brute"}. The \code{"cb"} should generally be much
    faster. If locations do not have identical coordinate values on
    the axis used for the nearest neighbor ordering then \code{"cb"} 
    and \code{"brute"} should produce identical neighbor sets. 
    However, if there are identical coordinate values on the axis used 
    for nearest neighbor ordering, then \code{"cb"} and \code{"brute"} 
    might produce different, but equally valid, neighbor sets, 
    e.g., if data are on a grid. }

  \item{std.by.sp}{a logical value indicating whether the covariates are standardized
    separately for each species within the corresponding range for each species (\code{TRUE}) 
    or not (\code{FALSE}). Note that if \code{range.ind} is specified in \code{data.list}, 
    this will result in the covariates being standardized differently for each species
    based on the sites where \code{range.ind == 1} for that given species. If \code{range.ind} is not specified 
    and \code{std.by.sp = TRUE}, this will simply be equivalent to standardizing 
    the covariates across all locations prior to fitting the model. Note that the covariates
    in \code{occ.formula} should still be standardized across all locations. This can be done
    either outside the function, or can be done by specifying \code{scale()} in the model formula 
    around the continuous covariates.}

  \item{n.factors}{the number of factors to use in the spatial factor model approach. 
    Note this corresponds to the number of factors used for each spatially-varying
    coefficient that is estimated in the model.
    Typically, the number of factors is set to be small (e.g., 4-5) relative to the 
    total number of species in the community, which will lead to substantial 
    decreases in computation time. However, the value can be anywhere 
    between 1 and N (the number of species in the community).}

  \item{n.batch}{the number of MCMC batches  in each chain to run for the Adaptive MCMC 
    sampler. See Roberts and Rosenthal (2009) for details.}
  
  \item{batch.length}{the length of each MCMC batch to run for the Adaptive 
    MCMC sampler. See Roberts and Rosenthal (2009) for details.}
  
  \item{accept.rate}{target acceptance rate for Adaptive MCMC. Defaul is 
    0.43. See Roberts and Rosenthal (2009) for details.}
  
  \item{n.omp.threads}{a positive integer indicating
   the number of threads to use for SMP parallel processing within chains. The package must
   be compiled for OpenMP support. For most Intel-based machines, we
   recommend setting \code{n.omp.threads} up to the number of
   hyperthreaded cores. Note, \code{n.omp.threads} > 1 might not
   work on some systems.}
 
  \item{verbose}{if \code{TRUE}, messages about data preparation, 
    model specification, and progress of the sampler are printed to the screen. 
    Otherwise, no messages are printed.}
 
  \item{n.report}{the interval to report Metropolis sampler acceptance and
    MCMC progress. Note this is specified in terms of batches and not overall
    samples for spatial models.}

  \item{n.burn}{the number of samples out of the total \code{n.samples} to 
    discard as burn-in for each chain. By default, the first 10\% of samples is discarded.}
  
  \item{n.thin}{the thinning interval for collection of MCMC samples. The
    thinning occurs after the \code{n.burn} samples are discarded. Default 
    value is set to 1.}

  \item{n.chains}{the number of chains to run in sequence.}

  \item{...}{currently no additional arguments}
}

\note{
 Some of the underlying code used for generating random numbers from the
 Polya-Gamma distribution is taken from the \pkg{pgdraw} package written
 by Daniel F. Schmidt and Enes Makalic. Their code implements Algorithm
 6 in PhD thesis of Jesse Bennett Windle (2013) \url{https://repositories.lib.utexas.edu/handle/2152/21842}.
}

\references{
  Doser, J. W., Finley, A. O., Saunders, S. P., Kery, M., Weed, A. S., &
  Zipkin, E. F. (2024A). Modeling complex species-environment
  relationships through spatially-varying coefficient occupancy models.
  Journal of Agricultural, Biological and Environmental Statistics.
  \doi{10.1007/s13253-023-00595-6}.

  Doser, J. W., Kery, M., Saunders, S. P., Finley, A. O., Bateman, B. L.,
  Grand, J., Reault, S., Weed, A. S., & Zipkin, E. F. (2024B). Guidelines
  for the use of spatially varying coefficients in species distribution
  models. Global Ecology and Biogeography, 33(4), e13814.
  \doi{10.1111/geb.13814}.

  Datta, A., S. Banerjee, A.O. Finley, and A.E. Gelfand. (2016)
  Hierarchical Nearest-Neighbor Gaussian process models for large
  geostatistical datasets. \emph{Journal of the American Statistical
    Association}, \doi{10.1080/01621459.2015.1044091}.

  Finley, A.O., A. Datta, B.D. Cook, D.C. Morton, H.E. Andersen, and
  S. Banerjee. (2019) Efficient algorithms for Bayesian Nearest Neighbor
  Gaussian Processes. \emph{Journal of Computational and Graphical
    Statistics}, \doi{10.1080/10618600.2018.1537924}.

  Polson, N.G., J.G. Scott, and J. Windle. (2013) Bayesian Inference for
  Logistic Models Using Polya-Gamma Latent Variables.
  \emph{Journal of the American Statistical Association}, 108:1339-1349.

  Roberts, G.O. and Rosenthal J.S. (2009) Examples  of adaptive MCMC.
  \emph{Journal of Computational and Graphical Statistics}, 18(2):349-367.

  Bates, Douglas, Martin Maechler, Ben Bolker, Steve Walker (2015).
  Fitting Linear Mixed-Effects Models Using lme4. Journal of
  Statistical Software, 67(1), 1-48. \doi{10.18637/jss.v067.i01}.
}

\author{
  Jeffrey W. Doser \email{doserjef@msu.edu}, \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\value{
  An object of class \code{svcMsPGOcc} that is a list comprised of: 

  \item{beta.comm.samples}{a \code{coda} object of posterior samples
    for the community level occurrence regression coefficients.}
  
  \item{alpha.comm.samples}{a \code{coda} object of posterior samples
    for the community level detection regression coefficients.}

  \item{tau.sq.beta.samples}{a \code{coda} object of posterior samples
    for the occurrence community variance parameters.}

  \item{tau.sq.alpha.samples}{a \code{coda} object of posterior samples
    for the detection community variance parameters.}

  \item{beta.samples}{a \code{coda} object of posterior samples
    for the species level occurrence regression coefficients.}

  \item{alpha.samples}{a \code{coda} object of posterior samples
    for the species level detection regression coefficients.}

  \item{theta.samples}{a \code{coda} object of posterior samples
    for the species level correlation parameters for each spatially-varying coefficient.}

  \item{lambda.samples}{a \code{coda} object of posterior samples
    for the latent spatial factor loadings for each spatially-varying coefficient.}

  \item{z.samples}{a three-dimensional array of posterior samples for
    the latent occurrence values for each species.}

  \item{psi.samples}{a three-dimensional array of posterior samples for
    the latent occupancy probability values for each species.}

  \item{w.samples}{a four-dimensional array of posterior samples for 
    the latent spatial random effects for each spatial factor within each 
    spatially-varying coefficient. Dimensions correspond to MCMC sample, 
    factor, site, and spatially-varying coefficient.}

  \item{sigma.sq.psi.samples}{a \code{coda} object of posterior samples
    for variances of random intercepts included in the occurrence portion
    of the model. Only included if random intercepts are specified in 
    \code{occ.formula}.}

  \item{sigma.sq.p.samples}{a \code{coda} object of posterior samples
    for variances of random intercpets included in the detection portion 
    of the model. Only included if random intercepts are specified in 
    \code{det.formula}.}

  \item{beta.star.samples}{a \code{coda} object of posterior samples
    for the occurrence random effects. Only included if random intercepts
    are specified in \code{occ.formula}.}

  \item{alpha.star.samples}{a \code{coda} object of posterior samples
    for the detection random effects. Only included if random intercepts 
    are specified in \code{det.formula}.}

  \item{like.samples}{a three-dimensional array of posterior samples
    for the likelihood value associated with each site and species. 
    Used for calculating WAIC.}

  \item{rhat}{a list of Gelman-Rubin diagnostic values for some of the model
    parameters.}

  \item{ESS}{a list of effective sample sizes for some of the model parameters.}

  \item{run.time}{MCMC sampler execution time reported using \code{proc.time()}.}

  The return object will include additional objects used for 
  subsequent prediction and/or model fit evaluation. Note that detection 
  probability estimated values are not included in the model object, but can
  be extracted using \code{fitted()}.
}

\examples{
set.seed(400)

# Simulate Data -----------------------------------------------------------
J.x <- 10
J.y <- 10 
J <- J.x * J.y
n.rep <- sample(5, size = J, replace = TRUE)
N <- 6
# Community-level covariate effects
# Occurrence
beta.mean <- c(0.2, -0.2, 0.3, -0.1, 0.4)
p.occ <- length(beta.mean)
tau.sq.beta <- c(0.6, 1.5, 0.4, 0.5, 0.3)
# Detection
alpha.mean <- c(0, 1.2, -0.5)
tau.sq.alpha <- c(1, 0.5, 1.3)
p.det <- length(alpha.mean)
# Random effects
psi.RE <- list(levels = 15, 
               sigma.sq.psi = 0.7)
p.RE <- list(levels = 20, 
             sigma.sq.p = 0.5)
# Draw species-level effects from community means.
beta <- matrix(NA, nrow = N, ncol = p.occ)
alpha <- matrix(NA, nrow = N, ncol = p.det)
for (i in 1:p.occ) {
  beta[, i] <- rnorm(N, beta.mean[i], sqrt(tau.sq.beta[i]))
}
for (i in 1:p.det) {
  alpha[, i] <- rnorm(N, alpha.mean[i], sqrt(tau.sq.alpha[i]))
}
# Number of spatial factors for each SVC
n.factors <- 2
# The intercept and first two covariates have spatially-varying effects
svc.cols <- c(1, 2, 3)
p.svc <- length(svc.cols)
q.p.svc <- n.factors * p.svc
# Spatial decay parameters
phi <- runif(q.p.svc, 3 / 0.9, 3 / 0.1)
# A length N vector indicating the proportion of simulated locations
# that are within the range for a given species.
range.probs <- runif(N, 0.4, 1)
factor.model <- TRUE
cov.model <- 'spherical'
sp <- TRUE

dat <- simMsOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, N = N, beta = beta, alpha = alpha,
                psi.RE = psi.RE, p.RE = p.RE, phi = phi, sp = sp, svc.cols = svc.cols,
                cov.model = cov.model, n.factors = n.factors, 
                factor.model = factor.model, range.probs = range.probs)

y <- dat$y
X <- dat$X
X.re <- dat$X.re
X.p <- dat$X.p
X.p.re <- dat$X.p.re
coords <- dat$coords
range.ind <- dat$range.ind

# Prep data for spOccupancy -----------------------------------------------
# Occurrence covariates
occ.covs <- cbind(X, X.re)
colnames(occ.covs) <- c('int', 'occ.cov.1', 'occ.cov.2', 'occ.cov.3', 
                        'occ.cov.4', 'occ.factor.1')
# Detection covariates
det.covs <- list(det.cov.1 = X.p[, , 2], 
                 det.cov.2 = X.p[, , 3], 
                 det.factor.1 = X.p.re[, , 1]) 
# Data list
data.list <- list(y = y, coords = coords, occ.covs = occ.covs, 
                  det.covs = det.covs, range.ind = range.ind)
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
                   alpha.comm.normal = list(mean = 0, var = 2.72), 
                   tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
                   tau.sq.alpha.ig = list(a = 0.1, b = 0.1), 
                   phi.unif = list(a = 3 / 1, b = 3 / .1)) 
inits.list <- list(alpha.comm = 0, 
                   beta.comm = 0, 
                   beta = 0, 
                   alpha = 0,
                   tau.sq.beta = 1, 
                   tau.sq.alpha = 1, 
                   z = apply(y, c(1, 2), max, na.rm = TRUE)) 
# Tuning
tuning.list <- list(phi = 1)

# Number of batches
n.batch <- 2
# Batch length
batch.length <- 25
n.burn <- 0
n.thin <- 1
n.samples <- n.batch * batch.length

# Note that this is just a test case and more iterations/chains may need to 
# be run to ensure convergence.
out <- svcMsPGOcc(occ.formula = ~ occ.cov.1 + occ.cov.2 + occ.cov.3 + 
                                  occ.cov.4 + (1 | occ.factor.1),
                  det.formula = ~ det.cov.1 + det.cov.2 + (1 | det.factor.1),
                  data = data.list,
                  inits = inits.list,
                  n.batch = n.batch,
                  n.factors = n.factors,
                  batch.length = batch.length,
                  std.by.sp = TRUE,
                  accept.rate = 0.43,
                  priors = prior.list,
                  svc.cols = svc.cols,
                  cov.model = "spherical",
                  tuning = tuning.list,
                  n.omp.threads = 1,
                  verbose = TRUE,
                  NNGP = TRUE,
                  n.neighbors = 5,
                  search.type = 'cb',
                  n.report = 10,
                  n.burn = n.burn,
                  n.thin = n.thin,
                  n.chains = 1) 

summary(out)
}
