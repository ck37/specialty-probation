load_all_libraries = function() {
  # Hide all the extra messages displayed when loading these libraries.
  suppressMessages({
    library(foreach)      # For parallel looping
    library(doMC)         # For multicore
    library(doSNOW)       # For multinode
    library(RhpcBLASctl)  # Accurate physical core detection
    library(SuperLearner)
    library(xgboost)      # For GBM
#    library(DSA)          # Install this package from Mark's software webpage
    library(earth)        # For mars
    library(tmle)
    library(glmnet)       # For ridge/lasso/glm
    library(polspline)    # For polymars
    library(MASS)         # For stepAIC
    library(rpart)        # For rpartPrune
    library(e1071)        # For SVM
    library(randomForest)
    library(gam)
    library(xtable)
    library(arm)          # For bayesglm
  })
}

# Review distributions of variables that match a string.
review_by_string = function(data, string, ignore_case = T) {
  cols = grep(string, names(data), ignore.case=ignore_case)
  print(names(data)[cols])
  for (col in cols) {
    cat(names(data)[col], ":\n")
    print(summary(data[, col]))
  }
  return(invisible(names(data)[cols]))
}

# Create a new dataset with just the desired fields, which are also renamed.
clean_dataset = function(data, keep_fields) {
  # Create a blank data frame with the correct number of rows.
  clean_data = data.frame(matrix(nrow=nrow(data)))
  # Remove the weird default column that was created.
  clean_data[, 1] = NULL
  dim(clean_data)

  cat("Saving fields... ")
  keep_fields
  names(keep_fields)
  for (name in names(keep_fields)) {
    cat(name, " ")
    clean_data[, keep_fields[[name]]] = data[, name]
  }
  cat("\n")
  cat("New dimensions:", paste(dim(clean_data), collapse=" "), "\n")
  cat("New variables:", paste(names(clean_data), collapse=", "))
  clean_data
}

impute_missing_values = function(data, add_indicators = T, skip_vars = c()) {
  # Loop over each feature.
  missing_indicators = NULL
  for (i in 1:ncol(data)){
    nas = sum(is.na(data[, i]))
    # Nothing to impute, continue to next column.
    if (nas == 0 || names(data)[i] %in% skip_vars) {
      #cat("Skipping", i, "\n")
      next
    }
    #cat("Found", nas, "missing values for", colnames(data)[i], "\n")


    missing_indicator = matrix(1*is.na(data[, i]))
    new_name = paste0("miss_", colnames(data)[i])
    colnames(missing_indicator) = new_name
    next_col = ncol(missing_indicators) + 1
    if (is.null(missing_indicators)) {
      missing_indicators = matrix(missing_indicator, nrow=nrow(data), ncol=1)
      colnames(missing_indicators) = new_name
    } else {
      missing_indicators = cbind(missing_indicators, missing_indicator)
    }
    #new_colnames = c(colnames(missing_indicators), )
    #cat("New colnames:", new_colnames, "\n")
    #colnames(missing_indicators) = new_colnames

    if (class(data[, i]) == "factor") {
      # Impute factors to the mode.
      # Choose the first mode in case of ties.
      data[is.na(data[, i]), i] = Mode(data[, i])[1]
    } else {
      # Impute numeric values to the median.
      data[is.na(data[, i]), i] = median(data[, i], na.rm=T)
    }
  }
  if (add_indicators) {
    #cat("New indicator names:\n")
    data = cbind(data, missing_indicators)
  }
  data
}

# Via http://stackoverflow.com/questions/2547402/is-there-a-built-in-function-for-finding-the-mode
Mode = function(x) {
  ux = unique(x)
  tab <- tabulate(match(x, ux))
  # This will return multiple modes if they exist.
  ux[tab == max(tab)]
}

estimate_effect = function(Y, A, W,
       family = "binomial", cv_folds = 10, digits = 5,
       sl_lib = c("SL.glmnet", "SL.step", "SL.glm.interaction"),
       parallel = NULL, cluster = NULL, crossvalidate = F, outer_cv_folds = 10) {

  ##########
  # Setup parallel SL if we can.
  ##########
  if (is.null(parallel)) {
    sl_fn = SuperLearner
    if (crossvalidate) {
      cv_sl_fn = function(...) {
        CV.SuperLearner(..., V = outer_cv_folds)
      }
    }
  } else if (parallel == "multicore") {
    cat("Running SL via multicore\n")
    sl_fn = function(...) {
      mcSuperLearner(...)
    }
    if (crossvalidate) {
      cv_sl_fn = function(...) {
        CV.SuperLearner(..., V = outer_cv_folds, parallel = parallel)
      }
    }
  } else if (parallel == "snow") {
    cat("Running SL via snow\n")
    sl_fn = function(...) {
      snowSuperLearner(cluster, ...)
    }
    if (crossvalidate) {
      cv_sl_fn = function(...) {
        CV.SuperLearner(..., V = outer_cv_folds, parallel = cluster)
      }
    }
  }

  # Number of observations.
  n = nrow(W)

  # Allow custom # of CV folds.
  cv_ctrl = SuperLearner.CV.control(V = cv_folds)

  # Stacked dataframe with A = a, A = 1, and A = 0.
  x = data.frame(W, A=A)
  x1 = x0 = x
  x1$A = 1
  x0$A = 0
  newdata = rbind(x, x1, x0)


  qinit = sl_fn(Y=Y, X=x, newX=newdata, cvControl = cv_ctrl, SL.library=sl_lib, family=family)
  QbarAW = qinit$SL.predict[1:n]
  Qbar1W = qinit$SL.predict[(n+1):(2*n)]
  Qbar0W = qinit$SL.predict[(2*n+1):(3*n)]

  # Confirm that qinit is not all NA.
  stopifnot(!all(is.na(qinit$SL.predict)))

  if (crossvalidate) {
    # No newdata argument in this call.
    qinit_cv = cv_sl_fn(Y=Y, X=x, cvControl = cv_ctrl, SL.library=sl_lib, family=family)
  }

  #######
  # Substitution estimator (g-computation).
  #######

  psihat_ss = mean(Qbar1W - Qbar0W)
  cat("Psihat SS:", round(psihat_ss, digits), "\n")

  # Stop early if psihat_ss is NA.
  stopifnot(!is.na(psihat_ss))

  gHatSL = sl_fn(Y=A, X=W, SL.library=sl_lib, cvControl = cv_ctrl, family="binomial")
  gHat1W = gHatSL$SL.predict
  gHat0W = 1 - gHat1W
  # Check if gHat1W is all NAs. If so, the model fitting failed.
  stopifnot(mean(is.na(gHat1W)) != 1)

  gHatAW = rep(NA, n)
  gHatAW[A == 1] = gHat1W[A == 1]
  gHatAW[A == 0] = gHat0W[A == 0]

  if (crossvalidate) {
    gHatSL_cv = cv_sl_fn(Y=A, X=W, SL.library=sl_lib, cvControl = cv_ctrl, family="binomial")
  }

  #########
  # IPTW estimator (not stabilized)
  #########

  # Classic version:
  psihat_iptw = mean((A == 1)*Y / gHatAW) - mean((A == 0)*Y / gHatAW)
  cat("Psihat IPTW:", round(psihat_iptw, digits), "\n")
  # Confirm that psihat_iptw is not NA.
  stopifnot(!is.na(psihat_iptw))

  # Stabilized version:
  wgt = 1 / gHatAW
  psihat_iptw_ht = mean((A == 1) * wgt * Y) / mean((A == 1) * wgt) -
    mean((A == 0) * wgt * Y) / mean((A == 0) * wgt)

  cat("Psihat IPTW HT:", round(psihat_iptw_ht, digits), "\n")


  #########
  # Clever covariate calculations.
  #########
  # Set clever covariate to -1/gHatAW or 1/gHatAW
  h_aw = (2*A - 1)/gHatAW

  h_1w = 1/gHat1W
  # Note: make sure to specify -1 here:
  h_0w = -1/gHat0W

  logitUpdate = glm(Y ~ -1 + offset(qlogis(QbarAW)) + h_aw, family=family)
  eps = coef(logitUpdate)
  QbarAW_star = plogis(qlogis(QbarAW) + eps * h_aw)
  Qbar1W_star = plogis(qlogis(Qbar1W) + eps * h_1w)
  Qbar0W_star = plogis(qlogis(Qbar0W) + eps * h_0w)

  ##########
  # TMLE estimate
  ##########
  psihat_tmle = mean(Qbar1W_star - Qbar0W_star)
  cat("Psihat TMLE:", round(psihat_tmle, digits), "\n")

  ###########
  # Efficient inference (Lab 6)
  ###########

  # Calculate inference curve.
  ic = h_aw*(Y - QbarAW_star) + Qbar1W_star - Qbar0W_star - psihat_tmle

  # SE of the influence curve.
  ic_se = sqrt(sum(var(ic)) / n)

  # Confidence interval of the estimate.
  ci = psihat_tmle + c(-1, 1) * ic_se * 1.96

  # P-value
  tmle_p = 2 * pnorm(-abs(psihat_tmle / ic_se), lower.tail=T)

  #########
  # Finalization.

  cat("Max g-weight:", round(max(wgt), 1), "\n")

  # Return the results.
  results = list(qinit = qinit, ghat = gHatSL, influence_curve = ic,
                 psihat_ss = psihat_ss, psihat_iptw = psihat_iptw,
                 psihat_iptw_ht = psihat_iptw_ht, psihat_tmle = psihat_tmle,
                 tmle_se = ic_se, tmle_ci = ci, tmle_p = tmle_p, weights=wgt)

  if (crossvalidate) {
    results = c(results, list(qinit_cv = qinit_cv, ghat_cv = gHatSL_cv))
  }

  class(results) = "estimate_effect"
  results
}

# Custom result printing.
# TODO: make prettier.
print.estimate_effect = function(obj, digits = 4) {
  cat("Simple substitution estimate:", round(obj$psihat_ss, digits), "\n")
  cat("IPTW estimate:", round(obj$psihat_iptw, digits), "\n")
  cat("IPTW HT estimate:", round(obj$psihat_iptw_ht, digits), "\n")
  cat("TMLE estimate:", round(obj$psihat_tmle, digits), "\n")
  cat("TMLE SE:", paste(round(obj$tmle_se, digits)), "\n")
  cat("TMLE CI:", paste(round(obj$tmle_ci, digits)), "\n")
  cat("TMLE p:", round(obj$tmle_p, digits+2), "\n")
}

# Setup parallel processing, either multinode or multicore.
# By default it uses a multinode cluster if available, otherwise sets up multicore via doMC.
# Libraries required: parallel, doSNOW, doMC, RhpcBLASctl, foreach
setup_parallelism = function(conf = NULL, type="either", allow_multinode = T,
                             machine_list = Sys.getenv("SLURM_NODELIST"),
                             cpus_per_node = as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE")),
                             outfile = "") {
  # Indicator for multi-node parallelism.
  multinode = F

  # Check if we have multiple machine access.
  if (allow_multinode) {
    machines = strsplit(machine_list, ",")[[1]]
    if (length(machines) > 1) {
      cat("Have multi-node access for parallelism with", length(machines), "machines:", machines, "\n")
      # NOTE: this may be a bad config if the nodes have different core counts.
      cores = rep(machines, each = cpus_per_node)
      cat("Multi-node cores enabled:", cores, "\n")
      multinode = T
    }
  }

  if (!multinode) {
    # Count of physical cores, unlike parallel:detectCores() which is logical cores (threads).
    cores = RhpcBLASctl::get_num_cores()
    cat("Local physical cores detected:", cores, "\n")

    if (exists("conf") && !is.null(conf) && "num_cores" %in% names(conf)) {
      cores = conf$num_cores
      cat("Using", cores, " local cores due to conf settings.\n")
    }
  }

  if (multinode || type %in% c("cluster", "doParallel")) {
    # Outfile = "" allows output from within foreach to be displayed while in RStudio.
    # TODO: figure out how to suppress the output from makeCluster()
    #capture.output({ cl = parallel::makeCluster(cores, outfile = outfile) })
    cat("Starting multinode cluster with cores:", cores, "\n")
    capture.output({ cl = parallel::makeCluster(cores, type="PSOCK", outfile = outfile) })
    # doParallel supports multicore and multinode parallelism, but may require
    # explicitly exporting functions and libraries across the cluster.
    registerDoSNOW(cl)
    setDefaultCluster(cl)
    parallel_type = "snow"
  } else {
    # doMC only supports multicore parallelism, not multi-node.
    registerDoMC(cores)
    cl = NA
    # TODO: is this needed since we've already done registerDoMC()?
    options(mc.cores = cores)
    parallel_type = "multicore"
  }

  # Make sure the BLAS is not competing with the SL parallelism.
  omp_threads = omp_get_max_threads()
  if (!is.null(omp_threads) && omp_threads > 1) {
    omp_set_num_threads(1)
  }

  # TODO: need to figure out difference between get_max_threads and get_num_procs.
  # They are not always both consistently set to 1 (i.e. on Benten).
  omp_threads = omp_get_num_procs()
  # If omp_get_num_procs() returns NULL we can safely plan on using 1 thread.
  omp_threads = ifelse(is.null(omp_threads), 1, omp_threads)
  cat("Our BLAS is setup for", blas_get_num_procs(), "threads and OMP is", omp_threads, "threads.\n")
  #cat("Multicore parallel is setup to use", getOption("mc.cores"), "cores.\n")

  cat("doPar backend registered:", getDoParName(), "\n")
  cat("Workers enabled:", getDoParWorkers(), "\n")
  # Return invisibily so that NA is not printed.
  return(invisible(cl))
}

# Stop the cluster if parallel:makeCluster() was used, but nothing needed if doMC was used.
stop_cluster = function(cluster_obj) {
  # Check if this cluster was created using parallel:makeCluster
  if (inherits(cluster_obj, "cluster")) {
    stopCluster(cluster_obj)
  } else {
    #cat("No cluster shutdown required.\n")
    invisible()
  }
}

### From Courtney's create_SL_library.R
# TODO: adjust configuration so that DSA is not as incredibly slow.
# CK: we are disabling substitution and deletion to try to speed up the algorithm.
# We are also reducign maxize from 2 * ncol(X) to 0.5 * ncol(X).
# May also want to add rank.cutoffs.
SL.DSA <- function(Y, X, newX, family, obsWeights, maxsize = 5, maxorderint = 2, maxsumofpow = 2, Dmove = T, Smove = T, vfold = 5, ...) {
  require('DSA')
  dsaweights <- matrix(obsWeights, nrow = (vfold +1), ncol = nrow(X), byrow = TRUE)
  fit.DSA <- DSA(Y ~ 1, data = data.frame(Y, X), family = family, maxsize = maxsize, maxorderint = maxorderint, maxsumofpow = maxsumofpow, Dmove = Dmove, Smove = Smove, vfold = vfold, weights = dsaweights)
  pred <- predict(fit.DSA, newdata = newX)
  if(family$family == "binomial") { pred <- 1 / (1 + exp(-pred))}
  fit <- list(object = fit.DSA)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.DSA")
  return(out)
}

#
predict.SL.DSA <- function(object, newdata, family, ...) {
  require('DSA')
  pred <- predict(object = object$object, newdata = newdata)
  if (family$family == "binomial") {
    pred <- 1 / (1 + exp(-pred))
  }
  return(pred)
}


create.SL.glmnet <- function(alpha = c(0, 0.25, 0.50, 0.75, 1)) {
  # TODO: don't use global vars here.
  for(mm in seq(length(alpha))){
    eval(parse(text = paste('SL.glmnet.', alpha[mm], '<- function(..., alpha = ', alpha[mm], ') SL.glmnet(..., alpha = alpha)', sep = '')), envir = .GlobalEnv)
  }
  invisible(TRUE)
}

SL.xgboost = function(Y, X, newX, family, obsWeights, id, ntrees = 1000,
                      max_depth=4, shrinkage=0.1, minobspernode=10,
                      early.stop.round=NULL, ...) {
  # TODO: Convert to an xgboost compatible data matrix, using the sample weights.
  #xgmat = xgb.DMatrix(as.matrix(X), label=Y, weight = obsWeights)
  # We are not using sample weights currently:
  xgmat = xgb.DMatrix(data=as.matrix(X), label=Y)

  # TODO: support early stopping, which requires a watchlist. See ?xgb.train

  if (family$family == "gaussian") {
    # We restrict xgboost to 1 thread so that SL can control parallelization.
    model = xgboost(data=xgmat, objective="reg:linear", nround = ntrees, max_depth = max_depth, minchildweight = minobspernode, eta = shrinkage, early.stop.round=early.stop.round, verbose=0, nthread = 1)
  }
  if (family$family == "binomial") {
    # We restrict xgboost to 1 thread so that SL can control parallelization.
    model = xgboost(data=xgmat, objective="binary:logistic", nround = ntrees, max_depth = max_depth, minchildweight = minobspernode, eta = shrinkage, early.stop.round = early.stop.round, verbose=0, nthread = 1)
  }
  if (family$family == "multinomial") {
    # TODO: test this.
    model = xgboost(data=xgmat, objective="multi:softmax", nround = ntrees, max_depth = max_depth, minchildweight = minobspernode, eta = shrinkage, num_class=length(unique(Y)), early.stop.round = early.stop.round, verbose=0, nthread = 1)
  }
  pred = predict(model, newdata=data.matrix(newX))
  fit = vector("list", length = 0)
  class(fit) = c("SL.xgboost")
  out = list(pred = pred, fit = fit)
  return(out)
}

create.SL.xgboost = function(tune = list(ntrees = c(1000), max_depth = c(4), shrinkage = c(0.1), minobspernode = c(10)), detailed_names = F) {
  # Create all combinations of hyperparameters, for grid-like search.
  tuneGrid = expand.grid(tune, stringsAsFactors=F)

  names = rep("", nrow(tuneGrid))

  name_prefix = "SL.xgb"
  for (i in seq(nrow(tuneGrid))) {
    g = tuneGrid[i,]
    # TODO: update this to not use the global environment by default.
    if (detailed_names) {
      name = paste(name_prefix, g$ntrees, g$max_depth, g$shrinkage, g$minobspernode, sep=".")
    } else {
      name = paste(name_prefix, i, sep=".")
    }
    names[i] = name
    eval(parse(text = paste0(name, "= function(..., ntrees = ", g$ntrees, ", max_depth = ", g$max_depth, ", shrinkage=", g$shrinkage, ", minobspernode=", g$minobspernode, ") SL.xgboost(..., ntrees = ntrees, max_depth = max_depth, shrinkage=shrinkage, minobspernode=minobspernode)")), envir = .GlobalEnv)
  }
  results = list(grid = tuneGrid, names = names)
  invisible(results)
}


create_SL_lib = function(num_cols = NULL, xgb = T, rf = T, dsa = F, glmnet = T, gam=T, glmnet_size = 11, detailed_names = F) {

  glmnet_libs = c()
  if (glmnet) {
    alpha_params = seq(0, 1, length.out=glmnet_size)
    # TODO: don't use global vars here.
    create.SL.glmnet(alpha = alpha_params)
    glmnet_libs = paste0("SL.glmnet.", alpha_params)
    cat("Glmnet:", length(glmnet_libs), "configurations. Alphas:", paste(alpha_params, collapse=", "), "\n")
  }

  # Create xgboost models.
  xgb_libs = c()
  xgb_grid = NA
  if (xgb) {

    # Slower, ideal configuration search (intended for servers).
    xgb_tune = list(ntrees = c(50, 200, 500, 1000, 3000), max_depth = c(1, 2, 3), shrinkage = c(0.01, 0.1, 0.2, 0.4), minobspernode = c(10))

    # Faster, less ideal configuration search (intended for laptops):
    # BUT disable for now - we have so few observations that we can use the server version.
    if (F & get_num_cores() < 8) {
      xgb_tune = list(ntrees = c(1000, 2000), max_depth = c(1, 2), shrinkage = c(0.1, 0.2), minobspernode = c(10))
    }

    # TODO: don't use global vars here.
    xgb_results = create.SL.xgboost(xgb_tune, detailed_names = detailed_names)
    xgb_grid = xgb_results$grid
    xgb_libs = xgb_results$names

    cat("XGBoost:", length(xgb_libs), "configurations.\n")
    print(xgb_grid)
  }

  rf_libs = c()
  rf_grid = NA
  if (rf) {
    if (!is.null(num_cols)) {
      # Much better is to send in how many columns are in the dataset.
      rf_tune = list(mtry = unique(round(exp(log(num_cols)*exp(c(-0.96, -0.71, -0.48, -0.4, -0.29, -0.2))))), nodesize = 1)
    } else {
      rf_tune = list(mtry = c(1, 5, 10), nodesize = 1)
    }

    rf_models = create.SL.randomForest(rf_tune, detailed_names = detailed_names)

    rf_libs = rf_models$names
    rf_grid = rf_models$grid

    cat("Random Forest:", length(rf_libs), "configurations.\n")
    print(rf_grid)
  }

  gam_libs = c()
  if (gam) {
    gam_degrees = c(2, 3, 4)
    gam_models = create.SL.gam(gam_degrees)
    gam_libs = gam_models$names
  }

  # TODO: see if we want to tweak the hyperparameters of any of these singular models.
  # Remove "SL.polymars", for now. (CK 5/4/16)
  lib = c(glmnet_libs, xgb_libs, rf_libs, "SL.svm", "SL.bayesglm",
          "SL.stepAIC", "SL.earth", "SL.rpartPrune", gam_libs)

  if (dsa) {
    # WARNING: super duper slow :/
    lib = append(lib, "SL.DSA")
  }

  cat("Total models:", length(lib), "\n")
  results = list(lib = lib, xgb_grid = xgb_grid, xgb_libs = xgb_libs,
                 glmnet_libs = glmnet_libs, rf_grid = rf_grid, rf_libs = rf_libs)
  results
}

create.SL.randomForest <- function(tune = list(mtry = c(1, 5, 10), nodesize = c(1, 5, 10)), detailed_names = F) {
  tuneGrid <- expand.grid(tune, stringsAsFactors = FALSE)
  names = rep("", nrow(tuneGrid))
  name_prefix = "SL.rf"
  for (mm in seq(nrow(tuneGrid))) {
    g = tuneGrid[mm,]
    if (detailed_names) {
      name = paste(name_prefix, tuneGrid[mm, 1], tuneGrid[mm, 2], sep=".")
    } else {
      name = paste(name_prefix, mm, sep=".")
    }
    names[mm] = name
    eval(parse(file = "", text = paste0(name, "<- function(..., mtry = ", g$mtry, ", nodesize = ", g$nodesize, ") SL.randomForest(..., mtry = mtry, nodesize = nodesize)")), envir = .GlobalEnv)
  }
  results = list(grid = tuneGrid, names = names)
  invisible(results)
}

# Review meta-weights from a CV.SuperLearner
cvSL_review_weights = function(cv_sl) {
  meta_weights = coef(cv_sl)
  means = colMeans(meta_weights)
  sds = apply(meta_weights, MARGIN=2, FUN=function(col) { sd(col) })
  mins = apply(meta_weights, MARGIN=2, FUN=function(col) { min(col) })
  maxs = apply(meta_weights, MARGIN=2, FUN=function(col) { max(col) })
# Combine the stats into a single matrix.
  sl_stats = cbind("mean(weight)"=means, "sd(weight)"=sds, "min(weight)"=mins, "max(weight)"=maxs)
  sl_stats
}

# Modified from SuperLearnerExtra:
# Creates gam wrappers in the global environment with different degrees.
# The default value for deg.gam in SL.gam is 2
create.SL.gam <- function(degree = c(3, 4)) {
  names = rep("", length(degree))
  for(mm in seq(length(degree))){
    name = paste0('SL.gam.', degree[mm])
    names[mm] = name
    eval(parse(text = paste0(name, '<- function(..., deg.gam = ', degree[mm], ') SL.gam(..., deg.gam = deg.gam)')), envir = .GlobalEnv)
  }
  results = list(names = names, degree = degree)
  invisible(results)
}

# Remove constant columns in the dataframe (sd = 0).
remove_constant_columns = function(data) {
  # Remove columns that are constant (sd = 0). (ignore warning about NA)
  suppressWarnings({
    col_sds = apply(data, MARGIN=2, FUN=sd)
  })
  # We need the !is.na() so columns with NA for sd (e.g. character cols) aren't included.
  bad_cols = col_sds == 0 & !is.na(col_sds)
  # TODO: handle situation where dataframe has no column names.
  bad_col_names = colnames(data)[bad_cols]
  if (sum(bad_cols) > 0) {
    cat("Constant columns:", bad_col_names, "\n")
    data = data[, !colnames(data) %in% bad_col_names]
  }
  data
}
