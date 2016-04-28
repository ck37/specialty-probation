load_all_libraries = function() {
  # Hide all the extra messages displayed when loading these libraries.
  suppressMessages({
    library(foreach)   # For multicore
    library(doMC)       # For multicore (not used at the moment)
    library(doSNOW) # For multicore
    library(RhpcBLASctl) # Accurate physical core detection
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
       family = "binomial",
       sl_lib = c("SL.glmnet", "SL.step", "SL.glm.interaction"),
       parallel = NULL,
       cluster = NULL) {
  # Recreate dataframe based on Y, A, W.
  # TODO: not sure that we even need this.
  data = data.frame(Y=Y, A=A, W)

  n = nrow(W)

  # Create x df without the outcome variable.
  # x = subset(data, select=-Y)
  # TODO: double-check this.
  x = data.frame(W, A=A)
  names(x)
  x1 = x0 = x
  x1$A = 1
  x0$A = 0
  newdata = rbind(x, x1, x0)
  dim(newdata)

  # Use parallel SL if we can.
  if (is.null(parallel)) {
    sl_fn = SuperLearner
  } else if (parallel == "multicore") {
    cat("Running SL via multicore\n")
    sl_fn = function(...) {
      mcSuperLearner(...)
    }
  } else if (parallel == "snow") {
    cat("Running SL via snow\n")
    sl_fn = function(...) {
      snowSuperLearner(cluster, ...)
    }
  }

  qinit = sl_fn(Y=Y, X=x, newX=newdata, SL.library=sl_lib, family=family)
  qinit
  QbarAW = qinit$SL.predict[1:n]
  Qbar1W = qinit$SL.predict[(n+1):(2*n)]
  Qbar0W = qinit$SL.predict[(2*n+1):(3*n)]
  # tail(cbind(data$A, QbarAW, Qbar1W, Qbar0W))

  psihat_ss = mean(Qbar1W - Qbar0W)

  gHatSL = sl_fn(Y=A, X=W, SL.library=sl_lib, family="binomial")
  gHat1W = gHatSL$SL.predict
  gHat0W = 1 - gHat1W

  gHatAW = rep(NA, n)
  gHatAW[A == 1] = gHat1W[A == 1]
  gHatAW[A == 0] = gHat0W[A == 0]

  psihat_iptw = mean((A == 1)*Y / gHatAW) - mean((A == 0)*Y / gHatAW)

  h_aw = (A == 1)/gHat1W - (A == 0)/gHat0W
  # Same as?
  # h_aw = (A == 1)/gHatAW - (A == 0)/gHatAW
  # Yup:
  # mean(h_aw == h_aw2)

  h_1w = 1/gHat1W
  # Note: make sure to specify -1 here:
  h_0w = -1/gHat0W

  logitUpdate = glm(Y ~ -1 + offset(qlogis(QbarAW)) + h_aw, family=family)
  eps = coef(logitUpdate)
  QbarAW_star = plogis(qlogis(QbarAW) + eps * h_aw)
  Qbar1W_star = plogis(qlogis(Qbar1W) + eps * h_1w)
  Qbar0W_star = plogis(qlogis(Qbar0W) + eps * h_0w)

  psihat_tmle = mean(Qbar1W_star - Qbar0W_star)

  # Inference (Lab 6)

  ic = h_aw*(Y - QbarAW_star) + Qbar1W_star - Qbar0W_star - psihat_tmle

  ic_se = sqrt(sum(var(ic)) / n)

  ci = psihat_tmle + c(-1, 1) * ic_se

  # 2 * pnorm(psihat_tmle / ic_se, lower.tail=F)
  tmle_p = 2 * pnorm(-abs(psihat_tmle / ic_se), lower.tail=T)

  # Return the results.
  results = list(qinit = qinit, ghat = gHatSL, influence_curve = ic,
                 psihat_ss = psihat_ss, psihat_iptw = psihat_iptw,
                 psihat_tmle = psihat_tmle,
                 tmle_se = ic_se, tmle_ci = ci, tmle_p = tmle_p)
  class(results) = "tmle"
  results
}

# Custom result printing.
# TODO: make prettier.
print.tmle = function(obj) {
  cat("Simple substitution estimate:", obj$psihat_ss, "\n")
  cat("IPTW estimate:", obj$psihat_iptw, "\n")
  cat("TMLE estimate:", obj$psihat_tmle, "\n")
  cat("TMLE CI:", paste(obj$tmle_ci), "\n")
  cat("TMLE p:", obj$tmle_p, "\n")
}

# Setup parallel processing, either multinode or multicore.
# By default it uses a multinode cluster if available, otherwise sets up multicore via doMC.
# Libraries required: parallel, doParallel, doMC, RhpcBLASctl, foreach
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
    registerDoParallel(cl)
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
