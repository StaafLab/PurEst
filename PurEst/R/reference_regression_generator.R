# reference_regression_generator

reference_regression_generator <- function(

  beta_values, # Matrix of beta values of CoGs of certain samples
  tumour_purities, # Purity values of the samples to be analysed
  set_seed = FALSE, # Boolean to specify if a seed should be set
  seed_num = 2000, # Seed to be used if set_seed == TRUE
  cores = 1, # Cores to be used to run the function
  regression_details = FALSE # Include the corrected betas based on the Staaf Aine method and the populations identified for each sample in each CpG in the output

) {

  # ===========================
  # CONFIGURING PARALLELIZATION
  # ===========================

  cat("\nUsing", cores,"core(s)\n\n")

  #Creating the cluster to run the process in parallel
  cl <- makeCluster(cores)
  registerDoParallel(cl)

  #Making sure that all packages have access to the flexmix package. Using invisible()
  #to avoid printing anything to the terminal.
  invisible(clusterEvalQ(cl, {library("flexmix")}))


  # ================================
  # DETERMINING VARIANCE OF EACH CPG
  # ================================

  #Create a vector with the variance of each cpg (row)
  cpg_variance <- apply(beta_values, 1, var)
  names(cpg_variance) <- rownames(beta_values)


  # =======================
  # CALCULATING REGRESSIONS
  # =======================

  cat("\nCalculating reference regressions\n\n")

  # Setting seed if necessary
  if (set_seed == TRUE) {

    #Adding seed to each row of the beta value dataframe
    betaRun <- cbind(seed=seed_num:seed_num+nrow(beta_values), beta_values)

  } else {

    betaRun <- beta_values

  }

  #Storing sample names
  betaNames <- colnames(beta_values)

  # Initializing progress bar and specifying options
  pbo <- pboptions(type = "txt", char="=", txt.width=80)

  #Running the analysis in parallel with a progress bar
  res <- pbapply(cl = cl, #ClusterS to run the process
                 MARGIN = 1, #Apply the function to the rows
                 FUN = adjustBeta, #Function to correct betas
                 purity=tumour_purities, #Purity values
                 snames=betaNames, #Sample names
                 seed=set_seed, #Specify if the seed has been added to the data or not
                 betaRun) #Beta values+the added seed

  # Stop clusters used in parallelization
  stopCluster(cl)


  # =================
  # GENERATING OUTPUT
  # =================

  cat("\nGenerating output...\n\n")


  # Generating output list. Short and extended versions
  if (regression_details == TRUE) {

    # Creating a list to add the results.
    result_list <- list(
      betas.original = do.call("rbind",lapply(res,function(x) x$y.orig)), #Original beta values
      betas.tumor = do.call("rbind",lapply(res,function(x) x$y.tum)), #Corrected tumor beta values
      betas.microenvironment = do.call("rbind",lapply(res,function(x) x$y.norm)), #Corrected microenvironment beta values
      cpg.populations =  do.call("rbind",lapply(res,function(x) x$groups)), #Methylation patterns (populations) of each CpG
      cpg.variance = cpg_variance, # Variance value of the beta values of each CpG
      reg.slopes = do.call("rbind",lapply(res,function(x) x$res.slopes)), #Slopes of the populations
      reg.intercepts = do.call("rbind",lapply(res,function(x) x$res.int)), #Intercepts of the populations
      reg.RSE = do.call("rbind",lapply(res,function(x) x$res.rse)), #Residual standard error
      reg.df = do.call("rbind",lapply(res,function(x) x$res.df)) #Degrees of freedom of the reversed regressions
    )

  } else {

    # Creating a list to add the results.
    result_list <- list(
      cpg.variance = cpg_variance, # Variance value of the beta values of each CpG
      reg.slopes = do.call("rbind",lapply(res,function(x) x$res.slopes)), #Slopes of the populations
      reg.intercepts = do.call("rbind",lapply(res,function(x) x$res.int)), #Intercepts of the populations
      reg.RSE = do.call("rbind",lapply(res,function(x) x$res.rse)), #Residual standard error
      reg.df = do.call("rbind",lapply(res,function(x) x$res.df)) #Degrees of freedom of the reversed regressions
    )

  }


  cat("\n\n=================\n")
  cat ("PROCESS FINISHED")
  cat("\n=================\n\n")


  # Returning output
  return(result_list)

}
