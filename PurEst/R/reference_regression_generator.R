# reference_regression_generator

reference_regression_generator <- function(

  beta_values, # Matrix of beta values of CoGs of certain samples
  tumour_purities, # Purity values of the samples to be analysed
  set_seed = FALSE, # Boolean to specify if a seed should be set
  seed = 2000, # Seed to be used if set_seed == TRUE
  cores = 1 # Cores to be used to run the function

) {

  # ===========================
  # CONFIGURING PARALLELIZATION
  # ===========================

  cat("\nUsing", cores,"core(s)\n\n")

  #Creating the cluster to run the process in parallel
  cl <- makeCluster(arguments$cores)
  registerDoParallel(cl)

  #Making sure that all packages have access to the flexmix package. Using invisible()
  #to avoid printing anything to the terminal.
  invisible(clusterEvalQ(cl, {library("flexmix")}))


}
