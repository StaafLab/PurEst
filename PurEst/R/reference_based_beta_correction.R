# Reference based beta correction

reference_based_beta_correction <- function(

  betas_to_correct,
  purities_samples_to_correct,
  only_certain_CpGs = FALSE,
  CpGs_to_correct_vector,
  refitting,
  reference_regressions,
  reference_betas,
  reference_purities,
  use_seed = TRUE,
  seed_num = 2000,
  cores = 1

) {


  # ===========================
  # CONFIGURING PARALLELIZATION
  # ===========================

  cat("\nUsing", cores,"core(s)\n\n")

  #Creating the cluster to run the process in parallel
  cl <- makeCluster(cores)


  # =============================================
  # CORRECT BETAS REFITTING REFERENCE REGRESSIONS
  # =============================================

  if (refitting == TRUE) {

    #Making sure that all cores have access to the flexmix package. Using invisible()
    #to avoid printing anything to the terminal
    registerDoParallel(cl)
    invisible(clusterEvalQ(cl, {library("flexmix")}))

    # PROCESSING PREDICTED PURITIES

    cat("\nPreprocessing the data...\n\n")

    # Getting predicted 1-purities from PurEst output
    predicted_1mPurities <- purities_samples_to_correct$`Estimated_1-Purities`

    # Removing samples with more than one estimates (if any)
    if (nrow(predicted_1mPurities[which(predicted_1mPurities[,2]!=1),])!=0) {

      #Calculate the number of samples to remove
      samples_to_remove <- nrow(predicted_1mPurities[which(predicted_1mPurities[,2]!=1),]) / 2

      #Print warining message
      cat("\n", samples_to_remove, "samples have more than one predicted purity. Samples removed from the beta correction.\n")

      #Filtering samples with more than one purity values
      predicted_1mPurities <- predicted_1mPurities[which(predicted_1mPurities[,2]==1),]
    }

    # Transforming the predicted_purities dataframe into a vector
    predicted_purities_vec <- 1 - as.numeric(predicted_1mPurities[,3]) # Using purity, no 1-Purity
    names(predicted_purities_vec) <- predicted_1mPurities[,1]


    # FILTERING CPGS

    cat("\nChecking CpGs to be corrected...\n\n")

    # Use only the specified CpGs if that option has been selected
    if (only_certain_CpGs) {

      #Keeping only CpGs of interest
      betas_to_correct <- betas_to_correct[CpGs_to_correct_vector,]

    }


    # Checking if the CpGs are included in the reference data
    if (sum(!(rownames(betas_to_correct) %in% rownames(reference_betas))) != 0) {

      # Printing warning message
      cat("\n",  sum(!(rownames(betas_to_correct) %in% rownames(reference_betas))), "CpG(s) is/are not included into the reference cohort, so it/they can not be corrected.\n\n")

      # Filtering not included CpGs
      betas_to_correct <- t[rownames(betas_to_correct) %in% rownames(reference_betas),]
    }

    # Remove CpGs from the cohort dataset that are not included into the data to correct to speed up the process.
    reference_betas <- reference_betas[rownames(reference_betas) %in% rownames(betas_to_correct),]

    #Sorting the cohort betas dataframe based on the rownames of betas_to_correct
    reference_betas <- reference_betas[rownames(betas_to_correct),]


    # MERGING REFERENCE AND PREDICTED VALUES TO REFIT THE REFERENCE REGRESSIONS

    # Creating a single purity vector
    purities <- c(reference_purities, predicted_purities_vec)

    # Creating a single betas dataframe
    betas <- cbind(reference_betas, betas_to_correct)

    #Removing sample purities not included into the beta dataset. It generates errors
    purities <- purities[colnames(betas)]


    # RUNNING BETA CORRECTION

    cat("\nCorrecting betas refitting the reference regressions...\n\n")

    #Acding seed if necessary
    if (use_seed) {

      #Adding seed to each row of the beta value dataframe
      betaRun <- cbind(seed=seed_num:seed_num+nrow(betas),betas)

    } else {

      betaRun <- betas

    }

    #Storing sample names
    betaNames <- colnames(betas)

    # Initializing progress bar and specifying options
    pbo <- pboptions(type = "txt", char="=", txt.width=80)

    #Running the analysis in parallel with a progress bar
    res <- pbapply(cl = cl, #ClusterS to run the process
                   MARGIN = 1, #Apply the function to the rows
                   FUN = adjustBeta, #Function to correct betas
                   purity=purities, #Purity values
                   snames=betaNames, #Sample names
                   seed=use_seed, #Specify if the seed has been added to the data or not
                   betaRun #Beta values+the added seed
    )

    # Stop clusters used in parallelization
    stopCluster(cl)


    # GENERATING RESULT LIST

    cat("\n\nGenerating output...\n\n")


    # Creating a list to add the results
    result_list <- list(
      betas.original = do.call("rbind",lapply(res,function(x) x$y.orig)), #Original beta values
      betas.tumor = do.call("rbind",lapply(res,function(x) x$y.tum)), #Corrected tumor beta values
      betas.microenvironment = do.call("rbind",lapply(res,function(x) x$y.norm)) #Corrected microenvironment beta values
    )

    # Creating a list to add the parameters of the correction regressions
    reg_list <- list(
      cpg.populations =  do.call("rbind",lapply(res,function(x) x$groups)), #Methylation patterns (populations) of each CpG
      reg.slopes = do.call("rbind",lapply(res,function(x) x$res.slopes)), #Slopes of the populations
      reg.intercepts = do.call("rbind",lapply(res,function(x) x$res.int)), #Intercepts of the populations
      reg.RSE = do.call("rbind",lapply(res,function(x) x$res.rse)), #Residual standard error
      reg.df = do.call("rbind",lapply(res,function(x) x$res.df)) #Degrees of freedom of the regressions
    )

    # Filtering results. Keeping only CpGs that were intended to be corrected
    lapply(names(result_list), function(n) {
      result_list[[n]][,names(predicted_purities_vec)]
    })



    cat("\n=================\n")
    cat ("PROCESS FINISHED")
    cat("\n=================\n")

    # Return output
    return(
      list(
        "Corrected_betas" = result_list,
        "Regression_parameters" = reg_list
      )
    )

  # =====================================================
  # CORRECT BETAS WITHOUT REFITTING REFERENCE REGRESSIONS
  # =====================================================


  } else {

    registerDoSNOW(cl)
    clusterExport(cl, c("identify_regression", "correcting_betas"))

    # PROCESSING PREDICTED PURITIES

    cat("\nPreprocessing the data...\n\n")

    # Getting predicted 1-purities from PurEst output
    predicted_1mPurities <- purities_samples_to_correct$`Estimated_1-Purities`


    # Removing samples with more than one estimates (if any)
    if (nrow(predicted_1mPurities[which(predicted_1mPurities[,2]!=1),])!=0) {

      #Calculate the number of samples to remove
      samples_to_remove <- nrow(predicted_1mPurities[which(predicted_1mPurities[,2]!=1),]) / 2

      #Print warining message
      cat("\n", samples_to_remove, "samples have more than one predicted purity. Samples removed from the beta correction.\n")

      #Filtering samples with more than one purity values
      predicted_1mPurities <- predicted_1mPurities[which(predicted_1mPurities[,2]==1),]
    }


    # Transforming the predicted_1mPurities dataframe into a vector
    predicted_1mPurities_vec <- predicted_1mPurities[,3]
    names(predicted_1mPurities_vec) <- predicted_1mPurities[,1]

  }


  # Reformatting reference regression parametersç
  my_slopes <- reference_regressions$reg.slopes
  my_intercepts <- reference_regressions$reg.intercepts
  my_RSE <- reference_regressions$reg.RSE
  my_df <- reference_regressions$reg.df


  # FILTERING CPGS

  cat("\nChecking CpGs...\n\n")

  # Use only the specified CpGs if that option has been selected
  if (only_certain_CpGs) {

    #Keeping only CpGs of interest
    betas_to_correct <- betas_to_correct[CpGs_to_correct_vector,]

  }


  # Checking if the CpGs are included in the reference regressions
  if (sum(!(rownames(betas_to_correct) %in% rownames(my_slopes))) != 0) {

    # Printing warning message
    cat("\n",  sum(!(rownames(betas_to_correct) %in% rownames(my_slopes))), "CpG(s) is/are not included into the refernce cohort, so it/they can not be corrected.\n\n")

    # Filtering not included CpGs
    betas_to_correct <- betas_to_correct[rownames(betas_to_correct) %in% rownames(my_slopes),]
  }

  # Remove CpGs from the regressions that are not included into the data to correct to speed up the process.
  my_slopes <- my_slopes[rownames(my_slopes) %in% rownames(betas_to_correct),]
  my_intercepts <- my_intercepts[rownames(my_intercepts) %in% rownames(betas_to_correct),]


  # CORRECTING BETAS BASED ON REFERENCE REGRESSIONS

  cat("\nCorrecting betas without refitting reference regressions...\n\n")


  # Configure progress bar
  p_bar <- txtProgressBar(min=0,
                          max=nrow(betas_to_correct),
                          style=3,
                          width=80)


  # Creating a function to follow the execution of the script
  progress <- function(n) setTxtProgressBar(p_bar, n)
  opts <- list(progress = progress)

  # Correcting betas through a parallelized for loop
  output <- foreach(cpg = rownames(betas_to_correct), .packages = "Kendall", .options.snow = opts) %dopar% {

    # ASSIGN REGRESSION TO THE DIFFERENT SAMPLES FOR EACH CPG

    identified_regressions <- identify_regression(
      vec_betas = as.numeric(betas_to_correct[cpg,]),
      vec_estimated_1mPurity = predicted_1mPurities_vec,
      vec_slopes = as.numeric(my_slopes[cpg,]),
      vec_intercepts = as.numeric(my_intercepts[cpg,])
    )


    # CORRECTING BETAS BASED ON THE IDENTIFIED REGRESSIONS

    return(
      list(

        "Tumour" =  correcting_betas(
          slopes_vec = as.numeric(identified_regressions$Slope),
          intercepts_vec = as.numeric(identified_regressions$Intercept),
          distances_vec = as.numeric(identified_regressions$Distance),
          to_correct = "Tumor"
        ),
        "Microenvironment" = correcting_betas(
          slopes_vec = as.numeric(identified_regressions$Slope),
          intercepts_vec = as.numeric(identified_regressions$Intercept),
          distances_vec = as.numeric(identified_regressions$Distance),
          to_correct = "Microenvironment"
        )
      )
    )

    # Stop clusters used in parallelization
    stopCluster(cl)

    # PROCESSING OUTPUT

    cat("\n\nGenerating output...\n\n")

    #Converting output list into dataframe of corrected tumor betas
    corrected_tumor <- as.data.frame(do.call(rbind, lapply(output, function(item) item[[1]])))
    colnames(corrected_tumor) <- colnames(betas_to_correct)
    rownames(corrected_tumor) <- rownames(betas_to_correct)

    #Converting output list into dataframe of corrected tumor betas
    corrected_microenvironment <- as.data.frame(do.call(rbind, lapply(output, function(item) item[[2]])))
    colnames(corrected_microenvironment) <- colnames(betas_to_correct)
    rownames(corrected_microenvironment) <- rownames(betas_to_correct)


    cat("\n\n=================\n")
    cat ("PROCESS FINISHED")
    cat("\n=================\n\n")

    return(
      list(
        "Corrected_tumour" = corrected_tumor,
        "Corrected_microenvironment"  = corrected_microenvironment
      )
    )

  }


}
