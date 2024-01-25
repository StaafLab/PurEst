#!/usr/bin/Rscript

## -SCRIPT'S NAME: purity_estimator.r
#
## - DESCRIPTION: 
#
#   This script predicts the purity of cancer samples based on the CpG's beta values and refernce
#   regressions that reflect the different methylation patterns of the different populations generated
#   with the ref_regression_calculator.r script.
# 
## - USED R PACKAGES:
#
#   *OPTPARSE. Parsing command line arguments
#   *PARALLEL. Parallelization of the script
#   *DOSNOW. Parallelization of the script
#   *FOREACH. Parallelization of the script
#   *KENDALL. Analyse progess through progress bar
#
## - USER DEFINED FUNCTIONS:
#   
#   * predicting_purity(). This function creates a dataframe containing the predicted 1-Purity intervals 
#                          from each CpG of each sample
#   * purity_coverage(). This function calculates a 1-Purity estimate per sample
#                        based on the intervals predicted for all the CpGs of each sample.
#
## - PROCEDURE:
#
#   1. Installing (if necessary) and loading packages, configuring command line arguments and sourcing
#      functions to be used.
#
#   2. Loading the required data. The R objects containing the parameters of the regressions and the beta
#      values which will be used for the analysis are loaded.
#
#   3. Filtering refrence regressions based on variance.
#
#   4. Configure parallelization . The number of cores to be used must be specified using the
#      corresponding command line flag.
#
#   5. Running the analysis per each sample. First the predicting_purity() function is used to predict
#      1-Purity intervals per each CpG of the sample analysed. Then, the purity_coverage() function is
#      used to calculate a single 1-Purity estimate and an interval per each sample based on all the 
#      predicted intervals calculated for each CpG. If the analysed CpG is not included into the reference 
#      regression dataset it will be ignored.
#
#   6. Storing the output data. A .rds file containing the CpGs used to estimate the purity of each sample can
#      be also obtained if line 331 is uncommented.
#
## - INPUT FILES:
#
#    -> Matrix stored as an R object with the SLOPES of the calculated regressions. All the parameters of
#       the regressions must be stored in the same directory. The name of the file must end with "reg.slopes.rds".
#
#    -> Matrix stored as an R object with the INTERCEPTS of the calculated regressions. All the parameters of
#       the regressions must be stored in the same directory. The name of the file must end with "reg.intercepts.rds".
#
#    -> Matrix stored as an R object with the RESIDUAL STANDARD ERRORS of the calculated regressions. All the parameters of
#       the regressions must be stored in the same directory. The name of the file must end with "reg.RSE.rds".
#
#    -> Matrix stored as an R object with the DEGREES OF FREEDOM of the calculated regressions. All the parameters of
#       the regressions must be stored in the same directory. The name of the file must end with "reg.df.rds".
#
#    -> Matrix stored as an R object with the BETA VALUES of the CpGs of the SAMPLES to analyse. The CpG Ids have to be 
#       the row names and the sample names the column names of the object.
#
#    -> Named vector stored as an R object with the VARIANCE of the CpGs used to calculate the reference regressions. 
#
## - OUTPUT FILES:
#
#    -> R object, whose name must be specified by the user, containing a list with the estimate 1-Purity values and intervals
#       predicted per each sample. If more than one estimates are obtained per samples (possible but extremely rare), more than
#       one 1-Pur. values and intervals will be generated.
#
#    -> TSV file, whose name must be specified by the user, containing the identified estimates (it verw exceptrional cases it could be different to 1) 
#       the estimated 1-Purity values and intervals predicted per each sample. If more than one estimates are obtained, one independent line per estimate
#       will be created in the tsv file.
#
## - USAGE:
#
#     The script must be run on the command line using the following flags. Keep in mind that the functions that the script calls
#     must be stored in the same directory of the main script.
#
#     """
#     Rscript path_to_script/purity_estimator.r -c [cores] -a [alpha] -r [threshold_rse] -s [threshold_slope] -p [percentage_to_interval] 
#     -v [variance_threshold] -d [path_to_regression_data] -b [path_to_betas] -o [output_filename] -l [output_location]
#     """
#     
#     *The function of the command line options is the following; 
#
#       -c: Number of cores to run the program
#       -a: The alpha value used to determine the prediction intervals from the regressions
#       -s: Minimum slope allowed per regression. The regressions with lower slopes will be ignored
#       -p: Percentage of the maximum coverage detected to include in estimated the 1-Purity interval
#       -v: CpG beta value variance cutoff to filter reference regressions.
#       -d: The directory containing the regression parameters must be entered here
#       -b: The path to the R object containing the betas to analyse must be entered here
#       -o: The name of the output R object containing the predicted values must be entered here.
#       -l: The name of the location output R object containing the predicted values must be entered here
#
## - VERSION: 1.0
#
## - DATE: 29/10/2023
#
## - AUTHOR: Iñaki Sasiain Casado
## - AFFILIATION: Johan Staaf lab @ Lund University / Oncology & Pathology

# ==================================
# INSTALLING AND/OR LOADING PACKAGES
# ==================================

#Install and load optparse to parse command line arguments
if(!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse", quietly = TRUE) }

suppressPackageStartupMessages(library(optparse))

#Install and load parallel to run the script in parallel
if(!requireNamespace("parallel", quietly = TRUE)) {
  install.packages("parallel", quietly = TRUE) }

suppressPackageStartupMessages(library(parallel))

#Install and load doSnow to run the script in parallel
if(!requireNamespace("doSNOW", quietly = TRUE)) {
  install.packages("doSNOW", quietly = TRUE) }

suppressPackageStartupMessages(library(doSNOW))

#Install and load foreach to run the script in parallel
if(!requireNamespace("foreach", quietly = TRUE)) {
  install.packages("foreach", quietly = TRUE) }

suppressPackageStartupMessages(library(foreach))

#Install and load Kendall to create a progress bar
if(!requireNamespace("Kendall", quietly = TRUE)) {
  install.packages("Kendall", quietly = TRUE) }

suppressPackageStartupMessages(library(Kendall))

# ==================================
# CONFIGURING COMMAND LINE ARGUMENTS
# ==================================

#Defining command line arguments
argument_list <- list(

  make_option(c("-c", "--cores"), type="integer", default=1,  
              help="Number of cores to be used to run the program [default %default]",
              metavar = "[number]"),
  
  make_option(c("-a", "--alpha"), type="double", default=0.7, 
              help="Alpha value to determine the prediction intervals of each CpG [default %default]",
              metavar= "[floating number]"),

  make_option(c("-s", "--min_slope"), type="double", default=0.25,
              help="Minimum slope allowed per CpG regression [default %default]", 
              metavar="[floating number]"),

  make_option(c("-v", "--variance_threshold"), type="numeric", default=0.05,
              help="Only the CpGs whose betas' variance are over this threshold will be used to determine the refrence regressions. Default [%default]",
              metavar = "[variance_threshold]"),

  make_option(c("-p", "--percentage_to_interval"), type="double", default=0.96,
              help="Percentage of the maximum coverage to include in the 1-Purity interval [default %default]",
              metavar="[floating number]"),

  make_option(c("-d", "--regression_data"), type="character",
              help="The directory containing the regression parameters must be entered here.",
              metavar="[directory]"),

  make_option(c("-b", "--betas_to_analyse"), type="character",
              help="The path to the R object contaoining the betas to analyse must be entered here.",
              metavar="[directory]"),

  make_option(c("-o", "--output_filename"), type="character", default="output",
              help="The name of the output R object containing the predicted values must be entered here. This name will also be used as the prefix of the file containing the cpgs used per sample. Default [%default]",
              metavar="[filename]"),

  make_option(c("-l", "--output_location"), type="character", default="./",
              help="The name of the location output R object containing the predicted values must be entered here. Default [%default]",
              metavar="[path_to_directory]")
              
)

#Parsing arguments
arguments <- parse_args(OptionParser(option_list=argument_list, 
                                    description="This program estimates the Purity values of samples based on the beta values of each sample's CpGs."))


# ==============================
# LOADING THE REQUIRED FUNCTIONS
# ==============================

# Getting the directory of the functions sourced on this script 
# (They must be located on the same directory than the main script)
dir <- commandArgs()[4]
dir <- gsub("--file=", "", dir)

#Getting the complete path of the functions to be called
fun1 <- gsub("purity_estimator.r", "predicting_purity.r", dir)
fun2 <- gsub("purity_estimator.r", "purity_coverage.r", dir)

#Sourcing the functions
source(fun1)
source(fun2)

# ========================
# PREPARING THE INPUT DATA
# ========================


#Reading the R objects containing the regression data as dataframes
my_slopes <- readRDS(list.files(arguments$regression_data, pattern="*reg.slopes.rds", full.names=TRUE))
my_intercepts <- readRDS(list.files(arguments$regression_data, pattern="*reg.intercepts.rds", full.names=TRUE))
my_RSE <- readRDS(list.files(arguments$regression_data, pattern="*reg.RSE.rds", full.names=TRUE))
my_df <- readRDS(list.files(arguments$regression_data, pattern="*reg.df.rds", full.names=TRUE))


#Reading the R objects containing variance of reference CpGs
my_CpG_variance <- readRDS(list.files(arguments$regression_data, pattern="*CpG_variance.rds", full.names=TRUE))

#Reading beta values
unadj_validation <- readRDS(arguments$betas_to_analyse)

#Create a list to append all the predicted purity intervals
list_of_predicted_intervals <- list()


# =================================================
# FILTERING REFERENCE REGRESSIONS BASED ON VARIANCE
# =================================================

#Generate a vector with the CpGs to filter
cpgs_to_keep <- names(my_CpG_variance[my_CpG_variance >= arguments$variance_threshold])

#Filtering regression objects
my_slopes <- my_slopes[cpgs_to_keep,]
my_intercepts <- my_intercepts[cpgs_to_keep,]
my_RSE <- my_RSE[cpgs_to_keep,]
my_df <- my_df[cpgs_to_keep,]

# ===========================
# CONFIGURING PARALLELIZATION
# ===========================

# Printing the number of cores to be used
cat("\nUsing", arguments$cores, "cores\n")

# Creatinhg clusters to run the script in  parallel
cl <- makeCluster(arguments$cores)

#Registering the clusters
registerDoSNOW(cl)

# ========================================
# RUNNING THE ANALYSIS WITH A PROGRESS BAR
# ========================================

# Printing command line message
cat("\nRunning the analysis...\n\n")

# Getting the names of the samples to analyse
samples <- colnames(unadj_validation)

# Defining the progress bar
p_bar <- txtProgressBar(min = 0, 
                        max = length(samples), 
                        style = 3)

# Creating a function to follow the execution of the script
progress <- function(n) setTxtProgressBar(p_bar, n)
opts <- list(progress = progress)

# Running the sourced functions in parallel for each sample. The execution level will be followed through a progress bar
out_list <- foreach(s = samples, .packages = "Kendall", .options.snow = opts) %dopar% {

  # Defining an empty matrix with the cpg ids as rownames to add the all the 1-Purity predicted intervals for all 
  # the CpGs of a sample
  interval_mat <- matrix(ncol=2, nrow=length(rownames(unadj_validation)))
  rownames(interval_mat) <- rownames(unadj_validation)
  

  # Predicting all the 1-Purity intervals for each CpG of each sample and append them to the empty interval_mat
  for (cpg in rownames(unadj_validation)) {

    # The following if statement will be used to take into account only cpgs included into the
    # refernce regression dataset
    if (cpg %in% rownames(my_slopes)) {

      interval_mat[cpg,] <- predicting_purity(beta=unadj_validation[cpg, s],
                                              slopes=my_slopes[cpg, ],
                                              intercepts=my_intercepts[cpg, ],
                                              RSE=my_RSE[cpg, ],
                                              degrees_of_freedom=my_df[cpg, ],
                                              arguments$min_slope,
                                              alpha=arguments$alpha)
    
    }
  }

  # Calculate the 1-Purity estimate and interval for the sample analysed.
  # The results with be shown in list named with the sample id
  list(name = s, 
       value = purity_value_per_sample(
                      pred_purity_confidence=interval_mat,
                      interval_threshold=100*(1-arguments$percentage_to_interval)),
        cpgs = rownames(na.omit(interval_mat))
       )
}

# Append the list defined for each sample to the list containing the predicted values for all the samples.
# The sample id is used to identify each element of the list
list_of_predicted_intervals <- setNames(lapply(out_list, function(x) x$value), sapply(out_list, function(x) x$name))
list_of_used_cpgs <- setNames(sapply(out_list, function(x) x$cpgs), sapply(out_list, function(x) x$name))


# =========================
# CREATING OUTPUT R OBJECTS
# =========================

# Stop the defined clusters
stopCluster(cl)

cat("\n\nSaving output files...\n")

# Save the list of predicted values as an R object
saveRDS(list_of_predicted_intervals, file=paste(arguments$output_location, arguments$output_filename, ".rds", sep=""))

# Save the list containing the cpgs used to estimate each purity value. Uncomment the following line to generate
# the a file containing the cpgs used for the purity estimation of each sample.
#saveRDS(list_of_used_cpgs, file=paste(arguments$output_location, arguments$output_filename, ".used_cpgs.rds", sep=""))


# ==========================
# CREATING OUTPUT TEXT FILES
# ==========================

# Create a vector with the column names of the output dataframe
cols <- c("#sample", "num_of_est", "Estimate_1-purity", "low_bound", "top_bound")

# Creating a dataframe with the columns below
output_tsv <- data.frame(matrix(nrow=0, ncol=length(cols)))

# Appending the values of the list of predicted intervals
for (sample in names(list_of_predicted_intervals)) {

  # A different row will be appended per each detected estimate per sample
  for (num in length(list_of_predicted_intervals[[sample]][["1-Pur_estimates"]])) {

    # Creating vector with the data to append
    row <- c(sample, 
             length(list_of_predicted_intervals[[sample]][["1-Pur_estimates"]]), # This will indicate the number of estimates detected
             list_of_predicted_intervals[[sample]][["1-Pur_estimates"]][num],
             list_of_predicted_intervals[[sample]][["interval(s)"]][[num]][1],
             list_of_predicted_intervals[[sample]][["interval(s)"]][[num]][2]
             )
    
    # Appending row to dataframe
    output_tsv <- rbind(output_tsv, row)
          
  }

}

#Setting column names
colnames(output_tsv) <- cols

# Saving text file (tsv)
write.table(output_tsv, 
            file=paste(arguments$output_filename, ".tsv", sep=""),
            row.name=FALSE,
            col.names=TRUE,
            quote=FALSE,
            sep="\t")


# Print message to show the end of the execution
cat("\n\n**********************\n")
cat("   PROCESS FINISHED\n")
cat("**********************\n")