#!/usr/bin/Rscript

## -SCRIPT'S NAME: ref_regression_calculator.r
#
## - DESCRIPTION: 
#
#   This script generates the reference regressions' parameters and variance of the used reference CpGs
#   required to be used in subsequent steps of the PurEst workflow.
# 
## - USED R PACKAGES:
#
#   *OPTPARSE. Parsing command line arguments
#   *DOPARALLEL. Parallelizing execution
#   *PARALLEL. Parallelizing script
#   *FLEXMIX. Flexible mixture modelling used in the adjustBeta() function.
#
## - USER DEFINED FUNCTIONS:
#   
#   *adjustBeta(). This function corrects beta values based on the sample purity generating
#    corrected betas for the actual cancer cells and the tumor microenvironment. It also outputs
#    the parameters of the regressions used for the correction
#
## - PROCEDURE:
#
#   1. Installing (if necessary) and loading packages, configuring command line arguments and sourcing 
#      user defined functions.
#
#   2. Configuring parallelization.
#
#   3. Running the adjustBeta() function for each CpG using a parallelized apply function.
#
#   4. Adding the produced results to a result list.
#
#   5. Saving each element of the result list as an independent R object. As the prupose of this script is 
#      to generate the reference data required in the PurEst workflow, the corrected tumour and microenvironment
#      betas following the original Staaf & Aine approach are not provided. In order to get them lines 209, 210 and 
#      211 should be uncommented.
#
## - INPUT FILES:
#
#    -> Beta values per sample should be previously saved in .rds format. The file should contain a matrix with Illumina 
#       CpG identifiers as row names and samples as column names.
#
#    -> Named vector stored as a .rds containing the purity values of the samples used for the analysis
#
#
## - OUTPUT FILES:
#
#
#    -> .rds file containing a matrix with the slopes of the regressions used for the purity estimation and the beta correction
#
#    -> .rds file containing a matrix with the intercepts of the regressions used for the purity estimation and the beta correction.
#
#    -> .rds file containing a matrix with the Residual Standard Error of the regressions used for the purity estimation and the beta correction.
#
#    -> .rds file containing a matrix with the degrees of freedom of the regressions used for the purity estimation and the beta correction.
#
#    -> .rds file containing a named vector with the variance of the betas of each CpG to be used for the purity estimation.
#
## - USAGE:
#
#     The script must be run on the command line using the following flags. 
#
#     """
#     Rscript path_to_script/new_purity_corrector.r -c [num_of_cores] -b [path_to_betas] 
#     -p [path_to_puritied] -o [path_to_save_output_files] -n [prefix_output_files]
#     """
#     
#     *The function of the command line options are the following; 
#
#       -c: Number of cores to be used to run the program. Default: 1.
#       -b: The path to the file with the beta values to be analysed must be entered here. The file must be an R object containing a matrix with the CpGs as rows and samples as columns.
#       -p: The path to the file with the purity values of the samples to be analysed must be entered here. The file must be an R object containing a named vector.
#       -o: The path to the location where the output files will be saved must be entered here. The output is an R object. Default: working directory.
#       -n: The prefix to be used to name the output files. Default: output.
#
## - VERSION: 1.0
#
## - DATE: 22/01/2024
#
## - AUTHOR: Mattias Aine  (mattias.aine@med.lu.se)
## - ADAPTED BY: Iñaki Sasiain
## - AFFILIATION: Johan Staaf lab @ Lund University / Oncology & Pathology

# =============================
# LOADING THE REQUIRED PACKAGES
# =============================

options(repos = "https://cran.r-project.org/")

if(!requireNamespace("doParallel", quietly = TRUE)) {
  install.packages("doParallel") }

suppressPackageStartupMessages(library(doParallel))

if(!requireNamespace("parallel", quietly = TRUE)) {
  install.packages("parallel") }

suppressPackageStartupMessages(library(parallel))

if(!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse") }

suppressPackageStartupMessages(library(optparse))

# ==================================
# CONFIGURING COMMAND LINE ARGUMENTS
# ==================================

argument_list <- list(

  make_option(c("-c", "--cores"), type="integer", default=1,  
              help="Number of cores to be used to run the program [default %default]",
              metavar = "[number]"),

  make_option(c("-b", "--input_betas"), type="character",  
              help="The path to the file with the beta values to be analysed must be entered here. The file must be an R object containing a dataframe with the CpGs as rows and samples as columns.",
              metavar = "[file path]"),

  make_option(c("-p", "--input_purity"), type="character",  
              help="The path to the file with the purity values of the samples to be analysed must be entered here. The file must be an R object containing a dictionary vector.",
              metavar = "[file path]"),

  make_option(c("-o", "--output"), type="character", default="./",
              help="The path to the location where the output files will be saved must be entered here. The output is an R object. Default [%default]",
              metavar = "[file path]"),

  make_option(c("-n", "--output_name"), type="character", default="output",
              help="The prefix to be used to name the output files. Default [%default]",
              metavar = "[file path]")

)

#Parsing command line arguments
arguments <- parse_args(OptionParser(option_list=argument_list, 
                                    description="This program corrects methylation beta values providing parameters of the regressions used for the correction."))


# =====================================
#   SOURCING THE CORRECT BETAS FUNCTION
# =====================================

#Getting path of the file containing the function to be sourced
dir <- commandArgs()[4]
dir <- gsub("--file=", "", dir)
dir <- gsub("ref_regression_calculator.r", "new_function_correctBetas.r", dir)

#Sourcing function
source(dir)

# ===========================
# CONFIGURING PARALLELIZATION
# ===========================

cat("\nUsing", arguments$cores,"cores\n\n")

#Creating the cluster to run the process in parallel
cl <- makeCluster(arguments$cores)  
registerDoParallel(cl)  

#Making sure that all packages have access to the flexmix package. Using invisible()
#to avoid printing anything to the terminal.
invisible(clusterEvalQ(cl, {library("flexmix")}))


# ====================================
# LOADING AND VARIANCE BASED FILTERING
# ====================================

#Loading the data
unadjusted_betas <- readRDS(arguments$input_beta)
purities <- readRDS(arguments$input_purity)

#Create a vector with the variance of each cpg (row)
cpg_variance <- apply(unadjusted_betas, 1, var)
names(cpg_variance) <- rownames(unadjusted_betas)

# ==================
# ANALYSING THE DATA
# ==================

#Adding seed to each row of the beta value dataframe
betaRun <- cbind(seed=1:nrow(unadjusted_betas),unadjusted_betas)

#Storing sample names
betaNames <- colnames(unadjusted_betas)

#Running the analysis in parallel
res <- parRapply(cl = cl, #ClusterS to run the process
                 betaRun, #Beta values+the added seed
                 adjustBeta, #Function to correct betas
                 purity=purities, #Purity values
                 snames=betaNames, #Sample names
                 seed=TRUE) #The seed has been added in the data


# ====================
# CREATING RESULT LIST
# ====================

# Creating a list to add the results. If the user wants to get the corrected betas based on the original Staaf & Aine
# approach and the methylation patterns detected per CpG the corresponding line should be uncommented.
result_list <- list(

  #betas.original = do.call("rbind",lapply(res,function(x) x$y.orig)), #Original beta values
  #betas.tumor = do.call("rbind",lapply(res,function(x) x$y.tum)), #Corrected tumor beta values
  #betas.microenvironment = do.call("rbind",lapply(res,function(x) x$y.norm)), #Corrected microenvironment beta values
  #cpg.populations =  do.call("rbind",lapply(res,function(x) x$groups)), #Methylation patterns (populations) of each CpG
  reg.slopes = do.call("rbind",lapply(res,function(x) x$res.slopes)), #Slopes of the populations
  reg.intercepts = do.call("rbind",lapply(res,function(x) x$res.int)), #Intercepts of the populations
  reg.RSE = do.call("rbind",lapply(res,function(x) x$res.rse)), #Residual standard error
  reg.df = do.call("rbind",lapply(res,function(x) x$res.df)) #Degrees of freedom of the reversed regressions

)

# =====================
# CREATING OUTPUT FILES
# =====================

#Creating output files per each dataframe of the output_list list
lapply(names(result_list), function(n) {
  saveRDS(result_list[[n]],file=paste(arguments$output, arguments$output_name,"_",n,".rds",sep=""))
})

#Generating output file with variace values of each CpG
saveRDS(cpg_variance, file=paste(arguments$output, arguments$output_name,"_CpG_variance.rds",sep=""))

# Stop clusters used in parallelization
stopCluster(cl)

#cat("\n=================\n")
#cat ("PROCESS FINISHED")
#cat("\n=================\n")