#!/usr/bin/Rscript

#####======================================================================#####
###  Function for correcting Illumina 450/850K beta values for tumor purity  ###
#####======================================================================#####

##Author: Mattias Aine  (mattias.aine@med.lu.se)
##Adapted by: Iñaki Sasiain
##Affiliation: Johan Staaf lab @ Lund University / Oncology & Pathology

# =============================
# LOADING THE REQUIRED PACKAGES
# =============================

#Loading Fleximix for clustering
if(!requireNamespace("flexmix", quietly = TRUE)) {
  install.packages("flexmix") }

if(! "flexmix" %in% names(sessionInfo()$otherPkgs) ) {
  library("flexmix") } 

# ==================================
# DEFINING THE adjustBeta() FUNCTION
# ==================================

## Input = betas (methylation), purity estimates and sample names
## Output = corrected betas and line parameters for the populations identified

adjustBeta <- function(methylation=NULL,purity=NULL,snames=NULL,nmax=3,nrep=3,seed=TRUE) {

    #If a seed is provided as the first element of the methylation vector (seed=TRUE) 
    #set seed for the cluster detremination and remove it from the methylation vector 
    if(seed) {
        set.seed(as.integer(methylation[1]))
        methylation<-methylation[-1]
    }

    #Defining variables for clustering and regression
    x <- as.numeric(purity)
    x2 <- 1-as.numeric(purity)
    y <- as.numeric(methylation)

    #Calculate global correlation between beta and purity
    gl.corr <- suppressWarnings(cor(x,y))
    gl.corr[is.na(gl.corr)] <- 0
    gl.corr <- round(gl.corr,3)

    #Add small gaussian noise to x (avoid errors when large number of zero samples)
    y2<-y+rnorm(length(y),mean=0,sd=.005)  

    #Modeling the methylation patterns (CpG populations)
    model <- stepFlexmix(y2 ~ x,k = 1:nmax, nrep = nrep,verbose = FALSE)
    model <- getModel(model, "BIC") 
    cl <- clusters(model)
    #Make sure clusters are numbered 1 to 3, odd cases exist where one pop has zero members from flexMix
    #can rename because flexmix object not used after this
    cl <- as.integer(factor(cl))



    ## CALCULATING REGESSIONS

    #Determining the b_vs_pur regressions for the clusters identified
    b_vs_pur <- lapply(1:nmax,function(z) { 
        if(z %in% cl) {
            lm(y[cl==z]~x[cl==z]) #Beta VS 1-Purity regression

        } else { NA } #If less than 3 populatins were detected NA is added

    })

    #Determining the b_vs_1mp regressions for the clusters identified
    b_vs_1mp <- lapply(1:nmax,function(z) { 
        if(z %in% cl) {
            lm(y[cl==z]~x2[cl==z]) #Beta VS 1-Purity regression
        } else { NA } #If less than 3 populatins were detected NA is added
    })


    ## DETERMINING THE FUNCTION'S OUTPUT FROM THE CALCULATED REGERSSIONS

    output_list <- list() #Creating a list to store the data

    #Adding the population to which each CpG belongs to the list
    output_list$groups <- cl    

    #Adding the number of population to which each CpG belongs to the list
    output_list$n.groups <- length(levels(factor(cl)))

    #Adding the correlation between betas and purity to the list
    output_list$glob.corr <- gl.corr

    #Adding corrected microenvironment betas to the list
    output_list$y.norm <- sapply(
        X = 1:nmax,
        FUN = function(z) {
                if (z %in% cl) {
                    n_vals <- coefficients(b_vs_pur[[z]])[1]+residuals(b_vs_pur[[z]])
                    names(n_vals)<-snames[cl==z]
                    n_vals

                } else {
                    NULL
                }
            })
    output_list$y.norm <- unlist(output_list$y.norm)[snames]

    #Adding corrected tumor betas to the list
    output_list$y.tum <- sapply(
        X = 1:nmax,
        FUN = function(z) {
                if (z %in% cl) {
                    n_vals <- coefficients(b_vs_1mp[[z]])[1]+residuals(b_vs_1mp[[z]])
                    names(n_vals)<-snames[cl==z]
                    n_vals

                } else {
                    NULL
                }
            })
    output_list$y.tum <- unlist(output_list$y.tum)[snames]


    ##in very rare instances flexmix calls 3 populations but one groups has zero members. -> has parameters for non-existant pop in output?!
    #output_list$res.int <- round(as.numeric(unlist(lapply(slot(model,"components"),function(z) slot(z[[1]],"parameters")$coef[1]))),3) 
    #output_list$res.slopes <- round(as.numeric(unlist(lapply(slot(model,"components"),function(z) slot(z[[1]],"parameters")$coef[2]))),3)

    #Getting line parameters from b_vs_1mp: INTERCEPTS
    output_list$res.int <- sapply(
        X = 1:nmax,
        FUN = function(z) {
                if (z %in% cl) {
                    int <- coefficients(b_vs_1mp[[z]])[1]
                    round(as.numeric(int),3)

                } else {
                    NA
                }
            },
            simplify = TRUE)

    #Getting line parameters from b_vs_1mp: SLOPES
    output_list$res.slopes <- sapply(
        X = 1:nmax,
        FUN = function(z) {
                if (z %in% cl) {
                    int <- coefficients(b_vs_1mp[[z]])[2]
                    round(as.numeric(int),3)

                } else {
                    NA
                }
            },
            simplify = TRUE)


    #Getting line parameters from b_vs_1mp: RESIDUAL STANDARD ERROR
    output_list$res.rse <- sapply(
        X = 1:nmax,
        FUN = function(z) {
                if (z %in% cl) {
                    RSE <- summary(b_vs_1mp[[z]])$sigma
                    round(as.numeric(RSE),6)

                } else {
                    NA
                }
            },
            simplify = TRUE)

    #Getting line parameters from b_vs_1mp: RESIDUAL STANDARD ERROR
    output_list$res.df <- sapply(
        X = 1:nmax,
        FUN = function(z) {
                if (z %in% cl) {
                    df.residual(b_vs_1mp[[z]])

                } else {
                    NA
                }
            },
            simplify = TRUE)

    #Cap the corrected betas to 0 and 1
    output_list$y.tum[output_list$y.tum > 1] <- 1
    output_list$y.tum[output_list$y.tum < 0] <- 0
    output_list$y.norm[output_list$y.norm > 1] <- 1
    output_list$y.norm[output_list$y.norm < 0] <- 0

    #Round the betas to three decimal places
    output_list$y.tum <- round(output_list$y.tum,3)
    output_list$y.norm <- round(output_list$y.norm,3)

    #Adding the rounded original betas to the list
    output_list$y.orig <- round(methylation,3)    

    return(output_list)
}