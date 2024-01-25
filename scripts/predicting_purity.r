#!/usr/bin/Rscript

#####======================================================================#####
###      Function to estimate purity intervals from each CpG of a sample     ###
#####======================================================================#####

##Author: Iñaki Sasiain Casado
##Affiliation: Johan Staaf lab @ Lund University / Oncology & Pathology


## Input = betas (methylation), reference regression parameters and alogorithm's parameters
## Output = One minus purity estimated of a CpG of a sample.

predicting_purity <- function(beta,slopes,intercepts,RSE,degrees_of_freedom,slope_threshold, alpha) {

    # CHECK IF ARGUMENTS ARE EMPTY
    if (any(sapply(list(beta, slopes, intercepts, RSE, degrees_of_freedom, slope_threshold, alpha), is.null))) {
        stop("One or more input arguments are empty.")
    }

    #Identifying the regressions whose residual standard error is over the threshold
    # and whose slope is too close to 0 to be informative
    to_be_ignored <- (slopes > -slope_threshold & slopes < slope_threshold) | is.na(slopes)

    # Removing the values that do not meet the conditions previously analysed.
    slopes <- slopes[!to_be_ignored]
    intercepts <- intercepts[!to_be_ignored]
    RSE <- RSE[!to_be_ignored]
    degrees_of_freedom <- degrees_of_freedom[!to_be_ignored]

    # Check that all the regressions are not uninformative and execute the following code 
    # to estimate the 1-Purity value. If they weren't NA would be assigned to the 1-Purity value,
    # so that CpG would not be used to estimate the final 1-Purity of the sample

    if (length(slopes)!=0) {

        #Creating a vector to store the regressions in which the beta value analysed is included using booleans
        matches <- ((beta - intercepts) / slopes >= 0) & ((beta - intercepts) / slopes <= 1)

        # If the length of the matches vector when matches==TRUE is one, a unique CpG methylation pattern (regression)
        # can be assigned to each beta value. If this value is different to one no signoficant association to a pattern can 
        # be made
        if (sum(matches)==1) {

            #Store the identified population (indexes of matches vector) into a variable
            identified_pop <- which(matches)

            #Predicting the 1-purity value from the inverse regression
            predicted_1_minus_p <- (beta - intercepts[identified_pop]) /slopes[identified_pop]

            # Get the prediction interval of the beta identified
            beta_int <- c(
               beta - qt((1-alpha)/2, degrees_of_freedom[identified_pop]) * RSE[identified_pop],
               beta + qt((1-alpha)/2, degrees_of_freedom[identified_pop]) * RSE[identified_pop]
             )
            

            #In some cases if the degrees of freedom are too low the result of qt(), the quantile of the Student's 
            #t-distribution can be undefined (NaN). Those cases are very uncommon,  but they can generate bugs in the
            #code. Therefore, if that happened, the one_minus_purity interval would be set to (NA, NA) in order to avoid
            #possible bugs

            if (anyNA(beta_int)) {

                one_minus_purity <- c(NA,NA)

            } else {

                #Getting the prediction interval of the 1 - purity value. Extrapolate from beta interval.
                one_minus_purity <- (beta_int - intercepts[identified_pop]) / slopes[identified_pop]
                one_minus_purity <- sort(one_minus_purity)
            
                one_minus_purity[one_minus_purity < 0] <- 0
                one_minus_purity[one_minus_purity > 1] <- 1

            }
            
        } else {
            #If the CpG can not be assigned to a single population (methylation pattern) the 1-purity value will not be assigned.
            one_minus_purity <- c(NA,NA)
        }
        
        #If all the regressions are uninformative NA value will be assigned to one_minus_purity.
        } else {
            one_minus_purity <- c(NA,NA)
        }

    #Setting the interval's endpoints to 0 or 1 if one of the limits lower or higher than those values.
    #The NAs are ignored
    if (!(is.na(one_minus_purity[1] & is.na(one_minus_purity[2])))) {
        one_minus_purity[1] <- max(0, one_minus_purity[1])
        one_minus_purity[2] <- min(1, one_minus_purity[2])
    }   
    
    #Returning the 1-purity vector
    return(one_minus_purity)
}