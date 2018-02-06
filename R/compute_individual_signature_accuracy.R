#' Compute Individual Signature Accuracy
#'
#' Each n-of-1 mutation signature method computes the exposure of each
#' signature (its contribution to the mutation burden). These numeric
#' outcomes are stored in the "exposure" column of the
#' calculated_exposures table in simulated data output. In order to
#' evaluate each method's accuracy at assigning exposure values, this method
#' uses various distance metrics to compare the true signatures which gave 
#' rise to a simulated mutation with the reconstructed signatures estimated 
#' by the signature method.
#' 
#' Unlike compute_quantitative_accuracy(), this function computes the individual
#' accuracy of each signature. This allows one to estimate how accurately each
#' individual signature is predicted with a given method.
#'
#' @param simulated_data        The data object output by simulate_samples()
#'                              with exposures calculated using one of the
#'                              signature exposure calculation methods.
#'
#' @return A data frame with each signature's accuracy
#'
#' @import dplyr
#' @export

compute_individual_signature_accuracy <- function(simulated_data) {
    calculated_exposures <- simulated_data$calculated_exposures %>%
        rename(calculated_exposure = exposure)

    merged <- simulated_data$true_exposures %>%
        rename(true_exposure = exposure) %>%
        left_join(calculated_exposures, by = 'signature') %>%
        replace_na(list(calculated_exposure = 0)) %>%
        mutate(
            error = calculated_exposure - true_exposure
        ) %>%
        select(signature, true_exposure, calculated_exposure, error)

    return(
        merged
    )
}


