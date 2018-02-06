SA_runner <- function(simulated_data) {
    reference_signatures <- simulated_data$reference_signatures
    catalog <- simulated_data$catalog
    signature_names <- simulated_data$signature_names
    n_mutations <- simulated_data$n_mutations

    exposures_output <- findSigExposures(
        M = catalog$count %>% as.matrix(),
        P = reference_signatures %>% reference_signatures_as_matrix(simulated_data$catalog),
        decomposition.method = decomposeSA
    ) %>%
        .$exposures

    exposures <- tibble(
        signature = rownames(exposures_output),
        exposure = (exposures_output %>% as.numeric) * n_mutations,
        signature_present = (exposures_output %>% as.numeric) > 0
    )

    return(list(
        method_name = 'SignatureEstimation_SA',
        exposures = exposures
    ))
}

