SignIT_runner <- function(simulated_data) {

    reference_signatures <- simulated_data$reference_signatures
    catalog <- simulated_data$catalog
    signature_names <- simulated_data$signature_names
    n_mutations <- simulated_data$n_mutations

    signit_cores <- getOption('signit_cores')
    if (is.null(signit_cores)) {
        signit_cores <- 8
    }
    exposures <- get_exposures(
        mutation_catalog = catalog,
        reference_signatures = reference_signatures,
        n_chains = 8,
        n_cores = signit_cores,
        n_iter = 200,
        n_adapt = 200
    )

	signature_present_table <- exposures$signatures_present

	exposures <- exposures %>%
        get_exposure_summary_table() %>%
		inner_join(signature_present_table, by = 'signature') %>%
        select(
            signature,
            exposure = mean_exposure,
            signature_present
        )

    return(list(
        method_name = 'SignIT',
        exposures = exposures
    ))
}


