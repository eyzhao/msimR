SignIT_runner <- function(simulated_data) {
    reference_signatures <- simulated_data$reference_signatures
    catalog <- simulated_data$catalog
    signature_names <- simulated_data$signature_names
    n_mutations <- simulated_data$n_mutations

    load('/projects/ezhao_prj/papers/SignIT-paper/analysis/scripts/stan_dso.RData')

    exposures <- get_exposures(
        mutation_catalog = catalog,
        reference_signatures = reference_signatures,
        stan_model = stan_dso,
        n_chains = 8,
        n_cores = 8,
        n_iter = 200,
        n_adapt = 200
    ) %>%
        get_exposure_summary_table() %>%
        mutate(
            signature_present = `2.5%` > (n_mutations / 2000)
        ) %>%
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


