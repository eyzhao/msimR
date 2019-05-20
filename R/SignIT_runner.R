SignIT_runner <- function(simulated_data) {

	signature_present <- function(exposure_output) {
	  exposure_output$exposure_chain %>%
		plyr::ddply('signature', function(z) {
		  normal_fit <- fitdistr(z$exposure, 'normal')
		  tribble(
			~mean, ~sd,
			normal_fit$estimate['mean'],
			normal_fit$estimate['sd']
		  ) %>%
			mutate(
			  lCI = mean - 1.96 * sd,
			  uCI = mean + 1.96 * sd
			)
		}) %>%
		mutate(
		  signature_present = lCI > 0
		)
	}

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
    )

	signature_present_table <- exposures %>% signature_present() %>% select(signature, signature_present)

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


