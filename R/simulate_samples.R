#' Simulate samples from reference signatures
#'
#' 

simulate_samples <- function(
    reference_signatures = wtsi_30_snv_signatures,
    n_mutations = 1000,
    n_active_signatures = 5,
    chosen_signatures = c(),
    exposure_sampler = uniform_exposure_sampler,
    perturbation_percent_deviation = 0
) {
    # Validate inputs

    sample_exposures <- exposure_sampler(n_active_signatures, n_mutations)

    stopifnot(
        'mutation_type' %in% names(reference_signatures),
        dim(reference_signatures)[2] > 1,
        sum(sample_exposures) - n_mutations < 1,
        all(sample_exposures > 0),
        all(chosen_signatures %in% names(reference_signatures))
    )

    # Perturb the reference signatures

    perturbed_reference_signatures <- reference_signatures %>%
        perturb_signatures(perturbation_percent_deviation)

    # Transfer exposure values to appropriate signatures

    signature_names = reference_signatures %>% 
        select(-mutation_type) %>%
        colnames()

    if (length(chosen_signatures) == 0) {
        chosen_signatures <- signature_names %>%
            sample(n_active_signatures, replace = FALSE)
    }

    exposure_table <- tibble(
        signature = chosen_signatures, 
        exposure = sample_exposures
    ) %>% bind_rows(tibble(
        signature = signature_names %>% .[! . %in% chosen_signatures],
        exposure = 0
    )) %>%
    mutate(signature = factor(signature, levels = signature_names)) %>%
    arrange(signature)

    # Simulate the Mutation Catalog

    catalog <- get_simulated_mutation_catalog(perturbed_reference_signatures, exposure_table)

    # Assemble the list object to return

    list(
        catalog = catalog,
        true_exposures = exposure_table,
        n_mutations = n_mutations,
        perturbed_reference_signatures = perturbed_reference_signatures,
        perturbation_percent_deviation = perturbation_percent_deviation,
        reference_signatures = reference_signatures,
        signature_names = signature_names,
        chosen_signatures = chosen_signatures
    )
}


