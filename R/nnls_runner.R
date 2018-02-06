nnls_runner <- function(simulated_data) {
    reference_signatures <- simulated_data$reference_signatures
    catalog <- simulated_data$catalog
    signature_names <- simulated_data$signature_names
    n_mutations <- simulated_data$n_mutations

    stopifnot(
        identical(reference_signatures$mutation_type, catalog$mutation_type)
    )

    signature_matrix <- reference_signatures %>% 
        select(-mutation_type) %>%
        as.matrix

    exposure_vector = nnls(signature_matrix, catalog$count)$x
    return(list(
        method_name = 'nnls',
        exposures = tibble(
            signature = signature_names,
            exposure = exposure_vector,
            signature_present = exposure_vector > n_mutations * 0.05
        )
    ))
}


