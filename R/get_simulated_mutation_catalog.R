get_simulated_mutation_catalog <- function(signatures, exposures) {
    signature_names <- signatures %>% select(-mutation_type) %>% colnames()

    stopifnot(
        signature_names == exposures$signature
    )

    signature_matrix <- signatures %>%
        select(-mutation_type) %>%
        as.matrix

    ideal_catalog <- signature_matrix %*% exposures$exposure
    probability_vector <- ideal_catalog / sum(ideal_catalog)

    simulated_catalog <- rmultinom(
        n = 1, 
        size = sum(ideal_catalog),
        prob = probability_vector
    )

    return(tibble(
        mutation_type = signatures$mutation_type,
        count = simulated_catalog %>% as.numeric
    ))
}
