perturb_signatures <- function(signatures, percent_deviation = 10) {
    signature_names <- signatures %>% select(-mutation_type) %>% colnames

    signatures %>% 
        gather(signature, proportion, -mutation_type) %>% 
        mutate(
            signature = factor(signature, levels = signature_names),
            deviated_proportion = sapply(
                proportion, 
                function(p) { 
                    rnorm(
                        n = 1, 
                        mean = p, 
                        sd = p * percent_deviation / 100
                    )
                }
            ),
            deviated_proportion = if_else(
                condition = deviated_proportion > 0, 
                true = deviated_proportion, 
                false = 0
            )
        ) %>% 
        select(-proportion) %>% 
        spread(signature, deviated_proportion)
}


