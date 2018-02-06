uniform_exposure_sampler <- function(n_signatures, n_mutations) {
    values = runif(
        n_signatures, 
        min = 10,
        max = 1000
    )

    n_mutations * (values / sum(values))
}
