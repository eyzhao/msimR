#' Compute Quantitative Accuracy Score
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
#' The metrics reported are L2 normalized Euclidean distance and cosine distance.
#' In L2 normalized Euclidean distance, the exposure vectors \eqn{\mathbf{x}, \mathbf{y}} are
#' normalized as in \eqn{\mathbf{x'} = frac{\mathbf{x}}{\lVert \mathbf{x} \rVert}} into
#' unit vectors. Then the Euclidean distance \eqn{\lVert \mathbf{y'} - \mathbf{x'} \rVert} is computed.
#'
#' The cosine distance is calculated as
#' \eqn{1 - \frac{\mathbf{x} \cdot \mathbf{y}}{\lVert \mathbf{x} \rVert \lVert \mathbf{y} \rVert}}.
#'
#' @param simulated_data        The data object output by simulate_samples()
#'                              with exposures calculated using one of the
#'                              signature exposure calculation methods.
#'
#' @return A list with named items that report the distance metrics.
#'
#' @import dplyr
#' @export

compute_quantitative_accuracy <- function(simulated_data) {
    calculated_exposures <- simulated_data$calculated_exposures %>%
        rename(calculated_exposure = exposure)

    merged <- simulated_data$true_exposures %>%
        rename(true_exposure = exposure) %>%
        left_join(calculated_exposures, by = 'signature') %>%
        replace_na(list(calculated_exposure = 0))

    euclidean_distance <- merged %>%
        mutate(
            calculated_exposure = calculated_exposure / sqrt(sum(calculated_exposure^2)),
            true_exposure = true_exposure / sqrt(sum(true_exposure^2))
        ) %>% mutate(
            square_difference = (calculated_exposure - true_exposure)^2
        ) %>%
        .$square_difference %>% sum %>% sqrt

    cosine_similarity <- function(a, b) { as.numeric(a %*% b / sqrt(sum(a^2) * sum(b^2))) }
    cosine_distance <- 1 - cosine_similarity(merged$calculated_exposure, merged$true_exposure)

    return(list(
        cosine_distance = cosine_distance, 
        euclidean_distance = euclidean_distance
    ))
}


