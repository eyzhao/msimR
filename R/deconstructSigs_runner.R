#' Computes exposures by the deconstructSigs method
#' 
#' deconstructSigs is a method for computing n-of-1 mutation signatures
#' reported in Genome Biology 
#' (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0893-4)
#' and implemented as a package on CRAN 
#' (https://cran.rstudio.com/web/packages/deconstructSigs/index.html).
#'
#' This method runs deconstructSigs on a simulated data object and
#' returns the associated signature exposures. Because deconstructSigs
#' employs its own feature selection, any signature with a non-zero
#' exposure value is considered to be a contributor to the mutation burden.
#' Thus, this method calls a signature present if there is at least
#' one mutation contributed by that signature.
#'
#' @param simulated_data        Simulated data object obtained by
#'                              \code{\link{simulate_signatures}}.
#'
#' @return A list containing a table of signature exposures and the name of the method used.
#'
#' @import deconstructSigs 
#' @import dplyr
#' @import tidyr
#'
#' @export

deconstructSigs_runner <- function(simulated_data) {
    reference_signatures <- simulated_data$reference_signatures
    catalog <- simulated_data$catalog
    signature_names <- simulated_data$signature_names
    n_mutations <- simulated_data$n_mutations

    deconstructSigs_reference_signatures <- reference_signatures %>% 
        gather(signature, proportion, -mutation_type) %>% 
        mutate(signature = factor(signature, levels = signature_names)) %>%
        spread(mutation_type, proportion) %>%
        as.data.frame %>% 
        `rownames<-`(NULL) %>%
        column_to_rownames('signature')

    deconstructSigs_catalog <- simulated_data$catalog %>% 
        mutate(count = count / sum(count)) %>%
        spread(mutation_type, count) %>% 
        as.data.frame %>%
        `rownames<-`('test')

    deconstructSigs_exposures <- whichSignatures(
        tumor.ref = deconstructSigs_catalog,
        signatures.ref = deconstructSigs_reference_signatures
    ) %>%
        .$weights %>% 
        gather(signature, exposure) %>%
        filter(signature %in% rownames(deconstructSigs_reference_signatures))

    exposures <- tibble(signature = signature_names) %>%
        left_join(deconstructSigs_exposures, by = 'signature') %>%
        replace_na(list(exposure = 0)) %>%
        mutate(
            exposure = exposure * n_mutations,
            signature_present = exposure > 1
        )

    return(
        list(
            method_name = 'deconstructSigs',
            exposures = exposures
        )
    )
}
