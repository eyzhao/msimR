#' Compute Qualitative Accuracy Score
#'
#' Each n-of-1 mutation signature method estimates which signatures are present
#' and which absent in each given patient/sample.
#' In order to evaluate each method's accuracy at assigning presence/absence
#' of signatures, this method computes the sensitivity, specificity, ppv, npv,
#' and overall accuracy.
#'
#' @param simulated_data        The data object output by simulate_samples()
#'                              with exposures calculated using one of the
#'                              signature exposure calculation methods.
#'
#' @return A list with named items that report the distance metrics.
#'
#' @import dplyr
#' @export

compute_qualitative_accuracy <- function(simulated_data) {
    signature_present_table <- simulated_data$true_exposures %>% 
      mutate(signature_present = exposure > 0) %>% 
      inner_join(
        simulated_data$calculated_exposures %>% 
        select(signature, calculated_signature_present = signature_present)
      )

    classification_table <- signature_present_table %>%
      mutate(
        metric = paste('accuracy_class', signature, sep='|'),
        value = case_when(
          signature_present & calculated_signature_present ~ 1, # true positive
          ! signature_present & ! calculated_signature_present ~ 2, # true negative
          ! signature_present & calculated_signature_present ~ -1, # false positive
          signature_present & ! calculated_signature_present ~ -2, # false negative
        )
      ) %>%
      select(metric, value)

    metrics <- signature_present_table %>%
      summarise(
        sensitivity = sum(calculated_signature_present & signature_present) / sum(signature_present), 
        specificity = sum(! calculated_signature_present & ! signature_present) / sum(! signature_present), 
        ppv = sum(calculated_signature_present & signature_present) / sum(calculated_signature_present), 
        npv = sum(! calculated_signature_present & ! signature_present) / sum(! calculated_signature_present), 
        accuracy = (sum(signature_present & calculated_signature_present) + sum(! signature_present & ! calculated_signature_present))/ n()
      ) %>%
      as.list()

    metrics$classification_table <- classification_table

    return(metrics)
}


