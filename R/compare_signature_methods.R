compare_signature_methods <- function(
    function_list,
    simulation_args = list(),
    full_output = TRUE,
    run_in_parallel = FALSE,
    n_cores = NULL
) {
    simulated_data <- do.call(
        'simulate_samples', 
        args = simulation_args
    )

    simulation_results <- plyr::llply(function_list, function(exposure_calculation_function) {
        start_time <- Sys.time()
        exposure_data <- exposure_calculation_function(simulated_data)
        simulated_data$method <- exposure_data$method_name
        simulated_data$calculated_exposures <- exposure_data$exposures
        end_time <- Sys.time()

        simulated_data$quantitative_accuracy <- compute_quantitative_accuracy(simulated_data)
        simulated_data$signature_errors <- compute_individual_signature_accuracy(simulated_data)
        simulated_data$run_time <- as.double(end_time - start_time, units = 'secs')

        message(sprintf('Completed %s. Runtime:', simulated_data$method))
        message(end_time - start_time)

        if (full_output) {
            return(simulated_data)
        } else {
            summary_table <- tibble(
               method = simulated_data$method,
               n_mutations = simulated_data$n_mutations,
               n_signatures = length(simulated_data$chosen_signatures),
               perturbation_percent_deviation = simulated_data$perturbation_percent_deviation,
               euclidean_distance = simulated_data$quantitative_accuracy$euclidean_distance,
               cosine_distance = simulated_data$quantitative_accuracy$cosine_distance,
               run_time = simulated_data$run_time
            ) %>% gather(metric, value, -method)

            signature_error_table <- simulated_data$signature_errors %>%
                gather(type, value, -signature) %>%
                unite(metric, type, signature, sep='|') %>%
                mutate(
                    method = simulated_data$method
                ) %>%
                select(method, metric, value)

            return(bind_rows(summary_table, signature_error_table))
        }
    })

    if (full_output) {
        return(simulation_results)
    } else {
        return(do.call('bind_rows', simulation_results))
    }
}
