test_signature_method <- function(
    exposure_calculation_function,
    n_simulations = 1000, 
    simulation_args = list(),
    full_output = TRUE,
    run_in_parallel = FALSE,
    n_cores = NULL
) {
    if (run_in_parallel) {
        if (is.null(n_cores)) {
            registerDoParallel()
            message(sprintf('Running in parallel with maximum number of cores: %s.', getDoParWorkers()))
        } else if (n_cores <= 1) {
            warning("Can't run in parallel with only one core. Running sequentially instead.")
            run_in_parallel = FALSE
        } else if (n_cores > 1) {
            registerDoParallel(n_cores)
            message(sprintf('Running in parallel with %s cores.', getDoParWorkers()))
        } else {
            stop('Error in parallel arguments provided to test_signature_method')
        }
    }

    simulation_results <- plyr::llply(1:n_simulations, function(simulation_index) {
        simulated_data <- do.call(
            'simulate_samples', 
            args = simulation_args
        )

        start_time <- Sys.time()
        exposure_data <- exposure_calculation_function(simulated_data)
        simulated_data$method <- exposure_data$method_name
        simulated_data$calculated_exposures <- exposure_data$exposures
        end_time <- Sys.time()

        simulated_data$quantitative_accuracy <- compute_quantitative_accuracy(simulated_data)
        simulated_data$run_time <- as.double(end_time - start_time, units = 'secs')

        if (full_output) {
            return(simulated_data)
        } else {
            return(tibble(
               method = simulated_data$method,
               n_mutations = simulated_data$n_mutations,
               n_signatures = length(simulated_data$chosen_signatures),
               included_signatures = paste(simulated_data$chosen_signatures, collapse = ','),
               perturbation_percent_deviation = simulated_data$perturbation_percent_deviation,
               euclidean_distance = simulated_data$quantitative_accuracy$euclidean_distance,
               cosine_distance = simulated_data$quantitative_accuracy$cosine_distance,
               run_time = simulated_data$run_time
            ))
        }
    }, .parallel = run_in_parallel)

    if (full_output) {
        return(simulation_results)
    } else {
        return(do.call('bind_rows', simulation_results))
    }
}
