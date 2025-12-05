#' Load Parameters from an Excel File
#'
#' Reads an Excel file containing the parameters required for base
#' editing detection and standardizes the column names. The first
#' seven columns are expected to be in a specific order.
#'
#' @param path A \code{character} string specifying the readable
#'   file path to the Excel parameters spreadsheet (e.g.,
#'   \code{"params.xlsx"}).
#' @return A \code{data.frame} where the first seven column names
#'   are set to \code{"sample_name"}, \code{"sample_file"},
#'   \code{"ctrl_file"}, \code{"motif"}, \code{"motif_fwd"},
#'   \code{"wt"}, and \code{"edit"}.
#' @importFrom readxl read_excel
#' @export
#' @examples
#' # Assuming "path/to/params.xlsx" exists and has the correct format.
#' # params_df <- load_parameters_file("path/to/params.xlsx")
load_parameters_file <- function(path){
    params <- readxl::read_excel(path)
    fixed_names <- c("sample_name", "sample_file", "ctrl_file", 
                     "motif", "motif_fwd", "wt", "edit")
    colnames(params)[seq_along(fixed_names)] <- fixed_names
    params
}

#' Load an Example Parameters Table
#'
#' Loads the pre-packaged example parameters table provided with the
#' \code{MultiEditR} package, automatically setting the correct
#' paths for the example sequence files.
#'
#' @return A \code{data.frame} containing the example parameters
#'   with absolute paths to the sample and control files in the
#'   package's \code{extdata} directory.
#' @export
#' @examples
#' example_params <- load_example_params()
load_example_params <- function(){
    ext_path <- system.file("extdata", package = "MultiEditR")
    params <- load_parameters_file(file.path(ext_path, "parameters.xlsx"))
    params$sample_file <- file.path(ext_path, params$sample_file)
    params$ctrl_file <- file.path(ext_path, params$ctrl_file)
    params
}

#' Save an Example Parameters Excel Spreadsheet
#'
#' Loads the example parameters table using \code{load_example_params}
#' and writes it to a specified file path as an Excel spreadsheet.
#' This is useful for users to modify a template parameter file.
#'
#' @param path A \code{character} string specifying the writable
#'   file path where the Excel spreadsheet should be saved (e.g.,
#'   \code{"my_params.xlsx"}).
#' @return The function **invisibly returns the path** to the created Excel file.
#'   Its main effect is writing the file to disk at the specified path.
#' @importFrom writexl write_xlsx
#' @export
#' @examples
#' # Save the example parameters to a temporary file:
#' # save_example_params(file.path(tempdir(), "template_params.xlsx"))
save_example_params <- function(path){
    params <- load_example_params()
    writexl::write_xlsx(params, path)
}

#' Run Base Edit Detection Over Multiple Samples in Batch
#'
#' Executes the base edit detection process (\code{detect_edits}) for
#' multiple samples defined in a parameter table or Excel file.
#' This function leverages \code{BiocParallel} for parallel processing
#' to speed up the analysis of many samples.
#'
#' @param params A \code{data.frame} or a \code{character} string
#'   specifying the path to an Excel file (\code{.xlsx}) containing
#'   the input parameters. The data frame must include the following
#'   columns: \code{sample_name}, \code{sample_file}, \code{ctrl_file},
#'   \code{motif}, \code{motif_fwd}, \code{wt}, and \code{edit}.
#'   Optional columns include \code{p_value} and \code{phred_cutoff}.
#'   See \code{load_example_params} for the required format.
#' @param BPPARAM A \code{BiocParallel::bpparam()} object specifying
#'   the parallel execution environment (e.g., number of cores).
#'   Default is the object returned by \code{BiocParallel::bpparam()}.
#' @return A \code{list} of \code{MultiEditR} objects, one for each
#'   row in the \code{params} table. If an analysis fails for a
#'   specific sample, the corresponding list element will contain
#'   an error message and the flag \code{completed = FALSE}.
#' @import BiocParallel
#' @importFrom  methods is
#' @export
#' @examples
#' # 1. Create a parameter file (if one does not exist)
#' param_path <- file.path(tempdir(), "batch_params.xlsx")
#' # save_example_params(param_path)
#'
#' # 2. Load the parameters data frame
#' # params_df <- load_parameters_file(param_path)
#'
#' # 3. Run the batch detection (use default BPPARAM for safety)
#' # results_list <- detect_edits_batch(params = params_df)
#'
#' # View results for the first sample:
#' # results(results_list[[1]])
detect_edits_batch <- function(params = NULL, 
                               BPPARAM=BiocParallel::bpparam()) {
    if (is.null(params)) {
        stop("detect_edits_batch requires a parameters data.frame")
    }
    
    if (length(class(params)) == 1 && is(params, "character")){
        message("params is a string, assuming it is a path",
                "to an xlsx sheet containing the parameters.",
                "Attempting to load. ")
        params <- load_parameters_file(params)
    }
    
    fits <- BiocParallel::bplapply(seq_len(nrow(params)),
                                   FUN = function(i) { tryCatch( {
                                       fit <- detect_edits(sample_file = params$sample_file[i],
                                                           ctrl_file = params$ctrl_file[i],
                                                           motif = params$motif[i],
                                                           wt = params$wt[i],
                                                           edit = params$edit[i],
                                                           motif_fwd = ifelse(is.null(params$motif_fwd[i]), TRUE, params$motif_fwd[i]),
                                                           p_value = ifelse(is.null(params$p_value[i]), 0.01, params$p_value[i]),
                                                           phred_cutoff = ifelse(is.null(params$phred_cutoff[i]), 0.001, params$phred_cutoff[i])
                                       )
                                       fit$sample_data$sample_name <- params$sample_name[i]
                                       fit$statistical_parameters$sample_name <- params$sample_name[i]
                                       fit$sample_name <- params$sample_name[i]
                                       fit$completed <- TRUE
                                       return(fit)
                                   },
                                   error = function(e){
                                       list(sample_name = params$sample_name[i],
                                            completed = FALSE,
                                            error = e)
                                   })
                                   }, BPPARAM = BPPARAM)
    return(fits)
}  

#' Get a Single Data Frame Containing All Batch Test Results
#'
#' Consolidates the primary results table (\code{sample_data}) from
#' a list of successful \code{MultiEditR} objects (the output of
#' \code{detect_edits_batch}) into a single, comprehensive \code{data.frame}.
#' Samples that failed the analysis are automatically excluded.
#'
#' @param fits The \code{list} of \code{MultiEditR} objects returned
#'   by \code{detect_edits_batch}.
#' @return A consolidated \code{data.frame} containing all key
#'   results across all successful samples, including the sample
#'   name, motif details, editing significance, p-values, and base
#'   percentages. The columns are renamed for clarity:
#'   \itemize{
#'     \item \code{target_position_in_motif}: Position relative to the motif.
#'     \item \code{position_in_sample_trace}: Position in the raw sample trace.
#'     \item \code{sample_max_base}: The primary called base in the sample.
#'     \item \code{sample_secondary_base}: The secondary called base in the sample.
#'     \item \code{expected_base}: The expected wild-type base at the site.
#'   }
#' @importFrom dplyr mutate select tibble bind_rows
#' @importFrom plyr ldply 
#' @export
#' @examples
#' # Requires a list of fit objects from a prior batch analysis:
#' # results_list <- detect_edits_batch(params_df)
#' # batch_table <- get_batch_results_table(results_list)
get_batch_results_table <- function(fits){
    # local definitions
    index <- target_base <- max_base <- NULL
    sample_name <- sample_file <- NULL
    sample_secondary_call <- sample_secondary_base <- NULL
    motif <- passed_trimming <- target_position_in_motif <- NULL
    expected_base <- ctrl_max_base <- sample_max_base <- NULL
    edit_sig <- edit_pvalue <- edit_padjust <- NULL
    A_perc <- C_perc <- G_perc <- T_perc <- NULL
    
    # toss fits which failed
    fits <- fits[vapply(fits, FUN = "[[", "completed", FUN.VALUE = logical(1))]
    
    tbl <- list()
    for (n in seq_along(fits)) {
        tbl[[n]] <- fits[[n]][[1]] |>
                         dplyr::mutate(target_position_in_motif = target_base) |>
                         dplyr::mutate(position_in_sample_trace = index) |>
                         dplyr::mutate(sample_max_base = max_base) |>
                         dplyr::mutate(sample_secondary_base = sample_secondary_call) |>
                         dplyr::select(sample_name, passed_trimming, target_position_in_motif, motif, 
                                       expected_base, ctrl_max_base, sample_max_base, sample_secondary_base,
                                       edit_sig, edit_pvalue, edit_padjust, A_perc, C_perc, G_perc, T_perc, sample_file)
    }
    tbl <- dplyr::bind_rows(tbl)
    return(tbl)
}

#' Get a Single Data Frame Containing All Batch Critical Statistics
#'
#' Consolidates the **Filliben correlation coefficients** from the
#' statistical noise modeling (\code{statistical_parameters}) for
#' all successful samples processed in a batch. The Filliben
#' coefficient measures the linearity of the residuals' QQ-plot
#' and indicates the goodness-of-fit of the ZAGA model to the
#' background noise.
#'
#' @param fits The \code{list} of \code{MultiEditR} objects returned
#'   by \code{detect_edits_batch}. Only samples flagged as
#'   \code{completed = TRUE} are included.
#' @return A consolidated \code{data.frame} where each row
#'   corresponds to a successful sample, and the columns show the
#'   Filliben correlation coefficient for each base's ZAGA model
#'   (e.g., \code{A fillibens coef.}, \code{C fillibens coef.}, etc.).
#' @importFrom dplyr mutate select bind_rows
#' @importFrom tidyr pivot_wider
#' @importFrom plyr ldply 
#' @export
#' @examples
#' # Requires a list of fit objects from a prior batch analysis:
#' # results_list <- detect_edits_batch(params_df)
#' # stats_table <- get_batch_stats_table(results_list)
get_batch_stats_table <- function(fits) {
    # local definitions
    base <- sample_name <- fillibens <- NULL
    
    # toss fits which failed
    fits <- fits[vapply(fits, FUN = "[[", "completed", FUN.VALUE = logical(1))]
    
    tbl <- list()
    for (n in seq_along(fits)) {
        tbl <- fits[[n]][[2]] |>
                         dplyr::mutate(base = paste0(base, " fillibens coef.")) |>
                         dplyr::select(sample_name, base, fillibens) |>
                         tidyr::pivot_wider(names_from = base, values_from = fillibens)
    }
    tbl <- bind_rows(tbl)
    
    return(tbl)
}

#' Create an HTML Batch Report for MultiEditR Results
#'
#' Generates a comprehensive HTML report summarizing the results of a
#' batch base editing analysis performed by \code{detect_edits_batch}.
#' The report uses an R Markdown template to include summary tables,
#' statistics, and visualizations similar to the package's Shiny app.
#'
#' @param batch_results The \code{list} of \code{MultiEditR} objects
#'   returned by \code{detect_edits_batch}.
#' @param params The \code{data.frame} of input parameters that was
#'   originally provided to \code{detect_edits_batch}.
#' @param path A \code{character} string specifying the file path and
#'   name for the output HTML report. Defaults to
#'   \code{"./multiEditR_batch_results.html"}.
#' @export
#' @return The function invisibly returns the file path of the
#'   generated HTML report. It primarily creates the file on disk.
#' @importFrom rmarkdown render
#' @examples
#' # Requires a list of fit objects and the original parameters:
#' # results_list <- detect_edits_batch(params_df)
#' # report_path <- file.path(tempdir(), "my_analysis_report.html")
#' # create_MultiEditR_report(results_list, params_df, path = report_path)
create_MultiEditR_report <- function(batch_results, params, 
                                     path = "./multiEditR_batch_results.html") {
    
    template <- system.file(package = "MultiEditR", "batch_report_template.Rmd")
    
    message("Rendering document, this could take a couple of minutes. ",
            "When it is done, open the html in a web browser.")
    
    rmarkdown::render(template,
                      params = list(params.tbl = params,
                                    results.list = batch_results),
                      output_dir = dirname(path),
                      output_file = basename(path), 
                      envir = new.env())
}
