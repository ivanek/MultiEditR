#' Detect Base Editing Events in Sanger Sequencing Data
#'
#' The primary function for detecting and quantifying base editing events
#' by aligning a sample sequence to a control sequence, identifying
#' target motif locations, and using a Zero-Adjusted Gamma (ZAGA)
#' distribution model to establish a null distribution for background
#' noise, followed by p-value calculation and adjustment for edits.
#'
#' The algorithm performs the following key steps:
#' 1. Load sample and control sequences, applying reverse complement
#'    to the control if it improves alignment score with the sample.
#' 2. Perform global pairwise alignment to map control positions onto
#'    sample positions.
#' 3. Use Mott's algorithm to determine and flag low-quality regions
#'    for trimming in the sample.
#' 4. Find all perfect matches of the target \code{motif} in the control
#'    sequence.
#' 5. Use non-motif and non-trimmed sample data to fit the ZAGA
#'    distribution for base-specific noise (the "null distribution").
#' 6. Calculate BH-adjusted p-values for all wild-type (\code{wt}) to
#'    edited (\code{edit}) base transitions at motif locations.
#'
#' @param sample_file A \code{character} string path to the sample
#'   Sanger sequence file (\code{.ab1}).
#' @param ctrl_file A \code{character} string path to the control
#'   sequence file (\code{.ab1} or \code{.fa}).
#' @param motif A \code{character} string holding the motif of interest.
#' @param motif_fwd A \code{logical}. If \code{TRUE} (default), the
#'   motif is matched as-is; if \code{FALSE}, its reverse complement
#'   is used for matching against the control.
#' @param wt A \code{character} string specifying the wild-type base
#'   (e.g., "C") to be tested for editing in motif locations.
#' @param edit A \code{character} string specifying the predicted
#'   edited base (e.g., "T" for C-to-T).
#' @param phred_cutoff A \code{numeric} cutoff for Mott trimming.
#'   Lower values are more stringent. Default is \code{0.00001}.
#' @param p_value A \code{numeric} cutoff for Benjamini-Hochberg
#'   adjusted significance. Default is \code{0.01}.
#' @return A \code{list} of class \code{multieditR} containing the
#'   full analysis results, sequences, and intermediate data,
#'   including \code{sample_data} (the primary results table),
#'   \code{statistical_parameters} (ZAGA fits), and sequencing
#'   objects.
#' @importFrom readr read_lines
#' @importFrom stringr str_split_1
#' @importFrom dplyr group_by ungroup select left_join distinct case_when summarize pull n rename mutate arrange filter inner_join ungroup add_rownames if_else
#' @importFrom sangerseqR readsangerseq makeBaseCalls
#' @importFrom Biostrings countPattern matchPattern 
#' @importFrom pwalign pairwiseAlignment alignedSubject alignedPattern
#' @importFrom tibble rownames_to_column tibble as_tibble
#' @export
#' @examples
#' # Use package-supplied example files
#' sample_file <- system.file("extdata", "RP272_cdna_wt.ab1",
#'   package = "MultiEditR")
#' ctrl_file <- system.file("extdata", "RP272_cdna_ko.ab1",
#'   package = "MultiEditR")
#' motif <- "AGTAGCTGGGATTACAGATG"
#' fit <- detect_edits(sample_file, ctrl_file,
#'   motif = motif, motif_fwd = TRUE,
#'   wt = "A", edit = "G",
#'   phred_cutoff = 0.001, p_value = 0.05)
#'
#' # Access primary results table:
#' tbl <- results(fit)
#' # writexl::write_xlsx(tbl, file.path(tempdir(), "my_results.xlsx")) 
#'
#' # Plot visualizations:
#' plot_sample_chromatogram(fit)
#' plot_raw_sample(fit) 
detect_edits <- function(sample_file, ctrl_file, motif, motif_fwd, wt, edit,
                         phred_cutoff = 0.00001, p_value = 0.01) {
    # local definitions
    A_area <- C_area <- G_area <- T_area <- NULL
    A_perc <- C_perc <- G_perc <- T_perc <- NULL
    position <- trimmed <- everything <- edit_sig <- NULL
    expected_motif <- passed_trimming <- ctrl_max_base <- NULL
    ctrl_post_aligned_index <- expected_base <- max_base <- NULL
    edit_pvalue <- edit_padjust <- NULL
    sample_primary_call <- sample_secondary_call <- control_primary_call <- NULL
    index <- is_trimmed <- Tot.Area <- max_base_height <- NULL
    sig <- base <- perc <- tally <- motif_id <- motif_found <- NULL
    
    # get the sequences
    sample_sanger <- sangerseqR::readsangerseq(sample_file)
    sample_seq <- sangerseqR::makeBaseCalls(sample_sanger) |> 
        sangerseqR::primarySeq() |>
        as.character()
    secondary_seq <- sangerseqR::makeBaseCalls(sample_sanger) |> 
        sangerseqR::secondarySeq() |>
        as.character()
    
    # get the control sequence and put it into "ctrl_seq"
    if (is_file_ab1(ctrl_file)){
        ctrl_sanger <- sangerseqR::readsangerseq(ctrl_file)
        base_calls <- sangerseqR::makeBaseCalls(ctrl_sanger)
        ctrl_seq <- base_calls |> 
            sangerseqR::primarySeq() |> 
            as.character()
    }else{
        ctrl_sanger <- NA
        fasta_lines <- read_lines(ctrl_file)
        ctrl_seq <- paste0(fasta_lines[2:length(fasta_lines)], collapse = "")
    }
    
    # figure out if the control should be rev-com
    ctrl_is_revcom <- FALSE
    if (is_revcom_ctrl_better(sample_seq, ctrl_seq)){
        message("Control sequence aligns better to sample sequence",
                "when rev-com.", "Applying revcom to control sequence.")
        ctrl_seq <- revcom(ctrl_seq)
        ctrl_is_revcom <- TRUE
    }
    
    # revcom motif if necessary
    motif_orig <- motif
    if (!motif_fwd){
        motif <- revcom(motif)
    }
    
    # make sure at least one motif is findable in the control
    n_motif_alignments_to_ctrl <-Biostrings::countPattern(pattern = Biostrings::DNAString(motif), 
                                                          subject = Biostrings::DNAString(ctrl_seq),
                                                          max.mismatch = 1)
    if (n_motif_alignments_to_ctrl == 0){
        stop("Motif not found in control sequence while allowing 1 mismatch.",
             "Are you sure you have motif_fwd correct?")
    }
    
    
    #### this will hold our main sequence table
    sample_df <- make_sample_df(sample_sanger) |>
        dplyr::rename(raw_sample_position = "position")
    
    
    # apply the control sequence to the sample sequence. 
    control_alignment <- pwalign::pairwiseAlignment(pattern = DNAString(ctrl_seq), 
                                                    subject = DNAString(sample_seq))
    
    aligned_sample <- as.character(pwalign::alignedSubject(control_alignment))
    aligned_sample <- stringr::str_split_1(aligned_sample, "")
    aligned_control <- as.character(pwalign::alignedPattern(control_alignment))
    aligned_control <- stringr::str_split_1(aligned_control, "")
    
    # figure out the raw sample position
    raw_sample_position <- cumsum(aligned_sample != "-")
    raw_control_position <- cumsum(aligned_control != "-")
    
    # Create a tibble with the alignment and relative positions
    alignment_df <- tibble(raw_sample_position,
                           sample_primary_call = aligned_sample, 
                           raw_control_position,
                           control_primary_call = aligned_control)
    
    # bind the sample_df, which holds quality scores and percentages
    alignment_df <- dplyr::left_join(alignment_df, sample_df, by="raw_sample_position")
    
    # also bind the secondary call in the sample
    secondary_call <- tibble(sample_secondary_call = stringr::str_split_1(secondary_seq, ""),
                             raw_sample_position = seq_len(nchar(secondary_seq)))
    alignment_df <- dplyr::left_join(alignment_df, secondary_call, by="raw_sample_position")
    
    # find out what should be trimmed using Mott's algo, which Mitch implemented
    trim_points <- get_trim_points(sample_file, cutoff = phred_cutoff)
    
    if (trim_points[1] == -1){
        warning("Low quality scores.",
                "Trimming would remove most of sequence.",
                "Skipping trimming.")
        start_pos <- -1
        end_pos <- Inf
    }else{
        start_pos <- trim_points[1]
        end_pos <- trim_points[2]    
    }
    
    alignment_df$trimmed <- dplyr::case_when(alignment_df$raw_sample_position < start_pos ~ TRUE,
                                             alignment_df$raw_sample_position >= end_pos ~ TRUE,
                                             .default = FALSE)
    
    # find all alignments of the motif to the control sequence, now we allow zero
    # mismatches because if the motif is tiny, we don't want clutter
    
    motif_alignments_to_ctrl <- Biostrings::matchPattern(pattern = DNAString(motif), 
                                                         subject = DNAString(ctrl_seq), 
                                                         max.mismatch = 0)
    motif_alignments <- motif_alignments_to_ctrl@ranges |> 
        as.data.frame() |>
        as_tibble()
    
    message("Total of ", nrow(motif_alignments),
            " alignment(s) of motif to control sequence found.")
    
    # add in the motif positions. -1 for non-motif, otherwise counting from 1 up for each detection
    alignment_df$motif_found <- -1
    
    for (i in seq_len(nrow(motif_alignments))) {
        mstart <- motif_alignments$start[i]
        mend <- motif_alignments$end[i]
        mrange <- seq(from = mstart, to = mend, by = 1)
        alignment_df <- alignment_df |>
            mutate(motif_found = if_else(raw_control_position %in% mrange, i, motif_found),
                   motif = motif_found)
    }  
    
    motif_part_of_sample <- alignment_df |>
        filter(motif_found != -1)
    
    ## lets check if any of the motif is being suggested for trimming, and warn if so
    if (any(motif_part_of_sample$trimmed)){
        warning("Part of the sample sanger sequence overlapping the motif",
                "is low quality and flagged for trimming.",
                "It is not being trimmed during edit detection,",
                "but this may affect results")
    }
    
    # this is all just renames to be compatible with make_ZAGA_df
    motif_part_of_sample <- motif_part_of_sample |>
        mutate(expected_motif = control_primary_call,
               ctrl_max_base = expected_motif ,
               index = raw_sample_position,
               ctrl_index = raw_control_position)
    # this should be used for generating the NULL; only keep the trimmed part
    nonmotif_part_of_sample <- alignment_df |>
        dplyr::filter(motif_found == -1) |>
        dplyr::filter(!trimmed)
    
    # pass it to the ZAGA function
    zaga_parameters <- make_ZAGA_df(nonmotif_part_of_sample, p_adjust = p_value) |>
        dplyr::mutate(sample_file = sample_file)
    
    output_stats <- calculate_edit_pvalue(motif_part_of_sample,
                                          zaga_parameters, wt, edit, p_value)
    
    # rearrange results to match previous version
    sample_data <- output_stats |>
        dplyr::left_join(alignment_df |> select(-motif),
                         by = c("raw_sample_position", "sample_primary_call", 
                                "raw_control_position", "control_primary_call", 
                                "A", "C", "G", "T", "max_base", 
                                "A_area", "C_area", "G_area", "T_area", 
                                "A_perc", "C_perc", "G_perc", "T_perc", "sample_secondary_call", 
                                "trimmed", "motif_found"))  |>
        dplyr::select(-motif) |>
        dplyr::mutate(motif = motif) |>
        tibble::rownames_to_column(var = "target_base") |> 
        dplyr::mutate(target_base = as.numeric(target_base)) |>
        dplyr::filter(!is.na(edit_sig)) |>
        dplyr::mutate(ctrl_max_base = expected_motif) |>
        dplyr::mutate(sample_file = sample_file) |>
        dplyr::mutate(expected_base = expected_motif) |>
        dplyr::mutate(sample_file = sample_file) |>
        dplyr::mutate(passed_trimming = !trimmed) |>
        dplyr::select(passed_trimming, target_base, motif, ctrl_max_base, 
                      expected_base, max_base, sample_secondary_call,
                      A_perc, C_perc, G_perc, T_perc, 
                      edit_pvalue, edit_padjust, edit_sig, index, sample_file)
    
    if (!motif_fwd){
        target_base <- as.numeric(sample_data$target_base)
        sample_data$target_base <- nchar(motif) - (target_base - 1)
    }
    
    raw_sample_df <- alignment_df |>
        dplyr::mutate(Tot.Area = A_area + C_area + G_area + T_area) |>
        dplyr::group_by(raw_sample_position) |>
        dplyr::mutate(max_base_height = max(A_area, C_area, G_area, T_area)) |> 
        dplyr::mutate(index = raw_sample_position) |>
        dplyr::mutate(is_trimmed = trimmed) |>
        ungroup() |>
        select(raw_sample_position, raw_control_position, sample_primary_call,
               control_primary_call, sample_secondary_call, motif, is_trimmed, 
               A_area, C_area, G_area, T_area, Tot.Area, 
               A_perc, C_perc, G_perc, T_perc, 
               index, max_base, max_base_height) |>
        ungroup()
    
    output_sample_alt <- output_stats |>
        dplyr::left_join(alignment_df, by=c("raw_sample_position", "sample_primary_call", 
                                            "raw_control_position", "control_primary_call", 
                                            "A", "C", "G", "T", "max_base", 
                                            "A_area", "C_area", "G_area", "T_area", 
                                            "A_perc", "C_perc", "G_perc", "T_perc", "sample_secondary_call", 
                                            "trimmed", "motif_found", "motif")) |>
        dplyr::mutate(perc = 100 * .data[[paste0(edit, "_perc")]]) |>
        dplyr::mutate(index = raw_sample_position) |>
        dplyr::mutate(base = edit) |>
        dplyr::mutate(sig = ifelse(edit_sig, "Significant", "Non-significant")) |>
        dplyr::group_by(sig) |>
        dplyr::mutate(tally = dplyr::n()) |>
        dplyr::select(index, base, perc, sig, tally) |>
        dplyr::filter(!is.na(sig))
    
    motif_positions <- motif_part_of_sample |>
        dplyr::mutate(ctrl_post_aligned_index = index) |>
        dplyr::mutate(motif_id = 1) |> # relic from when > 1 motif could be found
        dplyr::select(ctrl_post_aligned_index, motif_id)
    
    output <- list(
        "sample_data" = sample_data,
        "statistical_parameters" = zaga_parameters,
        "sample_sanger" = sample_sanger,
        "sample_fastq" = sample_seq,
        "ctrl_sanger" = ctrl_sanger,
        "ctrl_fastq" = ctrl_seq,
        "ctrl_is_revcom" = ctrl_is_revcom, 
        "motif" = motif_orig,
        "motif_fwd" = motif_fwd,
        "expected_change" = edit,
        "intermediate_data" = list("raw_sample_df"= raw_sample_df,
                                   "sample_alt"= motif_part_of_sample |>
                                       filter(expected_motif == wt) |>
                                       mutate(ctrl_file = ctrl_file) |>
                                       mutate(sample_file = sample_file),
                                   "output_sample_alt" = output_sample_alt,
                                   "motif_positions" = motif_part_of_sample$raw_sample_position)
    )
    class(output) <- "multieditR"
    
    return(output)
}
