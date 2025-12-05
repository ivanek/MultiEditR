# MultiEditR 

MultiEditR detects and quantifies base edits from Sanger sequencing traces using motif-guided statistical detection of relative trace heights.
The package loads sample Sanger sequence(s), control sequences or FASTA files, and identifies the presence and magnitude of edits within specific motifs.
Algorithms were developed by Mitch Kluesner and implemented by Jeremy Chacón.

## Installation

```{r}
Install the development version from GitHub:
install.packages("remotes")
remotes::install_github("ivanek/MultiEditR")
```
## Overview

MultiEditR provides:
 - Detection and quantification of base edits within user-defined motifs.
 - Support for Sanger .ab1 files or FASTA control sequences.
 - Automatic reverse-complement checking and correction for control sequences.
 - Statistical modeling of trace height data to call edits.
 - Visualizations: chromatograms, editing barplots, raw trace signals.
 - A batch mode to process many samples and generate an HTML report.

## Example

Below is a minimal example using two bundled Sanger traces.
We test for "A" → "G" edits within a specified motif.

```{r}
library(MultiEditR)

sample_file <- system.file("extdata", "RP272_cdna_wt.ab1", package = "MultiEditR")
ctrl_file   <- system.file("extdata", "RP272_cdna_ko.ab1", package = "MultiEditR")

motif <- "AGTAGCTGGGATTACAGATG"
wt    <- "A"
edit  <- "G"

fit <- detect_edits(
  sample_file  = sample_file,
  ctrl_file    = ctrl_file,
  p_value      = 0.0001,
  phred_cutoff = 0.0001,
  motif        = motif,
  motif_fwd    = TRUE,
  wt           = wt,
  edit         = edit
)
```

### Working with results

```{r}
# Retrieve full results table
tbl <- results(fit)

# Save the table (optional)
tmpDir <- tempdir()
writexl::write_xlsx(tbl, file.path(tmpDir, "my_results.xlsx"))
```

### Visualization

```{r}
plot_sample_chromatogram(fit)
plot_control_chromatogram(fit)
plot_raw_sample(fit)
plot_editing_barplot(fit)
```

## Batch Mode

MultiEditR can run many samples at once using a parameter spreadsheet, then create a combined HTML report.

### Load example parameters

```{r}
params <- load_example_params()
head(params)
```

### Process all samples

```{r}
fits <- detect_edits_batch(params)
```

### Summaries across all samples

```{r}
data.tbl  <- get_batch_results_table(fits)
stats.tbl <- get_batch_stats_table(fits)
```

### Inspect a single fit

```{r}
fit1 <- fits[[1]]
results(fit1)
plot_sample_chromatogram(fit1)
```

### Generate an HTML report

```{r}
create_MultiEditR_report(fits, params, file.path(tmpDir, "my_html_report.html"))
```

## Citation

If you use MultiEditR in published work, please cite the associated repositories and manuscripts.

### Original repositories: 
- [jeremymchacon/MultiEditR_Package](https://github.com/jeremymchacon/multiEditR_Package)
- [MoriarityLab/multiEditR.pckg](https://github.com/MoriarityLab/multiEditR.pckg)

### Original manuscript:
Kluesner MG et al. MultiEditR: The first tool for the detection and 
quantification of RNA editing from Sanger sequencing demonstrates comparable 
fidelity to RNA-seq. Mol Ther Nucleic Acids. 2021 Jul 21;25:515-523. 
doi: [10.1016/j.omtn.2021.07.008](https://doi.org/10.1016/j.omtn.2021.07.008). 
PMID: [34589274](https://pubmed.ncbi.nlm.nih.gov/34589274/); 
PMCID: [PMC8463291](https://pmc.ncbi.nlm.nih.gov/articles/PMC8463291/).
