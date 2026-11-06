#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
if (length(args) < 3L || length(args) > 6L) {
  stop("Usage: Rscript plotData.R INPUT_RDATA CHROM:START-END OUTPUT.pdf [CELLS=800] [PET_COUNT=3] [VC_SCORE=0.001]")
}
script <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())])
source(file.path(dirname(normalizePath(script)), "source.r"))
input <- normalizePath(args[1], mustWork = TRUE)
region <- args[2]
output <- args[3]
cells <- if (length(args) >= 4L) as.integer(args[4]) else 800L
pet.count <- if (length(args) >= 5L) as.numeric(args[5]) else 3
vc.score <- if (length(args) >= 6L) as.numeric(args[6]) else 0.001
if (is.na(cells) || cells < 1 || any(!is.finite(c(pet.count, vc.score)))) stop("Invalid plotting parameters")
if (!grepl("^.+:[0-9]+-[0-9]+$", region)) stop("Region must use CHROM:START-END")
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
format <- tolower(tools::file_ext(output))
if (!format %in% c("pdf", "png")) stop("Output extension must be .pdf or .png")
files <- basename(input)
data_list <- loadData(file = files, path = paste0(dirname(input), "/"))
run(region, file = files, direction = FALSE, cell_number = cells,
    PETcount = pet.count, VCscore = vc.score, type = "PET", condition = "PETDNARNA",
    sample_num = -1, clean = FALSE, showPETnum = TRUE, outfmt = format, output = output)
