# Explicit input paths for ChAIR single-cell and ensemble tracks.
readViewerManifest <- function(manifest, output) {
  manifest <- normalizePath(manifest, mustWork = TRUE)
  paths <- read.delim(manifest, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("track", "path") %in% names(paths)) || anyDuplicated(paths$track)) stop("Manifest requires unique track and path columns")
  required <- c("genes", "rna", "atac", "peaks", "pet", "rna_rpkm", "atac_rpkm", "loops")
  if (!all(required %in% paths$track)) stop("Missing tracks: ", paste(setdiff(required, paths$track), collapse = ", "))
  files <- setNames(paths$path, paths$track)
  files <- vapply(files, function(path) {
    if (!grepl("^(/|[A-Za-z]:)", path)) path <- file.path(dirname(manifest), path)
    normalizePath(path, mustWork = TRUE)
  }, character(1))
  read.track <- function(key, columns, optional = FALSE) {
    text.columns <- c("chrom", "chrom1", "chrom2", "barcode", "name", "transcription", "gene", "strand", "strand1", "strand2", "type", "itemRgb")
    empty <- as.data.frame(setNames(lapply(columns, function(x) if (x %in% text.columns) character(0) else numeric(0)), columns))
    if (!key %in% names(files)) {
      if (optional) return(empty)
      stop("Missing track: ", key)
    }
    connection <- if (grepl("\\.gz$", files[[key]], ignore.case = TRUE)) gzfile(files[[key]], "rt") else file(files[[key]], "rt")
    on.exit(close(connection))
    first <- readLines(connection, n = 1)
    if (!length(first)) return(empty)
    close(connection)
    connection <- if (grepl("\\.gz$", files[[key]], ignore.case = TRUE)) gzfile(files[[key]], "rt") else file(files[[key]], "rt")
    result <- read.table(connection, header = FALSE, stringsAsFactors = FALSE, sep = "\t", quote = "", comment.char = "")
    if (ncol(result) != length(columns)) stop(key, " requires ", length(columns), " columns; found ", ncol(result))
    names(result) <- columns
    result
  }
  interval <- c("chrom", "start", "end", "barcode", "count")
  coverage <- c("chrom", "start", "end", "value")
  loop7 <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "count")
  loop10 <- c(loop7, "cov1", "cov2", "score")
  gene <- read.track("genes", c("chrom", "start", "end", "transcription", "gene", "score", "strand", "type"))
  gene$ID <- paste(gene$gene, gene$transcription, sep = "-")
  rna <- read.track("rna", interval)
  atac <- read.track("atac", interval)
  pet <- read.track("pet", c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "barcode", "name", "score", "strand1", "strand2"))
  pet <- pet[pet$chrom1 == pet$chrom2 & pet$end2 - pet$start1 >= 1000, , drop = FALSE]
  loops <- read.track("loops", loop10)
  loops <- loops[loops$chrom1 == loops$chrom2 & loops$count > 0, , drop = FALSE]
  rna.cov <- read.track("rna_rpkm", coverage)
  atac.cov <- read.track("atac_rpkm", coverage)
  empty.loop <- loops[FALSE, loop7, drop = FALSE]
  chromstat <- read.track("chromstat", c("chrom", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb"), optional = TRUE)
  split.chrom <- function(table, key = "chrom") {
    chromosomes <- unique(c(gene$chrom, rna$chrom, atac$chrom, pet$chrom1, loops$chrom1))
    result <- split(table, table[[key]])
    for (chrom in setdiff(chromosomes, names(result))) result[[chrom]] <- table[FALSE, , drop = FALSE]
    result
  }
  tracks <- list(genelist = split.chrom(gene), scATACfragmentlist = split.chrom(atac),
    scATACpeaklist = split.chrom(read.track("peaks", c("chrom", "start", "end"))),
    scATACARCCOVlist = split.chrom(atac.cov), scATACARCCOVrpkmlist = split.chrom(atac.cov),
    scATACPETlist = split.chrom(pet, "chrom1"), scATAClooplist = split.chrom(empty.loop, "chrom1"),
    scATAClooptwopeaklist = split.chrom(empty.loop, "chrom1"), scATACloopgivenanchorlist = split.chrom(loops, "chrom1"),
    scATACGEXCOVlist = split.chrom(rna.cov), scATACGEXCOVrpkmlist = split.chrom(rna.cov),
    scATACGEXdupCOVlist = split.chrom(rna.cov), scATACGEXdupCOVrpkmlist = split.chrom(rna.cov),
    scRNAfragmentlist = split.chrom(rna), scRNAdupfragmentlist = split.chrom(rna), chromstatlist = split.chrom(chromstat))
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  save(list = names(tracks), envir = list2env(tracks), file = output)
  message("Saved ", normalizePath(output), "; RNA/ATAC/PET cells: ", length(unique(rna$barcode)), "/", length(unique(atac$barcode)), "/", length(unique(pet$barcode)))
}
