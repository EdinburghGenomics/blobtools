#!/usr/bin/env Rscript

# Script copied from qc_tools_python/scripts makes CSV tables suitable for QC reports.
#  Usage: parseBlobTable.R {output.csv} {params.pct_limit} {input.stats.txt} ...

cl <- commandArgs(trailingOnly = TRUE)

outfile <- cl[1]
percent_limit <- as.numeric(cl[2])
files <- cl[3:length(cl)]

# Get the main table from the file

read.blob_table <- function(file, comment.char = "##", ...) {
  lines <- readLines(file, warn = FALSE)

  # Get the input file names and derive names from the path

  bamlines <- lines[grep('## bam[0-9]', lines)]
  bamnos <- sub('.*(bam[0-9]+).*', '\\1', bamlines)
  bamnames <- sub('.*\\/([^\\/]+__[^\\/]+)\\/.*', '\\1', bamlines)

  # Now read the data lines

  lines <- lines[grep(paste0("^", comment.char), lines, invert = TRUE)]
  lines <- sub("^# ", '', lines)
  data <- read.delim(text = paste(lines, collapse = "\n"), sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

  # Substitute the bam number for names

  for (i in 1:length(bamnos)){
    colnames(data) <- gsub(paste0('^', bamnos[i], '_'), paste0(bamnames[i], '_'), colnames(data))
  }

  data
}

tables <- lapply(files, function(infile){

  blob_data <- read.blob_table(infile)

  # Get the actual read numbers so we can calculate a more precise percentage
  # than Blobtools provides. Then we can see how close to 0 some of these 0
  # percentages are.

  blob_data_reads <- blob_data[, grep('.*read_map$', colnames(blob_data)), drop = FALSE]
  blob_data_reads <- blob_data_reads[, grep('covsum', colnames(blob_data_reads), invert = TRUE), drop = FALSE]
  colnames(blob_data_reads) <- sub('_read_map', '', colnames(blob_data_reads))
  blob_data_reads <- data.frame(apply(blob_data_reads, 2, function(x) as.numeric(gsub(',', '', x))), row.names = rownames(blob_data_reads), check.names = FALSE)

  # Blobtools' percentage values are relative to the origninal read counts
  # (mapped and unmapped). If we want to re-calculate the percentages we need to
  # re-infer that original read cout. We can use the 'all' percent to do that.
  # It won't be completely right due to rounding error in the 'all' percentage
  # itself. But it should be close enough. Of course we might find it more
  # useful in future to have percentages relative to the mapped reads.

  blob_data_percent <- blob_data[, grep('.*read_map_p$', colnames(blob_data)), drop = FALSE]
  blob_data_percent <- blob_data_percent[, grep('covsum', colnames(blob_data_percent), invert = TRUE), drop = FALSE]
  colnames(blob_data_percent) <- sub('_read_map_p', '', colnames(blob_data_percent))
  blob_data_percent <- data.frame(apply(blob_data_percent, 2, function(x) as.numeric(sub('%', '', x))), row.names = rownames(blob_data_percent), check.names = FALSE)

  # Do the percent re-calculation

  totals <- blob_data_reads['all', ] * (1/(blob_data_percent['all',] / 100))
  new_blobdata_percent <- blob_data_percent
  used_names <- rownames(blob_data_percent)[! rownames(blob_data_percent) %in% c('all', 'other')]

  if (ncol(blob_data_percent) == 1){
    new_blobdata_percent[used_names,] <- ((blob_data_reads[used_names, , drop = FALSE] / totals)*100)[used_names,, drop = FALSE]
  }else{
    new_blobdata_percent[used_names,] <- do.call(rbind, apply(blob_data_reads[used_names, , drop = FALSE], 1, function(x) (x/ totals)*100))
  }

  # Take out the 'all' and 'no hit'. Note that we leave this until the end so that
  # the 'apply' above doesn't drop all our dimensions when there's a single row.

  new_blobdata_percent[! rownames(blob_data) %in% c('all', 'no-hit', 'undef'), ,drop=FALSE]
})

# Merge the table from the different files, which can have different rows

all_rows <- unique(unlist(lapply(tables, rownames)))
merged_table <- t(do.call(cbind, lapply(tables, function(x) x[all_rows,, drop = FALSE])))
merged_table[is.na(merged_table)] <- 0
colnames(merged_table) <- all_rows

# Remove anything at less than 1% in all samples

max_percent <- signif(max(merged_table), 4)

merged_table <- merged_table[, apply(merged_table, 2, function(x) any(x >= percent_limit)) , drop = FALSE]

if (ncol(merged_table) > 0){

  # Sort by most abundant

  merged_table <- merged_table[, order(colSums(merged_table), decreasing = TRUE), drop = FALSE]

  # Prepare for output

  final <- cbind(rownames(merged_table), round(merged_table, 2))
  colnames(final)[1] <- 'Library ID'

  # Write those outputs

  write.csv(final, file = outfile, row.names = FALSE, quote = FALSE)
}else{
  taxlevel <- unlist(strsplit(basename(files[1]), '\\.'))[2]

  fileConn<-file(outfile)
  writeLines(paste0('No ', taxlevel, ' is represented by at least ', percent_limit, '%', ' of reads (max ', max_percent, '%)'), fileConn)
  close(fileConn)
}
