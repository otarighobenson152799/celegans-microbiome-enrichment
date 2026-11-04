suppressPackageStartupMessages({
  library(dada2); library(microbiome); library(yaml)
})
args <- commandArgs(trailingOnly = TRUE)
CFG  <- read_yaml(sub("--config=", "", args[grep("--config=", args)]))
RAW  <- file.path(CFG$data_raw, "fastq")
PROC <- CFG$data_processed
FILT <- "/tmp/dada2_filtered_compost_fonly"
dir.create(PROC, showWarnings = FALSE, recursive = TRUE)
dir.create(FILT, showWarnings = FALSE, recursive = TRUE)

fqs <- sort(list.files(RAW, pattern = "_1\\.fastq\\.gz$", full.names = TRUE))
snames <- sub("_1\\.fastq\\.gz$", "", basename(fqs))
names(fqs) <- snames
filtF <- file.path(FILT, paste0(snames, "_F.fastq.gz"))
names(filtF) <- snames

cat(sprintf("[compost-fonly] %d samples\n", length(fqs)))

if (all(file.exists(filtF))) {
  cat("[compost-fonly] filtered files exist, skipping filter+trim\n")
} else {
  cat("[compost-fonly] filter+trim (truncLen=0)...\n")
  out <- filterAndTrim(fqs, filtF, truncLen = 0, maxN = 0, maxEE = 2,
                       truncQ = 2, rm.phix = TRUE, compress = TRUE,
                       multithread = FALSE)
  print(out)
}

cat("[compost-fonly] learnErrors (nbases=1e7, single-thread)...\n")
errF <- learnErrors(filtF, nbases = 1e7, multithread = FALSE)

cat("[compost-fonly] dada (F-only, multithread=TRUE)...\n")
dF <- dada(filtF, err = errF, multithread = TRUE)

seqtab <- makeSequenceTable(dF)
saveRDS(seqtab, file.path(PROC, "seqtab.rds"))
cat("[compost-fonly] SAVED seqtab.rds (raw, pre-chimera)\n")

cat(sprintf("[compost-fonly] raw: %d samples x %d ASVs, total reads %d\n",
            nrow(seqtab), ncol(seqtab), as.integer(sum(seqtab))))
cat(sprintf("[compost-fonly] depth range: %d - %d (median %d)\n",
            as.integer(min(rowSums(seqtab))),
            as.integer(max(rowSums(seqtab))),
            as.integer(median(rowSums(seqtab)))))

cat("[compost-fonly] removeBimeraDenovo...\n")
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = TRUE, verbose = TRUE)
saveRDS(seqtab_nochim, file.path(PROC, "seqtab_nochim.rds"))
cat("[compost-fonly] SAVED seqtab_nochim.rds\n")

cat(sprintf("[compost-fonly] after chimera: %d samples x %d ASVs\n",
            nrow(seqtab_nochim), ncol(seqtab_nochim)))
cat(sprintf("[compost-fonly] depth range: %d - %d (median %d)\n",
            as.integer(min(rowSums(seqtab_nochim))),
            as.integer(max(rowSums(seqtab_nochim))),
            as.integer(median(rowSums(seqtab_nochim)))))

meta <- read.delim(file.path(CFG$data_raw, "ena_metadata.tsv"),
                   stringsAsFactors = FALSE)
alias_lookup <- setNames(meta$sample_alias, meta$run_accession)
srr_ids <- sub("_F\\.fastq\\.gz$", "", rownames(seqtab_nochim))
new_names <- alias_lookup[srr_ids]
if (any(is.na(new_names))) {
  cat("[compost-fonly] WARNING: no alias for:",
      paste(srr_ids[is.na(new_names)], collapse=","), "\n")
  new_names[is.na(new_names)] <- srr_ids[is.na(new_names)]
}
rownames(seqtab_nochim) <- new_names
cat("[compost-fonly] sample names:", paste(rownames(seqtab_nochim), collapse=","), "\n")

rel    <- sweep(seqtab_nochim, 1, rowSums(seqtab_nochim), "/")
pseudo <- rel + 1e-6
clr    <- log(pseudo) - rowMeans(log(pseudo))

write.table(t(clr), file.path(PROC, "clr_matrix.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

samples <- data.frame(sample = rownames(seqtab_nochim),
                      depth  = as.integer(rowSums(seqtab_nochim)))
write.table(samples, file.path(PROC, "samples.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

taxa <- data.frame(sequence = colnames(seqtab_nochim))
write.table(taxa, file.path(PROC, "taxa.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("[compost-fonly] WROTE clr_matrix.tsv, samples.tsv, taxa.tsv\n")
cat("[compost-fonly] DONE — next: classify_compost.R\n")
