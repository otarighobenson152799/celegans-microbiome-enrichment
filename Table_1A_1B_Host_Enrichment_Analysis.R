suppressPackageStartupMessages({ library(dplyr) })

PROC <- "/mnt/d/ce_microbiome/data_compost/processed"
RES  <- "/mnt/d/ce_microbiome/results_compost"
dir.create(RES, showWarnings = FALSE, recursive = TRUE)

clr <- read.delim(file.path(PROC, "clr_matrix.tsv"),
                  row.names = 1, check.names = FALSE)
tax <- read.delim(file.path(PROC, "taxa_classified.tsv"), stringsAsFactors = FALSE)

cat(sprintf("[enrich] CLR: %d ASVs x %d samples\n", nrow(clr), ncol(clr)))
cat(sprintf("[enrich] taxonomy rows: %d (columns: %s)\n",
            nrow(tax), paste(colnames(tax), collapse=",")))

if (!("Genus" %in% colnames(tax)) || !("sequence" %in% colnames(tax))) {
  stop("taxa_classified.tsv needs 'sequence' and 'Genus' columns")
}
seq2genus <- setNames(tax$Genus, tax$sequence)

samples <- colnames(clr)
info <- data.frame(
  sample = samples,
  env    = sub("_.*$", "", samples),
  type   = sub("^[^_]+_([WS])[0-9]+$", "\\1", samples),
  stringsAsFactors = FALSE
)
cat("\n[enrich] Sample design:\n")
print(table(info$env, info$type))

rel_proxy <- exp(clr)
rel <- sweep(rel_proxy, 2, colSums(rel_proxy), "/")

enrich_list <- list()
for (env in sort(unique(info$env))) {
  w <- info$sample[info$env == env & info$type == "W"]
  s <- info$sample[info$env == env & info$type == "S"]
  if (length(w) == 0 || length(s) == 0) next
  worm_mean <- rowMeans(rel[, w, drop = FALSE])
  soil_mean <- rowMeans(rel[, s, drop = FALSE])
  eps <- 1e-6
  lfc <- log2((worm_mean + eps) / (soil_mean + eps))
  enrich_list[[env]] <- lfc
  cat(sprintf("[enrich] %s: %d W / %d S | lfc range %.2f to %.2f\n",
              env, length(w), length(s), min(lfc), max(lfc)))
}

enrich_mat <- do.call(cbind, enrich_list)

genus <- seq2genus[rownames(enrich_mat)]
keep  <- !is.na(genus) & nzchar(genus)
enrich_mat <- enrich_mat[keep, , drop = FALSE]
genus      <- genus[keep]
cat(sprintf("\n[enrich] %d/%d ASVs have genus assignment\n",
            sum(keep), length(keep)))

df <- data.frame(genus = genus, enrich_mat,
                 check.names = FALSE, stringsAsFactors = FALSE)
agg <- aggregate(. ~ genus, data = df, FUN = mean, na.rm = TRUE)
rownames(agg) <- agg$genus
agg <- agg[, -1, drop = FALSE]

write.table(agg, file.path(RES, "enrichment_by_genus.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

enriched_in_all <- apply(agg, 1, function(x) all(x > 1))
conserved <- agg[enriched_in_all, , drop = FALSE]
if (nrow(conserved) > 0) {
  conserved <- conserved[order(-rowMeans(conserved)), , drop = FALSE]
}
write.table(conserved, file.path(RES, "conserved_signature.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

cat(sprintf("\n[enrich] CONSERVED (>1 lfc in all %d envs): %d genera\n",
            ncol(agg), nrow(conserved)))
if (nrow(conserved) > 0) print(round(conserved, 2))

cat("\n[enrich] TOP 15 genera by MEAN log2FC (open-mind check):\n")
ord <- order(-rowMeans(agg))
top15 <- agg[ord[1:min(15, nrow(agg))], , drop = FALSE]
print(round(top15, 2))
cat("\n[enrich] Wrote enrichment_by_genus.tsv + conserved_signature.tsv\n")
