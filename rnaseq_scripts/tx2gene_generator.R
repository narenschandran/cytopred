
cdna <- readLines("Homo_sapiens.GRCh38.cdna.all.fa")

subset_lines <- cdna[grepl("ENS", cdna)]
split_lines <- strsplit(subset_lines, "[>/ ]")

tx2gene <- do.call("rbind", lapply(split_lines, function(sl) {
    c(sl[2], strsplit(sl[grepl("gene_symbol", sl)], ":")[[1]][2])
}))

colnames(tx2gene) <- c("transcript_id", "gene_id")

write.csv(tx2gene, "tx2gene.csv", row.names=F)
