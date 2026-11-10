args <- commandArgs(trailingOnly=TRUE)
repo <- normalizePath(if(length(args))args[[1]] else getwd(),mustWork=TRUE)
out <- file.path(repo,'.local','checks');dir.create(out,recursive=TRUE,showWarnings=FALSE)
files <- list.files(file.path(repo,'analysis'), pattern='\\.R$', recursive=TRUE, full.names=TRUE)
for (p in files) parse(file=p)
source(file.path(repo,'analysis/bulk_rnaseq/04_collect_rsem_matrices.R'))
temp <- tempfile('rsem_fixture_');dir.create(temp)
input <- file.path(temp,'input');dir.create(input)
a <- data.frame(gene_id=c('g1','g2'), expected_count=c(10,20), TPM=c(1,2))
b <- data.frame(gene_id=c('g2','g1'), expected_count=c(300,400), TPM=c(30,40))
write.table(a,file.path(input,'A.genes.results'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(b,file.path(input,'B.genes.results'),sep='\t',quote=FALSE,row.names=FALSE)
x <- collect_rsem_matrices(input,file.path(temp,'output'))
stopifnot(identical(rownames(x$expected_counts),c('g1','g2')),
          identical(as.numeric(x$expected_counts$A),c(10,20)),
          identical(as.numeric(x$expected_counts$B),c(400,300)),
          identical(as.numeric(x$tpm$A),c(1,2)),
          identical(as.numeric(x$tpm$B),c(40,30)))
b$gene_id <- c('g1','g3')
write.table(b,file.path(input,'B.genes.results'),sep='\t',quote=FALSE,row.names=FALSE)
bad <- try(collect_rsem_matrices(input,file.path(temp,'bad_output')),silent=TRUE)
stopifnot(inherits(bad,'try-error'))
writeLines(c(paste('R source files parsed:',length(files)),
 'PASS: expected_count and TPM collected consistently for both samples.',
 'PASS: reordered genes aligned by gene_id.',
 'PASS: mismatched gene sets rejected.',
 'No biological analysis was rerun.'),file.path(out,'repository_R_checks.txt'))
cat('PASS:',length(files),'R files parsed; RSEM fixtures passed.\n')
