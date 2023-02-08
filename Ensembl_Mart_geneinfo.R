
mart_geneinfo <- function() {
  
  library(biomaRt)
  
  #import Biomart annotations
  mart=useMart('ENSEMBL_MART_ENSEMBL',
               dataset='hsapiens_gene_ensembl', 
               host="https://feb2014.archive.ensembl.org",
               path="/biomart/martservice")
  
  # retrieve features
  bm<-getBM(attributes=c('chromosome_name',
                         'start_position',
                         'end_position',
                         'strand',
                         'ensembl_gene_id',
                         'ensembl_transcript_id',
                         'hgnc_symbol',
                         'external_transcript_id',
                         'external_gene_id',
                         'gene_biotype',
                         'transcript_biotype',
                         'transcript_start'),
            mart = mart)
  # change chr names
  bm$chr <- paste0("chr",bm$chromosome_name)
  bm$chromosome_name <- NULL
  # change strand symbols
  bm$strand[which(bm$strand == -1)] <- "-"
  bm$strand[which(bm$strand == 1)] <- "+"
  # add an ENSG:ENST column to be able to query the EP output files
  bm$ENSG_ENST <- paste0(bm$ensembl_gene_id,":", bm$ensembl_transcript_id)
  # dataframe for gene of interest
  return(bm <- as_tibble(bm))
  detach("package:biomaRt", unload=TRUE)
  rm(mart)
}