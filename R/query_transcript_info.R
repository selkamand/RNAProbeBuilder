
query_transcript_sequence <- function(transcript_id, database = c('GRCh38_EnsDb_v86', "hg19_ucsc")){
  # assertions::assert_subset(database, c('GRCh38_EnsDb_v86', 'hg19_ucsc'))
  #
  # if(database == "hg19_ucsc"){
  #   rlang::check_required(TxDb.Hsapiens.UCSC.hg19.knownGene)
  #   txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  #   df_res <- GenomicFeatures::transcripts(txdb, filter(list(tx_id = c(transcript_id))), columns = c("tx_id","tx_name"))
  #   return(df_res)
  # }
  # else if(database == "GRCh38_EnsDb_v86"){
  #  edb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
  #  dna <- ensembldb:::getGenomeTwoBitFile(edb)
  # }
}


default_biomart <- function(){
  biomaRt::useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
}

#' Load Specific Ensembl Biomart
#'
#' @param GRCh reference assembly ('38' or '37')
#'
#' @return By default, returns a gene database from the latest ensembl hg38/hg37 biomart
#' @export
load_biomart <- function(GRCh = c("38", "37")){
  GRCh <- rlang::arg_match(GRCh)
  biomaRt::useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh = GRCh)
}

query_transcript_sequence_scalar_api <- function(transcript_id, ensembl = default_biomart()){
  assertions::assert_string(transcript_id)
  cdna = biomaRt::getSequence(id = transcript_id, type = "ensembl_transcript_id_version", mart = ensembl, seqType = "cdna")
  #utr3prime = biomaRt::getSequence(id = transcript_id, type = "ensembl_transcript_id", mart = ensembl, seqType = "3utr")

  seq_cdna = cdna[[1]]

  if(length(seq_cdna) == 0) seq <- NA_character_

  return(seq_cdna)
}

# query_utr_5prime_size_scalar_api <- function(transcript_id, ensembl = default_biomart()){
#   assertions::assert_string(transcript_id)
#   utr5prime <- biomaRt::getSequence(id = transcript_id, type = "ensembl_transcript_id", mart = ensembl, seqType = "5utr")
#   utr5length <- nchar(utr5prime[[1]])
#   return(utr5length)
# }
#
# query_utr_5prime_size_api <- function(transcript_id, ensembl = default_biomart()){
#   assertions::assert_character_vector(transcript_id)
#   lengths = purrr::map_dbl(transcript_id, function(id) { query_utr_5prime_size_scalar_api(id, ensembl = ensembl)})
#   return(lengths)
# }
#
query_transcript_sequence_api <- function(transcript_id, ensembl = default_biomart()){
  assertions::assert_character(transcript_id)

  sequences = purrr::map_chr(transcript_id, function(id) { query_transcript_sequence_scalar_api(id, ensembl = ensembl)})
  names(sequences) <- transcript_id

  return(sequences)
}


#' #' Extract HGVS from vep-annotated VCF
#' #'
#' #' @param vcf  path to VCF
#' #' @param transcript_type which transcript type to use - will filter for HGVSc starting with 'ENST' for ensembl and 'NM' for refst
#' #' @return
#' #' @export
#' #'
#' probes_vcf_extract_hgvs <- function(vcf, transcript_type = c("ensembl", 'refseq'), exclude_noncoding_rna = TRUE, minimal_columns=TRUE){
#'   transcript_type <- rlang::arg_match(transcript_type)
#'
#'   vcf <- vcfR::read.vcfR(vcf, verbose = FALSE)
#'
#'   browser()
#'   df <- data.frame(
#'     VarID = vcfR::extract.info(vcf, element = "VCF"),
#'     HGVSg = vcfR::extract.info(vcf, element = "HGVSg"),
#'     HGVSc = vcfR::extract.info(vcf, element = "HGVSc"),
#'     SPDI = vcfR::extract.info(vcf, element = "SPDI")
#'   )
#'
#'   df <- tidyr::separate_rows(df, HGVSc, sep = ",")
#'
#'   if(transcript_type == "ensembl"){
#'     df <- dplyr::filter(df, startsWith(HGVSc, prefix = 'ENST'))
#'   }
#'   else if(transcript_type == "refseq"){
#'     df <- dplyr::filter(df, startsWith(HGVSc, prefix = 'NM'))
#'   }
#'
#'   if(exclude_noncoding_rna){
#'     df <- df[!grepl(x = df[['HGVSc']], pattern = ":n\\."), ]
#'   }
#'
#'   if(minimal_columns){
#'     df <- df[c('VarID', 'HGVSg', 'HGVSc')]
#'   }
#'   return(df)
#' }



