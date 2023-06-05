
query_transcript_sequence <- function(transcript_id, database = c('GRCh38_EnsDb_v86', "hg19_ucsc"){
  assertions::assert_subset(database, c('GRCh38_EnsDb_v86', 'hg19_ucsc'))

  if(database == "hg19_ucsc"){
    rlang::check_required(TxDb.Hsapiens.UCSC.hg19.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    df_res <- GenomicFeatures::transcripts(txdb, filter(list(tx_id = c(transcript_id))), columns = c("tx_id","tx_name"))
    return(df_res)
  }
  else if(database == "GRCh38_EnsDb_v86"){
   edb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
   dna <- ensembldb:::getGenomeTwoBitFile(edb)
  }



}

query_transcript_length <- function(){

}

query_transcript_sequence <- function(transcript_id){

}
