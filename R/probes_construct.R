probes_construct <- function(mutations, bp_upstream = 20, bp_downstream = 20, probe_type = c('mRNA','mRNA_no_U', 'cDNA')){

  # Assertions
  assertions::assert_no_duplicates(mutations)
  assertions::assert_character_vector(mutations)
  assertions::assert_number(bp_upstream)
  assertions::assert_number(bp_downstream)
  assertions::assert_subset(probe_type, c('mRNA','mRNA_no_U', 'cDNA'))

  probe_type <- rlang::arg_match(probe_type)
  assertions::assert(probe_type != 'cDNA', msg = 'cDNA probe types are not yet supported')

  # Convert mutations to dataframe
  df_mut <- hgvs_to_dataframe(mutations)

  #Confirm mutation types are sensible (no indels


  # Add start-end coords for probe
  df_mut[['p_start']]  <- df_mut[['position']] - bp_upstream
  df_mut[['p_end']]  <- df_mut[['position']] + bp_downstream


  # Retrieve wild type sequence

  # Compute Failure / Warning modes (too early / late in isoform)

  # Fail 1: start < 1
  # Fail 2: end > length of isoform
  # Fail 3: low sequence complexity

  return(df_mut)
}



hgvs_to_dataframe <- function(mutations, must_be_rna = TRUE){

  # Build Table
  df_mutations <- data.frame(raw = mutations)
  df_mutations[['transcript']] = sub(x = mutations, pattern = "(^.*?):.*$", "\\1")
  df_mutations[['position']] = as.numeric(sub(x = mutations, pattern = "^.*?:c\\.([0-9]+).*$", "\\1"))
  df_mutations[['type']] = sub(x = mutations, pattern = "^.*?:([cgpn])\\..*$", "\\1")
  df_mutations[['ref']] = sub(x = mutations, pattern = "^.*([ACGT]+)>.*$", "\\1")
  df_mutations[['alt']] = sub(x = mutations, pattern = "^.*>(.*)$", "\\1")

  # Test HGVS info extraction
  indices_where_extraction_failed <- unique(c(
    which(df_mutations[['transcript']] == df_mutations[['raw']]),
    which(df_mutations[['position']] == df_mutations[['raw']]),
    which(df_mutations[['type']] == df_mutations[['raw']]),
    which(df_mutations[['ref']] == df_mutations[['raw']]),
    which(df_mutations[['alt']] == df_mutations[['raw']])
  ))

  assertions::assert(
    length(indices_where_extraction_failed) == 0,
    msg = 'Improperly formatted HGVS mutations. Please check the following mutations: {df_mutations[["raw"]][indices_where_extraction_failed]}'
  )

  # Confirm mutation is in mRNA notation
  if(must_be_rna){
  assertions::assert(all(df_mutations[['type']] == "c"), msg = "Non mRNA based mutations detected. Please convert the following mutations to RNA space then retry: {df_mutations[['raw']][df_mutations[['type']] != 'c']}")
  }

  return(df_mutations)
}
