
#' Parse Vep annotated text file
#'
#' Read Mutations to build probes from in vep text File
#'
#' @param txt path to vep annotated text files
#' @param transcript_type Whether to use ensembl or refseq transcripts
#' @param exclude_noncoding_rna Should non-coding RNA mutations be excluded
#' @param minimal_columns Return only the minimal set of columns required
#'
#' @return data.frame describing each vep mutation
#' @export
probes_read_vep_txt <- function(txt, transcript_type = c("ensembl", 'refseq'), exclude_noncoding_rna = TRUE, minimal_columns=TRUE){
  transcript_type <- rlang::arg_match(transcript_type)

  df_vep <- read.csv(txt, header = TRUE, sep ="\t")
  df_vep <- dplyr::as_tibble(df_vep)

  #browser()
  df <- df_vep |>
    dplyr::select(VarID = `X.Uploaded_variation`, Gene, Location, Allele, HGVSc, Transcript = Feature, Feature_type, cDNA_position, CDS_position)

  # Fix VarID if '.'
  df[['VarID']] <- ifelse(df[['VarID']] == ".", yes = paste0(df[['Location']],':',df[['Allele']]), no = df[['VarID']])

  # Drop non-transcript features
  df <- df |>
    dplyr::filter(Feature_type == "Transcript")

  # Select only one transcript database (ensembl/gencode vs refseq)
  if(transcript_type == "ensembl"){
    df <- dplyr::filter(df, startsWith(HGVSc, prefix = 'ENST'))
  }
  else if(transcript_type == "refseq"){
    df <- dplyr::filter(df, startsWith(HGVSc, prefix = 'NM'))
  }

  # Exclude non-coding RNA
  if(exclude_noncoding_rna){
    df <- df[!grepl(x = df[['HGVSc']], pattern = ":n\\."), ]
  }

  # Exclude anything Noncoding
  df <- df[df[['cDNA_position']] != "-", ]

  # Fix cDNA position (split into multiple)
  df <- df |>
    dplyr::mutate(
      cDNA_mut_start = sub(x=cDNA_position, "([0-9]+)-?.*", "\\1") |> as.numeric(),
      cDNA_mut_end = sub(x=cDNA_position, "[0-9]+-([0-9]+)$", "\\1") |> as.numeric(),
      cDNA_mut_length = cDNA_mut_end - cDNA_mut_start + 1
      )

  if(minimal_columns){
    df <- df[c('VarID', 'Gene', 'HGVSc', 'Transcript', 'cDNA_position', 'cDNA_mut_start', 'cDNA_mut_end', 'cDNA_mut_length')]
  }
  return(df)
}


#' Construct Probes
#'
#' @param mutations (data.frame)
#' @param probe_type which type of probe sequence to return (cDNA, mRNA, mRNA_no_U)
#' @param ensembl biomart from which to fetch transcript (cDNA)  sequences. Define using [load_biomart()]
#'
#' @return data.frame describing 1 probe per mutation-transcript isoform combination.
#' Try piping into probes_collapse_duplicates to cleanup
#' @export
#'
#' @details
#' There are many ways to build probes. Our approach is as follows.
#'
#' 1. Define a fixed probe_size
#' 2. Dynamically calculate the number of bases upstream and downstream based on total probe size.
#' Done separately for WT and Mutant since mutant seq may have insertions
#'
#' Note where mutation and probe size does not allow for identical sizes of upstream and downstream flanking region,
#' we always place the extra base on the upstream size
#'
#'Key formulas used:
#'
#'Ts = Total probe size = defined by probe_size
#'
#'L = mutation length (SNV = 1)
#'
#'Uwt = Upstream bases (wild type) = ceil((Ts - 1)/2)
#'
#'Dwt = Downstream bases (wild type) = floor((Ts - 1)/2)
#'
#'Umut = Upstream bases (mutant probe) = ceil((Ts-L)/2)
#'
#'Dmut = Downstream bases (mutant probe) = floor((Ts-L)/2)
#'
probes_construct <- function(df_mut, ensembl, probe_size = 51, probe_type = c('cDNA', 'mRNA','mRNA_no_U')){

  # Assertions
  assertions::assert_dataframe(df_mut)
  assertions::assert_names_include(df_mut, c('HGVSc', 'VarID', 'cDNA_position', 'Transcript', 'cDNA_mut_start', 'cDNA_mut_end', 'cDNA_mut_length'))
  assertions::assert_no_duplicates(df_mut[['HGVSc']])

  #assertions::assert_character_vector(mutations)
  assertions::assert_number(probe_size)
  #assertions::assert_number(bp_downstream)
  assertions::assert_subset(probe_type, c('mRNA','mRNA_no_U', 'cDNA'))

  probe_type <- rlang::arg_match(probe_type)
  assertions::assert(probe_type == 'cDNA', msg = 'only cDNA probe types are currently supported')


  # Convert mutations to dataframe
  df_mut <- df_mut |>
   dplyr::mutate(hgvs_to_dataframe(HGVSc))

  #Confirm mutation types are sensible (no indels
  # TODO: implement


  # Dropping mutations which are too late in the transcript
  #(requires transcript length)

  # Retrieve full cdna sequence sequence (WT)
  #df_mut[["transcript_noversion"]] <- sub(x= df_mut[["ranscript"]], pattern = "\\..*$", replacement = "")
  df_mut[["cDNA"]] <- query_transcript_sequence_api(transcript_id = df_mut[["Transcript"]], ensembl = ensembl)


  # Mutate Sequence
  df_mut <- df_mut |>
    dplyr::mutate(cDNA_mutant = mutate_seq(sequence = cDNA, start_position = cDNA_mut_start, end_position = cDNA_mut_end, ref = ref, alt = alt, mut_type = mut_type))


  # Calculate how many upstream and downstream bases there should be based on total probe
  # size and mutation length
  df_mut['bp_upstream_wt'] = ceiling((probe_size - 1)/2)
  df_mut['bp_downstream_wt'] = floor((probe_size - 1)/2)
  df_mut['bp_upstream_mut'] = ceiling((probe_size - df_mut[['cDNA_mut_length']])/2)
  df_mut['bp_downstream_mut'] = floor((probe_size - df_mut[['cDNA_mut_length']])/2)

  # Add start-end coords for probe
  df_mut[['probe_wt_start']]  <- df_mut[['cDNA_mut_start']] - df_mut[['bp_upstream_wt']]
  df_mut[['probe_wt_end']]  <- df_mut[['cDNA_mut_start']] + df_mut[['bp_downstream_wt']]
  df_mut[["probe_wt_length_theoretical"]] <- df_mut[['probe_wt_end']] - df_mut[['probe_wt_start']] + 1

  # Add start-end coords for mutant (dels might need recentering)
  df_mut[['probe_mut_start']]  <- df_mut[['cDNA_mut_start']] - df_mut[['bp_upstream_mut']]
  df_mut[['probe_mut_end']]  <- df_mut[['cDNA_mut_end']] + df_mut[['bp_downstream_mut']]
  df_mut[["probe_mut_length_theoretical"]] <- df_mut[['probe_mut_end']] - df_mut[['probe_mut_start']] + 1

  # Subset cdna sequence to construct probe sequences
  df_mut[["probe_wt_seq"]] <- substr(df_mut[["cDNA"]], start = df_mut[["probe_wt_start"]], df_mut[["probe_wt_end"]])
  df_mut[["probe_mut_seq"]] <- substr(df_mut[["cDNA_mutant"]], start = df_mut[["probe_mut_start"]], df_mut[["probe_mut_end"]])

  # Calculate Probe Lengths
  df_mut[["probe_wt_length"]] <- nchar(df_mut[['probe_wt_seq']])
  df_mut[["probe_mut_length"]] <- nchar(df_mut[['probe_mut_seq']])

  df_mut[['probe_wt_expected_size']] <- df_mut[["probe_wt_length"]] == probe_size
  df_mut[['probe_mut_expected_size']] <- df_mut[["probe_mut_length"]] == probe_size

  if(!all(df_mut[['probe_wt_expected_size']])){
    df_unexpected_size <- df_mut[!df_mut[['probe_wt_expected_size']], c('VarID', 'Transcript', 'cDNA_position', 'ref', 'alt', 'probe_wt_seq')]
    cli::cli_rule('Warning: Bad probe size (WT)')
    cli::cli_alert_warning('The following mutations have an unexpected wt probe size. This usually means variant is too close to the start/end of the transcript.')
    print(df_unexpected_size)
    cli::cli_rule()
  }

  if(!all(df_mut[['probe_mut_expected_size']])){
    df_unexpected_size <- df_mut[!df_mut[['probe_mut_expected_size']], c('VarID', 'Transcript', 'cDNA_position', 'ref', 'alt', 'probe_wt_seq')]
    cli::cli_rule('Warning: Bad probe size (Mutant)')
    cli::cli_alert_warning('The following mutations have an unexpected mutant probe size. This usually means variant is too close to the start/end of the transcript.')
    print(df_unexpected_size)
    cli::cli_rule()
  }



  return(df_mut)

  # assertions::assert(all(df_mut[['wt_probe_expected_size']]), msg = paste0("There was a wt probe of unexpected size: {df_mut[!df_mut[['wt_probe_expected_size']]]}"))
  # assertions::assert(all(df_mut[['mut_probe_expected_size']]))


  # Compute Failure / Warning modes (too early / late in isoform)

  # Fail 1: start < 1
  # Fail 2: end > length of isoform
  # Fail 3: low sequence complexity

  # Convert to tibble
  df_mut <- dplyr::tibble(df_mut)

  return(df_mut)
}



hgvs_to_dataframe <- function(mutations, must_be_cdna = TRUE, exclude_raw = TRUE, minimal = TRUE){

  # Build Table
  df_mutations <- data.frame(raw = mutations)
  df_mutations[['transcript']] = sub(x = mutations, pattern = "(^.*?):.*$", "\\1")
  df_mutations[['position']] = as.numeric(sub(x = mutations, pattern = "^.*?:[cgpnr]\\.\\*?([-0-9]+).*$", "\\1"))
  df_mutations[['type']] = sub(x = mutations, pattern = "^.*?:([cgpnr])\\.\\*?.*$", "\\1")
  df_mutations[['ref']] = sub(x = mutations, pattern = "^.*([ACGTacgt]+)>.*$", "\\1")
  df_mutations[['alt']] = sub(x = mutations, pattern = "^.*>(.*)$", "\\1")


  # Deal with Duplicates
  df_mutations[['ref']] <- ifelse(df_mutations[['ref']] == mutations & grepl(x=mutations, pattern = "[0-9]+dup"), yes = "original", df_mutations[['ref']])
  df_mutations[['alt']] <- ifelse(df_mutations[['alt']] == mutations & grepl(x=mutations, pattern = "[0-9]+dup"), yes = "dup", df_mutations[['alt']])

  df_mutations[['mut_type']] <- dplyr::case_when(
    grepl(x=mutations, pattern = "[0-9]+dup") ~ "DUP",
    nchar(df_mutations[['ref']]) == nchar(df_mutations[['alt']]) ~ "SNV",
    .default = "unknown"
    )


  # Replace custom parsing with hgvsParseR::parseHGVS(strings = 'c.260_264+48dup')

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
  if(must_be_cdna){
  assertions::assert(all(df_mutations[['type']] %in% c("c")), msg = "Non cdna based mutations detected. Please convert the following mutations to RNA space then retry: {df_mutations[['raw']][df_mutations[['type']] != 'c']}")
  }

  if(exclude_raw){
    df_mutations <- df_mutations[names(df_mutations) != 'raw']
  }

  if(minimal){
     df_mutations <- df_mutations[,c('type', 'ref', 'alt', 'mut_type')]
  }
  return(df_mutations)
}

mutate_seq <- function(mut_type,sequence, start_position, end_position, ref, alt){
  vapply(seq_along(mut_type), FUN = function(i){
    if(mut_type[i] == "SNV")
      mutate_seq_snv(sequence[i], start_position[i], ref[i], alt[i])
    else if(mut_type[i] == "DUP")
      mutate_seq_dup(sequence[i], start_position[i])
    else
      return('not_run')
      },FUN.VALUE = character(1))
}

mutate_seq_snv <- function(sequence, position, ref, alt){
  vapply(X = seq_along(sequence), FUN = function(i){ mutate_seq_snv_scalar(sequence[i], position[i], ref[i], alt[i]) }, FUN.VALUE = character(1))
}

mutate_seq_snv_scalar <- function(sequence, position, ref, alt){
# sequence = Biostrings::DNAString(sequence)

 # Assert Ref Match
 ref_base_at_index <- substr(sequence, position, position)
 assertions::assert(ref_base_at_index == ref, msg = "Base {position} in reference sequence is a '{ref_base_at_index}', not '{ref}': {sequence} ")

 mutated <- sequence
 substr(mutated, start = position, stop = position) <- alt
 #mutated <- Biostrings::replaceLetterAt(x = sequence,at = position,letter = alt)
 return(mutated)
}

mutate_seq_dup <- function(sequence, position){
  vapply(X = seq_along(sequence), FUN = function(i){ mutate_seq_dup_scalar(sequence[i], position[i])}, FUN.VALUE = character(1))

}
mutate_seq_dup_scalar <- function(sequence, position){
  ref = substr(sequence, start = position, stop = position)
  firstpart = substr(sequence, start = 1, stop = position)
  lastpart = substr(sequence, start = position + 1, stop = nchar(sequence))
  mutated <- paste0(firstpart, ref,lastpart)
  #sequence = Biostrings::DNAString(sequence)
  #mutated <- Biostrings::replaceLetterAt(x = sequence,at = position:(position+1),letter = rep(sequence[position], times = 2))
  return(mutated)
}


#' Collapse duplicates
#'
#'
#' @param df_mut output from [probes_construct()]
#'
#' @return Duplicate probes collapsed into single probes (data.frame)
#' @export
#'
probes_collapse_duplicates <- function(df_mut){
  # Remove Redundant (identical) probes
  df_mut <- df_mut |>
    dplyr::group_by(probe_mut_seq) |>
    dplyr::summarise(
      #VarId = paste0(unique(VarID), collapse = "; "),
      VarID = paste0(unique(VarID), collapse = "; "),
      nVariants = dplyr::n_distinct(VarID),
      Transcripts = paste0(unique(Transcript), collapse = "; "),
      dplyr::across(-Transcript, .fns = function(vec){vec[1]})
    )
    #dplyr::mutate(ID = paste0(VarID, '|', Transcripts))

  if(any(df_mut[["nVariants"]] > 1)){
    cli::cli_rule('Warning: odd')
    cli::cli_alert_warning('Multiple different original variant IDs produce identical probe sequences')
    df_redundant_variants <- df_mut[df_mut[['nVariants']] > 1, c('VarID', 'Transcript', 'cDNA_position', 'ref', 'alt', 'probe_wt_seq')]
  }

  return(df_mut)
}

probes_longform_with_ids <- function(df_mut){
  df_mut |> tidyr::pivot_longer(c(probe_wt_seq, probe_mut_seq), names_to = "probe_class", values_to = "probe_sequence") |>
    dplyr::mutate(
      probe_class = dplyr::case_when(
        probe_class == "probe_wt_seq" ~ "WT",
        probe_class == "probe_mut_seq" ~ "MUTANT",
        .default = "PROBLEM!!"),
      ID = paste(VarID, probe_class, Transcripts, sep = " | ")
    )
}

probes_write_output <- function(df_mut, outdir, prefix = "probes"){
  df_mut <- probes_longform_with_ids(df_mut)

  # Create Dir
  if(!dir.exists(outdir)) {
    dir.create(outdir)
  }

  # Create Fasta File
  path_fasta = paste0(outdir, "/", prefix, '.fasta')
  cli::cli_alert_info('Writing to {path_fasta}')
  file.create(path_fasta)

  sequences <- df_mut[['probe_sequence']]
  names(sequences) <- df_mut[['ID']]

  purrr::walk(seq_along(sequences), .f = function(i){
    seq = sequences[i]
    seqname = names(sequences)[i]

    write(paste0('> ', seqname), file = path_fasta, append = TRUE)
    write(seq, file = path_fasta, append = TRUE)
  })


  # Create Output Table
  path_info = paste0(outdir, "/", prefix, '.csv')
  cli::cli_alert_info('Writing to {path_info}')

  write.csv(df_mut, file = path_info, row.names = FALSE)



}
