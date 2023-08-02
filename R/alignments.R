#
# probes_annotate_wt_vs_mut_plots_msa <- function(df_mut){
#   df_mut[["msa"]] <- purrr::map(seq_len(nrow(df_mut)), function(i) {generate_pairwise_msa_alignments_scalar(wt_seq = df_mut[['probe_wt_seq']][i], mut_seq = df_mut[['probe_mut_seq']][i]) })
#   return(df_mut)
# }
# generate_pairwise_msa_alignments_scalar <- function(wt_seq, mut_seq){
#   assertions::assert_string(wt_seq)
#   assertions::assert_string(mut_seq)
#   gg <- ggmsa::ggmsa(msa = c(Wildtype= wt_seq, Mutant=mut_seq), color = "Chemistry_NT", by_conservation = TRUE)
# }
#
# generate_pairwise_alignments_scalar <- function(){
#   alignment <- Biostrings::pairwiseAlignment(pattern = "ACTGCA", subject = "ACTTGCA")
#   pattern = Biostrings::alignedPattern(alignment)
#   subject = Biostrings::alignedSubject(alignment)
#   seq <- c(pattern, subject)
#   DECIPHER::BrowseSeqs(seq)
# }


gg_visualise_alignment <- function(seq1, seq2, seq1names = "wt", seq2names = "mutant", titles = NULL, interactive = FALSE){
  gglist <- lapply(seq_along(seq1), FUN = function(i){
    if(is.null(titles))
      t = NULL
    else
      t = titles[i]

    gg = gg_visualise_alignment_scalar(seq1 = seq1[i], seq2 = seq2[i], seq1name = seq1names[min(i, length(seq1names))], seq2name = seq2names[min(i, length(seq2names))], title = t, interactive = interactive)
    })
}


probe_generate_alignment_html <- function(gglist, prefix, median_size = 51, return_interactive_patch = FALSE){

  ggpatch <- patchwork::wrap_plots(gglist, ncol = 1)

  gginterpatch <- ggiraph::ggiraph(ggobj = ggpatch, height_svg = 1.5 * length(gglist), width_svg = median_size * 0.15)

  if(return_interactive_patch){
    return(gginterpatch)
  }

  htmlwidgets::saveWidget(widget = gginterpatch, file = paste0(prefix, ".alignments.html"), selfcontained = TRUE)
}


calculate_position <- function(sequence) {

  assertions::assert_character(sequence)
  position <- seq_along(sequence)  # Initialize position vector

  for (i in seq_along(sequence)) {
    if (sequence[i] == '-') {
      position[i:length(position)] <- position[i:length(position)] - 1
      position[i] <- NA  # Set position to NA for '-' elements
    }
  }

  return(position)
}

gg_visualise_alignment_scalar <- function(seq1, seq2, seq1name="subject", seq2name = "pattern", interactive = TRUE, return_df = FALSE, title = NULL, width_svg = NULL, height_svg = NULL){

  alignment <- Biostrings::pairwiseAlignment(subject=seq1, pattern = seq2)
  pattern <- Biostrings::alignedPattern(alignment)
  subject <- Biostrings::alignedSubject(alignment)

  df <- data.frame(
    subject = unlist(strsplit(as.character(pattern), split = "")),
    pattern = unlist(strsplit(as.character(subject), split = ""))
  )

  df[["position"]] <- seq_len(nrow(df))
  df[['difference']] <- paste0(df[['pattern']], ' > ', df[['subject']])
  df[["data_id"]] <- df[["position"]]

  df['different'] <- df[['subject']] != df[['pattern']]

  #browser()

  # Pivot Longer
  df <- df |>
    tidyr::pivot_longer(c(pattern, subject), names_to = "seqname", values_to = "sequence") |>
    dplyr::mutate(
      seqname = dplyr::case_when(
        seqname == "pattern" ~ seq1name,
        seqname == "subject" ~ seq2name,
        .default = "ERROR"
      )
    )

  # Position in Fasta (Excludes - characters)
  df <- df |>
    dplyr::mutate(position_in_fasta = calculate_position(sequence), .by = seqname)


  df[["tooltip"]] <- paste0(
    "Position (Alignment): ", df[['position']], "<br>",
    "Position (Fasta): ", df[['position_in_fasta']], "<br>",
    "Base: ", df[["sequence"]], "<br>",
    ifelse(df[["different"]], yes = paste0('Difference: ', df[['difference']], "<br>"), no = "")
  )


  if(return_df)
    return(df)

  if(is.null(width_svg)){
    width_svg = 0.15*nrow(df)
  }
  if(is.null(height_svg)){
    if(is.null(title)) height_svg <- 1.4
    else height_svg = 2
  }

  sequence_colours <- c(
    'T'="#66C2A5",
    '-' = "#FC8D62",
    'C'="#8DA0CB",
    'A'="#E78AC3",
    'U'="#A6D854",
    'G'="#FFD92F"
  )


  gg <- ggplot2::ggplot(df, ggplot2::aes(x=position,alpha = different, y = seqname, label = sequence, fill = sequence, tooltip = tooltip, data_id = data_id)) +
    ggiraph::geom_tile_interactive(show.legend = FALSE, ggplot2::aes(alpha = different, color = different, linetype = different), linewidth = 1) +
    ggplot2::geom_text(fontface = "bold", position = "identity", show.legend = FALSE, alpha= 0.8) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text = ggplot2::element_text(size = 11), axis.title = ggplot2::element_text(size = 11, face = "bold"), plot.title = ggplot2::element_text(hjust = 0.5,face = "bold", size = 13)) +
    ggplot2::scale_color_manual(values = c("TRUE"="red", "FALSE" = "white")) +
    ggplot2::scale_fill_manual(values = sequence_colours) +
    ggplot2::scale_alpha_manual(values = c("TRUE"=1, "FALSE" = 0.3)) +
    ggplot2::scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "blank")) +
    ggplot2::ggtitle(title) +
    ggplot2::xlab("Position") +
    ggplot2::ylab(NULL)

  if(interactive){
    ggiraph = ggiraph::ggiraph(ggobj=gg, width_svg = width_svg, height_svg = height_svg)
    return(ggiraph)
  }

  return(gg)
}
