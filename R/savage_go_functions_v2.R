#' Get Entrez gene IDs from mouse IDs
#'
#' @param gene_list A character vector of mouse gene symbols
#' @return A data frame with the original gene names and the associated Entrez IDs
#' @importFrom org.Mm.eg.db org.Mm.eg.db
#' @importFrom AnnotationDbi mapIds
#' @importFrom tidyr drop_na
#' @importFrom dplyr distinct
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @export
get_entrez <- function(gene_list) {
  gene_df <- data.frame(genes = gene_list)

  gene_df$entrezgene_id <- AnnotationDbi::mapIds(
    org.Mm.eg.db::org.Mm.eg.db,
    keys = gene_df$genes,
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  )

  gene_df <- gene_df %>%
    tidyr::drop_na(.data$entrezgene_id) %>%
    dplyr::distinct(.data$entrezgene_id, .keep_all = TRUE)

  return(gene_df)
}

#' Get human gene names and Entrez IDs from mouse gene names
#'
#' @param gene_list A character vector of mouse gene symbols
#' @return A data frame with the original mouse gene names and the associated human gene names and Entrez IDs
#' @importFrom biomaRt useEnsembl getBM
#' @importFrom dplyr left_join filter distinct
#' @importFrom tidyr drop_na
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @export
get_human_id <- function(gene_list) {
  gene_df <- data.frame(genes = gene_list)

  ensembl <- biomaRt::useEnsembl(biomart = "genes",
                                 dataset = "mmusculus_gene_ensembl")

  mouse_to_human <- biomaRt::getBM(
    attributes = c("external_gene_name", "hsapiens_homolog_associated_gene_name"),
    values = gene_df$genes,
    mart = ensembl
  )

  gene_df <- gene_df %>%
    dplyr::left_join(mouse_to_human, by = c("genes" = "external_gene_name")) %>%
    tidyr::drop_na(.data$hsapiens_homolog_associated_gene_name) %>%
    dplyr::filter(.data$hsapiens_homolog_associated_gene_name != "") %>%
    dplyr::distinct(.data$hsapiens_homolog_associated_gene_name, .keep_all = TRUE)

  hsa_ensembl <- biomaRt::useEnsembl(biomart = "genes",
                                     dataset = "hsapiens_gene_ensembl")

  hsa_name_to_entrez <- biomaRt::getBM(
    attributes = c("external_gene_name", "entrezgene_id"),
    values = gene_df$hsapiens_homolog_associated_gene_name,
    mart = hsa_ensembl
  )

  gene_df <- gene_df %>%
    dplyr::left_join(hsa_name_to_entrez, by = c("hsapiens_homolog_associated_gene_name" = "external_gene_name")) %>%
    tidyr::drop_na(.data$entrezgene_id) %>%
    dplyr::filter(.data$entrezgene_id != "") %>%
    dplyr::distinct(.data$entrezgene_id, .keep_all = TRUE)

  return(gene_df)
}

#' Run gene ontology on a list of mouse gene names
#'
#' @param gene_list A list of mouse gene names
#' @param gene_df A data frame with your mouse gene names and Entrez IDs
#' @param go_type The Gene Ontology type to run ("BP", "CC", or "MF")
#' @return An enrichResult object
#' @importFrom clusterProfiler enrichGO
#' @importFrom org.Mm.eg.db org.Mm.eg.db
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @export
go_gene_list <- function(gene_list, gene_df, go_type) {
  gene_df_subset <- dplyr::filter(gene_df, .data$genes %in% gene_list)

  go_ora <- clusterProfiler::enrichGO(
    gene = as.character(gene_df_subset$entrezgene_id),
    OrgDb = org.Mm.eg.db::org.Mm.eg.db,
    universe = as.character(gene_df$entrezgene_id),
    ont = go_type,
    readable = TRUE
  )
  return(go_ora)
}

#' Removes redundancy in GO terms using rrvgo
#'
#' @param go_ora An enrichResult object
#' @param go_type The Gene Ontology type ("BP", "CC", or "MF")
#' @return A data frame of reduced GO terms
#' @importFrom GOSemSim godata mgoSim
#' @importFrom rrvgo reduceSimMatrix
#' @importFrom stats setNames
#' @importFrom org.Mm.eg.db org.Mm.eg.db
#' @export
go_reduce <- function(go_ora, go_type) {
  # Note: GOSemSim often requires the actual OrgDb object
  sem_data <- GOSemSim::godata(annoDb = org.Mm.eg.db::org.Mm.eg.db, ont = go_type)

  simMatrix <- outer(go_ora$ID, go_ora$ID,
                     Vectorize(function(x, y) {
                       GOSemSim::mgoSim(x, y, semData = sem_data, measure = "Wang")
                     }))

  rownames(simMatrix) <- go_ora$ID
  colnames(simMatrix) <- go_ora$ID

  # Safer decimal conversion without eval(parse)
  decimal_values <- sapply(go_ora$GeneRatio, function(x) {
    parts <- as.numeric(unlist(strsplit(x, "/")))
    parts[1] / parts[2]
  })

  scores <- stats::setNames(decimal_values, go_ora$ID)

  reducedTerms <- rrvgo::reduceSimMatrix(
    simMatrix,
    scores = scores,
    threshold = 0.7,
    orgdb = org.Mm.eg.db::org.Mm.eg.db
  )
  return(reducedTerms)
}

#' Creates a plot of the reduced GO terms
#'
#' @param go_ora An enrichResult object
#' @param reducedTerms A data frame from go_reduce()
#' @param go_type The Gene Ontology type
#' @param plot_type Either "treemap" or "dotplot"
#' @param numTerms The number of GO terms to display
#' @return A ggplot2 object
#' @importFrom rrvgo treemapPlot
#' @importFrom dplyr arrange distinct filter
#' @importFrom ggplot2 ggplot geom_point aes scale_size theme_bw scale_fill_gradientn guides guide_colorbar scale_y_discrete xlab ylab theme element_text
#' @importFrom stringr str_wrap
#' @importFrom rlang .data
#' @export
go_reduced_plot <- function(go_ora, reducedTerms, go_type, plot_type = "dotplot", numTerms = 5) {

  if (plot_type == "treemap") {
    return(rrvgo::treemapPlot(reducedTerms))
  }

  if (plot_type == "dotplot") {
    sorted_terms <- dplyr::arrange(reducedTerms, desc(.data$score))
    unique_terms <- dplyr::distinct(sorted_terms, .data$parentTerm, .keep_all = TRUE)

    # Correcting factor levels for plotting
    unique_terms$parentTerm <- factor(unique_terms$parentTerm, levels = rev(unique_terms$parentTerm))

    go_current.df <- data.frame(go_ora)
    go_ora_reduced <- dplyr::filter(go_current.df, .data$ID %in% unique_terms$go)

    # Calculate ratio and sort
    go_ora_reduced$ratio_num <- sapply(go_ora_reduced$GeneRatio, function(x) {
      parts <- as.numeric(unlist(strsplit(x, "/")))
      parts[1] / parts[2]
    })

    go_ora_reduced <- dplyr::arrange(go_ora_reduced, desc(.data$ratio_num))

    if (nrow(go_ora_reduced) > numTerms) {
      go_ora_reduced <- go_ora_reduced[1:numTerms, ]
    }

    go_ora_reduced$Description <- factor(go_ora_reduced$Description, levels = rev(go_ora_reduced$Description))

    plot_out <- ggplot2::ggplot(data = go_ora_reduced) +
      ggplot2::geom_point(ggplot2::aes(x = .data$ratio_num,
                                       y = .data$Description,
                                       fill = .data$p.adjust,
                                       size = .data$Count),
                          shape = 21) +
      ggplot2::scale_size(limits = c(0, NA), range = c(0, 8)) +
      ggplot2::theme_bw() +
      ggplot2::scale_fill_gradientn(colors = c("red", "purple", "blue")) +
      ggplot2::guides(fill = ggplot2::guide_colorbar(reverse = TRUE)) +
      ggplot2::scale_y_discrete(labels = function(y) stringr::str_wrap(y, width = 30)) +
      ggplot2::xlab("Gene Ratio") +
      ggplot2::ylab("GO Term") +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 20),
                     axis.title = ggplot2::element_text(size = 18),
                     axis.text.x = ggplot2::element_text(size = 10),
                     axis.text.y = ggplot2::element_text(size = 18))

    return(plot_out)
  }
}
