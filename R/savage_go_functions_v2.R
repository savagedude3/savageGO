# Declare dependent packages in console, but not in this script itself
# These should be written to the Imports section of the DESCRIPTION file for the package
# use_package("clusterProfiler")
# use_package("dplyr")
# use_package("ggplot2")
# use_package("stringr")
# use_package("tidyr")
# use_package("ggplot2")
# use_package("AnnotationDbi")
# use_package("biomaRt")
# use_package("org.Mm.eg.db")
# use_package("enrichplot")
# use_package("RColorBrewer")
# use_package("rrvgo")
# use_package("stats")
# use_package("GOSemSim")

get_entrez <- function(gene_list){

  gene_df <- data.frame(gene_list)

  colnames(gene_df) <- c("genes")

  gene_df$entrezgene_id <- AnnotationDbi::mapIds(org.Mm.eg.db,
                                                 keys = gene_df$genes,
                                                 column = "ENTREZID",
                                                 keytype = "SYMBOL",
                                                 multiVals = "first")

  # remove rows with NAs in the column `entrezgene_id`
  gene_df <- gene_df %>%
    tidyr::drop_na(entrezgene_id)

  # drop duplicates
  gene_df <- gene_df %>%
    dplyr::distinct(entrezgene_id, .keep_all = TRUE)

  return(gene_df)
}

#TODO: get human entrez ids after entering mouse ids or gene names
get_human_id <- function(gene_list){

  gene_df <- data.frame(gene_list)

  colnames(gene_df) <- c("genes")

  ensembl <- biomaRt::useEnsembl(biomart = "genes",
                                 dataset = "mmusculus_gene_ensembl")

  mouse_to_human <- biomaRt::getBM(attributes = c("external_gene_name", "hsapiens_homolog_associated_gene_name"),
                                   values = gene_df$genes,
                                   mart = ensembl)

  # add this new info to gene_df
  gene_df <- gene_df %>%
    dplyr::left_join(mouse_to_human, by = c("genes" = "external_gene_name"))

  # remove rows with NAs in the column `hsapiens_homolog_associated_gene_name`
  gene_df <- gene_df %>%
    tidyr::drop_na(hsapiens_homolog_associated_gene_name)

  # remove rows with "" in the column `hsapiens_homolog_associated_gene_name`
  gene_df <- gene_df %>%
    dplyr::filter(hsapiens_homolog_associated_gene_name != "")

  # drop duplicates
  gene_df <- gene_df %>%
    dplyr::distinct(hsapiens_homolog_associated_gene_name, .keep_all = TRUE)

  #convert human ENSEMBL IDs to human Entrez IDs
  hsa_ensembl <- biomaRt::useEnsembl(biomart = "genes",
                                     dataset = "hsapiens_gene_ensembl")

  hsa_name_to_entrez <- biomaRt::getBM(attributes = c("external_gene_name", "entrezgene_id"),
                                       values = gene_df$hsapiens_homolog_associated_gene_name,
                                       mart = hsa_ensembl)

  # add this new info to gene_df
  gene_df <- gene_df %>%
    dplyr::left_join(hsa_name_to_entrez, by = c("hsapiens_homolog_associated_gene_name" = "external_gene_name"))

  # remove rows with NAs in the column `entrezgene_id`
  gene_df <- gene_df %>%
    tidyr::drop_na(entrezgene_id)

  # remove rows with "" in the column `entrezgene_id`
  gene_df <- gene_df %>%
    dplyr::filter(entrezgene_id != "")

  # drop duplicates
  gene_df <- gene_df %>%
    dplyr::distinct(entrezgene_id, .keep_all = TRUE)

  return(gene_df)
}


go_gene_list <- function(gene_list, gene_df, go_type){

  gene_df_subset <- dplyr::filter(gene_df, genes %in% gene_list)

  go_ora <- clusterProfiler::enrichGO(gene = as.character(gene_df_subset$entrezgene_id),
                                      OrgDb = org.Mm.eg.db,
                                      universe = as.character(gene_df$entrezgene_id),
                                      ont = go_type,
                                      readable = TRUE) # maps gene IDs to gene names
  #head(go_ora)
  return(go_ora)
}

#removes redundancy using the rrvgo package.
go_reduce <- function(go_ora, go_type){

  sem_data <- GOSemSim::godata(annoDb = "org.Mm.eg.db", ont= go_type)

  #using GOSemSim::mgoSim to make similarity matrix since the rrvgo::calculateSimMatrix function was dropping some terms
  simMatrix <- outer(go_ora$ID, go_ora$ID,
                     Vectorize(function(x, y) GOSemSim::mgoSim(x, y, semData=sem_data, measure="Wang")))

  rownames(simMatrix) <- go_ora$ID
  colnames(simMatrix) <- go_ora$ID

  #use the gene ratio from the data in that GO term to determine which should be kept
  # convert gene ratios from a ratio to a decimal
  fraction_list <- go_ora$GeneRatio
  decimal_values <- unname(sapply(fraction_list, function(x) eval(parse(text = x))))
  scores <- stats::setNames(decimal_values, go_ora$ID)

  reducedTerms <- rrvgo::reduceSimMatrix(simMatrix,
                                         scores = scores,
                                         threshold=0.7,
                                         orgdb="org.Mm.eg.db")
  return(reducedTerms)
}

go_reduced_plot <- function(go_ora, reducedTerms, go_type, plot_type = "dotplot", numTerms = 5){

  if(plot_type == "treemap"){
    rrvgo::treemapPlot(reducedTerms)
  }

  if(plot_type == "dotplot"){

    #first get only the parent terms
    sorted_terms <- dplyr::arrange(reducedTerms, desc(score))
    unique_terms <- dplyr::distinct(sorted_terms, parentTerm, .keep_all = TRUE)
    unique_terms$parentTerm <- factor(unique_terms$parentTerm, levels = rev(unique_terms$parentTerm))

    #get adjusted p value from the go_ora object
    go_current.df <- data.frame(go_ora)
    rownames(go_current.df) <- NULL

    go_ora_reduced <- dplyr::filter(go_current.df, ID %in% unique_terms$go)

    #calculate ratio and sort terms
    go_ora_reduced$ratio_num <- sapply(go_ora_reduced$GeneRatio, function(x) eval(parse(text = x)))
    go_ora_reduced <- dplyr::arrange(go_ora_reduced, desc(ratio_num))
    go_ora_reduced$Description <- factor(go_ora_reduced$Description, levels = rev(go_ora_reduced$Description))

    #max numTerms GO terms to plot
    if(length(go_ora_reduced$ID) > numTerms){
      go_ora_reduced <- go_ora_reduced[1:numTerms,]
    }


    plot_out <- ggplot2::ggplot(data = go_ora_reduced) + ggplot2::geom_point(aes(x = ratio_num,
                                                                                 y = Description,
                                                                                 fill = p.adjust,
                                                                                 size = Count),
                                                                             shape = 21) +
      ggplot2::scale_size(limits = c(0, NA), range = c(0,8)) +
      ggplot2::theme_bw() + ggplot2::scale_fill_gradientn(colors = c("red", "purple", "blue")) +
      ggplot2::guides(fill = guide_colorbar(reverse = TRUE)) +
      ggplot2::scale_y_discrete(labels = function(y) str_wrap(y, width = 30)) +
      ggplot2::xlab("Gene Ratio") + ggplot2::ylab("GO Term") +
      ggplot2::theme(plot.title = element_text(size = 20),
                     axis.title = element_text(size = 18),
                     axis.text.x = element_text(size = 10),
                     axis.text.y = element_text(size = 18))

    return(plot_out)

  }
}
