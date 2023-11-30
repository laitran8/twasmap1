bartesting1 <- function(phenotypes, genebass_files, twas_files){
  library(dplyr)
  library(ggplot2)
  
  # Ensure the lists of phenotypes and files are of the same length
  if (length(phenotypes) != length(genebass_files) | length(phenotypes) != length(twas_files)) {
    stop("Number of phenotypes and files do not match!")
  }
  
  results <- list()
  
  for (i in 1:length(phenotypes)) {
    phenotype <- phenotypes[i]
    genebass_file <- genebass_files[i]
    twas_file <- twas_files[i]
    
    # Read Genebass data
    gene_data <- read.csv(file = genebass_file)
    
    # Immediately check the structure of the gene_data to ensure the column is present
    print(names(gene_data))
    if (!"P.Value.Burden" %in% names(gene_data)) {
      stop("The 'P.Value.Burden' column does not exist in the gene_data.")
    }
    
    # Remove NA values before correction
    gene_data <- gene_data[!is.na(gene_data$P.Value.Burden), ]
    gene_data <- gene_data[gene_data$P.Value.Burden < 6.7e-7,]
    
    # Apply Bonferroni correction to the burden p-values
    gene_data$Bonferroni.P.Burden <- p.adjust(gene_data$P.Value.Burden, method = "bonferroni")
    
    # Read and filter TWAS data
    twas_data <- data.table::fread(twas_file, data.table = FALSE)
    twas_gtex_lasso <- twas_data %>%
      filter(stringr::str_detect(PANEL, 'GTEx'), MODEL == 'lasso')
    twas_gtex_lasso$Bonferroni.P <- p.adjust(twas_gtex_lasso$TWAS.P, method = "bonferroni")
    
    # Identify G_TWAS
    G_TWAS <- unique(twas_gtex_lasso[twas_gtex_lasso$Bonferroni.P < 0.05, "ID"])
    
    # Define G_burden based on the corrected p-values
    G_burden <- gene_data$Gene.Name
    
    # Find intersection of G_TWAS and G_burden
    G_common <- intersect(G_TWAS, G_burden)
    
    # Store results
    results[[phenotype]] <- list(G_TWAS = G_TWAS, G_burden = G_burden, G_common = G_common)
    
    # Plot results
    gene_count_df <- data.frame(
      Trait = phenotype,
      Gburden = length(G_burden),
      GTWAS = length(G_TWAS),
      Common = length(G_common)
    )
    
    # Reshape the data for plotting
    gene_count_df_long <- tidyr::gather(gene_count_df, "Set", "Count", -Trait)
    
    # Create a bar plot showing the number of genes in G_burden, G_TWAS, and common for the selected trait
    plot <- ggplot(data = gene_count_df_long, aes(x = Trait, y = Count, fill = Set)) +
      geom_bar(stat = "identity", position = "dodge") +
      geom_text(aes(label = Count), vjust = -0.5, position = position_dodge(0.9)) +
      ggtitle(paste("Number of genes in G_burden, G_TWAS and common for", phenotype)) +
      xlab("Trait") +
      ylab("Number of genes") +
      theme_minimal() +
      theme(legend.title = element_blank())
    print(plot)
    # Save the plot
    #ggsave(paste0(phenotype, "F_bar.png"), plot = plot, width = 10, height = 8)
  }
}


# Define input parameters
phenotypes <- "Height"
genebass_files <- "geneheight-burden-results-exomes_pLoF_continuous-50-both_sexes--irnt_2023_05_26_12_18_16.csv"
twas_files <- "50Height/50.dat"

# Call the function and store results
results <- bartesting1(phenotypes, genebass_files, twas_files)

# Optionally, print or write out the results
print(results)