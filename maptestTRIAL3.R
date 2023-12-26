maptesting1 <- function(phenotypes, genebass_files, twas_files) {
  library(dplyr)
  library(ggplot2)
  library(data.table) # For fread function
  library(purrr) # For map functions
  
  results <- list()
  
  for (i in 1:length(phenotypes)) {
    phenotype <- phenotypes[i]
    genebass_file <- genebass_files[i]
    twas_file <- twas_files[i]
    
    # Read in data
    gene_data <- read.csv(file = genebass_file)
    twas_data <- fread(twas_file, data.table = FALSE)
    
    # Filter Genebass data and get significant genes, their positions, and chromosome
    gene_data_filtered <- gene_data %>%
      filter(!is.na(P.Value.Burden), P.Value.Burden < 6.7e-7) %>%
      mutate(Bonferroni.P.Burden = p.adjust(P.Value.Burden, method = "bonferroni")) #remove they adjusted already
    
    # Before proceeding, make sure there is data to work with
    if(nrow(gene_data_filtered) == 0) {
      stop("Filtered gene data is empty. Check your filtering criteria.")
    }
    
    # Continue with the data processing...
    G_burden <- gene_data_filtered %>%
      select(Gene.Name, Chrom, Position) # Select the gene names, chromosome, and positions
    
    # Filter TWAS data and get significant genes, their positions, and chromosome
    twas_data_filtered <- twas_data %>%
      filter(stringr::str_detect(PANEL, 'GTEx'), MODEL == 'lasso') %>%
      mutate(Bonferroni.P = p.adjust(TWAS.P, method = "bonferroni"))
    
    G_TWAS <- twas_data_filtered %>%
      filter(Bonferroni.P < 0.05) %>%
      select(ID, CHR, P0, P1)
    
    # Find intersection of significant G_TWAS and G_burden
    G_common <- intersect(G_TWAS$ID, G_burden$Gene.Name)
    
    # Nearby genes analysis within various distances
    proximity_distances <- c(50000, 100000, 150000, 200000) # distances in base pairs
    nearby_lists <- list()
    for (distance in proximity_distances) {
      nearby_name <- paste0("Nearby_TWAS_", distance/1000, "kb")
      gene_data_filtered[[nearby_name]] <- map_lgl(seq_along(G_burden$Position), ~{    #check every Genebass gene for proximity to any TWAS gene. It returns a logical vector (TRUE or FALSE) indicating whether each Genebass gene passes the proximity condition.
        current_chrom <- G_burden$Chrom[.x]             #store the chromosome and position of the current Genebass gene being checked    #Generates a sequence along the Position vector of G_burden, essentially providing indices to iterate over
        current_position <- G_burden$Position[.x]  # store the position of the current Genebass gene being checked
        any(
          G_TWAS$CHR == current_chrom &
            G_TWAS$P0 - distance <= current_position &
            current_position <= G_TWAS$P1 + distance
        )
      })
      
      # Filter based on the calculated proximities and store in a list for later use
      nearby_lists[[nearby_name]] <- gene_data_filtered %>%
        filter(!!sym(nearby_name)) %>%
        select(Gene.Name, Chrom, Position, !!sym(nearby_name)) # Select only relevant columns for output
    }
    # Group genes into loci based on proximity and chromosome concordance
    # Position has no NAs and is finite
    if (any(is.na(gene_data_filtered$Position)) || length(gene_data_filtered$Position) == 0) {
      stop("Positions are either NA or empty.")
    }
    
    loci_groups <- gene_data_filtered %>%
      mutate(Locus = cut(Position,
                         breaks = seq(from = min(Position, na.rm = TRUE),
                                      to = max(Position, na.rm = TRUE),
                                      by = 500000),
                         labels = FALSE)) %>%
      group_by(Chrom, Locus) %>%
      summarise(Genes_in_Locus = list(Gene.Name), .groups = 'drop')
    
    # Store results for this phenotype
    results[[phenotype]] <- list(
      G_TWAS = G_TWAS,
      G_burden = G_burden,
      common_genes = G_common,
      nearby_lists = nearby_lists, # This should be populated from your existing nearby genes analysis code
      loci_groups = loci_groups
    )
  } # End of for loop
  
  return(results)
}


# Define input parameters
phenotypes <- "Height"
genebass_files <- "geneheight-burden-results-exomes_pLoF_continuous-50-both_sexes--irnt_2023_05_26_12_18_16.csv"
twas_files <- "50Height/50.dat"

# Call the function and store results
results <- maptesting1(phenotypes, genebass_files, twas_files)
print(results)
