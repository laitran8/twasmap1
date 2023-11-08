function (phenotypes, genebass_files, twas_files) 
{
    library(dplyr)
    library(ggplot2)
    if (length(phenotypes) != length(genebass_files) | length(phenotypes) != 
        length(twas_files)) {
        stop("Number of phenotypes and files do not match!")
    }
    for (i in 1:length(phenotypes)) {
        phenotype <- phenotypes[i]
        genebass_file <- genebass_files[i]
        twas_file <- twas_files[i]
        gene_data <- read.csv(file = genebass_file)
        twas_data <- data.table::fread(twas_file, data.table = FALSE)
        twas_gtex_lasso <- twas_data %>% filter(stringr::str_detect(PANEL, 
            "GTEx"), MODEL == "lasso")
        twas_gtex_lasso$Bonferroni.P <- p.adjust(twas_gtex_lasso$TWAS.P, 
            method = "bonferroni")
        GTWAS <- unique(twas_gtex_lasso[twas_gtex_lasso$Bonferroni.P < 
            0.05, "ID"])
        Gburden <- gene_data %>% filter(Gene.Name %in% GTWAS)
        Gburden$P.Value.Burden <- p.adjust(Gburden$P.Value.Burden, 
            method = "bonferroni")
        Gburden <- Gburden[Gburden$P.Value.Burden < 0.05, ]
        num_Gburden <- nrow(Gburden)
        num_GTWAS <- length(GTWAS)
        gene_burden_strict <- gene_data[gene_data$P.Value.Burden < 
            6.7e-07, ]
        gene_burden_strict <- gene_burden_strict[!is.na(gene_burden_strict$P.Value.Burden), 
            ]
        twas_lasso_strict <- twas_gtex_lasso[p.adjust(twas_gtex_lasso$TWAS.P, 
            method = "bonferroni") < 0.05, ]
        twas_genes_strict <- unique(twas_lasso_strict$ID)
        common_genes <- intersect(gene_burden_strict$Gene.Name, 
            twas_genes_strict)
        num_common <- length(common_genes)
        gene_count_df <- data.frame(Trait = phenotype, Gburden = num_Gburden, 
            GTWAS = num_GTWAS, Common = num_common)
        gene_count_df_long <- tidyr::gather(gene_count_df, "Set", 
            "Count", -Trait)
        plot <- ggplot(data = gene_count_df_long, aes(x = Trait, 
            y = Count, fill = Set)) + geom_bar(stat = "identity", 
            position = "dodge") + geom_text(aes(label = Count), 
            vjust = -0.5, position = position_dodge(0.9)) + ggtitle(paste("Number of genes in Gburden, GTWAS and common for", 
            phenotype)) + xlab("Trait") + ylab("Number of genes") + 
            theme_bw() + theme(plot.title = element_text(size = 16, 
            face = "bold"), legend.title = element_blank())
        filename <- paste0(phenotype, "_bar.png")
        ggsave(filename, plot = plot, width = 10, height = 8)
    }
}
