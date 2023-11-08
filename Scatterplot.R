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
        gene_data <- gene_data[!is.na(gene_data$P.Value.Burden), 
            ]
        twas_data <- data.table::fread(twas_file, data.table = FALSE)
        twas_data <- twas_data[!is.na(twas_data$TWAS.P), ]
        twas_gtex <- twas_data %>% filter(stringr::str_detect(twas_data$PANEL, 
            "GTEx"))
        twas_lasso <- twas_gtex[twas_gtex$MODEL == "lasso", ]
        all_genes <- merge(gene_data, twas_lasso, by.x = "Gene.Name", 
            by.y = "ID")
        all_genes <- all_genes %>% distinct(Gene.Name, .keep_all = TRUE)
        plot <- ggplot(data = all_genes, aes(x = -log10(P.Value.Burden), 
            y = -log10(TWAS.P))) + geom_point(color = "blue", 
            alpha = 0.5, size = 0.7) + ggtitle(paste("All genes for", 
            phenotype)) + xlab("-log10(Gene Bass Burden P-value)") + 
            ylab("-log10(TWAS P-value)") + geom_vline(xintercept = -log10(6.7e-07), 
            linetype = "dashed", color = "red") + geom_hline(yintercept = -log10(0.05), 
            linetype = "dashed", color = "red") + theme_bw() + 
            theme(plot.title = element_text(size = 16, face = "bold"))
        filename <- paste0(phenotype, "_scatter.png")
        ggsave(filename, plot = plot)
    }
}
