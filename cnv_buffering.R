# change R working directory to the project folder if required
#setwd("/path/to/folder")

# load R packages
library(tidyverse)

# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


# read table with multi-omics data (protein, mRNA and CNV)
multi_omics <- read_tsv("./multi_omics_data.txt.gz")

## simplify column names
multi_omics <- multi_omics %>%
  rename(protein = prot_log2FC, rna = rna_log2CPM, cnv = cnv_gistic2)

# this is how you could do it with base R:
#colnames(multi_omics)[3] <- "protein"
#colnames(multi_omics)[4] <- "rna"
#colnames(multi_omics)[5] <- "cnv"
# not so cleaner...

# read table with metadata (cancer type, experimental study, etc.)
metadata <- read_tsv("./samples_metadata.txt")

# read table containing the samples with protein, mRNA and CNV data
samples <- read_tsv("./samples_prot_rna_cnv.txt", col_names = F)

## add cancer type
samples <- samples %>%
  rename(sample = X1) %>%
  inner_join(metadata[, c("sample", "cancer")], by = "sample")


# make a simple barplot with the number of samples by cancer type
## geom_bar() counts the number of cases at each x position by default, making the height of the bar proportional to that number
## therefore we do not provide an y variable (it is computed internally by ggplot2)
samples_barplot <- ggplot(data = samples, mapping = aes(x = cancer, fill = cancer)) +
  geom_bar() +
  theme_classic()
samples_barplot

# in case you want to save the plot to a pdf file uncomment the following line
#ggsave(filename = file.path(".", "samples_barplot.pdf"), plot = samples_barplot, height = 4, width = 4)



# CNV buffering analysis


# set up a simple correlation function
# df must be a data.frame or tibble
# x and y the variables to select and correlate from df
# uses pearson correlation by default
corr <- function(df, x, y, method = "pearson"){
  a = df[[x]]
  b = df[[y]]
  cor <- broom::tidy(cor.test(a, b, method = method))
  
  stats = c("corr", "pval")
  cols = rep(list(c(x, y)), 2)
  nms <- map2_chr(.x = stats, .y = cols, .f = ~ str_c(c(.x, .y), collapse = "_"))
  
  res <- cor %>%
    select(estimate, p.value) %>%
    rename_all(.funs = ~ nms)
  
  return(res)
}


# CNV and mRNA correlation for each gene across samples
# after nesting the data by gene, we call the corr() function for each gene
cnv_rna_corr <- multi_omics %>%
  group_by(gene) %>%
  nest() %>%
  ungroup() %>%
  mutate(correlations = map(.x = data, .f = corr, x = "cnv", y = "rna")) %>%
  select(-data) %>%
  unnest() %>%
  select(-pval_cnv_rna)

# this is how you could do it with base R:
#cnv_rna_corr <- data.frame()
#for(g in unique(multi_omics$gene)){
#  temp_df <- multi_omics[multi_omics$gene == g, ]
#  corrs <- corr(temp_df, x = "cnv", y = "rna")
#  corrs <- cbind(g, corrs)
#  colnames(corrs)[1] <- "gene"
#  cnv_rna_corr <- rbind(cnv_rna_corr, corrs)
#}
# besides not looking so clean it is also much slower...


# CNV and protein correlation for each gene across samples
cnv_protein_corr <- multi_omics %>%
  group_by(gene) %>%
  nest() %>%
  ungroup() %>%
  mutate(correlations = map(.x = data, .f = corr, x = "cnv", y = "protein")) %>%
  select(-data) %>%
  unnest() %>%
  select(-pval_cnv_protein)


# merge both correlations
cnv_rna_protein_corr <- inner_join(cnv_rna_corr, cnv_protein_corr, by = "gene")


# calculate the attenuation potential as the correlation differences: cor(CNV, RNA) - cor(CNV, Protein)
# higher attenuation potentials may indicate genes whose CNVs are buffered at the protein level
cnv_rna_protein_corr <- cnv_rna_protein_corr %>%
  mutate(attenuation = corr_cnv_rna - corr_cnv_protein) %>%
  arrange(desc(attenuation))
cnv_rna_protein_corr



# example of the CNV attenuation at the protein level for the gene TOX4:
# component of the PTW/PP1 phosphatase complex (plays a role in the control of cell cycle progression, during the transition from mitosis into interphase)
# let's filter the multi-omics data by selecting the TOX4 gene
attenuation_example <- multi_omics %>%
  filter(gene == "TOX4")


# we will plot the RNA and protein measures across the different CNV levels:
# -2 = homozygous deletion; -1 = hemizygous deletion; 0 = heterozygotic/no change; 1 = gain; 2 = high level amplification

# CNV vs RNA
cnv_rna_example <- attenuation_example %>%
  ggplot(mapping = aes(x = cnv, y = rna, group = cnv, fill = as.character(cnv))) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 0.5) +
  scale_fill_discrete(guide = F) +
  theme_classic()
cnv_rna_example

# in case you want to save the plot to a pdf file uncomment the following line
#ggsave(filename = file.path(".", "cnv_rna_example.pdf"), plot = cnv_rna_example, height = 4, width = 4)


# CNV vs protein
cnv_protein_example <- attenuation_example %>%
  ggplot(mapping = aes(x = cnv, y = protein, group = cnv, fill = as.character(cnv))) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 0.5) +
  scale_fill_discrete(guide = F) +
  theme_classic()
cnv_protein_example

# in case you want to save the plot to a pdf file uncomment the following line
#ggsave(filename = file.path(".", "cnv_protein_example.pdf"), plot = cnv_protein_example, height = 4, width = 4)


# can you do now the same for a gene that is not attenuated?
# let's filter the data for positive attenuation values to exclude the genes that are more correlated with the protein than with the RNA
# they are not many and we think they are the result of sequencing errors or splicing isoforms (we might not be measuring the correct isoform at the RNA level)
cnv_rna_protein_corr %>%
  filter(attenuation > 0) %>%
  arrange(attenuation)

# try to reproduce the same plot as before for the PDCD5 gene (involved in cellular apoptosis)
# start by filtering the multi-omics dataset
# you can re-use the code



# make a scatterplot with the cor(CNV, RNA) in the x axis and the cor(CNV, Protein) in the y axis
# the genes with the highest attenuation at the protein level will be in the bottom-right corner
attenuation_plot1 <- ggplot(data = cnv_rna_protein_corr, mapping = aes(x = corr_cnv_rna, y = corr_cnv_protein)) +
  geom_point(size=1) +
  geom_abline(slope=1, intercept=0, linetype=2, color="black") +
  theme_classic() +
  coord_fixed() +
  scale_x_continuous(limits = c(-0.1, 0.8)) +
  scale_y_continuous(limits = c(-0.1, 0.8)) +
  labs(x = "corr(CNV, RNA)", y = "corr(CNV, Protein)")
attenuation_plot1

# in case you want to save the plot to a pdf file uncomment the following line
#ggsave(filename = file.path(".", "attenuation_plot1.pdf"), plot = attenuation_plot1, height = 4, width = 4)

# color the points by attenuation potential
attenuation_plot2 <- ggplot(data = cnv_rna_protein_corr, mapping = aes(x = corr_cnv_rna, y = corr_cnv_protein, color = attenuation)) +
  geom_point(size=1) +
  geom_abline(slope=1, intercept=0, linetype=2, color="black") +
  theme_classic() +
  coord_fixed() +
  scale_color_gradient(low = "#56B1F7", high = "#132B43") +
  scale_x_continuous(limits = c(-0.1, 0.8)) +
  scale_y_continuous(limits = c(-0.1, 0.8)) +
  labs(x = "corr(CNV, RNA)", y = "corr(CNV, Protein)")
attenuation_plot2

# in case you want to save the plot to a pdf file uncomment the following line
#ggsave(filename = file.path(".", "attenuation_plot2.pdf"), plot = attenuation_plot2, height = 4, width = 4)



# perform a GO enrichment analysis on the more attenuated genes (top 1000 genes)
# we will use a GO over-representation test based on the hypergeometric distribution:
# check this page for further details: https://yulab-smu.top/biomedical-knowledge-mining-book/enrichment-overview.html#ora-algorithm

# select the top 1000 attenuated genes
attenuated_genes <- cnv_rna_protein_corr %>% 
  top_n(n = 1000,  wt = attenuation)

# the attenuated genes are supposed to be enriched in protein complexes
# therefore we will select the Cellular Component (CC) ontology
enr_object <- clusterProfiler::enrichGO(gene = attenuated_genes$gene, universe = cnv_rna_protein_corr$gene, ont = "CC", OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "SYMBOL")


# extract the enrichment table from the enrichResult object
# transform it to a tibble
# select the significantly enriched categories (p.adjust < 0.05)
enr_table <- enr_object@result %>%
  as_tibble() %>%
  filter(p.adjust < 0.05)


# make a barplot with the top 30 enriched categories
# transform the adjusted P-values using the negative value of the logarithm base 10 (-log10):
# low values will have the highest values, and vice-versa
# reorder the categories using the -log10 adjusted P-values
categories <- enr_table %>%
  mutate(log10_p = -log10(p.adjust)) %>%
  top_n(n = 30, wt = log10_p) %>%
  mutate(Description = fct_reorder(.f=Description, .x=log10_p, .fun = function(x) x))

categories_barplot <- ggplot(data = categories, mapping = aes(x = Description, y = log10_p, fill = Count)) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  labs(x = "GO term", y = "Adjusted P-value (-log10)")
categories_barplot

# in case you want to save the plot to a pdf file uncomment the following line
#ggsave(filename = file.path(".", "categories_barplot.pdf"), plot = categories_barplot, height = 6, width = 6)













# Does the attenuation depend on the cancer type?
# To answer this question we will perform the CNV buffering analysis by cancer type
# To make it simpler, let's just compare the breast and ovarian cancers


# the protein abundances have missing values, which might make it impossible to calculate the CNV buffering of some genes in some tissues
# to make sure we will calculate the CNV buffering for the same set of genes along the cancer types,
# we will select the genes that have at least 10 samples with protein measures in each cancer type


# set up a function that will calculate the number of samples per cancer type
# after that it evaluates whether or not the number of samples is at least 10 in the 2 tissues
# the output will be TRUE if the last statement is true, and FALSE otherwise
select_genes <- function(df){
  grouped <- df %>%
    group_by(cancer) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    summarise(higher_10 = sum(n >= 10) == 2) %>%
    pull(higher_10)
  
  return(grouped)
}

# first we remove the missing values from the protein abundance
# then we add the cancer type using the table with the metadata
# after nesting the data by gene, we call the select_genes() function for each gene
multi_omics2 <- multi_omics %>%
  filter(!is.na(protein)) %>%
  inner_join(metadata[, c("sample", "cancer")], by = "sample") %>%
  filter(cancer != "COREAD") %>%
  group_by(gene) %>%
  nest() %>%
  ungroup() %>%
  mutate(genes_to_keep = map_lgl(.x = data, .f = select_genes)) %>%
  filter(genes_to_keep) %>%
  select(-genes_to_keep) %>%
  unnest()

# faster way of doing the same with less code (and no need for an external function)
# multi_omics2 <- multi_omics %>%
#   filter(!is.na(protein)) %>%
#   inner_join(metadata[, c("sample", "cancer")], by = "sample") %>%
#   filter(cancer != "COREAD") %>%
#   group_by(gene) %>%
#   filter(sum(table(cancer) >= 10) == 2) %>%
#   ungroup()


# CNV and mRNA correlation for each gene across samples for each cancer type
cnv_rna_corr2 <- multi_omics2 %>%
  group_by(cancer, gene) %>%
  nest() %>%
  ungroup() %>%
  mutate(correlations = map(.x = data, .f = corr, x = "cnv", y = "rna")) %>%
  select(-data) %>%
  unnest() %>%
  select(-pval_cnv_rna)


# CNV and protein correlation for each gene across samples for each cancer type
cnv_protein_corr2 <- multi_omics2 %>%
  group_by(cancer, gene) %>%
  nest() %>%
  ungroup() %>%
  mutate(correlations = map(.x = data, .f = corr, x = "cnv", y = "protein")) %>%
  select(-data) %>%
  unnest() %>%
  select(-pval_cnv_protein)


# merge both correlations
cnv_rna_protein_corr2 <- inner_join(cnv_rna_corr2, cnv_protein_corr2, by = c("cancer", "gene"))


# calculate attenuation potential as before
cnv_rna_protein_corr2 <- cnv_rna_protein_corr2 %>%
  mutate(attenuation = corr_cnv_rna - corr_cnv_protein) %>%
  arrange(desc(attenuation))


# make a scatterplot as before
# however, this time we will do one plot for each cancer type using the ggplot2 function facet_wrap()
attenuation_plot3 <- ggplot(data = cnv_rna_protein_corr2, mapping = aes(x = corr_cnv_rna, y = corr_cnv_protein, color = attenuation)) +
  geom_point(size=1) +
  geom_abline(slope=1, intercept=0, linetype=2, color="black") +
  theme_classic() +
  theme(
    panel.spacing = unit(1, "cm"),
    strip.background = element_blank(),
    legend.position = "bottom") +
  coord_fixed() +
  facet_wrap(facets = vars(cancer), nrow = 1) +
  scale_color_gradient(low = "#56B1F7", high = "#132B43") +
  scale_x_continuous(limits = c(-0.2, 1)) +
  scale_y_continuous(limits = c(-0.2, 1)) +
  labs(x = "corr(CNV, RNA)", y = "corr(CNV, Protein)")
attenuation_plot3

# in case you want to save the plot to a pdf file uncomment the following line
#ggsave(filename = file.path(".", "attenuation_plot3.pdf"), plot = attenuation_plot3, height = 4, width = 8)



# we will compare the gene attenuation scores between the two cancer types
# first let's make a new table with the attenuation scores in two separate columns
attenuation_scores <- cnv_rna_protein_corr2 %>%
  arrange(gene) %>%
  select(gene, cancer, attenuation) %>%
  pivot_wider(names_from = "cancer", values_from = "attenuation")


# make a scatterplot with the attenuation scores from ovarian cancer in the x axis and from breast cancer in the y axis
# add also a Pearson correlation coefficient using the function stat_cor() from the ggpubr package
library(ggpubr)
brca_ov_attenuation <- ggplot(data = attenuation_scores, mapping = aes(x = OV, y = BRCA)) +
  geom_point() +
  geom_smooth(method = "lm") +
  stat_cor() +
  theme_classic()
brca_ov_attenuation

# in case you want to save the plot to a pdf file uncomment the following line
#ggsave(filename = file.path(".", "brca_ov_attenuation.pdf"), plot = brca_ov_attenuation, height = 4, width = 4)

