library(readxl)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(data.table)
library(janitor)
library(tidyverse)
library(rstatix)
library(patchwork)

# Read in HD modifier variant summary 
hd_mod_gwas_hits <- fread("data/hd_modifier_gwas_summary.csv")

# Add a delays versus hastens column
hd_mod_gwas_hits$modifier_type <- ifelse(hd_mod_gwas_hits$years>0,"delays","hastens")

# Find robust modifier genes 
# (i.e., those not only  supported by a single rare variant)
rare_var_hd_mods <- hd_mod_gwas_hits %>%
  filter(single_rare_variant=="Yes") %>%
  .$gene_symbol

# Add whether DNA repair gene or HTT
dna_repair_hits <- c("PMS1", "MLH1", "MSH3", "PMS2","FAN1","LIG1")
hd_mod_gwas_hits$pathway <- ifelse(hd_mod_gwas_hits$candidate_mod_gene %in% dna_repair_hits, "DNA repair", "Other")
hd_mod_gwas_hits$pathway <- ifelse(hd_mod_gwas_hits$candidate_mod_gene %like% "HTT", "HTT", hd_mod_gwas_hits$pathway)
hd_mod_gwas_hits$pathway <- ifelse(hd_mod_gwas_hits$single_rare_variant== "Yes", "Rare variant", hd_mod_gwas_hits$pathway)

# Plot years and -log10pvalue
ggplot(hd_mod_gwas_hits, aes(x=abs(years), y=maf, col=pathway, shape=modifier_type)) +
  stat_smooth(method = lm, formula = y ~ log(x), se=FALSE, fullrange=TRUE, aes(shape = NULL, col=NULL), colour="#E69F00", linetype=2, size=0.5) +
  geom_point(size=5) + theme_minimal() +
  ylab("Minor allele frequency (%)") + xlab("Effect size (years)") +
  geom_text_repel(label=hd_mod_gwas_hits$candidate_mod_gene, fontface="italic", max.overlaps = Inf) +
  theme(legend.position="bottom",legend.title=element_blank()) +
  scale_shape_manual(values=c(24,25)) +
  ylim(0, max(hd_mod_gwas_hits$maf)) +
  xlim(0, max(abs(hd_mod_gwas_hits$years))) +
  scale_color_manual(values = c("#8967AC", "#D71920", "#0067B4","#C0C0C0"))
ggsave("results/hd_aoo_gwas_maf_eff.pdf", height=5.36, width=7.58, units='in')

# Look at exome data
# Remove missense only analysis (i.e., no deleterious filtering)
hd_exome_stats <- read_excel("data/pmid_35379994_supp_tables.xlsx",sheet = "table_s4", na = "NA") 

# Remove non-synonymous only analyses
hd_exome_stats <- hd_exome_stats %>%
  filter(Gene %in% hd_mod_gwas_hits$candidate_mod_gene,
         Analysis!="NS")

# Find significant genes
sig_exome_genes <- hd_exome_stats %>% 
  filter(p_value<0.05) %>%
  pull(Gene) %>%
  unique()

# Load in exome sequencing allele frequency information
hd_exome_direction <- read_excel("data/pmid_35379994_supp_tables.xlsx",sheet = "table_s5", na = "NA") 

# Remove only missense analysis (i.e., no deleterious filtering)
hd_exome_direction <- hd_exome_direction %>%
  filter(Analysis != "NS")

# Plot with facet by analysis
hd_exome_direction %>%
  filter(Gene %in% sig_exome_genes) %>%
  ggplot(., aes(x = Gene, y = `Allele frequency`, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~ Analysis) +
  labs(x = "Gene", y = "Allele frequency") +
  scale_fill_manual(values = c("#b73779", "#5ec962")) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        strip.text = element_text(size=12, face ="bold"),
        axis.text.y = element_text(size=12, face="bold"),
        axis.text.x = element_text(angle = 45, hjust = 1,  size=12, colour = "#8967AC", face="bold.italic"))
ggsave("results/hd_mod_exome_sig_alleles_analysis.pdf", height=5.36*0.8, width=5.36, units='in')

# Look at OMIM information
omim_genemap <- fread("data/genemap2.txt") %>%
  clean_names()
# Get stats
# Entries
omim_genemap %>%
  filter(phenotypes!="") %>%
  count()
# Filter for HD modifiers
omim_genemap_hd <- omim_genemap %>% filter(approved_gene_symbol %in% hd_mod_gwas_hits$gene_symbol) %>%
  filter(phenotypes!="")
omim_genemap_hd_select <- omim_genemap_hd %>% select(approved_gene_symbol, phenotypes)

# Read in clinvar info
clinvar <- fread("data/gene_specific_summary.txt", 
                 na.strings="-")
colnames(clinvar)[1] <- "gene"
clinvar$Other_alleles <- clinvar$Total_alleles - clinvar$Alleles_reported_Pathogenic_Likely_pathogenic
# - indicates no entry, therefore set to zero
clinvar <- clinvar %>% 
  mutate_if(is.integer, ~replace(., is.na(.), 0))

# Get HD modifiers
clinvar_hd <- clinvar %>% 
  filter(gene %in% hd_mod_gwas_hits$gene_symbol) %>%
  filter(gene!="HTT") %>%
  select(gene, Alleles_reported_Pathogenic_Likely_pathogenic, Other_alleles) %>%
  rename("Pathogenic/likely pathogenic alleles" = Alleles_reported_Pathogenic_Likely_pathogenic,
         "Other alleles" = Other_alleles)
clinvar_hd_robust <- clinvar_hd %>%
  filter(!(gene %in% rare_var_hd_mods))

# Get stats
clinvar %>%
  summarize(Total_submissions_sum = sum(Total_submissions),
            Total_alleles_sum = sum(Total_alleles),
            plp_sum = sum(Alleles_reported_Pathogenic_Likely_pathogenic),
            plp_mean = mean(Alleles_reported_Pathogenic_Likely_pathogenic),
            Unique_genes_count = n_distinct(gene))

clinvar_hd_melt <- melt(clinvar_hd, id.vars = "gene")
clinvar_hd_melt_robust <- clinvar_hd_melt %>%
  filter(!(gene %in% rare_var_hd_mods))

# Total number of alleles
sum(clinvar_hd_melt$value)
sum(clinvar_hd_melt_robust$value)

# Number of Path/Likely path
sum(clinvar_hd_robust$`Pathogenic/likely pathogenic alleles`)
sum(clinvar_hd_robust$`Pathogenic/likely pathogenic alleles`)/sum(clinvar_hd_melt_robust$value)*100

# Plot clinvar entries
clinvar_hd_melt %>%
  filter(!(gene %in% rare_var_hd_mods)) %>%
  ggplot(., aes(x=reorder(gene, -value), y=value, fill=variable))+
  geom_bar(position = "stack", stat = "identity") + theme_minimal() +
  xlab("Gene") + ylab("Alleles") +
  theme(axis.text.y = element_text(size=10, face="bold"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=12, face="italic"),
        legend.position="bottom", legend.title = element_blank()) +
  scale_fill_manual(values=c("#000000", "#808080"))
ggsave("results//hd_mod_clinvar_robust.pdf", height=5.36, width=7.58, units='in')

# Assess L2G infomation for modifier genes from Open Targets Genetics
l2g_combined <-  fread("data/hd_modifier_otg_04_02_23.tsv")
# Only look at robust signals
rare_var_hd_mods <- c("PBX1", "SYT9", "TMEM119", "GSG1L", "ALPK2")
l2g_combined_robust <- l2g_combined %>%
  filter(!(symbol %in% rare_var_hd_mods))

# Filter for select associations (remove count based, etc.)
l2g_combined_robust_select <-l2g_combined_robust %>% 
  filter(!(study.traitReported %like% " cell")) %>% 
  filter(!(study.traitReported %like% " count")) %>%
  filter(!(study.traitReported %like% " level")) %>%
  filter(!(study.traitReported %like% " fraction")) %>%
  filter(!(study.traitReported %like% " pressure")) %>%
  filter(!(study.traitReported %like% "Mean ")) %>%
  filter(!(study.traitReported %like% "Lung ")) %>%
  filter(!(study.traitReported %like% "Protein quant")) %>%
  filter(!(study.traitReported %like% "Height")) %>%
  filter(!(study.traitReported %like% " percent")) %>%
  filter(!(study.traitReported %like% "Waist")) %>%
  filter(!(study.traitReported %like% "Platelet")) %>%
  filter(!(study.traitReported %like% "crit")) %>%
  filter(!(study.traitReported %like% " density")) %>%
  filter(!is.na(study.pmid))

# Complete data set: number of signals
l2g_combined_robust %>% 
  filter(L2G>0.5) %>%
  nrow()

# Complete data set: signals per gene
l2g_combined_robust %>% 
  group_by(symbol) %>% 
  summarise(n = n()) %>%
  mutate(freq = n / sum(n) *100) %>%
  arrange(-n)

# Fine-mapped data set: number of signals
l2g_combined_robust %>% 
  filter(L2G>0.5) %>%
  nrow()

# Fine-mapped data set: signals per gene
l2g_combined_robust %>% 
  filter(L2G>0.5) %>%
  group_by(symbol) %>% 
  summarise(n = n()) %>%
  mutate(freq = n / sum(n) *100) %>%
  arrange(-n)


# Selected traits: signals per gene
l2g_combined_robust_select %>% 
  group_by(symbol) %>% 
  summarise(n = n()) %>%
  mutate(freq = n / sum(n) *100) %>%
  arrange(-n)

# Fine-mapped data set, selected traits: number of signals
l2g_combined_robust_select %>% 
  filter(L2G>0.5) %>% 
  nrow()

# Fine-mapped data set, selected traits: signals per gene
l2g_combined_robust_select %>% 
  filter(L2G>0.5) %>%
  group_by(symbol) %>% 
  summarise(n = n()) %>%
  mutate(freq = n / sum(n) *100) %>%
  arrange(-n)

# Add gene trait for plotting
l2g_combined_robust_select$gene_trait <- paste(l2g_combined_robust_select$symbol, l2g_combined_robust_select$study.traitReported, sep = " ")

# Add pathway information
dna_repair_hits <- c("PMS1", "MLH1", "MSH3", "PMS2","FAN1","LIG1")
l2g_combined_robust_select$pathway <- ifelse(l2g_combined_robust_select$symbol %in% dna_repair_hits, "DNA repair", "Other")

l2g_combined_robust_select %>% 
  filter(L2G>0.5) %>%
  nrow()

# Plot L2G results
l2g_combined_robust_select %>% 
  filter(L2G>0.5) %>%
  ggplot(., aes(x=L2G, y=reorder(gene_trait, L2G), col=pathway)) +
  geom_point(size=4) + theme_bw() + ylab("") +
  scale_x_continuous(limits = c(0, 1)) +
  geom_segment(aes(x=0, xend=L2G, 
                   y=gene_trait, yend=gene_trait)) +
  geom_vline(xintercept=0.5, linetype="dashed", color = "black") +
  theme(axis.text=element_text(size=10), legend.title=element_blank(),
        axis.title = element_text(size=12, face="bold"), legend.position="bottom") +
  scale_color_manual(values=c("#8967AC", "#0067B4")) 
ggsave("results/hd_aoo_otg_l2g_lolipop_robust.pdf", height=5.36, width=7.58*1.5, units='in')

# gnomAD constraint analysis
gnomad <- fread("data/gnomad.v2.1.1.lof_metrics.by_gene.txt")
hd_gnomad <- gnomad %>% filter(gene %in% hd_mod_gwas_hits$gene_symbol)
hd_gnomad$pathway <- ifelse(hd_gnomad$gene %in% dna_repair_hits, "DNA repair", "Other")
hd_gnomad$pathway <- ifelse(hd_gnomad$gene %like% "HTT", "HTT", hd_gnomad$pathway)
hd_gnomad$pathway <- ifelse(hd_gnomad$gene %in% rare_var_hd_mods, "Rare variant", hd_gnomad$pathway)

# Only look at robust signals
hd_gnomad_robust <- hd_gnomad %>%
  filter(!(gene %in% rare_var_hd_mods))
ggplot(hd_gnomad_robust, aes(x=-oe_lof, y=pLI, col=pathway)) + 
  geom_point() + theme_bw() + 
  geom_hline(yintercept=0.9,linetype=3) +
  scale_color_manual(values = c("#8967AC", "#D71920", "#0067B4", "#C0C0C0")) +
  labs(x="gnomAD -(observed/expected ratio LoF)", y="gnomAD pLI") +
  geom_text_repel(data=hd_gnomad_robust,aes(label=gene, fontface="bold.italic"), max.overlaps = Inf) +
  theme(legend.position="bottom", legend.title=element_blank(), 
        axis.text.y = element_text(size=12, face="bold"),
        axis.text.x = element_text(size=12, face="bold")) +
  theme(axis.text.x = element_text(size=10, face="bold"),
        axis.text.y = element_text(size=10, face="bold"))
ggsave("results/hd_mod_gnomad_constraint_robust.pdf", height=5.36, width=5.36, units='in')

# Make a master table that has information required for plotting
hd_mod_gene_info <- hd_mod_gwas_hits %>%
  select(gene_symbol, ensembl_id, 
         mod_gene_effect, pathway,
         hd_exome_support, cross_repeat_expansion,
         hpsc_msn_evidence) %>%
  distinct()

# Add min GWAS p-value by gene
hd_mod_gwas_min_p <- hd_mod_gwas_hits %>%
  group_by(gene_symbol) %>%
  summarise(min_gwas_p=min(continuous_p_value)) %>%
  distinct()

hd_mod_gene_info <- merge(hd_mod_gene_info, hd_mod_gwas_min_p, by="gene_symbol")
rm(hd_mod_gwas_min_p)

# Merge clinvar
hd_mod_gene_info <- merge(hd_mod_gene_info, clinvar_hd, by.x="gene_symbol", by.y="gene", all.x=T)

# Read in Druggable Genome Information (https://www.dgidb.org/, analyzed 7 June 2023)
dgib <- fread("data/dgidb_export_2023-06-07.tsv")
# Remove clinically actionable as this appears related to cancer
dgib <- dgib %>% 
  filter(category != "CLINICALLY ACTIONABLE")
table(dgib$category)
dgib_drug_genes <- dgib %>% 
  filter(category=="DRUGGABLE GENOME") %>%
  .$search_term
dgib_other_drug_genes <- dgib %>% 
  filter(category!="DRUGGABLE GENOME") %>%
  .$search_term
dgib_other_drug_genes <- unique(dgib_other_drug_genes)

# Add information to main data frame
hd_mod_gene_info$druggable_genome <- ifelse(hd_mod_gene_info$gene_symbol %in% dgib_drug_genes, "Yes", "No")
hd_mod_gene_info$druggable_other <- ifelse(hd_mod_gene_info$gene_symbol %in% dgib_other_drug_genes, "Yes", "No")
hd_mod_gene_info$druggable_any <- ifelse(hd_mod_gene_info$druggable_genome=="Yes"| hd_mod_gene_info$druggable_other=="Yes", "Yes", "No")

# Assess and add interaction data from Open Targets
# Subset high quality (mi score greater than 0.42)
# Threshold based on https://doi.org/10.1101/2023.02.07.23285407
hd_interactors <- fread("data/hd-mod-molecular-interactions-interactors.csv")

hd_int_high_q_count <- hd_interactors %>%
  filter(score>0.42) %>%
  group_by(gene) %>%
  summarise(int_0.42_count = n())
hd_int_count <- hd_interactors %>%
  group_by(gene) %>%
  summarise(int_all_count = n())

hd_mod_gene_info <- merge(hd_mod_gene_info, hd_int_high_q_count, by.x="gene_symbol", by.y="gene", all.x =TRUE)
hd_mod_gene_info <- merge(hd_mod_gene_info, hd_int_count, by.x="gene_symbol", by.y="gene", all.x =TRUE)
hd_mod_gene_info <- hd_mod_gene_info %>%
  mutate(int_0.42_count = replace_na(int_0.42_count, 0)) %>%
  mutate(int_all_count = replace_na(int_all_count, 0))
rm(hd_int_count, hd_int_high_q_count)

hd_mod_gene_info %>%
  filter(!(gene_symbol %in% rare_var_hd_mods)) %>%
ggplot(., aes(y=reorder(gene_symbol, int_0.42_count), x=int_0.42_count, fill=pathway))+
  geom_bar(position = "stack", stat = "identity") + theme_minimal() +
  ylab("Gene") + xlab("High quality interactions (Intact MI score > 0.42)") +
  geom_vline(xintercept=10, linetype=3) +
  theme(axis.text.x = element_text(size=10, face="bold"), 
        axis.text.y = element_text(size=12, face="bold.italic"),
        legend.position="bottom", legend.title = element_blank()) +
  scale_fill_manual(values=c( c("#8967AC", "#D71920", "#0067B4", "#C0C0C0"))) +
  geom_text(aes(label = int_0.42_count), 
            position = position_stack(vjust = 1), 
            hjust = -0.1, 
            size = 3.5, 
            fontface = "bold")
ggsave("results/hd_mod_interactors_robust.pdf", height=5.36, width=5.36, units='in')

# Read in protein atlas information
protein_atlas <- fread("data/proteinatlas.tsv", na.strings = "")
protein_atlas <- protein_atlas %>% 
  clean_names()
protein_atlas_select <- protein_atlas %>%
  dplyr::select(gene, rna_tissue_specificity, rna_single_cell_type_specificity)

hd_mod_gene_info <- merge(hd_mod_gene_info, protein_atlas_select, by.x="gene_symbol", by.y="gene", all.x=T)

# Plot RNA tissue specificity information for the different genes
# Count the number of genes in each rna_tissue_specificity category
data <- hd_mod_gene_info %>%
  group_by(rna_tissue_specificity) %>%
  mutate(gene_count = row_number()) %>%
  ungroup()

# Reorder the rna_tissue_specificity factor levels by the count of genes
ordered_levels <- data %>%
  count(rna_tissue_specificity) %>%
  arrange(desc(n)) %>%
  pull(rna_tissue_specificity)

# Order genes within each rna_tissue_specificity group by pathway
data <- data %>%
  mutate(rna_tissue_specificity = factor(rna_tissue_specificity, levels = ordered_levels)) %>%
  arrange(rna_tissue_specificity, pathway)

# Create a new column for the position of each gene within its group
data <- data %>%
  group_by(rna_tissue_specificity) %>%
  mutate(gene_position = row_number()) %>%
  ungroup()

# Plot data by enrichment type and pathway for robust genes
data %>%
  filter(!(gene_symbol %in% rare_var_hd_mods)) %>%
  ggplot(., aes(x = rna_tissue_specificity, y = gene_position, fill = pathway)) +
  geom_tile(color = "white") +
  geom_text(aes(label = gene_symbol), color = "white", fontface = "italic", size = 5) +
  labs(x = "RNA Tissue Specificity", y = "Gene Count", fill = "Pathway") +
  theme_minimal() +
  scale_fill_manual(values = c("#8967AC", "#D71920", "#0067B4", "#C0C0C0")) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, face = "bold"),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank())
ggsave("results/hd_mod_tissue_specificity_robust.pdf", height=7.58, width=7.58, units='in')
rm(data, ordered_levels)

# Look at gnomad data
gnomad_select <- gnomad %>%
  dplyr::select(gene, oe_lof, pLI)
hd_mod_gene_info <- merge(hd_mod_gene_info, gnomad_select, by.x="gene_symbol", by.y="gene", all.x=T)

# Add more columns for plotting
hd_mod_gene_info$genetic_evidence <- "Yes"
hd_mod_gene_info$knockdown_beneficial <- ifelse(hd_mod_gene_info$mod_gene_effect=="Hastens", "Yes", hd_mod_gene_info$mod_gene_effect)
hd_mod_gene_info$knockdown_beneficial <- ifelse(hd_mod_gene_info$mod_gene_effect=="Delays", "No", hd_mod_gene_info$knockdown_beneficial)
hd_mod_gene_info$knockdown_beneficial <- ifelse(hd_mod_gene_info$gene_symbol=="HTT", "Yes", hd_mod_gene_info$knockdown_beneficial)

# Add HTT as the effect
hd_mod_gene_info$mod_gene_effect <- ifelse(hd_mod_gene_info$gene_symbol=="HTT", "HTT", hd_mod_gene_info$mod_gene_effect)

# Base metrics on findings from https://www.nature.com/articles/s41588-024-01854-z
# https://www.nature.com/articles/s41588-024-01854-z
hd_mod_gene_info$pli_tolerance <- ifelse(hd_mod_gene_info$pLI<0.9, "Tolerant", "Intolerant")
hd_mod_gene_info$interacting_partners <- ifelse(hd_mod_gene_info$int_0.42_count<11, "Partners 0-10", "Partners 11 or more")

# Add whether the gene has been associated with a phenotypes in OMIM
hd_mod_gene_info$omim_phenotypes <- if_else(hd_mod_gene_info$gene_symbol %in% omim_genemap_hd$approved_gene_symbol, "omim_gene", "not_omim_gene") 

# Make a heatmap with key evidence
# Select key variables and reformat for plotting
hd_mod_gene_info_heatmap <- hd_mod_gene_info %>%
  select(gene_symbol, interacting_partners, 
         rna_tissue_specificity, pli_tolerance, 
         knockdown_beneficial, druggable_any, 
         genetic_evidence, omim_phenotypes,
         hd_exome_support, cross_repeat_expansion,
         hpsc_msn_evidence)

# Melt
hd_mod_gene_info_heatmap <- reshape2::melt(hd_mod_gene_info_heatmap, id="gene_symbol") %>%
  rename(feature = variable)

# Add a column with number of GWAS hits per locus
gwas_hits_loci <- hd_mod_gwas_hits %>%
  group_by(gene_symbol) %>%
  summarise(gwas_hits_locus=n())
hd_mod_gene_info_heatmap <- merge(hd_mod_gene_info_heatmap, gwas_hits_loci, by="gene_symbol")
rm(gwas_hits_loci)

# Perform replacements to group under favorable, unfavorable, etc.
hd_mod_gene_info_heatmap$value <- gsub("Unknown", NA, hd_mod_gene_info_heatmap$value)
hd_mod_gene_info_heatmap$value <- gsub("Not associated", 0, hd_mod_gene_info_heatmap$value)
hd_mod_gene_info_heatmap$value <- gsub("Group enriched", 1, hd_mod_gene_info_heatmap$value)
hd_mod_gene_info_heatmap$value <- gsub("Tissue enhanced", 1, hd_mod_gene_info_heatmap$value)
hd_mod_gene_info_heatmap$value <- gsub("Low tissue specificity", 0, hd_mod_gene_info_heatmap$value)
hd_mod_gene_info_heatmap$value <- gsub("Intolerant", 0, hd_mod_gene_info_heatmap$value)
hd_mod_gene_info_heatmap$value <- gsub("Tolerant", 1, hd_mod_gene_info_heatmap$value)
hd_mod_gene_info_heatmap$value <- gsub("No", 0, hd_mod_gene_info_heatmap$value)
hd_mod_gene_info_heatmap$value <- gsub("Yes", 1, hd_mod_gene_info_heatmap$value)
hd_mod_gene_info_heatmap$value <- gsub("Partners 0-10", 1, hd_mod_gene_info_heatmap$value)
hd_mod_gene_info_heatmap$value <- gsub("Partners 11 or more", 0, hd_mod_gene_info_heatmap$value)
hd_mod_gene_info_heatmap$value <- gsub("not_omim_gene", 1, hd_mod_gene_info_heatmap$value)
hd_mod_gene_info_heatmap$value <- gsub("omim_gene", 0, hd_mod_gene_info_heatmap$value)

# Replace unknown with 0.5 score for later plot ordering
hd_mod_gene_info_heatmap$value[is.na(hd_mod_gene_info_heatmap$value)] <- 0.5

# Add a plot value
hd_mod_gene_info_heatmap$plot_value <- gsub(1, "Favorable", hd_mod_gene_info_heatmap$value)
hd_mod_gene_info_heatmap$plot_value <- gsub(0.5, "Unresolved", hd_mod_gene_info_heatmap$plot_value)
hd_mod_gene_info_heatmap$plot_value <- gsub(0, "Unfavorable", hd_mod_gene_info_heatmap$plot_value)

# Re-level feature factor for plotting
gene_totals <- hd_mod_gene_info_heatmap %>%
  group_by(gene_symbol) %>%
  summarise(total_value = sum(as.numeric(as.character(value)), na.rm = TRUE))

desired_order <- c("interacting_partners", "rna_tissue_specificity", "pli_tolerance", 
                   "druggable_any", "knockdown_beneficial", "omim_phenotypes", 
                   "hpsc_msn_evidence", "hd_exome_support", "cross_repeat_expansion", "genetic_evidence")

# Re-level the feature variable with the desired order
hd_mod_gene_info_heatmap$feature <- factor(hd_mod_gene_info_heatmap$feature, levels = desired_order)
gene_totals <- gene_totals[order(gene_totals$total_value, decreasing = TRUE), ]
hd_mod_gene_info_heatmap$gene_symbol <- factor(hd_mod_gene_info_heatmap$gene_symbol, levels = gene_totals$gene_symbol)
rm(desired_order)
rm(gene_totals)

# Make a dotplot with GWAS hits 
single_rare_variant <- hd_mod_gwas_hits %>%
  filter(single_rare_variant=="Yes") %>%
  .$gene_symbol

# Only keep robust associations and remove single rare
hd_mod_dot_plot_robust <- hd_mod_gene_info_heatmap %>%
  filter(!(gene_symbol %in% rare_var_hd_mods)) %>%
  select(gene_symbol, gwas_hits_locus) %>%
  distinct() %>%
  mutate(is_single_rare_variant = gene_symbol %in% single_rare_variant) %>%
  ggplot(., aes(x = gene_symbol, y = gwas_hits_locus)) +
  theme_minimal() +
  ylim(0, max(hd_mod_gene_info_heatmap$gwas_hits_locus)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 2) +
  labs(y = "HD AOO GWAS hits", x=NULL) +
  theme(legend.title=element_blank(),
        axis.text.y = element_text(size=12, face="bold"),
        axis.text.x = element_blank())
hd_mod_dot_plot_robust

# Make two heatmaps
# 1) Human genetic evidence heatmap
human_genetic_features <- c("genetic_evidence", "omim_phenotypes", 
                            "hd_exome_support", "cross_repeat_expansion", 
                            "hpsc_msn_evidence")

hd_mod_human_heatmap_plot_robust <- hd_mod_gene_info_heatmap %>%
  filter(feature %in% human_genetic_features) %>%
  filter(!(gene_symbol %in% single_rare_variant)) %>%
  ggplot(., aes(x=gene_symbol, y=feature, fill=plot_value)) + 
  geom_tile(color="white", size=0.5) +
  coord_equal() +
  theme_minimal() +
  labs(x=NULL, y=NULL) +
  scale_fill_manual(values = c("#56B4E9", "#D55E00", "#7F7F7F")) +
  theme(legend.position="right", legend.title=element_blank(), 
        axis.text.y = element_text(size=12, face="bold"),
        axis.text.x = element_blank())
hd_mod_human_heatmap_plot_robust

# 2) Theoretical druggability heatmap
hd_mod_drug_heatmap_plot_robust <- hd_mod_gene_info_heatmap %>%
  filter(!(feature %in% human_genetic_features)) %>%
  filter(!(gene_symbol %in% single_rare_variant)) %>%
  ggplot(., aes(x=gene_symbol, y=feature, fill=plot_value)) + 
  geom_tile(color="white", size=0.5) +
  coord_equal() +
  theme_minimal() +
  labs(x=NULL, y=NULL) +
  scale_fill_manual(values = c("#21918c","#440154","#7F7F7F")) +
  theme(legend.position="right", legend.title=element_blank(), 
        axis.text.y = element_text(size=12, face="bold"),
        axis.text.x = element_text(size=12, face="bold.italic")) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))
hd_mod_drug_heatmap_plot_robust

# Combine the plots using patchwork (dotplot on top of heatmap)
hd_mod_dot_plot_robust/hd_mod_human_heatmap_plot_robust/hd_mod_drug_heatmap_plot_robust + 
  plot_layout(ncol = 1,  heights = c(1.5, 2, 2), 
              axes = "collect")
ggsave("results/hd_mod_druggability_hits_robust.pdf", height=5.36, width=7.58*0.8, units='in')