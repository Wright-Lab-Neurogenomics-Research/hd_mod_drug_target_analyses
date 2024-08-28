library(data.table)
library(ggplot2)
library(tidyverse)
library(readxl)
library(ggrepel)

# Read in HD top variants
hd_mod_gwas_hits <- fread("data/hd_modifier_gwas_summary.csv")
# Remove years etc. as this will be obtained from the GWAS summary stats
hd_mod_gwas_hits <- hd_mod_gwas_hits %>% select(-minor_allele, -maf, -years, -continuous_p_value)

# HD summary stats
gem_hd_9k <- fread('data/gem.euro.9k.summary.chr1-23.210125.txt')

# Rename all columns to have hd in front except SNP
# Use lower case for names
# Divide allele frequency by 100 to standardize with other data frames
gem_hd_9k <- gem_hd_9k %>%
  rename_with(., tolower) %>%
  rename(., hd_effect_allele = test_allele,
         hd_effect_allele_frequency = test_allele_frequency,
         hd_beta = beta,                      
         hd_se = se, 
         hd_p_value = "p-value") %>%
  mutate(hd_effect_allele_frequency = hd_effect_allele_frequency/100)

# Read in XDP summary stats (MAF cutoff 1%)
xdp_maf1 <- fread('data/XDP_summary_stats_MAF1.tsv')
# Perform column standardization as above with HD
xdp_maf1 <-   xdp_maf1 %>%
  rename_with(., tolower) %>%
  dplyr::select(-chromosome,-base_pair_position) %>%
  rename(., snp = variant_id,
         xdp_effect_allele = effect_allele,
         xpd_other_allele = other_allele,
         xdp_effect_allele_frequency = effect_allle_frequency,
         xdp_beta = beta,                      
         xdp_se = standard_error, 
         xdp_p_value = p_value)

# Merge GWAS together
repeat_expan_gwas <- merge(gem_hd_9k, xdp_maf1, by="snp", all.x=T)
repeat_expan_gwas_top <- merge(hd_mod_gwas_hits, repeat_expan_gwas, by="snp")

# Add whether DNA repair gene or HTT
dna_repair_hits <- c("PMS1", "MLH1", "MSH3", "PMS2","FAN1","LIG1")
repeat_expan_gwas_top$pathway <- ifelse(repeat_expan_gwas_top$candidate_mod_gene %in% dna_repair_hits, "DNA repair", "Other")
repeat_expan_gwas_top$pathway <- ifelse(repeat_expan_gwas_top$candidate_mod_gene %like% "HTT", "HTT", repeat_expan_gwas_top$pathway)


# Perform ar regression to calculate R2 and significance
hd_xdp_beta_reg <- summary(lm(repeat_expan_gwas_top$hd_beta ~ repeat_expan_gwas_top$xdp_beta))
hd_xdp_beta_reg
hd_xdp_beta_r2 <- round(hd_xdp_beta_reg$r.squared, digits=2)
hd_xdp_beta_p <- signif(hd_xdp_beta_reg$coefficients[2,4], digits=2)

# Nominally significant (to keep in line with HD exomes)
xdp_significant <- 0.05

# Plot effect sizes in HD and XDP GWAS (color by pathway)
repeat_expan_gwas_top %>%
  filter(complete.cases(xdp_p_value)) %>%
  ggplot(., aes(x=hd_beta, y=xdp_beta, col=pathway, shape=xdp_p_value < 0.05)) +
  geom_point() + theme_minimal() +
  stat_smooth(method = lm, formula = y ~ x, fullrange=TRUE, aes(shape = NULL, col=NULL), colour="#E69F00", fill="#fcd886", linetype=2, size=0.5) +
  geom_pointrange(aes(ymin=xdp_beta-xdp_se, ymax=xdp_beta+xdp_se)) +
  geom_pointrange(aes(xmin=hd_beta-hd_se, xmax=hd_beta+hd_se)) +
  scale_color_manual(values=c("#8967AC", "#0067B4")) +
  scale_shape_manual(values = c(1,19)) +
  geom_vline(xintercept = 0, size = 0.75, color = "grey") +
  geom_hline(yintercept = 0, size = 0.75, color = "grey") +
  xlab("HD  GWAS variant effect size (years)") + ylab("XDP GWAS variant effect size (years)") +
  theme(legend.position="bottom",legend.title=element_blank()) +
  ggtitle(paste0("Correlation between HD and XDP variant effect sizes\nR2 =", hd_xdp_beta_r2, ", P-value = ", hd_xdp_beta_p)) +
  geom_text_repel(aes(label=candidate_mod_gene, fontface="bold.italic"))
ggsave("results/hd_xdp_pathway.pdf", height=5.36, width=5.36, units='in')
