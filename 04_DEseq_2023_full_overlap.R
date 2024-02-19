library(tidyverse)

first_up <- read.delim(file="first_batch/trinity_full_LFC2upOrdered_results.txt")
first_up_genes <- rownames(first_up)

# 'TRINITY_DN158534_c0_g1' %in% first_up_genes

first_down <- read.delim(file="first_batch/trinity_full_LFC2downOrdered_results.txt")
first_down_genes <- rownames(first_down)


second_up <- read.delim(file="second_batch/trinity_full_LFC2upOrdered_results.txt")
second_up_genes <- rownames(second_up)

second_down <- read.delim(file="second_batch/trinity_full_LFC2downOrdered_results.txt")
second_down_genes <- rownames(second_down)

sum(first_up_genes %in% second_up_genes) # 4
# >TRINITY_DN181545_c0_g2_i2
# >TRINITY_DN181545_c0_g1_i1
# >TRINITY_DN158534_c0_g1_i1 ***
# >TRINITY_DN186731_c1_g3_i2 ***

sum(first_down_genes %in% second_down_genes) # 0

first_up_overlap <- first_up[which(first_up_genes %in% second_up_genes),] %>% 
    mutate(batch = 'first') %>% rownames_to_column('geneID')
second_up_overlap <- second_up[which(second_up_genes %in% first_up_genes),] %>% 
    mutate(batch = 'second') %>% rownames_to_column('geneID')
up_overlap <- rbind(first_up_overlap, second_up_overlap) %>% 
    arrange(geneID)
