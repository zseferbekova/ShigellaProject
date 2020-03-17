## Load dependencies
library(tidyverse)
library(ggplot2)
library(ggpubr)

## Read a file with statistics
stats <- data.frame(read.csv("../Data/total_stats.csv", sep=",", stringsAsFactors = FALSE))

## Calculate total genome size (chromosome + plasmid(s))
stats <- stats %>% 
  group_by(AssemblyID) %>%
  mutate(length_total.bp. = sum(length.bp.), 
         total_CDSs = sum(CDSs))

## Filter out plasmids and unnecessary columns
stats_for_pic <- stats %>% 
  ungroup() %>%
  filter(Type == 'Chromosome') %>%
  mutate(total_len_kb = round(length_total.bp./1000), 0) %>%
  select(organism = Organism, chrom_len_bp = length.bp., chrom_len_cds = CDSs, 
         total_len_kb, total_len_cds = total_CDSs,
         pathovar = Info,
  )

## Add a column with ISsaga counts (take the lowest number of IS)
stats_for_pic$is_num <- str_split(filter(stats, Type == "Chromosome")$ISTransposaseCount.ISsaga., "/", simplify = T)[,1]

## Add pathogenic group, non-pathogenic group
stats_for_pic$pathovar <- stats_for_pic$pathovar %>%
  replace(stats_for_pic$pathovar == "", "non-pathogenic\nE. coli")

path_ecoli <- stats_for_pic %>%
  filter(pathovar %in% c("EAEC","ExPEC","ETEC","STEC","EPEC","APEC","AIEC","EIEC","EHEC"))
path_ecoli$pathovar <- "pathogenic\nE. coli"

stats_for_pic$pathovar <- str_replace(stats_for_pic$pathovar, "shigellosis", "Shigella spp.")
stats_for_pic$is_num <- as.numeric(stats_for_pic$is_num)

all_ecoli <- stats_for_pic %>%
  filter(pathovar != "Shigella spp.") 
all_ecoli$pathovar <- "all E. coli"

stats_for_pic <- rbind(all_ecoli, path_ecoli, stats_for_pic)

## Set levels for a pretty plot and create a comparison list

## uncomment for 'STEC, ExPEC vs others' plot
# my_comparisons <- list( c("all E. coli", "STEC"),
#                         c("non-pathogenic\nE. coli", "STEC"),
#                         c("all E. coli", "ExPEC"),
#                         c("non-pathogenic\nE. coli", "ExPEC"),
#                         c("ExPEC", "STEC"))
# genome_sizes <- stats_for_pic %>%
#   filter(pathovar %in% c("all E. coli", "non-pathogenic\nE. coli", "STEC", "ExPEC"))
# genome_sizes$pathovar <- factor(genome_sizes$pathovar, levels=c("all E. coli", "non-pathogenic\nE. coli", "STEC", "ExPEC"))

## uncomment for 'Shigella vs others' plot
# my_comparisons <- list( c("Shigella spp.", "all E. coli"),
#                         c("Shigella spp.", "non-pathogenic\nE. coli"),
#                         c("Shigella spp.", "pathogenic\nE. coli"))
# genome_sizes <- stats_for_pic %>%
#   filter(pathovar %in% c("all E. coli", "non-pathogenic\nE. coli", "pathogenic\nE. coli", "Shigella spp."))
# genome_sizes$pathovar <- factor(genome_sizes$pathovar, levels=c("all E. coli", "non-pathogenic\nE. coli", "pathogenic\nE. coli", "Shigella spp."))

## uncomment for a plot with all pathovars
# genome_sizes <- stats_for_pic
# genome_sizes$pathovar <- factor(genome_sizes$pathovar, levels=c("non-pathogenic\nE. coli", "Shigella spp.", "STEC", "ExPEC", "ETEC", "APEC", "EPEC", "AIEC", "EAEC", "EHEC", "EIEC"))

## Boxplot
ggplot(data = genome_sizes, 
       aes(x = pathovar, y = total_len_kb)) +
  
  geom_boxplot(
    # custom outliers
    outlier.colour="red",
    outlier.fill="red",
    outlier.size=2,
    # custom boxes
    color="blue",
    fill="blue",
    alpha=0.2,
  ) +
  
  ## comment if plotting all pathovars
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", method.args = list(alternative = "less"), label = "p.signif",) +
  
  theme_classic() +
  theme(axis.text.x = element_text(angle = -15, hjust = 0.05)) +
  xlab("") +
  ylab("Genome size (kb)") 

ggsave("<Name>.png", device = "png", dpi = 720, width = 7, height = 4)
