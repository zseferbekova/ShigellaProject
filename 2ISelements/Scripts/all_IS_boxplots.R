## Import dependencies
library(ggplot2)
library(reshape2)
library(ggpubr)
library(tidyverse)

## Import total_stats.csv
stats = as.data.frame(read.csv("../Data/total_stats.csv", sep=",", stringsAsFactors = FALSE)) 

## Filter out plasmids and unnecessary columns
## Add a column with ISsaga counts (take the lowest number of IS)
stats_for_pic <- stats %>% 
  ungroup() %>%
  filter(Type == 'Chromosome') %>%
  mutate(Pathovar = Info,
         ISnumber = as.numeric(str_split(ISTransposaseCount.ISsaga., "/", simplify = T)[,1])) %>%
  select(Organism, Pathovar, ISnumber)

## Add pathogenic group, non-pathogenic group
stats_for_pic$Pathovar <- str_replace(stats_for_pic$Pathovar, "shigellosis", "Shigella spp.")
stats_for_pic$Pathovar <- stats_for_pic$Pathovar %>%
  replace(stats_for_pic$Pathovar == "", "non-pathogenic\nE. coli")

path_ecoli <- stats_for_pic %>%
  filter(Pathovar %in% c("EAEC","ExPEC","ETEC","STEC","EPEC","APEC","AIEC","EIEC","EHEC"))
path_ecoli$Pathovar <- "pathogenic\nE. coli"

all_ecoli <- stats_for_pic %>%
  filter(Pathovar != "Shigella spp.") 
all_ecoli$Pathovar <- "all E. coli"

stats_for_pic <- rbind(all_ecoli, path_ecoli, stats_for_pic)

## Set levels for a pretty plot and create a comparison list

## uncomment for 'STEC, ExPEC vs others' plot
# my_comparisons <- list( c("all E. coli", "STEC"),
#                         c("non-pathogenic\nE. coli", "STEC"),
#                         c("ExPEC", "all E. coli"),
#                         c("ExPEC", "non-pathogenic\nE. coli"),
#                         c("ExPEC", "STEC"))
# IS_numbers <- stats_for_pic %>%
#   filter(Pathovar %in% c("all E. coli", "non-pathogenic\nE. coli", "STEC", "ExPEC"))
# IS_numbers$Pathovar <- factor(IS_numbers$Pathovar, levels=c("all E. coli", "non-pathogenic\nE. coli", "STEC", "ExPEC"))

## uncomment for 'Shigella vs others' plot
# my_comparisons <- list( c("all E. coli", "Shigella spp."),
#                         c("non-pathogenic\nE. coli", "Shigella spp."),
#                         c("pathogenic\nE. coli", "Shigella spp."))
# IS_numbers <- stats_for_pic %>%
#   filter(Pathovar %in% c("all E. coli", "non-pathogenic\nE. coli", "pathogenic\nE. coli", "Shigella spp."))
# IS_numbers$Pathovar <- factor(IS_numbers$Pathovar, levels=c("all E. coli", "non-pathogenic\nE. coli", "pathogenic\nE. coli", "Shigella spp."))

## uncomment for a plot with all pathovars
# IS_numbers <- stats_for_pic %>%
#      filter(!(Pathovar %in% c("pathogenic\nE. coli", "all E. coli")))
# IS_numbers$Pathovar <- factor(IS_numbers$Pathovar, levels=c("non-pathogenic\nE. coli", "Shigella spp.", "STEC", "ExPEC", "ETEC", "APEC", "EPEC", "AIEC", "EAEC", "EHEC", "EIEC"))

## Boxplot
ggplot(data = IS_numbers, 
       aes(x = Pathovar, y = ISnumber)) +
  
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
  #theme(axis.text.x = element_text(angle = -15, hjust = 0.05)) +
  xlab("") +
  ylab("Number of IS elements") 

ggsave("<Name>.png", device = "png", dpi = 720, width = 4, height = 3)
