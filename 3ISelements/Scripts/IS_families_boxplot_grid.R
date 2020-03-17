## Import dependencies
library(tidyverse)
library(reshape2)
library(stringr)
library(ggpubr)

## Import ISsaga data for B1 phylogroup
IS_families <- as.data.frame(read.csv("../Data/B1_IS_elements_selected.csv", stringsAsFactors = F))
colnames(IS_families)[1] <- "strain"

## Import labels and map assembly names to strains
labels <- as.data.frame(read.csv("../Results/Trees/labels.txt", skip = 6, sep = "\t", row.names = 1, header=F, colClasses = c("character")))
IS_families$strain <- labels[IS_families$strain,]

## Melt to plot later
melted_f <- melt(IS_families,id.vars = "strain")
melted_f$species <- melted_f$strain %>%
  str_replace_all(c("Shigella"="S.", "Escherichia" = "E.")) %>%
  word(1, 2, " ")
melted_f$value <- as.numeric(melted_f$value)

## Change levels for a pretty plot
melted_f$variable <- factor(melted_f$variable, levels = c("IS1", "IS3", "IS4", "IS66", "IS21", "IS110", "IS91", "IS630"))

## Add comparisons
my_comparisons <- list( c("E. coli", "S. boydii"), 
                        c("E. coli", "S. flexneri"), 
                        c("E. coli", "S. sonnei") )

## Plot
ggplot(melted_f, aes(x=species, y=value)) +
  geom_boxplot(
    outlier.colour="red",
    outlier.fill="red",
    outlier.size=2,
    # custom boxes
    color="blue",
    fill="blue",
    alpha=0.2,
  ) +
  
  facet_wrap(~variable, scales = "free_y", ncol=4) +
  xlab("") +
  
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", ) +
  
  ylab("Number of IS elements") +
  theme_linedraw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    text = element_text(size = 14)
    )

## Save the plot
ggsave("<Name>.png", device = "png", dpi = 720, width = 14, height = 9)

