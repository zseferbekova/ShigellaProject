## Load dependencies
library(ape)
library(reshape2)
library(tidyverse)
library(vroom)
library(ggpubr)

## Import total_stats.csv and leave only assembly names and strain names
stats <- as.data.frame(read.csv("../Data/total_stats.csv"), stringsAsFactors = F)
mapped_names <- stats %>%
  filter(Type == 'Chromosome') %>%
  select(genome = AssemblyID, species = Organism)
rownames(mapped_names) <- mapped_names$genome

## Import the phylogenetic tree
tree <- read.tree('../Data/nt_tree_midpoint.nwk')

## Get a distance matrix from the tree
dist_matrix <- as.data.frame(cophenetic(tree), stringsAsFactors = F)

## Transform and melt the table to plot later
dist_matrix <- cbind(rownames(dist_matrix), dist_matrix)
colnames(dist_matrix)[1] <- 'genome'
dist_df <- dist_matrix %>%
  melt(id.vars="genome") %>%
  drop_na() # drop NAs of the upper triangle
colnames(dist_df) <- c('genome1', 'genome2', 'distance')

## Import a synteny blocks matrix for a cluster
df <- as.data.frame(vroom("../Results/Distance_matrices/phylogroup_<X>.txt", col_names = F, delim = '\t'))

## Transform the first column into valid names and change colnames
n <- as.integer(df[1,])
df <- df[-1,]
df <- as.data.frame(str_split_fixed(df, pattern = "\t", n = n), stringsAsFactors = F)
df[,1] <- paste0('GCA_',str_sub(df[,1], 1, 9), '.', str_sub(df[,1], 10, 10))
colnames(df)[1] <- 'genome'
colnames(df)[2:ncol(df)] <- df[,1][1:ncol(df)-1]
df[df==""] <- NA

## Melt the table to plot later
df <- df %>%
  melt(id.vars="genome") %>%
  drop_na()
colnames(df) <- c('genome1', 'genome2', 'num_blocks')

## Make a df with blocks number & distances
joint_df <- left_join(df, dist_df, by=c('genome1', 'genome2'))
joint_df$num_blocks <- as.numeric(joint_df$num_blocks)
joint_df$distance <- as.numeric(joint_df$distance)

## Add species column
joint_df$species1 <- mapped_names[joint_df$genome1,]$species
joint_df$species2 <- mapped_names[joint_df$genome2,]$species

## Use to plot all Shigella spp. vs E. coli
joint_df$group1 <- word(joint_df$species1, 1, 1, fixed(" "))
joint_df$group2 <- word(joint_df$species2, 1, 1, fixed(" "))

joint_df <- joint_df %>%
  mutate(group = ifelse((group1 == 'Escherichia') & (group2 == 'Escherichia'), '2 E.coli',
                        ifelse((group1 == 'Shigella') & (group2 == 'Shigella'), '2 Shigella', 'Shigella & E.coli')))

## Use to plot Shigella spp. vs each other
#joint_df$group1 <- word(joint_df$species1, 1, 2, fixed(" "))
#joint_df$group2 <- word(joint_df$species2, 1, 2, fixed(" "))

# joint_df <- joint_df %>%
#   filter(!((group1 == "Escherichia coli") | (group2 == "Escherichia coli"))) %>%
#   mutate(group = ifelse((group1 == 'Shigella flexneri') & (group2 == 'Shigella flexneri'), '2 S. flexneri',
#                         ifelse((group1 == 'Shigella boydii') & (group2 == 'Shigella boydii'), '2 S. boydii',
#                                ifelse((group1 == 'Shigella sonnei') & (group2 == 'Shigella sonnei'), '2 S. sonnei',
#                                     ifelse(((group1 == 'Shigella sonnei') & (group2 == 'Shigella boydii')) | ((group1 == 'Shigella boydii') & (group2 == 'Shigella sonnei')), 'S. sonnei vs S. boydii',
#                                            ifelse(((group1 == 'Shigella sonnei') & (group2 == 'Shigella flexneri')) | ((group1 == 'Shigella flexneri') & (group2 == 'Shigella sonnei')), 'S. sonnei vs S. flexneri',
#                                                   ifelse(((group1 == 'Shigella boydii') & (group2 == 'Shigella flexneri')) | ((group1 == 'Shigella flexneri') & (group2 == 'Shigella boydii')), 'S. boydii vs S. flexneri', '')))))))

## Plot phylogenetic distance vs number of syntenic blocks
ggplot(data = joint_df, 
            aes(x = num_blocks, y = distance)) +
  geom_point(aes(color = group), size = 3, shape=1) +
  #geom_smooth(method = 'lm', se = F, color = "grey") +
  theme_classic() + theme(text = element_text(size=15)) +
  labs(x = "Number of blocks", 
  y = "Phylogenetic distance",
  color = "2 genomes belong to")

## Save the plot
ggsave(file="<Name>1.png", device="png", dpi = 480, width = 10, height = 5)


## Select intersection for '2 E.coli' and 'Shigella & E.coli'
max_distance = max(joint_df[joint_df$group == '2 E.coli', 'distance'])
min_distance = min(joint_df[joint_df$group == 'Shigella & E.coli', 'distance'])

## Plot phylogenetic distance vs number of syntenic blocks
ggplot() +
  geom_point(data = filter(joint_df, (distance < min_distance) | (distance > max_distance), group != '2 Shigella'), 
             aes(x = num_blocks, y = distance), 
             size = 2, color='grey') +
  
  geom_point(data = filter(joint_df, distance >= min_distance, distance <= max_distance, group != '2 Shigella'), 
             aes(x = num_blocks, y = distance, color = group), size = 2, shape=1) +

  theme_classic() + theme(text = element_text(size=15)) +
  
  labs(x = "Number of blocks", 
       y = "Phylogenetic distance",
       color = "2 genomes belong to") +
  
  geom_hline(yintercept = min_distance, linetype = 'dashed') +
  geom_hline(yintercept = max_distance, linetype = 'dashed')

## Add a column with blocks/distance ratio
filtered_joint_df <- joint_df %>%
  filter(distance >= min_distance, distance <=max_distance) %>%
  mutate(ratio = num_blocks / distance)

## Plot boxplots and compare the ratio
ggplot(filter(filtered_joint_df, group %in% c("Shigella & E.coli", "2 E.coli")), 
       aes(x=group, y=ratio, color = group, fill = group)) +
  
  geom_boxplot(
    outlier.colour="red",
    outlier.fill="red",
    outlier.size=2,
    # custom boxes
    #color="blue",
    #fill="blue",
    alpha=0.2,
  ) +
  
  xlab("") +
  
  stat_compare_means(comparisons = list(c("Shigella & E.coli", "2 E.coli")), method='wilcox.test', method.args =list(alternative = "greater"), aes(label = paste0("p = ", ..p.format..))) +
  
  ylab("Blocks/distance ratio") +
  theme_linedraw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    text = element_text(size = 14)
  )

## Save the plot
ggsave("<Name>.png", device = "png", dpi = 720, width = 14, height = 9)

