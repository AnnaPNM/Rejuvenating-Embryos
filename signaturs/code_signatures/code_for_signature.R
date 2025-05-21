library(VennDiagram)
library(eulerr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(plotly)
library(sunburstR)
library(treemap)
library(treemapify)

setwd('~/Desktop/BI/ageing/')
# code for Venn diagrams -- many thanks to Anna P. 

senes_data <- read.csv('senescence_data.csv')
prot_data <- read.csv('proteostasis_signatures.csv')
auto_data <- read.csv('autophagy_signatures.csv')
gi_data <- read.csv("~/Desktop/BI/ageing/genomic_signatures_new.csv")
senes_genes <- senes_data$Gene
prot_genes <- prot_data$Gene
auto_genes <- auto_data$Gene
gi_genes <- gi_data$Gene

myCol <- c('darkslategray3', 'pink3','green4', '#ffffe0')

venn.diagram(
  x = list(senes_genes, prot_genes, auto_genes, gi_genes),
  category.names = c("Senescence", "Proteostasis", "Autophagy", "Genomic Instability"),
  filename = 'venn_signatures_4sets.png',
  

  imagetype = "png",
  height = 6, 
  width = 6,
  units = 'in',
  dpi = 500,
  
  # Настройки кругов
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Настройки чисел
  cex = 1.5,
  fontface = "bold",
  fontfamily = "sans",
  
  # Настройки подписей
  cat.cex = 1.2,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135, -45),
  cat.dist = c(0.05, 0.05, 0.05, 0.05),
  cat.fontfamily = "sans",
  margin = 0.1
)





fit <- euler(list(
  Senescence = unique(senes_genes),
  Proteostasis = unique(prot_genes),
  Autophagy = unique(auto_genes),
  Genomic_instability =unique(gi_genes)
))

p <- plot(fit,
     fills = myCol,
     quantities = TRUE,
     labels = TRUE,
     # main = "Common and specific signatures",
     alpha = 0.6,
     edges = F,
     legend = list(cex = 1.2))
ggsave('euler_4set.png', p, width = 8, height = 8, dpi = 500)


###### here I tried to summarise in visualisation
all_hallmarks <- read.csv('all_hallmarks.csv')

count_data <- all_hallmarks %>%
  group_by(Subprocess, Hallmark) %>%
  summarise(Gene_Count = n(), .groups = 'drop')

ggplot(count_data, aes(x = Subprocess, y = Gene_Count, fill = Hallmark)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = Gene_Count), position = position_dodge(width = 0.9), vjust = -0.5) +
  labs(title = "Количество генов по Subprocess и Hallmark",
       x = "Биологический процесс",
       y = "Количество генов",
       fill = "Hallmark") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

count_data <- all_hallmarks %>% 
  count(Subprocess, Hallmark, name = "Gene_Count")

hierarchy_data <- all_hallmarks %>%
  mutate(path = paste(Hallmark, Subprocess, sep = "-")) %>%
  count(path)

sunburst(data = hierarchy_data, 
         colors = list(range = c("#ffffe0", "#08306b"))) #color map is hard to ajust

## tree map
hierarchy_stats <- all_hallmarks %>%
  group_by(Hallmark, Subprocess) %>%
  summarise(Gene_Count = n_distinct(Gene), .groups = "drop") %>%
  filter(Gene_Count > 0)

max_size <- 400 
hierarchy_stats$Gene_Count <- pmin(hierarchy_stats$Gene_Count, max_size)

hierarchy_data_clean <- hierarchy_stats %>%
  filter(Gene_Count > 0) %>%
  mutate(Hallmark = factor(Hallmark, levels = c("Senescence", "Proteostasis", "Autophagy", "Genomic instability")))

hallmark_colors <- c(
  "Senescence" = "darkslategray3",
  "Proteostasis" = "pink3",
  "Autophagy" = "green4",
  "Genomic instability" = "#ffffe0"
)

fig <- ggplot(hierarchy_data_clean, aes(
  area = Gene_Count,
  fill = Hallmark,
  label = Subprocess,
  subgroup = Hallmark
)) +
  geom_treemap(color = "black") +
  geom_treemap_subgroup_border(color = "white", size = 2) +
  geom_treemap_subgroup_text(place = "centre", grow = TRUE, alpha = 0.5, colour = "black", fontface = "bold") +
  geom_treemap_text(colour = "black", place = "topleft", reflow = TRUE, size = 10) +
  scale_fill_manual(values = hallmark_colors) +
  theme(legend.position = "none") 



ggsave("treemap_horizontal.png", plot = fig, 
       width = 12,  
       height = 8,  
       units = "in",  
       dpi = 300)

