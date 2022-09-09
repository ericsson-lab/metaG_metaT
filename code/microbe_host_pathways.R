library(readxl)
library(tidyverse)
library(RColorBrewer)
library(plotly)
library(patchwork)
library(EnhancedVolcano)
library(cowplot)

df_25 <- read_csv("./data/kegg/KEGGREST_metaT_significant_pathways_top25.csv")

pathways_for_supp <- read_delim("data/select_sig_to_omit_pathways.txt",
                                delim = "\n")


process_data <- function(df){
  df %>% 
    mutate(pathways = gsub("'", "", pathways)) %>% 
    pivot_longer(-c(kegg, pathways),
                 names_to = "sample",
                 values_to = "count") %>% 
    filter(count != "NA") %>% 
    group_by(pathways, sample) %>% 
    summarize(count = sum(as.numeric(count))) %>% 
    mutate(gm = case_when(str_detect(string = sample,
                                     pattern = "GM1") ~ "GM1",
                          str_detect(string = sample,
                                     pattern = "GM4") ~ "GM4")) %>% 
    group_by(gm, pathways) 
}

data_top25 <- process_data(df_25)

length(unique(data_top25$pathways)) # 92 pathways

a <- data_top25 %>% 
  filter(!pathways %in% pathways_for_supp$pathway) %>% 
  ggplot( aes(x = sample, 
              y = pathways,
              fill = log10(count))) +
  geom_tile() +
  scale_fill_gradientn(colors=natparks.pals("Arches2"),
                       trans = "reverse",
                       name = "Log10\nCount", 
                       na.value = "#F3DAE4") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_bw() +
  
  theme(axis.text = element_text(face = "bold",
                                 color = "black"),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.spacing = unit(0.001, "lines"),
        strip.background  = element_rect(color = "black",
                                         fill = "#b6baba"),
        strip.text = element_text(face = "bold",
                                  color = "black", 
                                  size = 12),
        legend.title = element_text(face = "bold",
                                    color = "black",
                                    size = 10),
        legend.text = element_text(face = "bold",
                                   color = "black"),
        legend.position = "right")  +
  facet_wrap(~gm,
             scales = "free_x") +
  guides(fill = guide_colorbar(reverse = T,
                               ticks.colour = NA))


a
# ggsave("plots/select_pathways_heatmap.png",
#        heigh = 8,
#        width = 6,
#        units = c("in"),
#        bg = "white",
#        dpi = 600)


### Supplemental Figure 2
b <- data_top25 %>% 
  filter(pathways %in% pathways_for_supp$pathway) %>% 
  ggplot( aes(x = sample, 
              y = pathways,
              fill = log10(count))) +
  geom_tile() +
  scale_fill_gradientn(colors=natparks.pals("Arches2"),
                       trans = "reverse",
                       name = "Log10\nCount", 
                       na.value = "#F3DAE4") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_bw() +
  
  theme(axis.text = element_text(face = "bold",
                                 color = "black"),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.spacing = unit(0.001, "lines"),
        strip.background  = element_rect(color = "black",
                                         fill = "#b6baba"),
        strip.text = element_text(face = "bold",
                                  color = "black", 
                                  size = 12),
        legend.title = element_text(face = "bold",
                                    color = "black",
                                    size = 10),
        legend.text = element_text(face = "bold",
                                   color = "black"),
        legend.position = "right")  +
  facet_wrap(~gm,
             scales = "free_x") +
  guides(fill = guide_colorbar(reverse = T,
                               ticks.colour = NA))

b
# ggsave("plots/supplemental_pathways_heatmap.png",
#        heigh = 8,
#        width = 6,
#        units = c("in"),
#        bg = "white",
#        dpi = 600)  


a + b
# ggsave("plots/pathways_heatmap_combined.png",
#        heigh = 8,
#        width = 12,
#        units = c("in"),
#        bg = "white",
#        dpi = 600)  



### panel B kegg

kegg <- read_csv("data/volcano_plots/kegg_counts_metaT_DESeq2.csv")

kegg_volcano <- EnhancedVolcano(kegg,
                lab = kegg$kegg,
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-30, 30),
                ylim = c(-1,25),
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                pCutoff = 10e-3,
                FCcutoff = 2,
                colAlpha = 1, 
                labSize = 3) +  
  theme_set(theme_bw() +
              theme(strip.background = element_rect(colour = "black", fill = "white"), 
                    strip.text = element_text(size=16, family="Helvetica", face = "bold", colour = "black")) +
              theme(axis.text.x = element_text(size = 16, family = "Helvetica", face = "bold", colour = "black"),
                    axis.text.y = element_text(size = 16, family = "Helvetica", face = "bold", colour = "black"),
                    axis.title.y = element_text(size = 16, face = "bold", family = "Helvetica", colour = "black"),
                    axis.title.x=element_text(size = 16, face = "bold", family = "Helvetica", colour = "black"),
                    axis.ticks.x=element_blank()) +
              theme(legend.position = "none"))


# Pfam

pfam <- read_csv("data/volcano_plots/Pfam_counts_metaT_DESeq2.csv")

pfam_volcano <- EnhancedVolcano(pfam,
                lab = pfam$...1,
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-30, 30),
                ylim = c(-1,25),
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                pCutoff = 10e-3,
                FCcutoff = 2,
                colAlpha = 1, 
                labSize = 3) +  
  theme_set(theme_bw() +
              theme(strip.background = element_rect(colour = "black", fill = "white"), 
                    strip.text = element_text(size=16, family="Helvetica", face = "bold", colour = "black")) +
              theme(axis.text.x = element_text(size = 16, family = "Helvetica", face = "bold", colour = "black"),
                    axis.text.y = element_text(size = 16, family = "Helvetica", face = "bold", colour = "black"),
                    axis.title.y = element_text(size = 16, face = "bold", family = "Helvetica", colour = "black"),
                    axis.title.x=element_text(size = 16, face = "bold", family = "Helvetica", colour = "black"),
                    axis.ticks.x=element_blank()) +
              theme(legend.position = "none"))

cazy <- read_csv("data/volcano_plots/CAZymes_counts_metaT_DESeq2.csv")

cazy_volcano <- EnhancedVolcano(cazy,
                                lab = cazy$...1,
                                x = 'log2FoldChange',
                                y = 'pvalue',
                                xlim = c(-30, 30),
                                ylim = c(-1,25),
                                title = NULL,
                                subtitle = NULL,
                                caption = NULL,
                                pCutoff = 10e-3,
                                FCcutoff = 2,
                                colAlpha = 1, 
                                labSize = 3) +  
  theme_set(theme_bw() +
              theme(strip.background = element_rect(colour = "black", fill = "white"), 
                    strip.text = element_text(size=16, family="Helvetica", face = "bold", colour = "black")) +
              theme(axis.text.x = element_text(size = 16, family = "Helvetica", face = "bold", colour = "black"),
                    axis.text.y = element_text(size = 16, family = "Helvetica", face = "bold", colour = "black"),
                    axis.title.y = element_text(size = 16, face = "bold", family = "Helvetica", colour = "black"),
                    axis.title.x=element_text(size = 16, face = "bold", family = "Helvetica", colour = "black"),
                    axis.ticks.x=element_blank()) +
              theme(legend.position = "none"))
cazy_volcano


## Figure 5

a + (kegg_volcano/pfam_volcano/cazy_volcano)

volcano_plots<- plot_grid(kegg_volcano,
                          pfam_volcano,
                          cazy_volcano,
                          nrow = 3,
                          ncol = 1)
plot_grid(a,
          volcano_plots,
          ncol = 2,
          rel_widths = c(0.7, 0.3))

ggsave('plots/figure5.png',
       dpi = 600,
       units = c("in"),
       width = 9,
       height = 8,
       bg = "white")
