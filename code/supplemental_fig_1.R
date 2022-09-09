library(readxl)
library(tidyverse)
library(NatParksPalettes)
library(ggtext) 

df <- read_excel("data/MG_MT_results_summary.xlsx")


test <- cor.test(df$size.binner, df$genes)  # cor: 0.924, p < 2.2e-16
str(test)
test$p.value
test$estimate

df %>% 
  select(classification...16 ,size.binner, genes ) %>% 
  separate(classification...16, 
           into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), 
           sep = ";") %>% 
  select(Phylum, size.binner, genes) %>% 
  mutate(Phylum = case_when(str_detect(Phylum, "Firmicutes") ~ "Bacillota",
                            str_detect(Phylum, "Bacteroidota") ~ "Bacteroidota",
                            str_detect(Phylum, "Patescibacteria") ~ "Patescibacteria", 
                            str_detect(Phylum, "Proteobacteria") ~ "Pseudomonadota",
                            TRUE ~ "Unknown Phylum")) %>% 
  ggplot(aes(x = size.binner/1e6,y = genes, shape = Phylum, color = Phylum)) +
  geom_point(size = 4,
             alpha = 0.8) +
  scale_color_manual(values = natparks.pals("Arches2",7)) +
  scale_shape_manual(values = c(15, 16, 17, 18, 19, 0, 1)) +
  scale_y_continuous(labels = scales::comma) +
  theme_classic() +
  labs(y = "Gene Count",
       x = "Genome Size (Mb)",
       caption = paste0("R = ", round(test$estimate, 2), ", p < 0.001")) +

  theme(axis.title = element_text(face = "bold", 
                                  color = "black", 
                                  size = 14),
        axis.text = element_text(face = "bold", 
                                 color = "black", 
                                 size = 12),
        plot.caption = element_text(face = "bold", 
                                     color = "black", 
                                     size = 12),
        legend.position = c(0.15, 0.8),
        legend.title = element_text(face = "bold", 
                                    color = "black", 
                                    size = 14),
        legend.text = element_text(face = "bold", 
                                 color = "black", 
                                 size = 12)) 
ggsave("plots/sup_fig_1.png",
       width = 8, 
       height = 5,
       units = c("in"),
       bg = "white")
  

  