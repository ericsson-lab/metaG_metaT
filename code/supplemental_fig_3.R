library(tidyverse)
library(NatParksPalettes)
library(janitor)
library(patchwork)

merged_euks_abundance_all <- read_tsv(file = "data/eukaryotes/merged_eukulele_all_abundances.txt") %>% 
  clean_names()

taxon_long <- merged_euks_abundance_all %>% 
  pivot_longer(-full_classification,
               names_to = "sample",
               values_to = "counts") %>% 
  filter(str_detect(string = full_classification,
                    pattern = "Eukaryota;")) %>% 
  mutate(gm = case_when(str_detect(string = sample,
                                   pattern = "gm1") ~ "GM1",
                        str_detect(string = sample,
                                   pattern = "gm4") ~ "GM4")) %>% 
  separate(full_classification,
           into = c("Domain", "Supergroup", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";",
           fill = "warn") %>% 
  pivot_longer(-c(sample, counts, gm),
               names_to = "level",
               values_to = "taxon") %>% 
  mutate(taxon = str_replace_all(taxon, " ", ""))


rel_abund_table <- taxon_long %>% 
  drop_na() %>% 
  group_by(sample) %>% 
  mutate(rel_abund = counts/sum(counts))

phylum_table <- rel_abund_table %>% 
  filter(level == "Phylum") %>% 
  group_by(gm, sample, taxon) %>% 
  summarize(rel_abund = sum(rel_abund))

phylum_table %>% 
  group_by(gm, taxon) %>% 
  summarize(mean = mean(rel_abund)) %>% 
  arrange(desc(mean))

taxon_pool <- phylum_table %>% 
  group_by(taxon) %>% 
  summarise(pool = max(rel_abund) < 0.01, .groups = 'drop')


inner_join(phylum_table, taxon_pool, by = "taxon") %>% 
  mutate(taxon = if_else(pool, "Other", taxon)) %>% 
  group_by(sample, taxon, gm) %>% 
  summarise(rel_abund = sum(rel_abund)) %>% 
  ggplot(aes(x=sample, y = factor(taxon, level = rev(c("Ochrophyta",
                                                   "Haptophyta", 
                                                   "Dinoflagellata",
                                                   "Cryptophyta", 
                                                   "Ciliophora",
                                                   "Chlorophyta", 
                                                   "Other")))
             , 
             fill = rel_abund)) +
  geom_tile()  +
  scale_fill_gradientn(colors=natparks.pals("Arches2"),
                       trans = "reverse",
                       na.value = "#F3DAE4", 
                       name = "Relative\nAbundance",
                       labels = scales::percent, 
                       # limits = c(0, 10),
                       breaks = c(0, 0.025, 0.05, 0.075)) +
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



ggsave("plots/euk_rel_abund_phylum.png",
       width = 4,
       height = 3,
       units = c("in"),
       bg = "white")

