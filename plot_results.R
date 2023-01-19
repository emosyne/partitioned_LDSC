library(tidyverse)

(results <- data.table::fread("~/GoogleDrive/WORK/CF_PhD/NF_2HH/partitioned_LDSC/2023_01_19_NEURAL_100flank.results") %>% 
    pivot_longer(cols=!Category) %>% replace_na(list(value=1)))
# table(results$Category)

theme_set(theme_minimal())

pdf(file = paste0("/Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/partitioned_LDSC/", Sys.Date(), 
                  "_NEURAL_allEnh_100flank.pdf"), 
    width = 30, height = 10, pointsize = 24)
ggplot(aes(y=value, x = forcats::fct_inorder(Category), fill = Category), data= results)+ #[results$name=="Enrichment",]
  facet_grid(cols = vars(name), scales="free")+
  geom_bar(stat = "identity")+
  geom_text(aes(label = round(value,2), x = forcats::fct_inorder(Category), y = value), position = position_dodge(width = 0.2), hjust = -0.5, size=2)+
  # scale_fill_brewer(palette = "Set1") +
  coord_flip()+ theme(legend.position = "none")
dev.off()
