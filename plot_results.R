library(tidyverse)

(results <- data.table::fread("~/Google Drive/WORK/CF_PhD/NF_2HH/partitioned_LDSC/2022_11_10_baseline_allFANTOM_nonBRAINexp_anyFANTOMneuralEXP_PsychENCODEhighCONFforPFC_FANTOMsignESsignContactNeural_allBRAIN.tsv") %>% 
    pivot_longer(cols=!Category) %>% replace_na(list(value=1)))
# (results <- read_excel("~/Google Drive/WORK/CF_PhD/NF_2HH/partitioned_LDSC/results 5 Nov 22.xlsx") %>% 
#     as_tibble(.name_repair = "universal") %>% pivot_longer(cols=!Categories.) %>% 
#     mutate(Categories. = str_replace_all(Categories., "\\s+|:", "")))
table(results$Category)

theme_set(theme_minimal())

pdf(file = "/Users/eosimo/Google Drive/WORK/CF_PhD/NF_2HH/partitioned_LDSC/2022_11_10_baseline_allFANTOM_nonBRAINexp_anyFANTOMneuralEXP_PsychENCODEhighCONFforPFC_FANTOMsignESsignContactNeural_allBRAIN.pdf", 
    width = 30, height = 10, pointsize = 24)
ggplot(aes(y=value, x = forcats::fct_inorder(Category), fill = Category), data= results)+ #[results$name=="Enrichment",]
  facet_grid(cols = vars(name), scales="free")+
  geom_bar(stat = "identity")+
  geom_text(aes(label = round(value,3), x = forcats::fct_inorder(Category), y = value), position = position_dodge(width = 0.8), hjust = -0.5, size=2)+
  # scale_fill_brewer(palette = "Set1") +
  coord_flip()+ theme(legend.position = "none")
dev.off()
