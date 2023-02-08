library(tidyverse)

(results <- data.table::fread("~/GoogleDrive/WORK/CF_PhD/NF_2HH/partitioned_LDSC/2023_01_31_NEURAL_GRBorNot_100flank.results") %>% 
    replace_na(list(Enrichment_p=1)) )#%>% 
    #pivot_longer(cols=!Category) )
# table(results$Category)
results$Category <- sub("L2_0", "", results$Category)

(results <- results %>% select( c("Category","Enrichment","Enrichment_p")) %>% 
    filter(grepl("base|Coding|Intron|Promoter|Enhancer|UTR|GRB|Negative_34k|FANTOM|PsychENCODE|Brain", Category)) %>% 
    filter(!grepl("extend", Category))
    )

theme_set(theme_minimal())

pdf(file = paste0("/Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/partitioned_LDSC/", Sys.Date(), 
                  "_NEURAL_GRBorNot_100flank.pdf"), 
    width = 15, height = 10)
ggplot(aes(y=Enrichment, x = forcats::fct_rev(forcats::fct_inorder(Category)), fill = Enrichment_p), data= results)+
  #facet_grid(cols = vars(name), scales="free")+
  geom_bar(stat = "identity")+
  scale_fill_gradient(
    low = "navy",
    high = "yellow",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill"
  )+
  geom_text(show.legend=FALSE, aes(label = paste0("Enrichment = ",round(Enrichment,2),", p value = ",round(Enrichment_p,2), 
                               ifelse(test = Enrichment_p<0.05, yes = "*", no = ""),
                               ifelse(test = Enrichment_p<0.01, yes = "**", no = "")
  ), 
  x = forcats::fct_rev(forcats::fct_inorder(Category)), 
  y = 10, colour="red"), 
  #position = position_dodge(width = 1), 
  hjust = -0.5, size=4)+
  # scale_fill_brewer(palette = "Set1") +
  coord_flip()+ 
  theme(legend.position = "bottom", legend.key.width=unit(4,"cm"), text = element_text(size = 18))+
  xlab("LDSC-base Genomic partitions' enrichment towards schizophrenia GWAS") + ylab(NULL) 
dev.off()


##CARDIAC
(results <- data.table::fread("~/GoogleDrive/WORK/CF_PhD/NF_2HH/partitioned_LDSC/2023_01_31_CARDIAC_noFibro_GRBorNot_100flank.results") %>% 
    replace_na(list(Enrichment_p=1)) )#%>% 
#pivot_longer(cols=!Category) )
# table(results$Category)
results$Category <- sub("L2_0", "", results$Category)

(results <- results %>% select( c("Category","Enrichment","Enrichment_p")) %>% 
    filter(grepl("base|Coding|Intron|Promoter|Enhancer|UTR|GRB|FANTOM|Negative|Cardiac|Heart|eQTL", Category)) %>% 
    filter(!grepl("extend", Category))
)

theme_set(theme_minimal())

pdf(file = paste0("/Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/partitioned_LDSC/", Sys.Date(), 
                  "_CARDIAC_noFibro_GRBorNot_100flank.pdf"), 
    width = 15, height = 10)
ggplot(aes(y=Enrichment, x = forcats::fct_rev(forcats::fct_inorder(Category)), fill = Enrichment_p), data= results)+
  #facet_grid(cols = vars(name), scales="free")+
  geom_bar(stat = "identity")+
  scale_fill_gradient(
    low = "navy",
    high = "yellow",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill"
  )+
  geom_text(show.legend=FALSE, aes(label = paste0("Enrichment = ",round(Enrichment,2),", p value = ",round(Enrichment_p,2), 
                               ifelse(test = Enrichment_p<0.05, yes = "*", no = ""),
                               ifelse(test = Enrichment_p<0.01, yes = "**", no = "")
                               ), 
                x = forcats::fct_rev(forcats::fct_inorder(Category)), 
                y = 10, colour="red", ), 
            #position = position_dodge(width = 1), 
            hjust = -0.5, size=4)+
  # scale_fill_brewer(palette = "Set1") +
  coord_flip()+ 
  theme(legend.position = "bottom", legend.key.width=unit(4,"cm"), text = element_text(size = 18))+
  xlab("LDSC-base Genomic partitions' enrichment towards HCM GWAS") + ylab(NULL) 
dev.off()

