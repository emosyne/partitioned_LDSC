library(tidyverse)
library(gridExtra)
library(grid)

setwd("/Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/partitioned_LDSC")

underspace <- function(x,...){
  gsub('_',' ',x)
}


# ##SCZ
# (results <- data.table::fread("/Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/partitioned_LDSC/results/partitioned_LDSC/2023-07-19_NEURAL_fix_non-neural_post-viva.results") %>%
#     replace_na(list(Enrichment_p=1)) )#%>%
#     #pivot_longer(cols=!Category) )
# # table(results$Category)
# results$Category <- sub("L2_0", "", results$Category)
# 
# (results <- results %>% select( c(Partition = Category,Enrichment,Enrichment_p,Prop_SNPs = Prop._SNPs,Prop_heritability = Prop._h2 )) %>%
#     filter(grepl("base|Coding|Intron|Promoter|nhancer|UTR|GRB|Non-associated|FANTOM|PsychENCODE|Brain|eQTLs", Partition)) %>%
#     filter(!grepl("extend", Partition)) %>%
#     mutate(Partition = factor(Partition, ordered = T),
#            Enrichment = Enrichment-1,
#            Enrichment_p = p.adjust(method = "BH", Enrichment_p, n = nrow(results)))
#     )
# condition = "SCZ"
# 
# (plot = results[c(1:7,9:12),])
# plotname = "main"
# (plot = results[c(3,8,13:20),])
# plotname = "enhancer-based"



##CARDIAC
(results <- data.table::fread("/Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/partitioned_LDSC/results/partitioned_LDSC/2023-07-19_CARDIAC_fix_non-cardiac_post-viva.results") %>%
    replace_na(list(Enrichment_p=1)) )#%>%
#pivot_longer(cols=!Category) )
# table(results$Category)
results$Category <- sub("L2_0", "", results$Category)

(results <- results %>% select( c(Partition = Category,Enrichment,Enrichment_p,Prop_SNPs = Prop._SNPs,Prop_heritability = Prop._h2 )) %>%
    filter(grepl("base|Coding|Intron|Promoter|nhancer|UTR|GRB|Non-associated|FANTOM|PsychENCODE|Brain|eQTLs", Partition)) %>%
    filter(!grepl("extend", Partition)) %>%
    mutate(Partition = factor(Partition, ordered = T),
           Enrichment = Enrichment-1,
           Enrichment_p = p.adjust(method = "BH", Enrichment_p, n = nrow(results)))
)
condition = "HCM"
(plot = results[c(1:11),])
plotname = "main"
(plot = results[c(3,12:19),])
plotname = "enhancer-based"






theme_set(theme_minimal())
colours = values=MetBrewer::met.brewer("Johnson", 2)
(main = ggplot(aes(y = Enrichment, 
                   x = forcats::fct_rev(forcats::fct_inorder(underspace(Partition))), 
                   label = paste0("p = ",round(Enrichment_p,2), #"Enrichment = ",round(Enrichment,2),
                                  ifelse(test = Enrichment_p<0.05, yes = "*", no = ""),
                                  ifelse(test = Enrichment_p<0.01, yes = "**", no = "")),
                   fill = (Enrichment_p)), 
               data= plot)+
    geom_bar(stat = "identity")+
    scale_fill_gradientn(colours = c("red","gray","gray"), 
                         values = scales::rescale(c(0,0.06,1)), breaks = c(0,1),
                         guide = "colorbar", limits=c(0,1),name = "BH-adjusted p value",)+    
    scale_y_continuous(expand = expansion(mult = c(0, .2))) + 
    ggrepel::geom_text_repel(alpha=0.5, show.legend = F, min.segment.length = 0, hjust = -0.2, vjust = 0)+
    coord_flip()+ 
    theme(
      legend.position = "bottom",
      # legend.key.width=unit(4,"cm"), 
      text = element_text(size = 18),
      plot.margin = margin(t = 5.5, r = 0, b = 5.5, l = 5.5, "pt")
      # panel.grid.major.y = element_blank()
    )+
    xlab("Genomic partition") + ylab("Enrichment in GWAS signal") )

(main_grob = arrangeGrob(textGrob("A)", just = "left",
                                  gp = gpar(fontsize = 18, fontface = "bold", col="black")), 
                         textGrob(paste0("LDSC-based genomic partitions' enrichment towards ",condition, " GWAS"), 
                                  gp = gpar(fontsize = 18, fontface = "bold", col="navyblue")), 
                         main, 
                         layout_matrix=rbind(c(1,2),
                                             c(3,3)),
                         widths = c(0.1, 1), heights = c(0.1, 1)))




(long_res = plot %>% dplyr::select(Partition, Prop_SNPs, Prop_heritability) %>% pivot_longer(!Partition))
(side = ggplot(aes(x = value, 
                   y = forcats::fct_rev(forcats::fct_inorder(underspace(Partition))), 
                   group = underspace(name),
                   fill = underspace(name)),
               data= long_res)+
    #facet_grid(cols = vars(name), scales="free")+
    geom_bar(position="dodge", stat="identity")+
    scale_fill_manual(values=MetBrewer::met.brewer("Johnson", 2),
                      name = "")+
    theme(
      legend.position = "bottom", 
      text = element_text(size = 18),
      axis.text.y = element_blank(),
      plot.margin = margin(t = 5.5, r = 5.5, b = 47.5, l = 0, "pt"),
      # panel.grid.major.y = element_blank()
    )+
    xlab(NULL) + ylab(NULL) )
(side_grob = arrangeGrob(textGrob("B)", just = "left",
                                  gp = gpar(fontsize = 18, fontface = "bold", col="black")), 
                         textGrob(paste0("Proportion of SNPs in Partition (blue), and\nproportion of heritability for ", condition), 
                                  gp = gpar(fontsize = 12, fontface = "bold", col="navyblue", lineheight=0.8)), 
                         side, 
                         layout_matrix=rbind(c(1,2),
                                             c(3,3)),
                         widths = c(0.1, 1), heights = c(0.1, 1)))

tot_grob = arrangeGrob(
  textGrob(paste0("Proportion of heritability for ", condition, " for the ",plotname," genomic partitions"), 
           gp = gpar(fontsize = 22, fontface = "bold", col="navyblue")), 
  main_grob, side_grob, 
  layout_matrix=rbind(c(1,1),
                      c(2,3)),
  widths = c(1, 0.4), heights = c(0.1,1)
)
ggsave(
  filename = paste0("results/partitioned_LDSC/",Sys.Date(),"_",condition,"_",plotname,"_pcorr.pdf"), 
  tot_grob,  
  width = 10, height = 5, device = "pdf", scale = 1.5)

