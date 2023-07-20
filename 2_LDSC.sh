#https://bios25328.hakyimlab.org/post/2021/04/16/ld-score-regression/
#https://github.com/bulik/ldsc/wiki/Summary-Statistics-File-Format
#https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation
#mamut:

##PART ONE: GENERATE sumstats FILE FOR GWAS
#https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation

##### use ldsc v2 (py3) for munge_sumstats, use original version for LDSC. ####
# #MAMUT
# storage 
# cd ldsc_py3
# conda activate ldsc3

#IMPERIAL HPC
module load anaconda3/personal
source activate ldsc3
cd /rds/general/user/eosimo/home/lenhard_prs/ldsc_py3


zcat gwas/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz | tail -n +74 |\
    gzip > gwas/PGC3_SCZ_wave3.european.autosome.public.v3.tsv.gz

munge_sumstats.py \
    --sumstats gwas/PGC3_SCZ_wave3.european.autosome.public.v3.tsv \
    --snp ID --N-cas-col NCAS --N-con-col NCON --p PVAL \
    --merge-alleles data/w_hm3.snplist \
    --out gwas/formatted_hapmap/PGC3_SCZ_wave3_european_autosome_v3_py3

munge_sumstats.py \
    --sumstats gwas/hcm.gwama.sumstats_hg19_24Feb21.tsv \
    --signed-sumstats Z,0\
    --merge-alleles data/w_hm3.snplist \
    --out gwas/formatted_hapmap/HCM_Sean_hg19_24Feb21_european_autosome_v3_py3


# ONLY ONCE
#update frq files
cd data/1000G_EUR_Phase3_plink
for CHR in {1..22}
do
  plink --bfile 1000G.EUR.QC.$CHR --freq --out 1000G.EUR.QC.$CHR
done

# ONLY ONCE for modified baseline
#create baseline files matched to modern BIM and my annot files
for(i in c(1:22)){
  print(i)
  (baseline_annot<-data.table::fread(paste0("data/baseline/baseline.",i,".annot.gz")) %>% 
      dplyr::select(-CM)    )
  (bim<- data.table::fread(paste0("data/1000G_EUR_Phase3_plink/1000G.EUR.QC.",i,".bim"),
                           header = F,col.names = c("CHR","SNP" ,"CM", "BP" ,"A1","A2")) %>% 
      dplyr::select(-A1,-A2))
  
  (EPs_annot<-data.table::fread(paste0("data/EPs/ALL_BRAIN_EPs.",i,".annot.gz")))
  (EPs_annot <- cbind( bim[,c(1,4,2,3)], EPs_annot) %>% rename(EPs=ANNOT))
  
  (merge<-left_join(EPs_annot,baseline_annot, by=c("CHR","SNP" ,"BP")))
  sample_n(merge,10)
  
  # merge_match_bim<-merge[match(bim$SNP, merge$SNP),] 
  merge$base<-1
  merge[is.na(merge)] <- 0
  head(merge[,c(1:9)])
  
  nrow(bim)==nrow(merge)
  
  (base_merged <- merge %>%  select(-EPs))
  
  data.table::fwrite(
    base_merged, 
    file = paste0("data/baseline_Osimo/baseline.",i,".annot.gz"), sep = "\t"
  )
}
for CHR in {1..22}
do
    echo $CHR
    python ldsc.py\
      --l2\
      --bfile data/1000G_EUR_Phase3_plink/1000G.EUR.QC.$CHR\
      --ld-wind-cm 1\
      --annot data/baseline_Osimo/baseline.$CHR.annot.gz\
      --out data/baseline_Osimo/baseline.$CHR
done
#		--print-snps data/hapmap3_snps/hm.$CHR.snp\


# ONLY ONCE unpartitioned weights
for CHR in {1..22}
do
    echo $CHR
    python ldsc.py\
      --l2\
      --bfile data/1000G_EUR_Phase3_plink/1000G.EUR.QC.$CHR\
      --ld-wind-cm 1\
      --out data/weights_Osimo/weights.$CHR
done
#		--print-snps data/hapmap3_snps/hm.$CHR.snp\




#PART 2: CALCULATE Partitioned LD Scores
#https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial


# FILTER bed file WITH unique ENH
cd data/EPs/
for bedfile in *.bed
do
  echo $bedfile
  cat $bedfile | awk '!seen[$1,$2,$3]++' > uniq_$bedfile
done
# cat data/EPs/ALL_BRAIN_EPs.bed | awk -F'\t' 'BEGIN {OFS = FS} { gsub(/chr/,"", $1); print } ' | awk '!seen[$1,$2,$3]++' > data/EPs/ALL_BRAIN_EPs_nochr.bed



#Step 1: Creating an annot file for EPs
# and #Computing LD scores

## MAMUT
# storage 
# cd ldsc_py2.7
# conda activate ldsc

#IMPERIAL HPC

#!/bin/bash
#PBS -lselect=1:ncpus=2:mem=50gb
#PBS -lwalltime=2:0:0
#PBS -N ldsc

module load anaconda3/personal
source activate ldsc
cd /rds/general/user/eosimo/home/lenhard_prs/ldsc_py2.7

# cd data/EPs/
for CHR in {9..9}
do
  # for bedfile in uniq_*.bed
  # do
  # echo $bedfile
    echo $CHR
    #Creating an annot file for each chr
    python make_annot.py \
		  --bed-file data/EPs/2023-07-03_NON_Cardiac_significant.bed \
		  --bimfile data/1000G_EUR_Phase3_plink/1000G.EUR.QC.$CHR.bim \
		  --annot-file data/EPs/annot_files/NON_CARDIAC_21k_significant_EPs_Jul23.$CHR.annot.gz
  # done
    #Computing LD scores for each chr
    python ldsc.py\
      --l2\
      --thin-annot\
      --bfile data/1000G_EUR_Phase3_plink/1000G.EUR.QC.$CHR\
      --ld-wind-cm 1 \
      --annot data/EPs/annot_files/NON_CARDIAC_21k_significant_EPs_Jul23.$CHR.annot.gz\
      --out data/EPs/annot_files/NON_CARDIAC_21k_significant_EPs_Jul23.$CHR
done


# Step 2: Computing LD scores with an annot file.

#this is not needed if using thin-annot later:
##add cordinates from BIM file to annot files:
# library(tidyverse)
# for(i in c(1:22)){
#   # R merge all annot files
#   (temp = list.files(pattern=paste0("*.bed.",i,".annot.gz"), 
#                      path = "/mnt/storage/emanuele/ldsc_py2.7/data/EPs/annot_files/"))
  
#   for (j in 1:length(temp)) {
#     temp[j]
#     (col.name = str_split(string = temp[j], pattern = "\\.")[[1]][1])
#     # (x= assign(temp[j], read.csv(temp[j])))
#     (bim<- data.table::fread(paste0("/mnt/storage/emanuele/ldsc_py2.7/data/1000G_EUR_Phase3_plink/1000G.EUR.QC.",i,".bim"),
#                              header = F,col.names = c("CHR","SNP" ,"CM", "BP" ,"A1","A2")) %>% 
#        dplyr::select(-A1,-A2))
    
#     (EPs_annot<-data.table::fread(paste0("/mnt/storage/emanuele/ldsc_py2.7/data/EPs/annot_files/",temp[j])))
    
#     (EPs_annot2 <- cbind( bim[,c(1,4,2,3)], EPs_annot))
#     names(EPs_annot2)[names(EPs_annot2) == "ANNOT"] <- col.name
#     head(EPs_annot2)
    
#     data.table::fwrite(
#       EPs_annot2, 
#       file = paste0("/mnt/storage/emanuele/ldsc_py2.7/data/EPs/annot_files/",temp[j]), sep = "\t"
#     )
#   }
# }




# https://github.com/bulik/ldsc/wiki/Partitioned-Heritability
## 3- The following command will allow you to partition heritability:

#!/bin/bash
#PBS -lselect=1:ncpus=2:mem=50gb
#PBS -lwalltime=2:0:0
#PBS -N ldsc

module load anaconda3/personal
source activate ldsc
cd /rds/general/user/eosimo/home/lenhard_prs/ldsc_py2.7

python ldsc.py \
	--h2 gwas/formatted_hapmap/PGC3_SCZ_wave3_european_autosome_v3_py3.sumstats.gz \
	--overlap-annot \
  --w-ld-chr data/weights_Osimo/weights. \
	--frqfile-chr data/1000G_EUR_Phase3_plink/1000G.EUR.QC. \
  --ref-ld-chr data/baseline_Osimo/baseline.,\
data/EPs/annot_files/NEURAL_8k_GRB_Enhancers.,\
data/EPs/annot_files/NEURAL_14k_noGRB_significant_EPs_Feb23.,\
data/EPs/annot_files/21k_NEURAL_ENH_EXP_significant_ES_significant_contact_EPs.,\
data/EPs/annot_files/notNeural_20k_100flank.,\
data/EPs/annot_files/NON_NEURAL_8k_significant_EPs_Jul23.,\
data/EPs/annot_files/PsychENCODE_DER_03b_PFC_enhancers_18k_100flank.,\
data/EPs/annot_files/BRAIN_EP_eQTL_Jan23.,\
data/EPs/annot_files/Radina_GRBs_hg19_mm10.98.50.,\
data/EPs/annot_files/all_FANTOM5_hg19_enhancers. \
	--out stratified_LDSC_out/2023_07_17_NEURAL_GRBorNot_100flank_fixnotneural

#data/EPs/annot_files/2022-11-08_bed_34k_neg_enhs. \ #not working on 18 Jul 23 due to corrupt file in storage



python ldsc.py \
	--h2 gwas/formatted_hapmap/HCM_Sean_hg19_24Feb21_european_autosome_v3_py3.sumstats.gz \
	--w-ld-chr data/weights_Osimo/weights. \
	--overlap-annot \
	--frqfile-chr data/1000G_EUR_Phase3_plink/1000G.EUR.QC. \
  --ref-ld-chr data/baseline_Osimo/baseline.,\
data/EPs/annot_files/9k_CARDIAC_noFibro_Enhancers.,\
data/EPs/annot_files/CARDIAC_NoFibro_3k_GRB_Enhancers.,\
data/EPs/annot_files/CARDIAC_NoFibro_6k_noGRB_Enhancers.,\
data/EPs/annot_files/NON_CARDIAC_21k_significant_EPs_Jul23.,\
data/EPs/annot_files/Radina_GRBs_hg19_mm10.98.50.,\
data/EPs/annot_files/all_FANTOM5_hg19_enhancers.,\
data/EPs/annot_files/uniq_HEART_EP_eQTL_Nov22. \
	--out stratified_LDSC_out/2023_07_19_CARDIAC_noFibro_GRBorNot_100flank_fixnotcardiac

#data/EPs/annot_files/40k_notCARDIAC_Enhancers.,\
# data/EPs/annot_files/2022-11-08_bed_34k_neg_enhs. \

