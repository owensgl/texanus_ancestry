library(tidyverse)
library(cowplot)
library(PNWColors)

chr_lengths <- read_tsv("../wild_gwas_2018/Ha412HO.chrlengths.txt") %>%
  mutate(chr=gsub("Ha412HOChr","",chr))


dstat <- read_tsv("Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.ANN_CEN.ANN_STEX.DEB.Dstats.outliers.txt") %>%
  mutate(type="D",chr=gsub("Ha412HOChr","",chr)) %>%
  select(chr,start,end,type)
dstat %>% mutate(length = end-start) %>% summarise(sum_lengt=sum(length))
pcadmix <- read_tsv("Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.beagle.deb.pcadmixoutliers.txt") %>%
  rename(chr=Chr,start=start_bp,end=end_bp) %>%
  mutate(type="PCAdmix") %>% select(-deb_percent)

pcadmix %>% mutate(length = end-start) %>% summarise(sum_lengt=sum(length))

pdf("figures/Texanus_all_D_PCAdmix_dstatoutliers.pdf",height=6,width=6)
ggplot(chr_lengths) + geom_segment(aes(x=start/1000000,y=as.numeric(chr),xend=end/1000000,yend=as.numeric(chr))) +
  geom_point(data=rbind(dstat,pcadmix) %>%
               mutate(middle =   ((end - start)/2)+start,
                      chr=case_when(type=="D" ~ as.numeric(chr)-0.1,
                                    type=="PCAdmix" ~as.numeric(chr)+0.1)),
             aes(x=middle/1000000,y=as.numeric(chr),color=type)) +
  scale_y_continuous(breaks=seq(17)) +
  ylab("Chromosome") + xlab("Mbp") +
  scale_color_manual(values=c("#1d457f","#f57946"),name="Method") +
  ggtitle("Introgression outliers") +
  annotate("segment",x=159493853/1000000,xend=159493853/1000000,y=12.8,yend=13.2,color="#edd746",size=1) +
  annotate("segment",x=27209191/1000000,xend=27209191/1000000,y=14.8,yend=15.2,color="#edd746",size=1)
dev.off()

rbind(dstat,pcadmix) %>% arrange(chr,start) %>% View()
  mutate(middle =   ((end - start)/2)+start) %>%
  ggplot(.,aes(x=middle,y=chr))

#Overlap regions:
  # 13	159486826	159500880 159493853
  #Ha412HOChr13g0625071
  anthocyanidin 3-O-glucosyltransferase 2-like
  XP_022002422 
  # 15	27160973	27257409  27209191
