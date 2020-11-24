library(tidyverse)
library(ggthemes)

directory <- "/home/owens/working/texanus/WGS/fst"

ANNCEN_ANNTEX <- read_tsv(paste(directory,"/","Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.equalsize.ANNCEN_ANNTEX.fst.txt.gz",sep=""))
ANNCEN_ANNTEX$comparison <- "ANNCEN_ANNTEX"
ANNCEN_DEB <- read_tsv(paste(directory,"/","Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.equalsize.ANNCEN_DEB.fst.txt.gz",sep=""))
ANNCEN_DEB$comparison <- "ANNCEN_DEB"
ANNTEX_DEB <- read_tsv(paste(directory,"/","Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.equalsize.ANNTEX_DEB.fst.txt.gz",sep=""))
ANNTEX_DEB$comparison <- "ANNTEX_DEB"
ANNTEX_ARG <- read_tsv(paste(directory,"/","Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.equalsize.ANNTEX_ARG.fst.txt.gz",sep=""))
ANNTEX_ARG$comparison <- "ANNTEX_ARG"
ANNCEN_ARG <- read_tsv(paste(directory,"/","Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.equalsize.ANNCEN_ARG.fst.txt.gz",sep=""))
ANNCEN_ARG$comparison <- "ANNCEN_ARG"

genetic_map <- read_tsv("/home/owens/ref/Ha412HOv2.0-20181130.Nov22k22.geneticmap.extradivisions.txt") %>% 
  mutate(pos = pos-1) %>% rename(window=pos) %>%
  mutate(rate = lead(cM)-cM) %>% 
  mutate(rate = case_when(is.na(rate) ~ cM - lag(cM),
                          rate < 0 ~ cM - lag(cM),
                          TRUE ~ rate))


genetic_map_100k <- genetic_map %>% 
  mutate(row_id = toString(seq(0,9))) %>% 
  separate_rows(row_id) %>%
  mutate(row_id = as.numeric(row_id)) %>%
  mutate(window = window + (row_id*100000)) 

chr_lengths <- read_tsv("../wild_gwas_2018/Ha412HO.chrlengths.txt")

fst <- rbind(ANNCEN_ANNTEX, ANNCEN_DEB, ANNTEX_DEB,ANNCEN_ARG, ANNTEX_ARG)

#Remove inversion sites
non_inv_sites <- read_tsv(paste0(directory,"/../","Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.v3noinv.sites.txt"),col_names = c("chr","pos"))
fst <- fst %>% inner_join(non_inv_sites)
fst_cum <- chr_lengths %>%
  select(chr,end) %>%
  rename(chr_len = end) %>%
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(fst, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, pos) %>%
  mutate( poscum=pos+tot)

axisdf = fst_cum %>% group_by(chr) %>% summarize(center=( max(poscum) + min(poscum) ) / 2,
                                                 max= max(poscum),
                                                 min=min(poscum))
fst_cum %>%
  select(chr,poscum,FstNum,FstDenom,comparison) %>%
  filter(comparison == "ANNCEN_DEB" | comparison == "ANNTEX_DEB") %>%
  pivot_wider(names_from = comparison, values_from = c("FstNum","FstDenom")) %>%
  filter(!is.na(FstDenom_ANNCEN_DEB), !is.na(FstDenom_ANNTEX_DEB)) -> deb_snps

fst_cum %>%
  select(chr,poscum,FstNum,FstDenom,comparison) %>%
  filter(comparison == "ANNCEN_ARG" | comparison == "ANNTEX_ARG") %>%
  pivot_wider(names_from = comparison, values_from = c("FstNum","FstDenom")) %>%
  filter(!is.na(FstDenom_ANNCEN_ARG), !is.na(FstDenom_ANNTEX_ARG)) -> arg_snps

window_size = 1000000
pdf("figures/Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.equalsize.genomewidefst.pdf",width=16,height=6)
fst_cum %>%
  filter(comparison == "ANNCEN_ANNTEX") %>%
  mutate(window = floor(poscum/window_size)*window_size) %>% 
  group_by(comparison, chr,window) %>%
  summarize(window_fst = sum(FstNum)/sum(FstDenom),count=n()) %>%
  filter(count > 10) %>%
  filter(!is.na(window_fst)) %>%
  ggplot(aes(x=window,y=window_fst)) + 
  scale_x_continuous( label = gsub("Ha412HOChr","",axisdf$chr), breaks= axisdf$center ) +
  geom_vline(xintercept = axisdf$max[1:nrow(axisdf)-1],linetype="dotted") +
  #scale_linetype_manual(values = rep(c("solid", "dashed"), 22 )) +
  theme_bw() +
  theme(legend.position="bottom",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  geom_line(alpha=0.7) + 
  xlab("Chromosome") + ggtitle("H. annuus: Central - Texas") +
  ylab(expression(F[ST]))


arg_snps %>%
  mutate(window = floor(poscum/window_size)*window_size) %>% 
  group_by(chr,window) %>%
  summarize(ANNCEN_ARG = sum(FstNum_ANNCEN_ARG)/sum(FstDenom_ANNCEN_ARG),
            ANNTEX_ARG = sum(FstNum_ANNTEX_ARG)/sum(FstDenom_ANNTEX_ARG),
            count=n()) %>%
  filter(count > 10) %>%
  pivot_longer(c("ANNCEN_ARG","ANNTEX_ARG"),names_to="comparison",values_to="window_fst") %>%
  filter(!is.na(window_fst)) %>%
  ggplot(aes(x=window,y=window_fst,color=comparison)) + 
  scale_x_continuous( label = gsub("Ha412HOChr","",axisdf$chr), breaks= axisdf$center,
                      expand = c(0.01,0.01)) +
  geom_vline(xintercept = axisdf$max[1:nrow(axisdf)-1],linetype="dotted") +
  #scale_linetype_manual(values = rep(c("solid", "dashed"), 22 )) +
  theme_bw() +
  theme(legend.position="bottom",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  geom_line(alpha=0.7) + 
  scale_color_brewer(palette = "Set1",name="Comparison",labels=c("Central ANN", "Texas ANN")) +
  xlab("Chromosome") + ggtitle("H. annuus:H. argophyllus") +
  ylab(expression(F[ST]))

deb_snps %>%
  mutate(window = floor(poscum/window_size)*window_size) %>% 
  group_by(chr,window) %>%
  summarize(ANNCEN_DEB = sum(FstNum_ANNCEN_DEB)/sum(FstDenom_ANNCEN_DEB),
            ANNTEX_DEB = sum(FstNum_ANNTEX_DEB)/sum(FstDenom_ANNTEX_DEB),
            count=n()) %>%
  filter(count > 10) %>%
  pivot_longer(c("ANNCEN_DEB","ANNTEX_DEB"),names_to="comparison",values_to="window_fst") %>%
  filter(!is.na(window_fst)) %>%
  ggplot(aes(x=window,y=window_fst,color=comparison)) + 
  scale_x_continuous( label = gsub("Ha412HOChr","",axisdf$chr), breaks= axisdf$center,
                      expand = c(0.01,0.01)) +
  geom_vline(xintercept = axisdf$max[1:nrow(axisdf)-1],linetype="dotted") +
  #scale_linetype_manual(values = rep(c("solid", "dashed"), 22 )) +
  theme_bw() +
  theme(legend.position="bottom",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  geom_line(alpha=0.7) + 
  scale_color_brewer(palette = "Set1",name="Comparison",labels=c("Central ANN", "Texas ANN")) +
  xlab("Chromosome") + ggtitle("H. annuus:H. debilis") +
  ylab(expression(F[ST]))


deb_snps %>%
  mutate(window = floor(poscum/window_size)*window_size) %>% 
  group_by(chr,window) %>%
  summarize(ANNCEN_DEB = sum(FstNum_ANNCEN_DEB)/sum(FstDenom_ANNCEN_DEB),
            ANNTEX_DEB = sum(FstNum_ANNTEX_DEB)/sum(FstDenom_ANNTEX_DEB),
            count=n()) %>%
  filter(count > 10) %>%
  mutate(Fstdif = ( ANNTEX_DEB - ANNCEN_DEB)) %>%
  filter(!is.na(Fstdif)) %>% 
  ggplot() + 
  theme_few() + 
  geom_line(aes(x=window,y=Fstdif,color=Fstdif)) +
  scale_x_continuous( label = gsub("Ha412HOChr","",axisdf$chr), breaks= axisdf$center,
                      expand = c(0.01,0.01)) +
  geom_vline(xintercept = axisdf$max[1:nrow(axisdf)-1],linetype="dotted") +
  #scale_linetype_manual(values = rep(c("solid", "dashed"), 22 )) +
  theme_bw() +
  theme(legend.position="bottom",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  scale_color_gradient2(high="#E41A1C",mid="grey",low="#377EB8",
                        name="FST difference") +
  xlab("Chromosome") + ggtitle("H. debilis:H. annuus (Texas - Central)") +
  ylab(expression(paste(F[ST]," difference"))) +
  theme(legend.position = "none")
dev.off()


arg_snps %>%
  mutate(window = floor(poscum/window_size)*window_size) %>% 
  group_by(chr,window) %>%
  summarize(window_fst_cen = sum(FstNum_ANNCEN_ARG)/sum(FstDenom_ANNCEN_ARG),
            window_fst_tex = sum(FstNum_ANNTEX_ARG)/sum(FstDenom_ANNTEX_ARG),
            nloci=n()) %>%
  filter(nloci > 10)-> arg_windows 


t.test(arg_windows$window_fst_cen,arg_windows$window_fst_tex,paired=T)



deb_snps %>%
  mutate(ANNCEN_DEB = case_when(ANNCEN_DEB < 0 ~ 0,
                                TRUE ~ ANNCEN_DEB)) %>%
  mutate(ANNTEX_DEB = case_when(ANNTEX_DEB < 0 ~ 0,
                                TRUE ~ ANNTEX_DEB)) %>%
  mutate(fst_dif = ANNCEN_DEB - ANNTEX_DEB) %>%
  ggplot(.,aes(fst_dif)) + geom_histogram()



deb_snps %>%
  mutate(window = floor(poscum/window_size)*window_size) %>% 
  group_by(chr,window) %>%
  summarize(window_fst_cen = sum(FstNum_ANNCEN_DEB)/sum(FstDenom_ANNCEN_DEB),
            window_fst_tex = sum(FstNum_ANNTEX_DEB)/sum(FstDenom_ANNTEX_DEB),
            nloci=n()) %>%
  filter(nloci > 10)-> deb_windows 

deb_windows %>% 
  inner_join(genetic_map) %>%
  mutate(fst_dif = window_fst_cen - window_fst_tex) %>%
  ggplot(.,aes(x=rate,y=nloci)) + geom_point()  

t.test(deb_windows$window_fst_cen,deb_windows$window_fst_tex,paired=T)


##Make plot of fst windows
pdf("figures/Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.equalsize.fstwindows.pdf",height=6,width=6)
rbind(deb_windows %>%
        mutate(species = "Deb"),
      arg_windows %>%
        mutate(species = "Arg")) %>%
  ungroup() %>%
  select(-chr,-window,-nloci) %>%
  pivot_longer(-species, names_to = "comparison", values_to = "fst") %>%
  ggplot(.,aes(x=species,y=fst,fill=comparison)) + geom_boxplot() +
  scale_fill_brewer(palette = "Set1",labels=c("Central Ann", "South\nTexas Ann"),
                    name="Comparison") +
  ylab(expression(F[ST])) + xlab("Species") +
  theme_bw()
dev.off()

#Correlations
tmp1 <- deb_windows %>%
  mutate(fstdif = window_fst_tex - window_fst_cen) %>%
  select(chr,window,fstdif)


tmp2 <- fst_cum %>%
  filter(comparison == "ANNCEN_ANNTEX") %>%
  mutate(window = floor(poscum/window_size)*window_size) %>% 
  group_by(comparison, chr,window) %>%
  summarize(window_fst = sum(FstNum)/sum(FstDenom),count=n()) %>%
  filter(count > 10) %>%
  filter(!is.na(window_fst)) %>%
  ungroup() %>%
  select(chr,window,window_fst) %>%
  rename(ann_fst = window_fst)

pdf("figures/Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.equalsize.fstrelationship.pdf",height=6,width=6)

inner_join(tmp1, tmp2) %>%
  ggplot(.,aes(x=ann_fst,y=fstdif)) + 
  geom_point() + 
  geom_smooth(method="lm",color="light green") +
  theme_bw() + 
  ylab("H. debilis:H. annuus (Texas - Central)") +
  xlab("Central H. annuus:Texas H. annuus")
dev.off()

cor.test(inner_join(tmp1, tmp2) %>% pull(fstdif),
         inner_join(tmp1, tmp2) %>% pull(ann_fst))
