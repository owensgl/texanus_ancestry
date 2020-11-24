#plot FST between texanus and central annuus
library(tidyverse)
library(ggthemes)
library(PNWColors)
library(patchwork)
chosen_colors <- pnw_palette(name="Bay",n=3,type="discrete")

directory <- "/home/owens/working/texanus/WGS/fst"

ANNCEN_ANNTEX <- read_tsv(paste(directory,"/","Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.ANNCEN_ANNTEX.fst.txt.gz",sep=""))
ANNCEN_ANNTEX$comparison <- "ANNCEN_ANNTEX"
ANNCEN_DEB <- read_tsv(paste(directory,"/","Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.ANNCEN_DEB.fst.txt.gz",sep=""))
ANNCEN_DEB$comparison <- "ANNCEN_DEB"
ANNTEX_DEB <- read_tsv(paste(directory,"/","Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.ANNTEX_DEB.fst.txt.gz",sep=""))
ANNTEX_DEB$comparison <- "ANNTEX_DEB"
ANNTEX_ARG <- read_tsv(paste(directory,"/","Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.ANNTEX_ARG.fst.txt.gz",sep=""))
ANNTEX_ARG$comparison <- "ANNTEX_ARG"
ANNCEN_ARG <- read_tsv(paste(directory,"/","Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.ANNCEN_ARG.fst.txt.gz",sep=""))
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


#Load in inversion locations:
inv_locations <- read_tsv("../wild_gwas_2018/MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.inversions.regions.v1.txt") %>%
  mutate(full_id = paste0(chr,"_",mds)) %>%
  mutate(id = gsub("Ha412HOChr","",full_id)) %>%
  mutate(id = gsub("pos","",id)) %>%
  mutate(id = gsub("neg","",id)) %>%
  mutate(id = gsub("syn","",id)) %>%
  filter(species == "Annuus")

chr_lengths %>% 
  select(chr,end) %>%
  rename(chr_len = end) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(inv_locations, ., by=c("chr"="chr")) %>%
  # Add a cumulative position of each SNP
  arrange(chr, start) %>%
  mutate( startcum=start+tot,endcum=end+tot) %>%
  separate(chr, c("ref","chr_n"), "Chr") %>%
  mutate(chr_n = as.numeric(chr_n)) -> inv_locations

window_size = 1000000
#pdf("figures/Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.genomewidefst.pdf",width=16,height=6)
fst_plot1 <- fst_cum %>%
  filter(comparison == "ANNCEN_ANNTEX") %>%
  mutate(window = floor(poscum/window_size)*window_size) %>% 
  group_by(comparison, chr,window) %>%
  summarize(window_fst = sum(FstNum)/sum(FstDenom),count=n()) %>%
  filter(count > 10) %>%
  filter(!is.na(window_fst)) %>%
  ggplot(aes(x=window,y=window_fst)) + 
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
  xlab("Chromosome") + ggtitle("H. a. annuus:H. a. texanus") +
  ylab(expression(F[ST])) +
  geom_segment(data=inv_locations,
               aes(x=startcum,xend=endcum,
                   y=0,yend=0),
               size=2,color="#0f85a0",alpha=0.8) +
  theme(plot.title = element_text(face = "italic"))


fst_plot2 <- arg_snps %>%
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
  scale_color_manual(values=chosen_colors[c(1,3)],name="Comparison",labels=c("H. a. annuus", "H. a. texanus")) +
  xlab("Chromosome") + ggtitle("H. annuus:H. argophyllus") +
  ylab(expression(F[ST])) +
  theme(plot.title = element_text(face = "italic"))

fst_plot3 <- deb_snps %>%
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
  scale_color_manual(values=chosen_colors[c(1,3)],name="Comparison",labels=c("H. a. annuus", "H. a. texanus")) +
  xlab("Chromosome") + ggtitle("H. annuus:H. debilis") +
  ylab(expression(F[ST])) +
  theme(plot.title = element_text(face = "italic"),
        legend.position = "none")

genome_wide_plots <- (fst_plot1 / fst_plot3 / fst_plot2)
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
#dev.off()

arg_snps %>%
  mutate(window = floor(poscum/window_size)*window_size) %>% 
  group_by(chr,window) %>%
  summarize(window_fst_cen = sum(FstNum_ANNCEN_ARG)/sum(FstDenom_ANNCEN_ARG),
            window_fst_tex = sum(FstNum_ANNTEX_ARG)/sum(FstDenom_ANNTEX_ARG),
            nloci=n()) %>%
  filter(nloci > 10)-> arg_windows 


t.test(arg_windows %>% 
  inner_join(genetic_map) %>%
  mutate(cm_rank=100*(rank(rate)/length(rate))) %>%
  mutate(cm_quantile = floor(cm_rank/20)+1) %>%
  mutate(cm_quantile = case_when(cm_quantile > 5 ~ 5,
                                 TRUE ~ cm_quantile)) %>%
  mutate(fst_dif = window_fst_cen - window_fst_tex) %>%
  filter(cm_quantile == 5) %>% pull(window_fst_cen),
  arg_windows %>% 
    inner_join(genetic_map) %>%
    mutate(cm_rank=100*(rank(rate)/length(rate))) %>%
    mutate(cm_quantile = floor(cm_rank/20)+1) %>%
    mutate(cm_quantile = case_when(cm_quantile > 5 ~ 5,
                                   TRUE ~ cm_quantile)) %>%
    mutate(fst_dif = window_fst_cen - window_fst_tex) %>%
    filter(cm_quantile == 5) %>% pull(window_fst_tex))



rate_arg_plot <- arg_windows %>% 
  inner_join(genetic_map) %>%
  mutate(cm_rank=100*(rank(rate)/length(rate))) %>%
  mutate(cm_quantile = floor(cm_rank/20)+1) %>%
  mutate(cm_quantile = case_when(cm_quantile > 5 ~ 5,
                                 TRUE ~ cm_quantile)) %>%
  mutate(fst_dif = window_fst_cen - window_fst_tex) %>%
  pivot_longer(
    cols = starts_with("window_fst"),
    names_to = "comparison",
    values_to = "Fst",
  ) %>%
  ggplot(.,aes(x=as.factor(cm_quantile),y=Fst,fill=comparison)) + 
  #geom_jitter()  + 
  geom_boxplot() +
  scale_fill_manual(values=chosen_colors[c(1,3)],name="Comparison",labels=c("H. a. annuus", "H. a. texanus")) +
  ylab(expression(paste(F[ST]))) + xlab("Recombination rate percentile") +
  ggtitle("H. argophyllus") +
  theme_cowplot() +
  theme(legend.position = "none") +
  scale_x_discrete(labels=c("0-20", "20-40", "40-60","60-80","80-100")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.x = element_text(size=9)) 


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

rate_deb_plot <- deb_windows %>%
  inner_join(genetic_map) %>%
  mutate(cm_rank=100*(rank(rate)/length(rate))) %>%
  mutate(cm_quantile = floor(cm_rank/20)+1) %>%
  mutate(cm_quantile = case_when(cm_quantile > 5 ~ 5,
                                 TRUE ~ cm_quantile)) %>%
  mutate(fst_dif = window_fst_cen - window_fst_tex) %>%
  pivot_longer(
    cols = starts_with("window_fst"),
    names_to = "comparison",
    values_to = "Fst",
  ) %>%
  ggplot(.,aes(x=as.factor(cm_quantile),y=Fst,fill=comparison)) + 
  #geom_jitter()  + 
  geom_boxplot() +
  scale_fill_manual(values=chosen_colors[c(1,3)],name="Comparison",labels=c("H. a. annuus", "H. a. texanus")) +
  ylab(expression(paste(F[ST]))) + xlab("Recombination rate percentile") +
  ggtitle("H. debilis") +
  theme_cowplot() +
  scale_x_discrete(labels=c("0-20", "20-40", "40-60","60-80","80-100")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.x = element_text(size=9),
        legend.position="bottom") 
t.test(deb_windows$window_fst_cen,deb_windows$window_fst_tex,paired=T)


##Make plot of fst windows
#pdf("figures/Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.fstwindows.pdf",height=6,width=6)
boxplot_plot <- rbind(deb_windows %>%
  mutate(species = "Deb"),
  arg_windows %>%
    mutate(species = "Arg")) %>%
  ungroup() %>%
  select(-chr,-window,-nloci) %>%
  pivot_longer(-species, names_to = "comparison", values_to = "fst") %>%
  ggplot(.,aes(x=species,y=fst,fill=comparison)) + geom_boxplot() +
  scale_fill_manual(values=chosen_colors[c(1,3)],labels=c("H. a. annuus", "H. a. texanus"),
                    name="Comparison") +
  ylab(expression(F[ST])) + xlab("Species") +
  theme_cowplot() +
  theme(legend.position = "none")
dev.off()
  
#Correlations
tmp1.deb <- deb_windows %>%
  mutate(fstdif = window_fst_tex - window_fst_cen) %>%
  select(chr,window,fstdif)

tmp1.arg <- arg_windows %>%
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

#pdf("figures/Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.fstrelationship.pdf",height=6,width=6)

correlation_plot_deb <- inner_join(tmp1.deb, tmp2) %>%
  ggplot(.,aes(x=ann_fst,y=fstdif)) + 
  geom_point() + 
  geom_smooth(method="lm",color="grey") +
  #theme_bw() + 
  ylab("H. debilis:H. annuus\n(ssp. texanus - ssp. annuus)") +
  xlab("H. a. annuus:H. a. texanus")

correlation_plot_arg <- inner_join(tmp1.arg, tmp2) %>%
  ggplot(.,aes(x=ann_fst,y=fstdif)) + 
  geom_point() + 
  geom_smooth(method="lm",color="grey") +
  #theme_bw() + 
  ylab("H. argophyllu:H. annuus\n(ssp. texanus - ssp. annuus)") +
  xlab("H. a. annuus:H. a. texanus")

#dev.off()

cor.test(inner_join(tmp1.deb, tmp2) %>% pull(fstdif),
         inner_join(tmp1.deb, tmp2) %>% pull(ann_fst))

cor.test(inner_join(tmp1.arg, tmp2) %>% pull(fstdif),
         inner_join(tmp1.arg, tmp2) %>% pull(ann_fst))

#Master plot

other_plots <- boxplot_plot  | rate_deb_plot | rate_arg_plot

pdf("figures/Texanus_all_fst.v2.pdf",height=9,width=7)
(genome_wide_plots / other_plots)  + plot_annotation(tag_levels = 'A') +plot_layout(heights = c(1,1,1,3))
dev.off()
