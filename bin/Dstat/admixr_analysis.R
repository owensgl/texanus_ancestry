#Edit the ~/.Renviron file to give admixr the path to AdmixTools
library(admixr)
library(tidyverse)
library(PNWColors)
library(patchwork)
library(cowplot)
library(ggrepel)
library(ggforce)
five_colors <- pnw_palette(name="Bay",n=5,type="discrete")

#Run admixr for the full genome, for each texas sample
directory <- "/home/owens/working/texanus/WGS/"
prefix <- "Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO"
snps <- eigenstrat(paste(directory,"/",prefix,sep=""))
samples <- read_tsv(paste0(directory,"/",prefix,".ind"),col_names = c("sample","spacer","group"))
# 
# #Pick all texas individuals to run individually
# full_genome <- tibble(W=character(),X=character(),Y=character(),Z=character(),
#                       D=numeric(),stderr=numeric(),Zscore=numeric(),BABA=numeric(),
#                       ABBA=numeric(),nsps=numeric())
# for (i in 14:1){
#   j=(i*10)
#   tex_samples <- samples %>%
#     filter(group != "ANN_CEN") %>%
#     filter(grepl("ANN",group)) %>%
#     tail(j) %>%
#     head(10) %>%
#     pull(group)
#   tmp <- snps %>%
#     d(W = tex_samples, X = "ANN_CEN", Y = c("PET_FAL","DEB","ARG"), Z = "OUT") 
#   full_genome <- rbind(full_genome,tmp)
# }
# 
# 
arg_deb <- snps %>%
  d(W = "ANN_STEX", X = "ARG", Y = c("DEB"), Z = "OUT")
arg_deb_2 <- snps %>%
  d(W = "DEB", X = "PET_FAL", Y = c("ARG"), Z = "OUT")
arg_deb_3 <- snps %>%
  d(W = "DEB", X = "PET_PET", Y = c("ARG"), Z = "OUT")
arg_deb_4 <- snps %>%
  d(W = "DEB", X = "PET_PET", Y = c("ANN_CEN"), Z = "OUT")
arg_deb_5 <- snps %>%
  d(W = "DEB", X = "PET_PET", Y = c("ANN_STEX"), Z = "OUT")
# 
# 
# full_genome <- unique(full_genome)
# write_tsv(full_genome, "Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.Dstats.txt")
# 
# #Run admixr for the full genome, for each texas sample
# prefix <- "Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.v3noinv"
# snps <- eigenstrat(paste(directory,"/",prefix,sep=""))
# samples <- read_tsv(paste0(directory,"/",prefix,".ind"),col_names = c("sample","spacer","group"))
# 
# #Pick all texas individuals to run individually without inversions
# full_genome <- tibble(W=character(),X=character(),Y=character(),Z=character(),
#                       D=numeric(),stderr=numeric(),Zscore=numeric(),BABA=numeric(),
#                       ABBA=numeric(),nsps=numeric())
# for (i in 14:1){
#   j=(i*10)
#   tex_samples <- samples %>%
#     filter(group != "ANN_CEN") %>%
#     filter(grepl("ANN",group)) %>%
#     tail(j) %>%
#     head(10) %>%
#     pull(group)
#   tmp <- snps %>%
#     d(W = tex_samples, X = "ANN_CEN", Y = c("PET_FAL","DEB","ARG"), Z = "OUT") 
#   full_genome <- rbind(full_genome,tmp)
# }
# 
# full_genome <- unique(full_genome)
# write_tsv(full_genome, "Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.v3noinvDstats.txt")

###Plotting all Dstat scores
labels <- read_tsv("/home/owens/working/sample_info_apr_2018.tsv") %>%
  rename(sample = name) 
pop_loc <- read_tsv("pop_loc_allnum.txt") %>%
  rename(population = pop)
labels %>% inner_join(.,pop_loc) -> labels
labels %>%
  mutate(group = case_when(lat < 30 ~ "STex",
                           lat < 33 & long > -105 ~ "NTex",
                           TRUE ~ "Central")) -> labels

data <- read_tsv("Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.SVs.Dstat.txt")

data <- read_tsv("Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.v3noinvDstats.txt") %>%
  mutate(dataset = "no_SV") %>% rbind(.,data)

data <- read_tsv("Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.Dstats.txt") %>%
  mutate(dataset = "full_genome") %>% rbind(.,data)

data <- data %>%
  rename(sample = W) %>%
  inner_join(.,labels)

dstat_data <- data
#ANN_55, ANN_56, ANN_57
#ANN1363 <- most debilis introgressed

ind_d_unorganized <- data %>%
  filter(group == "STex") %>%
  filter(dataset=="full_genome") %>% 
  mutate(sig = case_when(abs(Zscore) > 1.96 ~ "p<=0.05",
                         TRUE ~ "p>0.05")) %>% 
  mutate(text_label = case_when(sample == "ANN1363" & Y == "ARG"  ~ "ANN1363",
                                #population == "ANN_55" ~ "ANN_55",
                                #population == "ANN_56" ~ "ANN_56",
                                #population == "ANN_57" ~ "ANN_57",
                                TRUE ~ "")) %>%
  ggplot(.,aes(x=fct_reorder(sample,D),y=D)) + 
  geom_point(aes(color=sig),width=0.1) +
  geom_linerange(aes(ymin=D-(stderr*1.96),ymax=D+(stderr*1.96),
                    color=sig)) +
  #theme_bw() +
  scale_color_brewer(palette = "Set1",name="Significance") +
  geom_hline(yintercept=0,linetype="dotted") +
  facet_wrap(~Y,nrow=3,scales="free_y") +
  theme(axis.text.x = element_blank()) +
  xlab("Sample") 

ind_d_organized <- data %>%
  filter(group == "STex" ) %>%
  mutate(Y = case_when(Y == "PET_FAL" ~ "PET",
                       TRUE ~ Y)) %>%
  filter(dataset=="full_genome") %>% 
  mutate(sig = case_when(abs(Zscore) > 1.96 ~ "p<=0.05",
                         TRUE ~ "p>0.05"),
         Y = fct_relevel(Y,c("DEB","ARG","PET"))) %>% 
  ggplot(.,aes(x=fct_reorder(sample,D),y=D)) + 
  theme_cowplot() +
  geom_point(aes(color=sig),width=0.1) +
  geom_linerange(aes(ymin=D-(stderr*1.96),ymax=D+(stderr*1.96),
                     color=sig)) +
  scale_color_manual(values=five_colors[c(1,5)],name="Significance") +
  geom_hline(yintercept=0,linetype="dotted") +
  facet_grid(Y~population,scales="free") +
  theme(axis.text.x = element_blank(),
        legend.position = "bottom",
        axis.text.y = element_text(size=6),
        strip.text.x = element_text(size = 8)) +
  xlab("Sample")  +
  geom_text_repel(data=data %>% filter(dataset=="full_genome",
                                       sample == "ANN1363",Y=="DEB"),
                  aes(label=sample),box.padding=0.3, nudge_y=2,size=2.5)

ind_d_correlation <- data %>%
  filter(group == "STex") %>%
  filter(dataset=="full_genome") %>%
  select(sample,Y,D,stderr) %>%
  pivot_wider(names_from = Y, values_from = c(D,stderr)) %>%
  mutate(text_label = case_when(sample == "ANN1363"  ~ "ANN1363",
                                TRUE ~ "")) %>%
  ggplot(.,aes(x=D_DEB,y=D_ARG)) + 
  geom_point(,width=0.1) +
  geom_linerange(aes(ymin=D_ARG-(stderr_ARG*1.96),ymax=D_ARG+(stderr_ARG*1.96))) +
  geom_errorbarh(aes(xmin=D_DEB-(stderr_DEB*1.96),xmax=D_DEB+(stderr_DEB*1.96),width=0)) +
  theme_cowplot() +
  geom_hline(yintercept = 0,linetype="dotted")+
  geom_vline(xintercept = 0,linetype="dotted") +
  geom_smooth(method="lm",color="dark grey") +
  ylab("H. argophyllus D") + xlab("H. debilis D") +
  geom_text_repel(aes(label=text_label),point.padding=1)


data %>%
  filter(group == "STex") %>%
  filter(dataset=="full_genome") %>%
  select(sample,Y,D,stderr) %>%
  pivot_wider(names_from = Y, values_from = c(D,stderr))-> tmp

summary(lm(D_DEB ~ D_ARG, tmp))
data %>%
  filter(group == "STex") %>%
  filter(dataset=="no_SV") %>%
  mutate(sig = case_when(abs(Zscore) > 1.96 ~ "p<=0.05",
                         TRUE ~ "p>0.05")) %>%
  ggplot(.,aes(x=fct_reorder(population,lat),y=D)) + geom_point(aes(color=sig)) +
  theme_bw() +
  scale_color_brewer(palette = "Set1",name="Significance") +
  geom_hline(yintercept=0,linetype="dotted") +
  facet_wrap(~Y)



data %>%
  filter(group == "STex") %>%
  filter(dataset=="no_SV") %>% 
  mutate(sig = case_when(abs(Zscore) > 1.96 ~ "p<=0.05",
                         TRUE ~ "p>0.05")) %>%
  ggplot(.,aes(x=fct_reorder(sample,D),y=D)) + 
  geom_point(aes(color=sig),width=0.1) +
  geom_linerange(aes(ymin=D-(stderr*1.96),ymax=D+(stderr*1.96),
                     color=sig)) +
  theme_bw() +
  scale_color_brewer(palette = "Set1",name="Significance") +
  geom_hline(yintercept=0,linetype="dotted") +
  facet_wrap(~Y,nrow=3,scales="free_y") +
  theme(axis.text.x = element_blank()) +
  xlab("Sample") +ylab("D without SV")

data %>%
  filter(group == "STex") %>%
  filter(Y != "PET_FAL") %>%
  filter(dataset=="full_genome" | dataset == "no_SV") %>%
  select(sample,X,Y,D,population,lat,long,group,dataset) %>%
  pivot_wider(names_from = dataset, values_from = D) %>%
  ggplot(.,aes(x=full_genome,y=no_SV))+ geom_point() + geom_smooth(method="lm") +
  facet_wrap(~Y,scales="free") + geom_abline(slope=1) +
  ylab("D without SV") + xlab("D full genome")

############################
#Load up genome wide D scores. 
############################
genetic_map <- read_tsv("/home/owens/ref/Ha412HOv2.0-20181130.Nov22k22.geneticmap.extradivisions.txt") %>% 
  mutate(rate = lead(cM)-cM) %>% 
  rename(start = pos) %>% mutate(end = start + 999999) %>%
  mutate(rate = case_when(is.na(rate) ~ cM - lag(cM),
                          rate < 0 ~ cM - lag(cM),
                          TRUE ~ rate))
chr_lengths <- read_tsv("../wild_gwas_2018/Ha412HO.chrlengths.txt")



dstats <- read_tsv("/home/owens/working/texanus/WGS/dstat/Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.ANN_CEN.ANN_STEX.ARG.Dstats.txt.gz")

window_snps <- 100
dstats %>%
  group_by(chr) %>%
  mutate(window = floor(row_number()/window_snps)) %>%
  group_by(chr,window) %>%
  summarize(middle= mean(pos),
            start= min(pos),end=max(pos),
            all_abba = sum(abba),
            all_baba = sum(baba),
            all_p3_abba = sum(p3_abba),
            all_p3_baba = sum(p3_baba),
            D = (sum(abba) - sum(baba)) / (sum(abba) + sum(baba)),
            fd = (sum(abba) - sum(baba)) / (sum(p3_abba) - sum(p3_baba))) %>% 
  mutate(fd = case_when(D < 0 ~ 0,
                        TRUE ~ fd)) %>%
  ungroup() %>%
  mutate(rank = percent_rank(fd),
         outlier = case_when(rank >0.998 ~ 1,
                             TRUE ~ 0)) -> dstats_windows_arg

big_window <- 2000
dstats %>%
  group_by(chr) %>%
  mutate(window = floor(row_number()/big_window)) %>%
  group_by(chr,window) %>%
  summarize(middle= mean(pos),
            start= min(pos),end=max(pos),
            all_abba = sum(abba),
            all_baba = sum(baba),
            all_p3_abba = sum(p3_abba),
            all_p3_baba = sum(p3_baba),
            D = (sum(abba) - sum(baba)) / (sum(abba) + sum(baba)),
            fd = (sum(abba) - sum(baba)) / (sum(p3_abba) - sum(p3_baba))) %>% 
  mutate(fd = case_when(D < 0 ~ 0,
                        TRUE ~ fd)) %>%
  ungroup() %>%
  mutate(rank = percent_rank(fd),
         outlier = case_when(rank >0.998 ~ 1,
                             TRUE ~ 0)) -> dstats_bigwindows_arg

#Add recombination rate to dstats_windows
dstats_windows_arg$cm_rate <- NA
for (i in 1:nrow(dstats_windows_arg)){
  print(i)
  chosen_start <- dstats_windows_arg$start[i]
  chosen_end <-dstats_windows_arg$end[i]
  chosen_chr <- dstats_windows_arg$chr[i]
  start_rate <- genetic_map %>% filter(chr == chosen_chr) %>% 
    filter(start <= chosen_start, end >= chosen_start) %>%
    pull(rate)
  end_rate <- genetic_map %>% filter(chr == chosen_chr) %>% 
    filter(start <= chosen_end, end >= chosen_end) %>%
    pull(rate)
  mean_rate <- mean(start_rate,end_rate)
  dstats_windows_arg$cm_rate[i] <- mean_rate
}

dstats_cumulative_arg <- chr_lengths %>%
  select(chr,end) %>%
  rename(chr_len = end) %>%
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(dstats_windows_arg, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, middle) %>%
  mutate( middle_cum=middle+tot)

dstats_curve_cumulative_arg <- chr_lengths %>%
  select(chr,end) %>%
  rename(chr_len = end) %>%
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(dstats_bigwindows_arg, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, middle) %>%
  mutate( middle_cum=middle+tot)

axisdf = dstats_cumulative_arg %>% group_by(chr) %>% summarize(center=( max(middle_cum) + min(middle_cum) ) / 2,
                                                           max= max(middle_cum),
                                                           min=min(middle_cum))

arg_genome_wide <- dstats_cumulative_arg %>%
  ggplot(.,aes(y=D,x=middle_cum,group=1)) + 
  geom_point(aes(color=as.factor(outlier)),alpha=0.4) +
  scale_color_manual(values=c("black","black")) +
  geom_line(data=dstats_curve_cumulative_arg,aes(x=middle_cum,y=D),color=five_colors[3],size=2,alpha=0.9) +
  geom_hline(yintercept=0,color="blue",linetype="dashed") +
  scale_x_continuous( label = gsub("Ha412HOChr","",axisdf$chr), breaks= axisdf$center,
                      expand = c(0.01,0.01)) +
  geom_vline(xintercept = axisdf$max[1:nrow(axisdf)-1],linetype="dotted") +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank()) +
  ggtitle("H. argophyllus") +
  xlab("Chromosome")


arg_correlation <- dstats_windows_arg %>%
  ggplot(.,aes(x=log10(cm_rate),y=D)) + geom_point() + geom_smooth(method="lm",color="grey") +
  ylab("D") + xlab("log10(cM/Mb)") + ggtitle("H. argophyllus")

###Now for deb


dstats <- read_tsv("/home/owens/working/texanus/WGS/dstat/Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.ANN_CEN.ANN_STEX.DEB.Dstats.txt.gz")

window_snps <- 100
dstats %>%
  group_by(chr) %>%
  mutate(window = floor(row_number()/window_snps)) %>%
  group_by(chr,window) %>%
  summarize(middle= mean(pos),
            start= min(pos),end=max(pos),
            all_abba = sum(abba),
            all_baba = sum(baba),
            all_p3_abba = sum(p3_abba),
            all_p3_baba = sum(p3_baba),
            D = (sum(abba) - sum(baba)) / (sum(abba) + sum(baba)),
            fd = (sum(abba) - sum(baba)) / (sum(p3_abba) - sum(p3_baba))) %>% 
  mutate(fd = case_when(D < 0 ~ 0,
                        TRUE ~ fd)) %>%
  ungroup() %>%
  mutate(rank = percent_rank(fd),
         outlier = case_when(rank >0.998 ~ 1,
                             TRUE ~ 0)) -> dstats_windows_deb

dstats %>%
  group_by(chr) %>%
  mutate(window = floor(row_number()/big_window)) %>%
  group_by(chr,window) %>%
  summarize(middle= mean(pos),
            start= min(pos),end=max(pos),
            all_abba = sum(abba),
            all_baba = sum(baba),
            all_p3_abba = sum(p3_abba),
            all_p3_baba = sum(p3_baba),
            D = (sum(abba) - sum(baba)) / (sum(abba) + sum(baba)),
            fd = (sum(abba) - sum(baba)) / (sum(p3_abba) - sum(p3_baba))) %>% 
  mutate(fd = case_when(D < 0 ~ 0,
                        TRUE ~ fd)) %>%
  ungroup() %>%
  mutate(rank = percent_rank(fd),
         outlier = case_when(rank >0.998 ~ 1,
                             TRUE ~ 0)) -> dstats_bigwindows_deb

#Add recombination rate to dstats_windows
dstats_windows_deb$cm_rate <- NA
for (i in 1:nrow(dstats_windows_deb)){
  print(i)
  chosen_start <- dstats_windows_deb$start[i]
  chosen_end <-dstats_windows_deb$end[i]
  chosen_chr <- dstats_windows_deb$chr[i]
  start_rate <- genetic_map %>% filter(chr == chosen_chr) %>% 
    filter(start <= chosen_start, end >= chosen_start) %>%
    pull(rate)
  end_rate <- genetic_map %>% filter(chr == chosen_chr) %>% 
    filter(start <= chosen_end, end >= chosen_end) %>%
    pull(rate)
  mean_rate <- mean(start_rate,end_rate)
  dstats_windows_deb$cm_rate[i] <- mean_rate
}

dstats_cumulative_deb <- chr_lengths %>%
  select(chr,end) %>%
  rename(chr_len = end) %>%
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(dstats_windows_deb, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, middle) %>%
  mutate( middle_cum=middle+tot)

dstats_curve_cumulative_deb <- chr_lengths %>%
  select(chr,end) %>%
  rename(chr_len = end) %>%
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(dstats_bigwindows_deb, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, middle) %>%
  mutate( middle_cum=middle+tot)

axisdf = dstats_cumulative_deb %>% group_by(chr) %>% summarize(center=( max(middle_cum) + min(middle_cum) ) / 2,
                                                               max= max(middle_cum),
                                                               min=min(middle_cum))

deb_genome_wide <- dstats_cumulative_deb %>%
  ggplot(.,aes(y=D,x=middle_cum,group=1)) + 
  geom_point(aes(color=as.factor(outlier)),alpha=0.4) +
  scale_color_manual(values=c("black","black")) +
  geom_line(data=dstats_curve_cumulative_deb,aes(x=middle_cum,y=D),color=five_colors[5],size=2,alpha=0.9) +
  geom_hline(yintercept=0,color="blue",linetype="dashed") +
  scale_x_continuous( label = gsub("Ha412HOChr","",axisdf$chr), breaks= axisdf$center,
                      expand = c(0.01,0.01)) +
  geom_vline(xintercept = axisdf$max[1:nrow(axisdf)-1],linetype="dotted") +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank()) +
  ggtitle("H. debilis") +
  xlab("Chromosome")


deb_correlation <- dstats_windows_deb %>%
  ggplot(.,aes(x=log10(cm_rate),y=D)) + geom_point() + geom_smooth(method="lm",color="grey") +
  ylab("D") + xlab("log10(cM/Mb)") + ggtitle("H. debilis")

summary(lm(D ~ log10(cm_rate), data=dstats_windows_deb))


################
#Now for ANN1363
################
dstats <- read_tsv("/home/owens/working/texanus/WGS/dstat/Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.ANN_CEN.ANN1363.DEB.Dstats.txt.gz")

window_snps <- 100
dstats %>%
  group_by(chr) %>%
  mutate(window = floor(row_number()/window_snps)) %>%
  group_by(chr,window) %>%
  summarize(middle= mean(pos),
            start= min(pos),end=max(pos),
            all_abba = sum(abba),
            all_baba = sum(baba),
            all_p3_abba = sum(p3_abba),
            all_p3_baba = sum(p3_baba),
            D = (sum(abba) - sum(baba)) / (sum(abba) + sum(baba)),
            fd = (sum(abba) - sum(baba)) / (sum(p3_abba) - sum(p3_baba))) %>% 
  mutate(fd = case_when(D < 0 ~ 0,
                        TRUE ~ fd)) %>%
  ungroup() %>%
  mutate(rank = percent_rank(fd),
         outlier = case_when(rank >0.998 ~ 1,
                             TRUE ~ 0)) -> dstats_windows_ann1363


dstats %>%
  group_by(chr) %>%
  mutate(window = floor(row_number()/big_window)) %>%
  group_by(chr,window) %>%
  summarize(middle= mean(pos),
            start= min(pos),end=max(pos),
            all_abba = sum(abba),
            all_baba = sum(baba),
            all_p3_abba = sum(p3_abba),
            all_p3_baba = sum(p3_baba),
            D = (sum(abba) - sum(baba)) / (sum(abba) + sum(baba)),
            fd = (sum(abba) - sum(baba)) / (sum(p3_abba) - sum(p3_baba))) %>% 
  mutate(fd = case_when(D < 0 ~ 0,
                        TRUE ~ fd)) %>%
  ungroup() %>%
  mutate(rank = percent_rank(fd),
         outlier = case_when(rank >0.998 ~ 1,
                             TRUE ~ 0)) -> dstats_bigwindows_ann1363


dstats_cumulative_ann1363 <- chr_lengths %>%
  select(chr,end) %>%
  rename(chr_len = end) %>%
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(dstats_windows_ann1363, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, middle) %>%
  mutate( middle_cum=middle+tot)

dstats_curve_cumulative_ann1363 <- chr_lengths %>%
  select(chr,end) %>%
  rename(chr_len = end) %>%
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(dstats_bigwindows_ann1363, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, middle) %>%
  mutate( middle_cum=middle+tot)

axisdf = dstats_cumulative_ann1363 %>% group_by(chr) %>% summarize(center=( max(middle_cum) + min(middle_cum) ) / 2,
                                                               max= max(middle_cum),
                                                               min=min(middle_cum))

ann1363_genome_wide <- dstats_cumulative_ann1363 %>%
  ggplot(.,aes(y=D,x=middle_cum,group=1)) + 
  geom_point(aes(color=as.factor(outlier)),alpha=0.4) +
  scale_color_manual(values=c("black","black")) +
  geom_line(data=dstats_curve_cumulative_ann1363,aes(x=middle_cum,y=D),color=five_colors[5],size=2,alpha=0.9) +
  geom_hline(yintercept=0,color="blue",linetype="dashed") +
  scale_x_continuous( label = gsub("Ha412HOChr","",axisdf$chr), breaks= axisdf$center,
                      expand = c(0.01,0.01)) +
  geom_vline(xintercept = axisdf$max[1:nrow(axisdf)-1],linetype="dotted") +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank()) +
  ggtitle("ANN1363 | H. debilis") +
  xlab("Chromosome")


####Load admixture results and plot them.
directory <- "/home/owens/working/texanus/WGS/admixture"

#pdf("Texanus.tranche90.snp.gwas.90.bi.0p2ld.admix.thin.admixture.pdf",height=15,width=15)
five_colors <- pnw_palette(name="Bay",n=5,type="discrete")
Q <- read_delim(paste(directory, "/Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.0p2ld.annargdeb.numchr.5.Q",sep=""),col_names=c("Q1","Q2","Q3","Q4","Q5"),
                delim=" ")
samples <- read_tsv(paste(directory, "/Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.admixr.annargdeb.groups.txt",sep=""),col_names = c("sample","group"))

data <- cbind(samples,Q)

vertical_lines <- c(data %>% filter(group== "ANN_CEN") %>% nrow(),
                    data %>% filter(group== "ANN_NTEX" | group == "ANN_CEN") %>% nrow(),
                    data %>% filter(group== "ANN_NTEX" | group == "ANN_CEN" | group== "ANN_STEX") %>% nrow(),
                    data %>% filter(group!= "DEB") %>% nrow())
pop_labels_positions <- c(data %>% filter(group== "ANN_CEN") %>% nrow()/2,
                          data %>% filter(group== "ANN_CEN") %>% nrow() + data %>% filter(group== "ANN_NTEX") %>% nrow()/2,
                          data %>% filter(group== "ANN_CEN") %>% nrow() + data %>% filter(group== "ANN_NTEX") %>% nrow() +
                            data %>% filter(group== "ANN_STEX") %>% nrow()/2,
                          data %>% filter(group== "ANN_CEN") %>% nrow() + data %>% filter(group== "ANN_NTEX") %>% nrow() +
                            data %>% filter(group== "ANN_STEX") %>% nrow() + data %>% filter(group== "ARG") %>% nrow()/2,
                          data %>% filter(group== "ANN_CEN") %>% nrow() + data %>% filter(group== "ANN_NTEX") %>% nrow() +
                            data %>% filter(group== "ANN_STEX") %>% nrow() + data %>% filter(group== "ARG") %>% nrow() +
                            data %>% filter(group== "DEB") %>% nrow()/2)
pop_labels <- c("ANN_NON-TEX","ANN_NTEX","ANN_STEX","ARG","DEB")
pop_labels_tibble <- tibble(pop=pop_labels,position=pop_labels_positions)

admixture_plot <- data %>% 
  gather(Q, value, Q1:Q5) %>%
  ggplot(aes(x = fct_reorder(sample, group), y = value, fill = factor(Q)))+
  geom_bar(stat = "identity", width = 1.1) +
  theme_cowplot() +
  theme(axis.text.x=element_blank()) +
  scale_fill_manual(values=five_colors[c(1,2,4,3,5)],name="Group") +
  #ggtitle("Admixture results") +
  geom_vline(xintercept = vertical_lines+0.5) +
  annotate("text",x=pop_labels_tibble$position[1],y=-0.07,label=pop_labels_tibble$pop[1],size=3) +
  annotate("text",x=pop_labels_tibble$position[2],y=-0.07,label=pop_labels_tibble$pop[2],size=3) +
  annotate("text",x=pop_labels_tibble$position[3],y=-0.07,label=pop_labels_tibble$pop[3],size=3) +
  annotate("text",x=pop_labels_tibble$position[4],y=-0.07,label=pop_labels_tibble$pop[4],size=3) +
  annotate("text",x=pop_labels_tibble$position[5]+1,y=-0.07,label=pop_labels_tibble$pop[5],size=3) +
  theme(legend.position="none",axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y = element_text(vjust = -5)) +
  coord_cartesian(ylim=c(-.11,1)) +
  ylab("Ancestry") 
  #ggtitle("ADMIXTURE")


####Load major structure results and plot them.
structure_directory <- "/home/owens/working/texanus/WGS/structure/clumpak/K=5/MajorCluster/CLUMPP.files"

five_colors <- pnw_palette(name="Bay",n=5,type="discrete")
Q <- read_delim(paste(structure_directory, "/ClumppIndFile.output",sep=""),col_names=c("A1","A2","A3","A4","A5",  "Q1","Q2","Q3","Q4","Q5"),
                delim=" ") %>% select(Q1:Q5)
samples <- read_tsv(paste(directory, "/Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.admixr.annargdeb.groups.txt",sep=""),col_names = c("sample","group"))

data <- cbind(samples,Q)


vertical_lines <- c(data %>% filter(group== "ANN_CEN") %>% nrow(),
                    data %>% filter(group== "ANN_NTEX" | group == "ANN_CEN") %>% nrow(),
                    data %>% filter(group== "ANN_NTEX" | group == "ANN_CEN" | group== "ANN_STEX") %>% nrow(),
                    data %>% filter(group!= "DEB") %>% nrow())
pop_labels_positions <- c(data %>% filter(group== "ANN_CEN") %>% nrow()/2,
                          data %>% filter(group== "ANN_CEN") %>% nrow() + data %>% filter(group== "ANN_NTEX") %>% nrow()/2,
                          data %>% filter(group== "ANN_CEN") %>% nrow() + data %>% filter(group== "ANN_NTEX") %>% nrow() +
                            data %>% filter(group== "ANN_STEX") %>% nrow()/2,
                          data %>% filter(group== "ANN_CEN") %>% nrow() + data %>% filter(group== "ANN_NTEX") %>% nrow() +
                            data %>% filter(group== "ANN_STEX") %>% nrow() + data %>% filter(group== "ARG") %>% nrow()/2,
                          data %>% filter(group== "ANN_CEN") %>% nrow() + data %>% filter(group== "ANN_NTEX") %>% nrow() +
                            data %>% filter(group== "ANN_STEX") %>% nrow() + data %>% filter(group== "ARG") %>% nrow() +
                            data %>% filter(group== "DEB") %>% nrow()/2)
pop_labels <- c("ANN_NON-TEX","ANN_NTEX","ANN_STEX","ARG","DEB")
pop_labels_tibble <- tibble(pop=pop_labels,position=pop_labels_positions)

structure_major <- data %>% 
  gather(Q, value, Q1:Q5) %>%
  ggplot(aes(x = fct_reorder(sample, group), y = value, fill = factor(Q)))+
  geom_bar(stat = "identity", width = 1.1) +
  theme_cowplot() +
  theme(axis.text.x=element_blank()) +
  scale_fill_manual(values=five_colors[c(2,4,3,5,1)],name="Group") +
  #ggtitle("Admixture results") +
  geom_vline(xintercept = vertical_lines+0.5) +
  annotate("text",x=pop_labels_tibble$position[1],y=-0.07,label=pop_labels_tibble$pop[1],size=3) +
  annotate("text",x=pop_labels_tibble$position[2],y=-0.07,label=pop_labels_tibble$pop[2],size=3) +
  annotate("text",x=pop_labels_tibble$position[3],y=-0.07,label=pop_labels_tibble$pop[3],size=3) +
  annotate("text",x=pop_labels_tibble$position[4],y=-0.07,label=pop_labels_tibble$pop[4],size=3) +
  annotate("text",x=pop_labels_tibble$position[5]+1,y=-0.07,label=pop_labels_tibble$pop[5],size=3) +
  theme(legend.position="none",axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y = element_text(vjust = -5)) +
  coord_cartesian(ylim=c(-.11,1)) +
  ylab("Ancestry") 
  #ggtitle("STRUCTURE")


####Load minor structure results and plot them.
structure_directory <- "/home/owens/working/texanus/WGS/structure/clumpak/K=5/MinorCluster1/CLUMPP.files"

five_colors <- pnw_palette(name="Bay",n=5,type="discrete")

Q <- read_delim(paste(structure_directory, "/ClumppIndFile.output",sep=""),col_names=c("A1","A2","A3","A4","A5",  "Q1","Q2","Q3","Q4","Q5"),
                delim=" ") %>% select(Q1:Q5)
samples <- read_tsv(paste(directory, "/Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.admixr.annargdeb.groups.txt",sep=""),col_names = c("sample","group"))

data <- cbind(samples,Q)


vertical_lines <- c(data %>% filter(group== "ANN_CEN") %>% nrow(),
                    data %>% filter(group== "ANN_NTEX" | group == "ANN_CEN") %>% nrow(),
                    data %>% filter(group== "ANN_NTEX" | group == "ANN_CEN" | group== "ANN_STEX") %>% nrow(),
                    data %>% filter(group!= "DEB") %>% nrow())
pop_labels_positions <- c(data %>% filter(group== "ANN_CEN") %>% nrow()/2,
                          data %>% filter(group== "ANN_CEN") %>% nrow() + data %>% filter(group== "ANN_NTEX") %>% nrow()/2,
                          data %>% filter(group== "ANN_CEN") %>% nrow() + data %>% filter(group== "ANN_NTEX") %>% nrow() +
                            data %>% filter(group== "ANN_STEX") %>% nrow()/2,
                          data %>% filter(group== "ANN_CEN") %>% nrow() + data %>% filter(group== "ANN_NTEX") %>% nrow() +
                            data %>% filter(group== "ANN_STEX") %>% nrow() + data %>% filter(group== "ARG") %>% nrow()/2,
                          data %>% filter(group== "ANN_CEN") %>% nrow() + data %>% filter(group== "ANN_NTEX") %>% nrow() +
                            data %>% filter(group== "ANN_STEX") %>% nrow() + data %>% filter(group== "ARG") %>% nrow() +
                            data %>% filter(group== "DEB") %>% nrow()/2)
pop_labels <- c("ANN_NON-TEX","ANN_NTEX","ANN_STEX","ARG","DEB")
pop_labels_tibble <- tibble(pop=pop_labels,position=pop_labels_positions)

structure_minor <- data %>% 
  gather(Q, value, Q1:Q5) %>%
  ggplot(aes(x = fct_reorder(sample, group), y = value, fill = factor(Q)))+
  geom_bar(stat = "identity", width = 1.1) +
  theme_cowplot() +
  theme(axis.text.x=element_blank()) +
  scale_fill_manual(values=five_colors[c(2,1,3,4,5)],name="Group") +
  #ggtitle("Admixture results") +
  geom_vline(xintercept = vertical_lines+0.5) +
  annotate("text",x=pop_labels_tibble$position[1],y=-0.07,label=pop_labels_tibble$pop[1],size=3) +
  annotate("text",x=pop_labels_tibble$position[2],y=-0.07,label=pop_labels_tibble$pop[2],size=3) +
  annotate("text",x=pop_labels_tibble$position[3],y=-0.07,label=pop_labels_tibble$pop[3],size=3) +
  annotate("text",x=pop_labels_tibble$position[4],y=-0.07,label=pop_labels_tibble$pop[4],size=3) +
  annotate("text",x=pop_labels_tibble$position[5]+1,y=-0.07,label=pop_labels_tibble$pop[5],size=3) +
  theme(legend.position="none",axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y = element_text(vjust = -5)) +
  coord_cartesian(ylim=c(-.11,1)) +
  ylab("Ancestry")  
  #ggtitle("STRUCTURE")


#######
#Make map plots of D
#######

#Mapping data
usa <- map_data('state')
states <- map_data("state")
target_state <- map_data('state')
lat_range <- c(25, 35)
long_range <- c(-105,-93)


dstat_data %>%
  filter(group == "STex" ) %>%
  filter(dataset=="full_genome") %>%
  filter(Y == "ARG") %>%
  group_by(population,lat,long) %>%
  summarise(mean_D = mean(D)) -> plotting_arg

dstat_data %>%
  filter(group == "STex" ) %>%
  filter(dataset=="full_genome") %>%
  filter(Y == "DEB") %>%
  group_by(population,lat,long) %>%
  summarise(mean_D = mean(D)) -> plotting_deb

map_plot <- ggplot(target_state, aes(long, lat)) +
  geom_map(map=target_state, aes(map_id=region), fill=NA, color="black") +
  coord_quickmap() +
  geom_point(data=plotting_arg,
             aes(x=long, y=lat, fill=mean_D),
             size=5,shape=21,color="white") +
  xlab("Longitude") + ylab("Latitude") +
  scale_x_continuous(limits = long_range, expand = c(0, 0)) +
  scale_y_continuous(limits = lat_range, expand = c(0, 0)) +
  scale_fill_viridis_c(name="Mean\nH. argophyllus D") +
  scale_color_viridis_c(option="inferno",name="mean_") +
  geom_label_repel(data=plotting_arg,aes(x=long,y=lat,label=population),
                  box.padding = 1,size=2.5) +
  theme_cowplot()

ggplot(target_state, aes(long, lat)) +
  geom_map(map=target_state, aes(map_id=region), fill=NA, color="black") +
  coord_quickmap() +
  geom_point(data=plotting_deb,
             aes(x=long, y=lat, color=mean_D),
             size=1) +
  scale_fill_brewer(name="Cluster",palette = "Set1") +theme_bw() +
  xlab("Longitude") + ylab("Latitude") +
  scale_x_continuous(limits = long_range, expand = c(0, 0)) +
  scale_y_continuous(limits = lat_range, expand = c(0, 0)) +
  scale_color_viridis_c(option="inferno",name="mean_") +
  theme_cowplot()




full_genome_plots <- deb_genome_wide / arg_genome_wide
structure_plots <- admixture_plot /structure_major / structure_minor

range_plot <- file.path("/home/owens/bin/texanus/figures/Species_ranges.png")
abba_plot <- file.path("/home/owens/bin/texanus/figures/ABBA-BABA.png")

range_plot <- ggplot() + 
  draw_image(range_plot) +
  theme_void()

abba_plot <- ggplot() + 
  draw_image(abba_plot) +
  theme_void()



pdf("figures/Texanus_all_admixture.v3.pdf",height=12,width=9,useDingbats = FALSE)
(ind_d_organized / (abba_plot | range_plot) / (map_plot | ind_d_correlation) / structure_plots ) +
  plot_annotation(tag_levels = 'A') +
  plot_layout(heights = c(2,2,2,2))
dev.off()

pdf("figures/Texanus_all_manhattan.v2.pdf",height=8,width=7)
full_genome_plots / ann1363_genome_wide
dev.off()

