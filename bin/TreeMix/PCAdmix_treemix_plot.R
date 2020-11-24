#PCAdmix and treemix trees:
#Plotting treemix likelihood tests
library(tidyverse)
library(ggpubr)
library(ggforce)
library(cowplot)
library(PNWColors)
library(patchwork)
library(ggrastr)

pallete_3 <- pnw_palette("Bay",3,type="discrete")
source("treemix_plotting_functions.R")
source("treemix_plotting_tidy_functions.R")

species <- c("arg","deb","petfal")
all_data <- tibble(species=character(),mig=numeric(),dataset=character(),likelihood=numeric())
#Add zero migration tests
pruned_zeromig <- read_delim("/home/owens/working/texanus/WGS/treemix/remapped_alltest/Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.0p2ld.treemix.llik",
                             ,delim=" ")[7] %>% pull()
full_zeromig <- read_delim("/home/owens/working/texanus/WGS/treemix/remapped_alltest/Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.treemix.llik",
                           ,delim=" ")[7] %>% pull()

for (chosen_species in species){
  
  
  directory <- paste0("/home/owens/working/texanus/WGS/treemix/",chosen_species,"_introgression_test")
  pruned <- list.files(directory, "Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.0p2ld.treemix")[
    grepl("llik", list.files(paste0("/home/owens/working/texanus/WGS/treemix/",chosen_species,"_introgression_test"), 
                             "Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.0p2ld.treemix"))]
  
  for (file in pruned){
    mig <- gsub("Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.0p2ld.treemix.","",file)
    mig <- gsub(".llik","",mig)
    likelihood <- read_delim(paste0(directory,"/",file),delim=" ")[7] %>% pull()
    tmp <- tibble(species=chosen_species,mig=mig,dataset="pruned",likelihood=likelihood)
    all_data <- rbind(all_data,tmp)
  }
  full <- list.files(directory, "Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.treemix")[
    grepl("llik", list.files(paste0("/home/owens/working/texanus/WGS/treemix/",chosen_species,"_introgression_test"), 
                             "Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.treemix"))]
  
  for (file in full){
    mig <- gsub("Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.treemix.","",file)
    mig <- gsub(".llik","",mig)
    likelihood <- read_delim(paste0(directory,"/",file),delim=" ")[7] %>% pull()
    tmp <- tibble(species=chosen_species,mig=mig,dataset="full",likelihood=likelihood)
    all_data <- rbind(all_data,tmp)
    
  }
  tmp <- tibble(species=chosen_species,mig=0,dataset="pruned",likelihood=pruned_zeromig)
  all_data <- rbind(all_data,tmp)
  
  tmp <- tibble(species=chosen_species,mig=0,dataset="full",likelihood=full_zeromig)
  all_data <- rbind(all_data,tmp)
  
  
}


treemix_test <- all_data %>%
  filter(dataset == "pruned") %>%
  filter(species != "petfal") %>%
  group_by(species,dataset) %>%
  mutate(max_likelihood = max(likelihood)) %>%
  mutate(rel_likelihood = max_likelihood - likelihood) %>%
  mutate(mig = as.numeric(mig)) %>%
  mutate(start = case_when(mig <= 0.01 ~ "yes",
                           TRUE ~ "no")) %>%
  ungroup() %>%
  mutate(species = case_when(species == "arg"~ "ARG",
                             species == "deb" ~ "DEB")) %>%
  ggplot(.,aes(x=mig,y=-rel_likelihood,color=species)) + 
  theme_cowplot() + 
  geom_line() + 
  facet_wrap(~species,scales="free_y") +
  ylab("Likelihood difference") + 
  xlab("Migration weight") +
  scale_color_manual(values=pallete_3[c(2,3)]) +
  theme(legend.position="none",
        axis.text.x = element_text(size=8))


#Trying to plot residuals
covariance <- read_delim(delim=" ", "/home/owens/working/texanus/WGS/treemix/remapped_alltest/Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.0p2ld.treemix.cov.gz",
                         col_names = c("species_1", "ANN_CEN", "ANN_NTEX", "ANN_STEX", "ARG", "DEB", "OUT", "PET_CAN", "PET_FAL", "PET_PET"),
                         skip=1) %>%
  pivot_longer(-species_1, names_to = "species_2", values_to = "covariance")
cov_se <- read_delim(delim=" ", "/home/owens/working/texanus/WGS/treemix/remapped_alltest/Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.0p2ld.treemix.covse.gz",
                     col_names = c("species_1", "ANN_CEN", "ANN_NTEX", "ANN_STEX", "ARG", "DEB", "OUT", "PET_CAN", "PET_FAL", "PET_PET"),
                     skip=1) %>%
  pivot_longer(-species_1, names_to = "species_2", values_to = "covse")

model <- read_delim(delim=" ", "/home/owens/working/texanus/WGS/treemix/remapped_alltest/Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.0p2ld.treemix.modelcov.gz") 
colnames(model) <- c("species_1",colnames(model)[1:9])
model <- model %>% pivot_longer(-species_1, names_to = "species_2", values_to = "modelcov")

population_order <- c("ARG","ANN_STEX","ANN_NTEX", "ANN_NON-TEX","PET_PET","PET_FAL","NIV","DEB","OUT")
treemix_matrix <- inner_join(covariance, cov_se) %>% inner_join(.,model) %>%
  mutate(ave_se = mean(covse)) %>%
  mutate(scaled_res = (covariance - modelcov)/ave_se) %>% 
  mutate(species_1 = case_when(species_1 == "PET_CAN" ~ "NIV",
                               species_1 == "ANN_CEN" ~ "ANN_NON-TEX",
                               TRUE ~ species_1),
         species_2 = case_when(species_2 == "PET_CAN" ~ "NIV",
                               species_2 == "ANN_CEN" ~ "ANN_NON-TEX",
                               TRUE ~ species_2)) %>%
  mutate(n1 = as.numeric(fct_relevel(species_1, rev(population_order)))) %>%
  mutate(n2 = as.numeric(fct_relevel(species_2, rev(population_order)))) %>%
  filter(n2 > n1) %>%

  ggplot(.,aes(x=fct_relevel(species_1, rev(population_order)),y=fct_relevel(species_2, (population_order)),fill=scaled_res)) + 
  theme_cowplot() +
  geom_tile() +
  scale_fill_gradient2(high="#dd4124",low="#00496f",mid="white",name="Scaled\nresidual fit") +
  theme(axis.text.x = element_text(size=8,angle = 60, hjust = 1),
        axis.text.y = element_text(size=8),
        legend.position = "left",
        text = element_text(size = 8)) +
  ylab("Population") + xlab("Species") 



# Plotted is the residual fit from the maximum likelihood tree in A. 
# We divided the residual covariance between each pair of populations  and  
# by the average standard error across all pairs. 
# We then plot in each cell  this scaled residual. 



#Load up treemix_plotting_tidy_functions.R
treemix_5 <- read_treemix("/home/owens/working/texanus/WGS/treemix/remapped_alltest/Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.0p2ld.treemix.5")
treemix_tree <- plot_treemix(treemix_5) + ylab("") + xlab("Drift parameter") +
  theme_void() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.text=element_blank(),
        text = element_text(size = 8),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_colour_gradient2("Migration\nweight", high = "red", mid = "yellow") +
  scale_x_continuous(expand = c(0, .01)) 



###########################
#PCAdmix plotting
###########################
sample_info <- read_tsv("sample_info_apr_2018.tsv") %>%
  rename(sample=name) %>%
  select(sample,population)
#Load control data
data_directory <- "/home/owens/working/texanus/WGS/PCAdmix_control/"
files <- list.files(data_directory,"bed")
data <- data_frame(filename = files) %>% # create a data frame
  # holding the file names
  mutate(file_contents = map(filename,          # read files into
                             ~ read_tsv(file.path(data_directory, .),
                                        col_types = cols(Chr = col_character()))) # a new data column
  )  
all_data_control <- unnest(data, cols = c(file_contents))
all_data_control %>%
  separate(filename,c("chr","info"),"_") %>% 
  separate(info,c("sample","strand")) %>% 
  select(-chr) %>%
  inner_join(sample_info) %>%
  rename(start_bp=`Start(bp)`,
         end_bp = `End(bp)`,
         start_cm= `Start(cM)`,
         end_cm=`End(cM)`) %>%
  mutate(length_bp = end_bp - start_bp,
         length_cm = end_cm - start_cm) -> all_data_control

all_data_control %>%
  filter(Confidence > 0.9) %>%
  group_by(sample,Vit) %>%
  summarize(total_bp = sum(length_bp),
            total_cm= sum(length_cm)) %>%
  group_by(sample) %>%
  mutate(percent_bp = total_bp / sum(total_bp),
         percent_cm = total_cm / sum(total_cm)) %>%
  filter(Vit != "ANN") %>%
  mutate(population = "Central" ) -> all_data_control_summarized

data_directory <- "/home/owens/working/texanus/WGS/PCAdmix/"
files <- list.files(data_directory,"bed")
data <- data_frame(filename = files) %>% # create a data frame
  # holding the file names
  mutate(file_contents = map(filename,          # read files into
                             ~ read_tsv(file.path(data_directory, .),
                                        col_types = cols(Chr = col_character()))) # a new data column
  )  
all_data <- unnest(data, cols = c(file_contents))

all_data %>%
  separate(filename,c("chr","info"),"_") %>% 
  separate(info,c("sample","strand")) %>% 
  select(-chr) %>%
  inner_join(sample_info) %>%
  rename(start_bp=`Start(bp)`,
         end_bp = `End(bp)`,
         start_cm= `Start(cM)`,
         end_cm=`End(cM)`) %>%
  mutate(length_bp = end_bp - start_bp,
         length_cm = end_cm - start_cm) -> all_data

pcadmix_boxplot <- all_data %>%
  filter(Confidence > 0.9) %>%
  group_by(sample,population,Vit) %>%
  summarize(total_bp = sum(length_bp),
            total_cm= sum(length_cm)) %>%
  group_by(sample,population) %>%
  mutate(percent_bp = total_bp / sum(total_bp),
         percent_cm = total_cm / sum(total_cm)) %>%
  filter(Vit != "ANN") %>%
  rbind(.,all_data_control_summarized) %>%
  ungroup() %>%
  mutate(population = case_when(population == "Central" ~ "NON-TEX", 
                                TRUE ~ population)) %>%
  ggplot(.,aes(x=population,y=percent_bp*100,fill=Vit)) + geom_boxplot() +
  theme_cowplot() +
  scale_fill_manual(values=pallete_3[2:3],name="Donor") +
  facet_wrap(~Vit,scale="free_y") +
  ylab("Percent introgression") + xlab("Population") +
  theme(axis.text.x=element_text(angle=60, hjust=1,size=6),
        legend.position = "none") +
  geom_hline(yintercept = 0,linetype="dotted") 



all_data %>% 
  group_by(Chr,sample,strand) %>%
  mutate(joined_start_bp = lag(end_bp),
         joined_start_cm = lag(end_cm),
         start_bp = case_when(is.na(joined_start_bp) ~ start_bp,
                              TRUE ~ joined_start_bp),
         start_cm = case_when(is.na(joined_start_cm) ~ start_cm,
                              TRUE ~ joined_start_cm)) %>% 
  group_by(Chr,start_cm,end_cm,start_bp,end_bp,Vit, .drop=F) %>%
  summarise(count=n()) %>% 
  group_by(Chr,start_cm,end_cm,start_bp,end_bp) %>%
  mutate(ANN_count = sum(count[Vit=="ANN"]),
         ANN_ARG_count = sum(count[Vit=="ARG"], count[Vit=="ANN"]),
         Total_count = sum(count)) %>%
  mutate(y_min = case_when(Vit == "ANN" ~ as.integer(0),
                           Vit == "ARG" ~ ANN_count,
                           Vit == "DEB" ~ ANN_ARG_count),
         y_max = case_when(Vit == "ANN" ~ ANN_count,
                           Vit == "ARG" ~ ANN_ARG_count,
                           Vit == "DEB" ~ Total_count)) ->  all_data_for_rect


pallete_3 <- pnw_palette("Bay",3,type="discrete")
chr_lengths <- read_tsv("../wild_gwas_2018/Ha412HO.chrlengths.txt") %>% mutate(Chr = gsub("Ha412HOChr","",chr)) %>%
  select(-chr)

#Make it cumulative so that it can be one long plot:
all_data_for_rect_cum <- chr_lengths %>%
  select(Chr,end) %>%
  rename(chr_len = end) %>%
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(all_data_for_rect, ., by=c("Chr"="Chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(Chr, start_bp) %>%
  mutate( start_bp_cum=start_bp+tot,end_bp_cum=end_bp+tot)

axisdf = all_data_for_rect_cum %>% group_by(Chr) %>% summarize(center=( max(end_bp_cum) + min(start_bp_cum) ) / 2,
                                                               max= max(end_bp_cum),
                                                               min=min(start_bp_cum))

pcadmix_genome <- all_data_for_rect_cum%>%
  ggplot(.,aes(xmin=start_bp_cum,xmax=end_bp_cum,ymin=y_min/218,ymax=y_max/218, fill=Vit)) + 
  geom_rect() +
  theme_cowplot() +
  scale_x_continuous( label = gsub("Ha412HOChr","",axisdf$Chr), breaks= axisdf$center,
                      expand = c(0.01,0.01)) +
  xlab("Chromosome") + ylab("Ancestry") +
  scale_fill_manual(values=pallete_3,name="Donor") +
  geom_vline(xintercept = axisdf$max[1:nrow(axisdf)-1],linetype="solid") +
  theme(legend.position="bottom",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) 

all_data_for_rect_cum%>%
  filter(Chr == "09") %>%
  filter(Vit == "ARG", count > 100) %>% View()
  ggplot(.,aes(xmin=start_bp,xmax=end_bp,ymin=y_min/218,ymax=y_max/218, fill=Vit)) + 
  geom_rect() +
  theme_cowplot() +
  xlab("Chromosome") + ylab("Ancestry") +
  scale_fill_manual(values=pallete_3,name="Donor") +
  theme(legend.position="bottom",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) 


################
#Plotting PCAdmix against recombination rate
################


#Stupid way of adding zero points
all_data_for_rect %>%
  mutate(bp_length = end_bp - start_bp,
         cm_length = end_cm - start_cm,
         rate = cm_length/bp_length) %>%
  ungroup() %>%
  select(Chr,start_bp,rate,Vit,count) -> all_data_rates

all_data_rates %>% select(Chr,start_bp,rate) -> regions

for (i in 1:nrow(regions)){
  print(i)
  chosen_chr <- regions$Chr[i]
  chosen_start_bp <- regions$start_bp[i]
  chosen_rate <- regions$rate[i]
  #Check if it exists in the rate file
  
  arg_present <- all_data_rates %>% filter(Chr == chosen_chr, start_bp == chosen_start_bp) %>%
    filter(Vit == "ARG") %>% nrow()
  deb_present <- all_data_rates %>% filter(Chr == chosen_chr, start_bp == chosen_start_bp) %>%
    filter(Vit == "DEB") %>% nrow()
  if (arg_present == 0){
    all_data_rates <- rbind(all_data_rates,tibble(Chr = chosen_chr,
                                                  start_bp = chosen_start_bp,
                                                  rate = chosen_rate,
                                                  Vit = "ARG",
                                                  count = 0))
  }
  if (deb_present == 0){
    all_data_rates <- rbind(all_data_rates,tibble(Chr = chosen_chr,
                                                  start_bp = chosen_start_bp,
                                                  rate = chosen_rate,
                                                  Vit = "DEB",
                                                  count = 0))
  }
}
all_data_rates %>%
  group_by(Chr,start_bp,rate) %>%
  tally()
pcadmix_recom_plot <- all_data_rates %>%
  filter(Vit != "ANN") %>%
  mutate(cm_rank=100*(rank(rate)/length(rate))) %>%
  mutate(cm_quantile = floor(cm_rank/20)+1) %>%
  mutate(cm_quantile = case_when(cm_quantile > 5 ~ 5,
                                 TRUE ~ cm_quantile)) %>%
  ggplot(.,aes(x=as.factor(cm_quantile),y=count/218,color=Vit)) +
  theme_cowplot() +
  geom_jitter(data=all_data_rates %>%
                filter(Vit != "ANN") %>%
                mutate(cm_rank=100*(rank(rate)/length(rate))) %>%
                mutate(cm_quantile = floor(cm_rank/20)+1) %>%
                mutate(cm_quantile = case_when(cm_quantile > 5 ~ 5,
                                               TRUE ~ cm_quantile)) %>%
                sample_frac(0.1), alpha=0.2) +
  geom_boxplot(outlier.shape=NA,fill=NA,color="black") + 
  facet_wrap(~Vit) +
  ylab("Introgression frequency") +
  scale_color_manual(values=pallete_3[2:3],name="Donor") +
  scale_x_discrete(labels=c("0-20", "20-40", "40-60","60-80","80-100")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.x = element_text(size=12),
        legend.position = "none") +
  xlab("Recombination rate percentile")  +
  coord_cartesian(ylim=c(0,0.4))



############################
#Dstats with recombination rate
############################

genetic_map <- read_tsv("/home/owens/ref/Ha412HOv2.0-20181130.Nov22k22.geneticmap.extradivisions.txt") %>% 
  mutate(rate = lead(cM)-cM) %>% 
  rename(start = pos) %>% mutate(end = start + 999999) %>%
  mutate(rate = case_when(is.na(rate) ~ cM - lag(cM),
                          rate < 0 ~ cM - lag(cM),
                          TRUE ~ rate))
window_snps <- 100

chr_lengths <- read_tsv("../wild_gwas_2018/Ha412HO.chrlengths.txt")

dstats_deb <- read_tsv("/home/owens/working/texanus/WGS/dstat/Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.ANN_CEN.ANN_STEX.DEB.Dstats.txt.gz")

dstats_deb %>%
  group_by(chr) %>%
  mutate(window = floor(row_number()/window_snps)) %>%
  group_by(chr,window) %>%
  summarize(middle= mean(pos),
            all_abba = sum(abba),
            all_baba = sum(baba),
            start= min(pos),end=max(pos),
            D = (sum(abba) - sum(baba)) / (sum(abba) + sum(baba)),
            fd = (sum(abba) - sum(baba)) / (sum(p3_abba) - sum(p3_baba))) %>% 
  mutate(fd = case_when(D < 0 ~ 0,
                        TRUE ~ fd)) %>%
  ungroup() %>%
  mutate(rank = percent_rank(fd),
         outlier = case_when(rank >0.98 ~ 1,
                             TRUE ~ 0)) -> dstats_deb_windows

#Add recombination rate to dstats_windows
dstats_deb_windows$cm_rate <- NA
for (i in 1:nrow(dstats_deb_windows)){
  print(i)
  chosen_start <- dstats_deb_windows$start[i]
  chosen_end <-dstats_deb_windows$end[i]
  chosen_chr <- dstats_deb_windows$chr[i]
  start_rate <- genetic_map %>% filter(chr == chosen_chr) %>% 
    filter(start <= chosen_start, end >= chosen_start) %>%
    pull(rate)
  end_rate <- genetic_map %>% filter(chr == chosen_chr) %>% 
    filter(start <= chosen_end, end >= chosen_end) %>%
    pull(rate)
  mean_rate <- mean(start_rate,end_rate)
  dstats_deb_windows$cm_rate[i] <- mean_rate
}
dstats_deb_windows$comparison <- "DEB"


dstats_arg <- read_tsv("/home/owens/working/texanus/WGS/dstat/Texanus.tranche90.snp.gwas.90.bi.remappedHa412HO.ANN_CEN.ANN_STEX.ARG.Dstats.txt.gz")

dstats_arg %>%
  group_by(chr) %>%
  mutate(window = floor(row_number()/window_snps)) %>%
  group_by(chr,window) %>%
  summarize(middle= mean(pos),
            all_abba = sum(abba),
            all_baba = sum(baba),
            start= min(pos),end=max(pos),
            D = (sum(abba) - sum(baba)) / (sum(abba) + sum(baba)),
            fd = (sum(abba) - sum(baba)) / (sum(p3_abba) - sum(p3_baba))) %>% 
  mutate(fd = case_when(D < 0 ~ 0,
                        TRUE ~ fd)) %>%
  ungroup() %>%
  mutate(rank = percent_rank(fd),
         outlier = case_when(rank >0.98 ~ 1,
                             TRUE ~ 0)) -> dstats_arg_windows

#Add recombination rate to dstats_windows
dstats_arg_windows$cm_rate <- NA
for (i in 1:nrow(dstats_arg_windows)){
  print(i)
  chosen_start <- dstats_arg_windows$start[i]
  chosen_end <-dstats_arg_windows$end[i]
  chosen_chr <- dstats_arg_windows$chr[i]
  start_rate <- genetic_map %>% filter(chr == chosen_chr) %>% 
    filter(start <= chosen_start, end >= chosen_start) %>%
    pull(rate)
  end_rate <- genetic_map %>% filter(chr == chosen_chr) %>% 
    filter(start <= chosen_end, end >= chosen_end) %>%
    pull(rate)
  mean_rate <- mean(start_rate,end_rate)
  dstats_arg_windows$cm_rate[i] <- mean_rate
}
dstats_arg_windows$comparison <- "ARG"

dstat_recom_plot <- rbind(dstats_arg_windows, dstats_deb_windows) %>%
  mutate(cm_rank=100*(rank(cm_rate)/length(cm_rate))) %>%
  mutate(cm_quantile = floor(cm_rank/20)+1) %>%
  ggplot(.,aes(x=as.factor(cm_quantile),y=D)) + 
  geom_jitter(data=rbind(dstats_arg_windows, dstats_deb_windows) %>%
                mutate(cm_rank=100*(rank(cm_rate)/length(cm_rate))) %>%
                mutate(cm_quantile = floor(cm_rank/20)+1) %>%
                sample_frac(0.5), aes(color=comparison),alpha=0.1) + 
  theme_cowplot() +
  geom_boxplot(outlier.shape = NA,fill=NA) + facet_wrap(~comparison) +
  scale_color_manual(values=pallete_3[2:3],name="Donor") +
  scale_x_discrete(labels=c("0-20", "20-40", "40-60","60-80","80-100")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.x = element_text(size=12),
        legend.position = "none") +
  xlab("Recombination rate percentile") +
  geom_hline(yintercept=0,color="black",linetype="dashed")


############################
#Put it together
############################

#Download and convert it to PNG
tree_plot <- file.path("/home/owens/bin/texanus/figures/tmp.png")

tree_plot <- ggplot() + 
  draw_image(tree_plot) +
  theme_void()

pdf("figures/Texanus_all_pcadmix_treemix.v3.pdf",height=12,width=7, useDingbats = FALSE)
(((tree_plot | treemix_matrix )/ (treemix_test | pcadmix_boxplot)  / pcadmix_genome ) /
  (pcadmix_recom_plot | dstat_recom_plot )) +plot_annotation(tag_levels = 'A')
dev.off()

