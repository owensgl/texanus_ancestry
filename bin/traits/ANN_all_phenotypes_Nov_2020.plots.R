#Plotting phenotypic measures of texanus
library(tidyverse)
library(grid)
library(gridExtra)
library(forcats)
library(missMDA)
library(FactoMineR)
library(ggpubr)
library(PNWColors)
library(patchwork)
labels <- read_tsv("/home/owens/working/sample_info_apr_2018.tsv") %>%
  rename(sample = name) 
expectations <- read_tsv("texanus_phenotype_expectations.txt")

#Compare specific phenotypes with different expectations
expectations %>% filter(!is.na(tex_expectations)) -> differing_phenotypes

pop_loc <- read_tsv("pop_loc_allnum.txt") %>%
  rename(population = pop)
labels %>% inner_join(.,pop_loc) -> labels
labels %>%
  mutate(group = case_when(lat < 30 ~ "STex",
                           lat < 33 & long > -105 ~ "NTex",
                           TRUE ~ "Non-Tex")) -> labels

phenotype <- read_tsv("../wild_gwas_2018/resources/annuus_all_phenotypes_May2019.txt") %>% 
  rename(sample = FID) %>% select(-IID) %>%
  select(-Guides_3_petals,-Guides_3_petals_mod,-Guides_individual_petal,-Guides_individual_petal_mod,
         -Dark_petal)

imputed_phenotype <- imputePCA(as.matrix(phenotype[,c(2:ncol(phenotype))]),method="Regularized")
pca <- PCA(imputed_phenotype,ncp=20)

pca_results <- as_tibble(pca$ind$coord)
colnames(pca_results) <- paste0("PC",1:20)
pca_results$sample <- phenotype$sample
pca_colors <- pnw_palette(name="Bay",n=3,type="discrete")
plot1 <- pca_results %>%
  inner_join(.,labels) %>%
  ggplot(.,aes(x=PC1,y=PC2,group=group,color=group)) + geom_point(aes(color=group),alpha=0.8) +
  theme_bw() +
  stat_ellipse() +
  scale_color_manual(values=pca_colors,name="Region",
                     labels=c("Non-Texas","North Texas","South Texas")) +
  theme(legend.position = "bottom")
plot2 <- pca_results %>%
  inner_join(.,labels) %>%
  ggplot(.,aes(x=PC1)) + geom_density(aes(fill=group),alpha=0.5) +
  theme_bw() +
  scale_fill_manual(values=pca_colors,name="Region",
                     labels=c("Non-Texas","North Texas","South Texas")) +
  ylab("Density") +
  theme(legend.position = "bottom")

highlight_phenotype_plots <- list()
for (i in c(3,5,10,14)){
  chosen_variable <- differing_phenotypes$variable[i]
  printed_variable <- gsub("_"," ",chosen_variable)
  if (printed_variable == "Distance_of_first_branching_from_ground"){
    printed_variable <- "Height_of_first_branch"
  }
  my_comparisons <- list( c("Non-Tex", "NTex"), c("Non-Tex", "STex") )
  max_y =  phenotype %>%
    rename(variable= chosen_variable) %>%
    select(variable,sample) %>% pull(variable) %>% max(na.rm=T) *1.2
  min_y = phenotype %>%
    rename(variable= chosen_variable) %>%
    select(variable,sample) %>% pull(variable) %>% min(na.rm=T) 
  p <- ggboxplot(  phenotype %>%
                     rename(variable= chosen_variable) %>%
                     select(variable,sample) %>%
                     inner_join(labels), x = "group", y = "variable",
                   color = "group", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
                   add = "jitter", shape = "group",outlier.shape = NA,
                   ylim=c(min_y,max_y)) +
    stat_compare_means(comparisons = my_comparisons,method="wilcox")+ # Add pairwise comparisons p-value
    ylab(paste0(printed_variable)) +
    scale_color_manual(values=pca_colors) +
    xlab("Region") +
    theme(legend.position = "none",
          text = element_text(size=10))
  highlight_phenotype_plots[[i]]  <- p
  
}
patch1 <- plot1 / plot2
patch2 <- (highlight_phenotype_plots[[3]] | highlight_phenotype_plots[[5]] ) / (highlight_phenotype_plots[[10]] | highlight_phenotype_plots[[14]] )
patch_together <- patch1 / patch2
  
pdf("figures/Texanus_phenotype_pca_v2.pdf",height=12,width=7, useDingbats = FALSE)
plot1 / plot2 / (highlight_phenotype_plots[[3]] | highlight_phenotype_plots[[5]] ) / 
  (highlight_phenotype_plots[[10]] | highlight_phenotype_plots[[14]] ) +
  plot_annotation(tag_levels = 'A')
dev.off()


#data.hcpc <- HCPC(pca, nb.clust=0, conso=0, min=2, max=10)



pca_cont <- as_tibble(pca$var$contrib[1:length(colnames(phenotype)[-1]),])
colnames(pca_cont) <- paste0("PC",1:20)

pca_cont$variable <- colnames(phenotype)[-1]

pca_cor <- as_tibble(pca$var$cor[1:length(colnames(phenotype)[-1]),])
colnames(pca_cor) <- paste0("PC",1:20)
pca_cor$variable <- colnames(phenotype)[-1]

pca_cor %>% select(variable,PC1) %>%
  mutate(relationship = case_when(PC1 > 0 ~ "pos",
                                  TRUE ~ "neg")) %>%
  select(-PC1) -> pca_cor

pca_cor <- pca_cor %>% inner_join(.,pca_cont %>% select(variable,PC1))

expectations <- read_tsv("texanus_phenotype_expectations.txt")

#Compare specific phenotypes with different expectations
expectations %>% filter(!is.na(tex_expectations)) -> differing_phenotypes

phenotype_plots <- list()
for (i in 1:nrow(differing_phenotypes)){
  chosen_variable <- differing_phenotypes$variable[i]
  printed_variable <- chosen_variable
  if (printed_variable == "Distance_of_first_branching_from_ground"){
    printed_variable <- "Height_of_first_branch"
  }
  my_comparisons <- list( c("Non-Tex", "NTex"), c("Non-Tex", "STex") )
  max_y =  phenotype %>%
    rename(variable= chosen_variable) %>%
    select(variable,sample) %>% pull(variable) %>% max(na.rm=T) *1.2
  min_y = phenotype %>%
    rename(variable= chosen_variable) %>%
    select(variable,sample) %>% pull(variable) %>% min(na.rm=T) 
  p <- ggboxplot(  phenotype %>%
                     rename(variable= chosen_variable) %>%
                     select(variable,sample) %>%
                     inner_join(labels), x = "group", y = "variable",
                 color = "group", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
                 add = "jitter", shape = "group",outlier.shape = NA,
                 ylim=c(min_y,max_y)) +
    stat_compare_means(comparisons = my_comparisons,method="wilcox")+ # Add pairwise comparisons p-value
    ylab(paste0(printed_variable," (",differing_phenotypes$tex_expectations[i],")")) +
    scale_color_manual(values=pca_colors) +
    xlab("Region") +
    theme(legend.position = "none",
    text = element_text(size=10))
  phenotype_plots[[i]]  <- p
    
}
pdf("figures/Texanus_all_phenotypes_v2.pdf",height=12,width=12,useDingbats = FALSE)
wrap_plots(phenotype_plots[[1]],phenotype_plots[[2]],phenotype_plots[[3]],
           phenotype_plots[[4]],phenotype_plots[[5]],phenotype_plots[[6]],
           phenotype_plots[[7]],phenotype_plots[[8]],phenotype_plots[[9]],
           phenotype_plots[[10]],phenotype_plots[[11]],phenotype_plots[[12]],
           phenotype_plots[[13]],phenotype_plots[[14]],phenotype_plots[[15]],
           phenotype_plots[[16]],phenotype_plots[[17]],phenotype_plots[[18]],
           phenotype_plots[[19]],phenotype_plots[[20]]) +plot_annotation(tag_levels = 'A')

dev.off()
