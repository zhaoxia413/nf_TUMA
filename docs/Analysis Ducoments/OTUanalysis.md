---
sort: 2
---

-   [1 Requires](#requires)
-   [2 OTU and cor-microbiome
    analysis](#otu-and-cor-microbiome-analysis)
    -   [2.1 Description](#description)
    -   [2.2 Rank-Abundance](#rank-abundance)
        -   [2.2.1 Code](#code)
        -   [2.2.2 Figures](#figures)
        -   [2.2.3 Effective sequences](#effective-sequences)
        -   [2.2.4 Taxon annotation](#taxon-annotation)
    -   [2.3 Core Microbiome](#core-microbiome)
        -   [2.3.1 Code](#code-1)
        -   [2.3.2 Figures](#figures-1)
    -   [2.4 Microbiome compostion](#microbiome-compostion)
        -   [2.4.1 Code](#code-2)
        -   [2.4.2 Figures](#figures-2)

[`Return`](./)

1 Requires
==========

<details>
<summary>
<font size=4>Requires</font>
</summary>

    library(tidyverse)
    library(ggthemes)
    library(ggsci)
    library(ggpubr)
    library(survminer)
    library(survival)
    library(survivalROC)
    library(reshape2)
    library(data.table)
    library(ggExtra)
    library(cowplot)
    library(pheatmap)
    library(scico)
    library(colorspace)
    library(RColorBrewer)
    library(lubridate)
    library(tableone)
    library(kableExtra)
    library(BiodiversityR)
    library(reactable)
    source("../R_function/colors.R")
    source("../R_function/surv_plot.R")
    source("../R_function/OTUanalysis.R")
    theme_set(theme_cowplot())
    "%ni%" <- Negate("%in%")
    options(stringsAsFactors = F)

</details>

2 OTU and cor-microbiome analysis
=================================

2.1 Description
---------------

*OTU（Operational Taxonomic Units*

In 16S metagenomics approaches, OTUs are cluster of similar sequence
variants of the 16S rDNA marker gene sequence. Each of these cluster is
intended to represent a taxonomic unit of a bacteria species or genus
depending on the sequence similarity threshold. Typically, OTU cluster
are defined by a 97% identity threshold of the 16S gene sequences to
distinguish bacteria at the genus level.

*Limited taxonomic resolution*

OTU resolution depends on the 16S approach which has some limits in
distinguishing at the species level

2.2 Rank-Abundance
------------------

A rank abundance curve or Whittaker plot is a chart used by ecologists
to display relative species abundance, a component of biodiversity. It
can also be used to visualize species richness and species evenness. It
overcomes the shortcomings of biodiversity indices that cannot display
the relative role different variables played in their calculation.

### 2.2.1 Code

    otu<-fread("../Data/Data/OTUtable_ori.csv",data.table = F)
    meta<-fread("../Data/Data/meta.csv",data.table = F)
    otu<-t(data.frame(row.names = otu$OTU,otu[,-1]))
    otu_relative <- otu / rowSums(otu)
    rank_dat <- data.frame()
    for (i in rownames(otu_relative)) {
      rank_dat_i <- data.frame(rankabundance(subset(otu_relative, 
                                                    rownames(otu_relative) == i), 
                                             digits = 6))[1:2]
      rank_dat_i$sample <- i
      rank_dat <- rbind(rank_dat, rank_dat_i)
    }
    rank_dat <- subset(rank_dat, abundance != 0)
    colnames(rank_dat)[3]="Samples"
    rank_dat<-merge(rank_dat,meta,by="Samples")
    p<-ggplot(rank_dat, aes(rank, log(abundance, 10), color = Site)) +
      geom_line(size=0.2) +
      scale_colour_manual(limits = c('Saliva','Stool'), values = c('darkblue','orange')) +
      labs(x = 'OTUs rank', y = 'Relative adundance (%)', color = NULL) +
      theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.key = element_rect(fill = 'transparent')) +
      scale_y_continuous(breaks = 0:-5, labels = c('100', '10', '1', '0.1', '0.01', '0.001'), limits = c(-5, 0))

### 2.2.2 Figures

    p

<img src="OTUanalysis_files/figure-markdown_strict/unnamed-chunk-3-1.png" width="50%" style="display: block; margin: auto;" />

### 2.2.3 Effective sequences

    ggplot(meta,aes(reorder(Samples,Sequences),Sequences,fill=Site))+
      geom_col()+
      scale_fill_jama()+
      facet_grid(~Site,scales = "free",space = "free")+
      theme(axis.text.x = element_blank())+
      xlab("Samples")

<img src="OTUanalysis_files/figure-markdown_strict/unnamed-chunk-4-1.png" width="50%" style="display: block; margin: auto;" />

### 2.2.4 Taxon annotation

    OTUstool<-fread("../Data/Data/OTUtabale_regaStool.csv",data.table = F)
    OTUsaliva<-fread("../Data/Data/OTUtabale_regaSaliva.csv",data.table = F)
    message("Taxon number on stool samples")

    ## Taxon number on stool samples

    for (i in 1:7) {

      print(paste0(colnames(OTUstool)[i],":",
                     length(levels(factor(OTUstool[,i])))))
    }

    ## [1] "Phylum:11"
    ## [1] "Class:17"
    ## [1] "Order:40"
    ## [1] "Family:71"
    ## [1] "Genus:188"
    ## [1] "Species:433"
    ## [1] "OTU:1260"

    message("Taxon number on saliva samples")

    ## Taxon number on saliva samples

    for (i in 1:7) {
      
      print(paste0(colnames(OTUsaliva)[i],":",
                     length(levels(factor(OTUsaliva[,i])))))
    }

    ## [1] "Phylum:11"
    ## [1] "Class:18"
    ## [1] "Order:45"
    ## [1] "Family:77"
    ## [1] "Genus:151"
    ## [1] "Species:362"
    ## [1] "OTU:722"

2.3 Core Microbiome
-------------------

One of the aims of the Human Microbiome Project when established in 2007
was to identify a human ‘core microbiome’, defined as a group of
microbial taxa or genes that are shared by all or most humans (Hamady &
Knight, 2009; Turnbaugh et al., 2007). These pioneering studies found
that a universal taxonomic core rarely exists across groups of humans,
even at the scale of the family unit (Yatsunenko et al., 2012), yet most
shared the same set of core microbial genes. This pattern appears to be
similar for most host species studied to date, with variable microbial
composition often underpinned by similar gene context across individuals
(Burke, Steinberg, Rusch, Kjelleberg, & Thomas, 2011; Louca et al.,
2018). Nevertheless, most ecologists aim to understand host–microbe
interactions at the organismal level, accounting for microbe functional
traits where possible; therefore, within this field, the core microbiome
has been largely applied to taxonomically defined microbial communities
with an aim to identify groups of microbes that are particularly
widespread across the host population (e.g. Ainsworth et al., 2015;
Grieneisen, Livermore, Alberts, Tung, & Archie, 2017; Hernandez‐Agreda,
Gates, & Ainsworth, 2017; Kembel et al., 2014; Lundberg et al., 2012;
Muletz Wolz, Yarwood, Campbell Grant, Fleischer, & Lips, 2018).

ref:Applying the core microbiome to understand host–microbe systems.J
Anim Ecol. 2020;89:1549–1558 (2020).
<a href="https://doi.org/10.1111/1365-2656.13229" class="uri">https://doi.org/10.1111/1365-2656.13229</a>

ref: Turnbaugh, P., Hamady, M., Yatsunenko, T. et al. A core gut
microbiome in obese and lean twins. Nature 457, 480–484 (2009).
<a href="https://doi.org/10.1038/nature07540" class="uri">https://doi.org/10.1038/nature07540</a>

Within an individual oral cavity, we found over 3600 unique sequences,
over 500 different OTUs or “species-level” phylotypes (sequences that
clustered at 3% genetic difference) and 88 - 104 higher taxa (genus or
more inclusive taxon). The predominant taxa belonged to Firmicutes
(genus Streptococcus, family Veillonellaceae, genus Granulicatella),
Proteobacteria (genus Neisseria, Haemophilus), Actinobacteria (genus
Corynebacterium, Rothia, Actinomyces), Bacteroidetes (genus Prevotella,
Capnocytophaga, Porphyromonas) and Fusobacteria (genus Fusobacterium).

ref:Zaura, E., Keijser, B.J., Huse, S.M. et al. Defining the healthy
“core microbiome” of oral microbial communities. BMC Microbiol 9, 259
(2009).
<a href="https://doi.org/10.1186/1471-2180-9-259" class="uri">https://doi.org/10.1186/1471-2180-9-259</a>

### 2.3.1 Code

    meta<-fread("../Data/Data/meta.csv",data.table = F)
    OTUstool<-fread("../Data/Data/OTUtabale_regaStool.csv",data.table = F)
    OTUsaliva<-fread("../Data/Data/OTUtabale_regaSaliva.csv",data.table = F)
    TaxonLevels<-colnames(OTUstool)[1:7]
    Stool<-OTUanalysis(OTUtable = OTUstool,
                        TaxonLevels = TaxonLevels,
                        topTaxonomyABvalue = c(0.01,0.01,0.1,0.2,0.2,0.2,0.2))

    ## Warning: `summarise_each_()` was deprecated in dplyr 0.7.0.
    ## Please use `across()` instead.

    Saliva<-OTUanalysis(OTUtable = OTUsaliva,
                      TaxonLevels = TaxonLevels,
                      topTaxonomyABvalue = c(0.01,0.01,0.05,0.05,0.05,0.05,0.05))
                      
    # Transform to compositional abundances
    ABstool<-Stool$TaxonomyAbundance$Genus
    ABsaliva<-Saliva$TaxonomyAbundance$Genus
    # Pick the core (>5% relative abundance in >50% of the samples)
    ABstool.core<-ABstool[apply(ABstool, 1, function(x){
      length(which(x==0))<34&length(which(x>0.05))!=0
    }),]%>%
      .[-grep("norank",rownames(.)),]%>%
      .[-grep("unclassified",rownames(.)),]

    mat1<-ABstool.core
    for (i in 1:nrow(mat1)) {
      mat1[i,][which(mat1[i,]>=0.05)]=2
      mat1[i,][which(mat1[i,]>0&mat1[i,]<0.05)]=1
    }

    pheatmap(mat1,border_color = NA,
                width = unit(cm,"5"),height = unit(cm,"5"),
             fontsize = 6,
             main = "Prevalence of coreMicrobiome (Stool)",
             legend_labels = c("0","<0.05",">0.05"),
             legend_breaks = c(0,1,2),
             show_colnames = F,
             color = scico(100,palette = "bilbao",direction = 1))

![](OTUanalysis_files/figure-markdown_strict/unnamed-chunk-6-1.png)

    ABsaliva.core<-ABsaliva[apply(ABsaliva, 1, function(x){
      length(which(x==0))<21&length(which(x>0.01))!=0
    }),]%>%
      .[-grep("norank",rownames(.)),]%>%
      .[-grep("unclassified",rownames(.)),]

    mat2<-ABsaliva.core
    for (i in 1:nrow(mat2)) {
      mat2[i,][which(mat2[i,]>=0.05)]=2
      mat2[i,][which(mat2[i,]>0&mat2[i,]<0.05)]=1
    }

    pheatmap(mat2,border_color = NA,
             width = unit(cm,"5"),height = unit(cm,"5"),
             fontsize = 6,
             main = "Prevalence of coreMicrobiome (Saliva)",
             legend_labels = c("0","<0.05",">0.05"),
             legend_breaks = c(0,1,2),
             show_colnames = F,
             color = scico(100,palette = "bilbao",direction = 1))

![](OTUanalysis_files/figure-markdown_strict/unnamed-chunk-6-2.png)

    ABstool.core<-data.frame(Samples=colnames(ABstool.core),t(ABstool.core))%>%
      mutate(.,Others=apply(.[,-1],1, function(x) 1-sum(x)))

    ABsaliva.core<-data.frame(Samples=colnames(ABsaliva.core),t(ABsaliva.core))%>%
      mutate(.,Others=apply(.[,-1],1, function(x) 1-sum(x)))

    coreM<-bind_rows(ABstool.core,ABsaliva.core)
    coreM<-merge(select(meta,c(Samples,patientID,Site,Time)),coreM,by="Samples")
    coreMp13<-subset(coreM,patientID=="Patient13")
    coreMp13<-coreMp13[,-c(1:2)]%>%group_by(Site,Time)%>%summarise_each(mean)
    coreMp13<-melt(coreMp13,id.vars = c("Site","Time"),value.name = "Abundance",
                 variable.name = "CoreMicrobiome")

    ## Warning in melt(coreMp13, id.vars = c("Site", "Time"), value.name =
    ## "Abundance", : The melt generic in data.table has been passed a grouped_df
    ## and will attempt to redirect to the relevant reshape2 method; please note that
    ## reshape2 is deprecated, and this redirection is now deprecated as well. To
    ## continue using melt methods from reshape2 while both libraries are attached,
    ## e.g. melt.list, you can prepend the namespace like reshape2::melt(coreMp13). In
    ## the next version, this warning will become an error.

    coreM1<-coreM[,-c(1:2)]%>%group_by(Site,Time)%>%summarise_each(mean)
    coreM1<-melt(coreM1,id.vars = c("Site","Time"),value.name = "Abundance",
                 variable.name = "CoreMicrobiome")

    ## Warning in melt(coreM1, id.vars = c("Site", "Time"), value.name = "Abundance", :
    ## The melt generic in data.table has been passed a grouped_df and will attempt
    ## to redirect to the relevant reshape2 method; please note that reshape2 is
    ## deprecated, and this redirection is now deprecated as well. To continue using
    ## melt methods from reshape2 while both libraries are attached, e.g. melt.list,
    ## you can prepend the namespace like reshape2::melt(coreM1). In the next version,
    ## this warning will become an error.

### 2.3.2 Figures

    p1<-ggplot(coreM1,aes(Time,Abundance,fill=CoreMicrobiome))+
      geom_area()+
      facet_wrap(~Site)+
      theme_minimal()+
      scale_fill_manual(name="coreMicrobiome",values=c(col31,col21,col16,col11))+
      theme(legend.text = element_text(size = 6),legend.box = "horizontal",
            plot.background = element_rect(colour = "black", 
                                           size = 1, linetype = 1,
                                           fill = "gray"),
            legend.margin = margin(0.1,unit="pt"),
            legend.key.size=unit(.1,"inches"),
            legend.text.align=0,
            legend.title=element_text(colour="black",size=8,face = "bold"),
            legend.direction ="vertical",
            legend.box.just="top",
            legend.spacing = unit(0.1,"cm"),
            legend.spacing.y = unit(0.1,"cm"),
            legend.spacing.x =unit(0.1,"cm"),
            legend.box.spacing = unit(0.1,"cm"),
            legend.justification=c(.4,.4),
            legend.position="top",
            plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"))+
      scale_x_continuous(n.breaks = 8,minor_breaks = NULL)+
      xlab("Cycles of treatment")
    p1

<img src="OTUanalysis_files/figure-markdown_strict/unnamed-chunk-7-1.png" width="80%" style="display: block; margin: auto;" />

    p2<-ggplot(coreMp13,aes(Time,Abundance,fill=CoreMicrobiome))+
      geom_area()+
      facet_wrap(~Site,scales = "free")+
      theme_minimal()+
      scale_fill_manual(name="coreMicrobiome (Patient13)",values=c(col31,col21,col16,col11))+
      theme(legend.text = element_text(size = 6),legend.box = "horizontal",
            plot.background = element_rect(colour = "black", 
                                           size = 1, linetype = 1,
                                           fill = "gray"),
            legend.margin = margin(0.1,unit="pt"),
            legend.key.size=unit(.1,"inches"),
            legend.text.align=0,
            legend.title=element_text(colour="black",size=8,face = "bold"),
            legend.direction ="vertical",
            legend.box.just="top",
            legend.spacing = unit(0.1,"cm"),
            legend.spacing.y = unit(0.1,"cm"),
            legend.spacing.x =unit(0.1,"cm"),
            legend.box.spacing = unit(0.1,"cm"),
            legend.justification=c(.4,.4),
            legend.position="top",
            plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"))+
      scale_x_continuous(n.breaks = 8,minor_breaks = NULL)+
      xlab("Cycles of treatment")
    p2

<img src="OTUanalysis_files/figure-markdown_strict/unnamed-chunk-7-2.png" width="80%" style="display: block; margin: auto;" />

2.4 Microbiome compostion
-------------------------

### 2.4.1 Code

    meta<-fread("../Data/Data/meta.csv",data.table = F)
    OTUstool<-fread("../Data/Data/OTUtabale_regaStool.csv",data.table = F)
    OTUsaliva<-fread("../Data/Data/OTUtabale_regaSaliva.csv",data.table = F)
    TaxonLevels<-colnames(OTUstool)[1:7]
    Stool<-OTUanalysis(OTUtable = OTUstool,
                        TaxonLevels = TaxonLevels,
                        topTaxonomyABvalue = c(0.01,0.01,0.1,0.2,0.2,0.2,0.2))

    Saliva<-OTUanalysis(OTUtable = OTUsaliva,
                      TaxonLevels = TaxonLevels,
                      topTaxonomyABvalue = c(0.01,0.01,0.05,0.05,0.05,0.05,0.05))

    stoolComp<-Stool$TaxonomyComposition%>%lapply(.,function(x){
      x<-mutate(x,Microname=rownames(x))
    })
    salivaComp<-Saliva$TaxonomyComposition%>%lapply(.,function(x){
      x<-mutate(x,Microname=rownames(x))
    })

    otu_abundance<-list()
    areaPlot<-list()
    for (i in seq_along(stoolComp)) {
      otu_abundance[[i]]<-merge(stoolComp[[i]],salivaComp[[i]],by="Microname",all = T)
      otu_abundance[[i]]<-reshape2::melt(otu_abundance[[i]],id.vars = "Microname",
                               variable.name = "Samples",value.name = "Abundance")
      otu_abundance[[i]]<-merge(meta,otu_abundance[[i]],by="Samples")
      otu_abundance[[i]]<-select(otu_abundance[[i]],c(Time,Site,Microname,Abundance))%>%
        group_by(Time,Site,Microname)%>%
        summarise_each(mean)
      otu_abundance[[i]]<-arrange(otu_abundance[[i]],desc(Abundance))
     ##plot
      areaPlot[[i]]<- ggplot(otu_abundance[[i]],aes(Time,Abundance,fill=Microname))+
        geom_area()+
        facet_wrap(~Site)+
        theme_minimal()+
        scale_fill_manual(name=names(stoolComp)[i],values=c(col31,col21,col16,col11))+
        theme(legend.text = element_text(size = 6),legend.box = "horizontal",
              plot.background = element_rect(colour = "black", 
                                             size = 1, linetype = 1,
                                             fill = "gray"),
              legend.margin = margin(0.1,unit="pt"),
              legend.key.size=unit(.1,"inches"),
              legend.text.align=0,
              legend.title=element_text(colour="black",size=8,face = "bold"),
              legend.direction ="vertical",
              legend.box.just="top",
              legend.spacing = unit(0.1,"cm"),
              legend.spacing.y = unit(0.1,"cm"),
              legend.spacing.x =unit(0.1,"cm"),
              legend.box.spacing = unit(0.1,"cm"),
              legend.justification=c(.4,.4),
              legend.position="top",
              plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"))+
        scale_x_continuous(n.breaks = 8,minor_breaks = NULL)+
        xlab("Cycles of treatment")
      names(otu_abundance)[i]=names(stoolComp)[i]
      names(areaPlot)[i]=names(stoolComp)[i]
    }

### 2.4.2 Figures

#### 2.4.2.1 Phylum compostion

<details>
<summary>
<font size=4>Figure</font>
</summary>

    print(areaPlot$Phylum)

<img src="OTUanalysis_files/figure-markdown_strict/unnamed-chunk-9-1.png" width="100%" style="display: block; margin: auto;" />
</details>

#### 2.4.2.2 Class compostion

<details>
<summary>
<font size=4>Figure</font>
</summary>

    print(areaPlot$Class)

<img src="OTUanalysis_files/figure-markdown_strict/unnamed-chunk-10-1.png" width="100%" style="display: block; margin: auto;" />
</details>

#### 2.4.2.3 Order compostion

<details>
<summary>
<font size=4>Figure</font>
</summary>

    print(areaPlot$Order)

<img src="OTUanalysis_files/figure-markdown_strict/unnamed-chunk-11-1.png" width="100%" style="display: block; margin: auto;" />
</details>

#### 2.4.2.4 Family compostion

<details>
<summary>
<font size=4>Figure</font>
</summary>

    print(areaPlot$Family)

<img src="OTUanalysis_files/figure-markdown_strict/unnamed-chunk-12-1.png" width="100%" style="display: block; margin: auto;" />
</details>

#### 2.4.2.5 Genus compostion

<details>
<summary>
<font size=4>Figure</font>
</summary>

    print(areaPlot$Genus)

<img src="OTUanalysis_files/figure-markdown_strict/unnamed-chunk-13-1.png" width="100%" style="display: block; margin: auto;" />
</details>

#### 2.4.2.6 Species compostion

<details>
<summary>
<font size=4>Figure</font>
</summary>

    print(areaPlot$Species)

<img src="OTUanalysis_files/figure-markdown_strict/unnamed-chunk-14-1.png" width="100%" style="display: block; margin: auto;" />
</details>

#### 2.4.2.7 OTU compostion

<details>
<summary>
<font size=4>Figure</font>
</summary>

    print(areaPlot$OTU)

<img src="OTUanalysis_files/figure-markdown_strict/unnamed-chunk-15-1.png" width="100%" style="display: block; margin: auto;" />
</details>
