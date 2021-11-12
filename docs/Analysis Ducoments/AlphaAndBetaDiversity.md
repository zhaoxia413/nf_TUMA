---
sort: 3
---

-   [1 Requires](#requires)
-   [2 Alpha and Beta Diversity
    Analysis](#alpha-and-beta-diversity-analysis)
    -   [2.1 Description](#description)
    -   [2.2 Alpha diversity](#alpha-diversity)
        -   [2.2.1 Relationship Alpha diversity and
            Response](#relationship-alpha-diversity-and-response)
        -   [2.2.2 Patients with both saliva and stool samples at
            baseline](#patients-with-both-saliva-and-stool-samples-at-baseline)
        -   [2.2.3 Paired baseline saliva and stool
            samples](#paired-baseline-saliva-and-stool-samples)
    -   [2.3 beta diversity](#beta-diversity)
        -   [2.3.1 meta](#meta)
        -   [2.3.2 anosim: Analysis of
            Similarities](#anosim-analysis-of-similarities)
        -   [2.3.3 adonis: Permutational Multivariate Analysis of
            Variance](#adonis-permutational-multivariate-analysis-of-variance)

[`Return`](./)

1 Requires
==========

<details>
<summary>
<font size=4>Requires</font>
</summary>

    library(data.table)
    library(tidyverse)
    library(ggthemes)
    library(ggsci)
    library(ggpubr)
    library(survminer)
    library(survival)
    library(survivalROC)
    library(reshape2)
    library(ggExtra)
    library(cowplot)
    library(GUniFrac)
    library(ade4)
    library(phangorn)
    library(fpc)
    library(vegan)
    source("../R_function/colors.R")
    "%ni%" <- Negate("%in%")
    theme_set(theme_few())
    options(stringsAsFactors = F)

</details>

2 Alpha and Beta Diversity Analysis
===================================

2.1 Description
---------------

**Ecosystem**

a community plus the physical environment that it occupies at a given
time

**Alpha diversity**

the diversity within a particular area or ecosystem; usually expressed
by the number of species (i.e., species richness) in that ecosystem

**Beta diversity**

a comparison of of diversity between ecosystems, usually measured as the
amount of species change between the ecosystems

**Gamma diversity**

a measure of the overall diversity within a large region.
Geographic-scale species diversity according to Hunter (2002: 448)

ref: Alpha, Beta, and Gamma Diversity. (2021, January 4). Retrieved
April 19, 2021, from
<a href="https://bio.libretexts.org/@go/page/17392" class="uri">https://bio.libretexts.org/@go/page/17392</a>

2.2 Alpha diversity
-------------------

### 2.2.1 Relationship Alpha diversity and Response

    meta<-fread("../Data/Data/meta.csv",data.table = F)
    alpha<-fread("../Data/Data/alpha_diversity.csv",data.table = F)[,c(1:4)]
    data<-merge(select(meta,c(Samples,Response,Cycle,Site)),alpha,by="Samples")
    data<-melt(data,id.vars = c("Samples","Response","Cycle","Site"),variable.name = "Index",value.name = "Value")
    data_BL<-subset(data,Cycle=="BL"&Response!="NE")
    data_treat<-subset(data,Cycle!="BL"&Response!="NE")

    data_BL$Group<-paste0(data_BL$Site,"_",data_BL$Response)
    data_BL$Index<-factor(data_BL$Index,levels = c("shannon","simpson","sobs"))

    data_treat$Group<-paste0(data_treat$Site,"_",data_treat$Response)
    data_treat$Index<-factor(data_treat$Index,levels = c("shannon","simpson","sobs"))

    p1<-ggplot(subset(data_BL,Site=="Saliva"),aes(Response,Value))+
      geom_dotplot(aes(fill=Response,color=Response),binaxis = "y",stackdir = "center",dotsize = .8)+
      stat_compare_means(ref.group = "NR",label = "p.signif")+
      stat_summary(fun.data = mean_sdl,fun.args = list(mult=1),geom = "errorbar",width=.1,size=0.5)+
      stat_summary(fun = median,geom = "crossbar",color="black",size=.4,width=.4)+
      facet_wrap(~Index,scales = "free")+ theme_few()+
      scale_color_manual(values = c("#FF0000FF","#00FFFFFF"))+
      scale_fill_manual(values = c("#FF0000FF","#00FFFFFF"))+
      ggtitle("Saliva (Baseline)")+
      theme_few(base_size = 8)+
      theme(axis.title.x = element_blank())
    p2<-ggplot(subset(data_BL,Site=="Stool"),aes(Response,Value))+
      geom_dotplot(aes(fill=Response,color=Response),binaxis = "y",stackdir = "center",dotsize = .8)+
      stat_compare_means(ref.group = "NR",label = "p.signif")+
      stat_summary(fun.data = mean_sdl,fun.args = list(mult=1),geom = "errorbar",width=.1,size=0.5)+
      stat_summary(fun = median,geom = "crossbar",color="black",size=.4,width=.4)+
      facet_wrap(~Index,scales = "free")+ theme_few()+
      scale_color_manual(values = c("#80FF00FF","#8000FFFF"))+
      scale_fill_manual(values = c("#80FF00FF","#8000FFFF"))+
      ggtitle("Stool (Baseline)")+
      theme_few(base_size = 8)+
      theme(axis.title.x = element_blank())
    p3<-ggplot(subset(data_treat,Site=="Saliva"),aes(Response,Value))+
      geom_dotplot(aes(fill=Response,color=Response),binaxis = "y",stackdir = "center",dotsize = .8)+
      stat_compare_means(ref.group = "NR",label = "p.signif")+
      stat_summary(fun.data = mean_sdl,fun.args = list(mult=1),geom = "errorbar",width=.1,size=0.5)+
      stat_summary(fun = median,geom = "crossbar",color="black",size=.4,width=.4)+
      facet_wrap(~Index,scales = "free")+ theme_few()+
      scale_color_jama()+
      scale_fill_jama()+
      ggtitle("Saliva (Treat)")+
      theme_few(base_size = 8)+
      theme(axis.title.x = element_blank())

    ##use the treated stool samples without duplication
    BL_treat<-fread("../Data/Data/16pt_BL_treat_pairs_stool.csv",data.table = F)
    BL_treat_alpha<-merge(select(BL_treat,c(Samples,Group,Response)),alpha,by="Samples")
    stool_treat<-subset(BL_treat_alpha,Group=="Treat")
    stool_treat1<-melt(stool_treat[,-2],id.vars = c("Samples","Response"),value.name = "Value",variable.name = "Index")
    stool_treat1$Index<-factor(stool_treat1$Index,levels = c("shannon","simpson","sobs"))
    p4<-ggplot(stool_treat1,aes(Response,Value))+
      geom_dotplot(aes(fill=Response,color=Response),binaxis = "y",stackdir = "center",dotsize = .8)+
      stat_summary(fun.data = mean_sdl,fun.args = list(mult=1),geom = "errorbar",width=.1,size=0.5)+
      stat_summary(fun = median,geom = "crossbar",color="black",size=.4,width=.4)+
      stat_compare_means(ref.group = "NR",label = "p.signif",method = "wilcox.test")+
      facet_wrap(~Index,scales = "free")+
      scale_color_manual(values = c("#FF0099FF","#FF9900FF"))+
      scale_fill_manual(values = c("#FF0099FF","#FF9900FF"))+
      ggtitle("Stool (Treat)")+
      theme_few(base_size = 8)+ 
      theme(axis.title.x = element_blank())

    plot_grid(p1,p2,p3,p4, labels = c("A","B","C","D"), ncol =2, nrow = 2)

    ## `stat_bindot()` using `bins = 30`. Pick better value with `binwidth`.
    ## `stat_bindot()` using `bins = 30`. Pick better value with `binwidth`.
    ## `stat_bindot()` using `bins = 30`. Pick better value with `binwidth`.
    ## `stat_bindot()` using `bins = 30`. Pick better value with `binwidth`.

<img src="AlphaAndBetaDiversity_files/figure-markdown_strict/unnamed-chunk-3-1.png" width="80%" style="display: block; margin: auto;" />

### 2.2.2 Patients with both saliva and stool samples at baseline

    paired_meta<-fread("../Data/Data/stool_saliva_BL_paried.csv",data.table = F)
    alpha_paired<-subset(alpha,Samples%in%paired_meta$Samples)
    alpha_paired<-merge(paired_meta,alpha_paired,by="Samples")%>%filter(.,Response!="NE")
    ggplot(alpha_paired,aes(Site,shannon,fill=Site))+
      geom_boxplot()+
      geom_line(aes(group=patientID))+
      geom_point(size=1)+
      theme_few(base_size = 8)+
      stat_compare_means(label = "p.signif")+
      facet_wrap(~Effect,scales = "free_x")+
      scale_fill_d3()

<img src="AlphaAndBetaDiversity_files/figure-markdown_strict/unnamed-chunk-4-1.png" width="30%" style="display: block; margin: auto;" />

### 2.2.3 Paired baseline saliva and stool samples

    load("../Data/Data/regaMicrcobiome.RData")
    phylum_saliva<-regaMicrobiome$SalivaMicrobiome$TaxonomyComposition$Phylum
    phylum_stool<-regaMicrobiome$StoolMicrobiome$TaxonomyComposition$Phylum
    saliva_meta<-subset(paired_meta,Site=="Saliva"&Response!="NE")
    stool_meta<-subset(paired_meta,Site=="Stool"&Response!="NE")
    phylum_saliva<-phylum_saliva[,which(colnames(phylum_saliva)%in%saliva_meta$Samples)]
    phylum_stool<-phylum_stool[,which(colnames(phylum_stool)%in%stool_meta$Samples)]
    phylum_saliva<-data.frame(Samples=colnames(phylum_saliva),t(phylum_saliva))
    phylum_stool<-data.frame(Samples=colnames(phylum_stool),t(phylum_stool))
    phylum_saliva<-merge(select(paired_meta,c(patientID,Response,Samples)),phylum_saliva,by="Samples")
    phylum_stool<-merge(select(paired_meta,c(patientID,Response,Samples)),phylum_stool,by="Samples")
    phylum_saliva<-melt(phylum_saliva[,-1],id.vars = c("patientID","Response"),variable.name = "Phylum",value.name = "Abundance")
    phylum_stool<-melt(phylum_stool[,-1],id.vars = c("patientID","Response"),variable.name = "Phylum",value.name = "Abundance")
    phylum_saliva$Site=rep("Saliva",171)
    phylum_stool$Site=rep("Stool",171)
    df_p<-bind_rows(phylum_saliva,phylum_stool)
    test<-subset(df_p,Phylum=="Firmicutes"&Site=="Stool")
    g<-test[order(test$Abundance,decreasing = TRUE),]
    df_p$patientID<-factor(df_p$patientID,levels=g$patientID)
    levels(factor(df_p$Phylum))
    df_p$Phylum<-factor(df_p$Phylum,levels = c("Firmicutes","Bacteroidetes","Proteobacteria",
                                               "Fusobacteriota","Actinobacteriota","Campilobacterota",
                                               "Patescibacteria","Spirochaetota","Desulfobacterota",
                                               "Synergistota","Verrucomicrobiota","Others"))

    p1<-ggplot(subset(df_p,Site=="Saliva"),aes(patientID,Abundance,fill=Phylum))+
      geom_bar(stat = "identity", width=1)+
      facet_grid(~Response,scales = "free",space = "free")+
      scale_fill_manual(values = col11)+
      theme_few(base_size = 8)+
      theme(axis.text.x = element_blank(),
            strip.text = element_text(face = "bold"),
            axis.title.x = element_blank(),
            legend.margin = margin(0.1,unit="pt"),
            legend.key.size=unit(.1,"inches"),
            legend.text.align=0,
            legend.title=element_text(colour="black",size=8,face = "bold"),
            legend.direction = "horizontal",
            legend.box.just="top",
            legend.spacing = unit(0.1,"cm"),
            legend.spacing.y = unit(0.1,"cm"),
            legend.spacing.x =unit(0.1,"cm"),
            legend.box.spacing = unit(0.1,"cm"),
            legend.justification=c(.4,.4),
            legend.position="top")+
      ggtitle("Patients with paired baseline saliva and stool samples")

    p2<-ggplot(subset(df_p,Site=="Stool"),aes(patientID,Abundance,fill=Phylum))+
      geom_bar(stat = "identity", width=1)+
      facet_grid(~Response,scales = "free",space = "free")+
      scale_fill_manual(values = col11,guide=FALSE)+
      theme_few(base_size = 8)+
      theme(axis.text.x = element_text(angle = 90,vjust = 1,hjust = 1,size=6),
            strip.text = element_blank(),
            axis.title.x = element_blank())

    ##  [1] "Actinobacteriota"  "Bacteroidetes"     "Campilobacterota" 
    ##  [4] "Firmicutes"        "Fusobacteriota"    "Patescibacteria"  
    ##  [7] "Proteobacteria"    "Spirochaetota"     "Others"           
    ## [10] "Desulfobacterota"  "Synergistota"      "Verrucomicrobiota"

    plot_grid(p1,p2, labels = c("A","B"), ncol =1, nrow = 2)

<img src="AlphaAndBetaDiversity_files/figure-markdown_strict/unnamed-chunk-6-1.png" width="60%" style="display: block; margin: auto;" />

2.3 beta diversity
------------------

### 2.3.1 meta

    meta_file<-subset(meta,Response!="NE"&Cycle=="BL")
    meta_file$Group<-paste0(meta_file$Site,"_",meta_file$Response)
    meta_file$newName<-paste0(meta_file$patientID,"_",meta_file$Cycle)

    meta_file_stool<-subset(meta_file,Site=="Stool")
    meta_file_saliva<-subset(meta_file,Site=="Saliva")
    meta_file_BL<-meta_file
    tree_file <- read.tree("../Data/Data/top500_OTU_phylo_tree.tre")
    rooted_tree <- midpoint(tree_file)
    otu<-fread("../Data/Data/top500_OTU_table_phylo_tree.csv",data.table = F)
    otu<-data.frame(row.names = otu$ID,otu[,-1])

    otu_stool<-otu[,which(colnames(otu)%in%meta_file_stool$Samples)]
    otu_saliva<-otu[,which(colnames(otu)%in%meta_file_saliva$Samples)]
    otu_BL<-otu[,which(colnames(otu)%in%meta_file_BL$Samples)]

    otu_stool<-data.frame(t(otu_stool[,order(names(otu_stool))]))
    otu_saliva<-data.frame(t(otu_saliva[,order(names(otu_saliva))]))
    otu_BL<-data.frame(t(otu_BL[,order(names(otu_BL))]))

    unifracs_stool<-GUniFrac(otu_stool, rooted_tree, alpha = c(0.0,0.5,1.0))$unifracs
    unifracs_saliva<-GUniFrac(otu_saliva, rooted_tree, alpha = c(0.0,0.5,1.0))$unifracs
    unifracs_BL<-GUniFrac(otu_BL, rooted_tree, alpha = c(0.0,0.5,1.0))$unifracs

    unifract_dist_stool <- unifracs_stool[, , "d_0.5"]# GUniFrac with alpha 0.5
    unifract_dist_saliva <- unifracs_saliva[, , "d_0.5"]
    unifract_dist_BL <- unifracs_BL[, , "d_0.5"]

    stool_fit<-hclust(as.dist(unifract_dist_stool),method = "ward.D2")
    saliva_fit<-hclust(as.dist(unifract_dist_saliva),method = "ward.D2")
    BL_fit<-hclust(as.dist(unifract_dist_BL),method = "ward.D2")

    stool_tree <- as.phylo(stool_fit)
    saliva_tree <- as.phylo(saliva_fit)
    BL_tree <- as.phylo(BL_fit)

    stool_groups<-as.factor(c(meta_file_stool$Group))
    saliva_groups<-as.factor(c(meta_file_saliva$Group))
    BL_groups<-as.factor(c(meta_file_BL$Group))

    plot_color_BL<-rainbow(length(levels(BL_groups)))[BL_groups]
    plot_color_stool<-rainbow(length(levels(stool_groups)))[stool_groups]
    plot_color_saliva<-rainbow(length(levels(saliva_groups)))[saliva_groups]

### 2.3.2 anosim: Analysis of Similarities

**Description** Analysis of similarities (ANOSIM) provides a way to test
statistically whether there is a significant difference between two or
more groups of sampling units.

    par(mfrow=c(3,1))
    plot(anosim(vegdist(otu_stool,method = "bray"), meta_file_stool$Group),xlab="",ylab="",col=c("gray","#8000FFFF","#80FF00FF"))
    plot(anosim(vegdist(otu_saliva,method = "bray"), meta_file_saliva$Group),xlab="",ylab="",col=c("gray","#00FFFFFF","#FF0000FF"))
    plot(anosim(vegdist(otu_BL,method = "bray"), meta_file_BL$Group),xlab="",ylab="",col=c("gray","#00FFFFFF","#FF0000FF","#8000FFFF","#80FF00FF"))

<img src="AlphaAndBetaDiversity_files/figure-markdown_strict/unnamed-chunk-8-1.png" width="60%" style="display: block; margin: auto;" />

### 2.3.3 adonis: Permutational Multivariate Analysis of Variance

\*\* Description\*\* Analysis of variance using distance matrices — for
partitioning distance matrices among sources of variation and fitting
linear models (e.g., factors, polynomial regression) to distance
matrices; uses a permutation test with pseudo-F ratios.

    adonis_stool<-adonis2(as.dist(unifract_dist_stool) ~ stool_groups)
    adonis_saliva<-adonis2(as.dist(unifract_dist_saliva) ~ saliva_groups)
    adonis_BL<-adonis2(as.dist(unifract_dist_BL) ~ BL_groups)

    # Calculate and display the NMDS plot (Non-metric Multidimensional Scaling plot)
    res_stool <- metaMDS(unifract_dist_stool,k = 2,distance = "bray")
    s.class(
      res_stool$points, col = unique(plot_color_stool), cpoint = 2, 
      fac = stool_groups,
      sub = paste("NMDS plot of OTU \n (p ",adonis_stool$`Pr(>F)`[1],", stress ",
                  round(res_stool$stress,3),")",sep="")
    )

<img src="AlphaAndBetaDiversity_files/figure-markdown_strict/unnamed-chunk-9-1.png" width="30%" style="display: block; margin: auto;" />

    res_saliva <- metaMDS(unifract_dist_saliva,k = 2,distance = "bray")
    s.class(
      res_saliva$points, col = unique(plot_color_saliva), cpoint = 2, 
      fac = saliva_groups,
      sub = paste("NMDS plot of OTU \n (p ",adonis_saliva$`Pr(>F)`[1],", stress ",round(res_saliva$stress,3),")",sep="")
    )

<img src="AlphaAndBetaDiversity_files/figure-markdown_strict/unnamed-chunk-9-2.png" width="30%" style="display: block; margin: auto;" />

    res_BL <- metaMDS(unifract_dist_BL,k = 2,distance = "bray")
    s.class(
      res_BL$points, col = unique(plot_color_BL), cpoint = 2, 
      fac = BL_groups,
      sub = paste("NMDS plot of OTU \n (p ",adonis_BL$`Pr(>F)`[1],", stress ",
                  round(res_BL$stress,3),")",sep="")
    )

<img src="AlphaAndBetaDiversity_files/figure-markdown_strict/unnamed-chunk-9-3.png" width="30%" style="display: block; margin: auto;" />

    ## Run 0 stress 0.1685622 
    ## Run 1 stress 0.1685622 
    ## ... Procrustes: rmse 1.470863e-05  max resid 6.422478e-05 
    ## ... Similar to previous best
    ## Run 2 stress 0.1740149 
    ## Run 3 stress 0.168407 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0129008  max resid 0.05542212 
    ## Run 4 stress 0.1685622 
    ## ... Procrustes: rmse 0.01289925  max resid 0.05536435 
    ## Run 5 stress 0.1987313 
    ## Run 6 stress 0.1685622 
    ## ... Procrustes: rmse 0.01291624  max resid 0.0554523 
    ## Run 7 stress 0.1668062 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0641624  max resid 0.3365322 
    ## Run 8 stress 0.1754553 
    ## Run 9 stress 0.1689635 
    ## Run 10 stress 0.2173753 
    ## Run 11 stress 0.1668062 
    ## ... Procrustes: rmse 9.950925e-06  max resid 4.231405e-05 
    ## ... Similar to previous best
    ## Run 12 stress 0.1689634 
    ## Run 13 stress 0.1913187 
    ## Run 14 stress 0.1668062 
    ## ... Procrustes: rmse 9.200012e-06  max resid 3.840477e-05 
    ## ... Similar to previous best
    ## Run 15 stress 0.1685622 
    ## Run 16 stress 0.1685622 
    ## Run 17 stress 0.168407 
    ## Run 18 stress 0.1668062 
    ## ... New best solution
    ## ... Procrustes: rmse 5.129414e-06  max resid 1.751809e-05 
    ## ... Similar to previous best
    ## Run 19 stress 0.1685622 
    ## Run 20 stress 0.1764919 
    ## *** Solution reached
    ## Run 0 stress 0.07988775 
    ## Run 1 stress 0.07988776 
    ## ... Procrustes: rmse 1.56881e-05  max resid 5.678516e-05 
    ## ... Similar to previous best
    ## Run 2 stress 0.08928076 
    ## Run 3 stress 0.1144174 
    ## Run 4 stress 0.09548344 
    ## Run 5 stress 0.09972515 
    ## Run 6 stress 0.08928077 
    ## Run 7 stress 0.09897446 
    ## Run 8 stress 0.1176597 
    ## Run 9 stress 0.08928076 
    ## Run 10 stress 0.1125236 
    ## Run 11 stress 0.08928077 
    ## Run 12 stress 0.09897451 
    ## Run 13 stress 0.09443043 
    ## Run 14 stress 0.07988776 
    ## ... Procrustes: rmse 1.826532e-05  max resid 6.424777e-05 
    ## ... Similar to previous best
    ## Run 15 stress 0.112533 
    ## Run 16 stress 0.08928076 
    ## Run 17 stress 0.09972511 
    ## Run 18 stress 0.0989745 
    ## Run 19 stress 0.09333011 
    ## Run 20 stress 0.09972529 
    ## *** Solution reached
    ## Run 0 stress 0.1175033 
    ## Run 1 stress 0.1536876 
    ## Run 2 stress 0.1175033 
    ## ... Procrustes: rmse 8.312975e-06  max resid 5.142803e-05 
    ## ... Similar to previous best
    ## Run 3 stress 0.1395504 
    ## Run 4 stress 0.1565893 
    ## Run 5 stress 0.1175033 
    ## ... Procrustes: rmse 7.229981e-06  max resid 3.982653e-05 
    ## ... Similar to previous best
    ## Run 6 stress 0.1175033 
    ## ... Procrustes: rmse 4.53565e-06  max resid 2.569473e-05 
    ## ... Similar to previous best
    ## Run 7 stress 0.1457774 
    ## Run 8 stress 0.1175033 
    ## ... Procrustes: rmse 6.942335e-06  max resid 4.602578e-05 
    ## ... Similar to previous best
    ## Run 9 stress 0.1426608 
    ## Run 10 stress 0.1616086 
    ## Run 11 stress 0.1542791 
    ## Run 12 stress 0.1516694 
    ## Run 13 stress 0.1175033 
    ## ... Procrustes: rmse 6.837063e-06  max resid 4.304509e-05 
    ## ... Similar to previous best
    ## Run 14 stress 0.1330585 
    ## Run 15 stress 0.1516694 
    ## Run 16 stress 0.1419744 
    ## Run 17 stress 0.1514029 
    ## Run 18 stress 0.1175033 
    ## ... Procrustes: rmse 1.856959e-05  max resid 8.940436e-05 
    ## ... Similar to previous best
    ## Run 19 stress 0.1233161 
    ## Run 20 stress 0.1233161 
    ## *** Solution reached
