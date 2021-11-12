---
sort: 6
---
-   [1 Requires](#requires)
-   [2 Patients with BL and treat stool and saliva
    sampls](#patients-with-bl-and-treat-stool-and-saliva-sampls)
    -   [2.1 Cluser Dendrograms](#cluser-dendrograms)
    -   [2.2 Phylum and Family
        composition](#phylum-and-family-composition)
-   [3 Autoencoder for the patients with completion
    data](#autoencoder-for-the-patients-with-completion-data)
    -   [3.1 K-means clustering](#k-means-clustering)
        -   [3.1.1 Choose k value](#choose-k-value)
        -   [3.1.2 Plot K-means clustering](#plot-k-means-clustering)
        -   [3.1.3 PFS survival based on Autoencoder
            data](#pfs-survival-based-on-autoencoder-data)
    -   [3.2 Representative features of the k-means
        clusters](#representative-features-of-the-k-means-clusters)
    -   [3.3 Table of representative OTUs for the k-means
        clusters](#table-of-representative-otus-for-the-k-means-clusters)
    -   [3.4 PCA plot based on the representative
        OTUs](#pca-plot-based-on-the-representative-otus)
-   [4 Reclassify the patients with incompletion data by the
    representative
    OTUs](#reclassify-the-patients-with-incompletion-data-by-the-representative-otus)
    -   [4.1 All the patients with baseline
        samples](#all-the-patients-with-baseline-samples)
        -   [4.1.1 Choose k value](#choose-k-value-1)
        -   [4.1.2 Plot K-means clustering](#plot-k-means-clustering-1)
    -   [4.2 All the patients with treated
        samples](#all-the-patients-with-treated-samples)
        -   [4.2.1 Choose k value](#choose-k-value-2)
        -   [4.2.2 Plot K-means clustering](#plot-k-means-clustering-2)
    -   [4.3 PFS survival of the Patients with incompletion
        data](#pfs-survival-of-the-patients-with-incompletion-data)

[`Return`](./)

1 Requires
==========

<details>
<summary>
<font size=4>Requires</font>
</summary>

    library(data.table)
    library(reshape2)
    library(ggthemes)
    library(ggsci)
    library(tidyverse)
    library(FactoMineR)
    library(corrplot)
    library(colortools)
    library(visibly)
    library(plotly)
    library(scico)
    library(factoextra)
    library(randomForest)
    library(ANN2)
    library(NeuralNetTools)
    library(ConsensusClusterPlus)
    library(survminer)
    library(survival)
    library(ggExtra)
    library(cowplot)
    library(corrplot)
    library(limma)
    source("../R_function/colors.R")
    theme_set(theme_cowplot())
    "%ni%" <- Negate("%in%")
    options(stringsAsFactors = F)

</details>

2 Patients with BL and treat stool and saliva sampls
====================================================

    paired_stool_saliva_TB<-fread("../Data/Data/paired_stool_saliva_TB_8patients.csv",data.table = F)
    clinical<-fread("../Data/Data/clinical_adver_41p.csv",data.table = F)[,-10]
    otu<-fread("../Data/Data/OTUtable_ori.csv",data.table = F)
    colnames(otu)[1]="OTUid"
    group<-paired_stool_saliva_TB[,c(1,2)]
    group<-group[-which(duplicated(group$patientID)),]
    group<-data.frame(row.names = group$patientID,Response=group$Response)
    paired_stool_saliva_TB<-paired_stool_saliva_TB[,-2]
    paired_stool_saliva_TB$Group<-paste0(paired_stool_saliva_TB$Group,"_",paired_stool_saliva_TB$Site)
    paired_stool_saliva_TB<-paired_stool_saliva_TB[,-4]
    OTU_sliva_stool_TB<-otu[,which(colnames(otu)%in%c("OTUid",paired_stool_saliva_TB$Samples))]
    OTU_sliva_stool_TB<-data.frame(row.names = OTU_sliva_stool_TB$OTUid,OTU_sliva_stool_TB[,-1])
    OTU_sliva_stool_TB<-data.frame(Samples=colnames(OTU_sliva_stool_TB),t(OTU_sliva_stool_TB))
    OTU_sliva_stool_TB<-merge(paired_stool_saliva_TB,OTU_sliva_stool_TB,by="Samples")
    levels(factor(OTU_sliva_stool_TB$Group))
    OTU_BL_Saliva<-subset(OTU_sliva_stool_TB,Group=="BL_Saliva")[,-c(1,3)]
    OTU_BL_Stool<-subset(OTU_sliva_stool_TB,Group=="BL_Stool")[,-c(1,3)]
    OTU_Treat_Saliva<-subset(OTU_sliva_stool_TB,Group=="Treat_Saliva")[,-c(1,3)]
    OTU_Treat_Stool<-subset(OTU_sliva_stool_TB,Group=="Treat_Stool")[,-c(1,3)]

    OTU_BL_Saliva_stat<-data.frame(OTUid=colnames(OTU_BL_Saliva)[-1],AV=apply(OTU_BL_Saliva[,-1], 2, mean))%>%
      filter(.,AV>5)
    OTU_BL_Stool_stat<-data.frame(OTUid=colnames(OTU_BL_Stool)[-1],AV=apply(OTU_BL_Stool[,-1], 2, mean))%>%
      filter(.,AV>5)
    OTU_Treat_Saliva_stat<-data.frame(OTUid=colnames(OTU_Treat_Saliva)[-1],AV=apply(OTU_Treat_Saliva[,-1], 2, mean))%>%
      filter(.,AV>5)
    OTU_Treat_Stool_stat<-data.frame(OTUid=colnames(OTU_Treat_Stool)[-1],AV=apply(OTU_Treat_Stool[,-1], 2, mean))%>%
      filter(.,AV>5)
    OTU_BL_Saliva<-OTU_BL_Saliva[,which(colnames(OTU_BL_Saliva)%in%c("patientID",levels(factor(OTU_BL_Saliva_stat$OTUid))))]
    OTU_BL_Stool<-OTU_BL_Stool[,which(colnames(OTU_BL_Stool)%in%c("patientID",levels(factor(OTU_BL_Stool_stat$OTUid))))]
    OTU_Treat_Saliva<-OTU_Treat_Saliva[,which(colnames(OTU_Treat_Saliva)%in%c("patientID",levels(factor(OTU_Treat_Saliva_stat$OTUid))))]
    OTU_Treat_Stool<-OTU_Treat_Stool[,which(colnames(OTU_Treat_Stool)%in%c("patientID",levels(factor(OTU_Treat_Stool_stat$OTUid))))]
    colnames(OTU_BL_Saliva)[-1]<-paste0(colnames(OTU_BL_Saliva)[-1],"_BL_Saliva")
    colnames(OTU_BL_Stool)[-1]<-paste0(colnames(OTU_BL_Stool)[-1],"_BL_Stool")
    colnames(OTU_Treat_Saliva)[-1]<-paste0(colnames(OTU_Treat_Saliva)[-1],"_Treat_Saliva")
    colnames(OTU_Treat_Stool)[-1]<-paste0(colnames(OTU_Treat_Stool)[-1],"_Treat_Stool")
    OTU_BL_Saliva_mat<-data.frame(row.names =OTU_BL_Saliva$patientID,OTU_BL_Saliva[,-1])%>%as.matrix()
    OTU_BL_Stool_mat<-data.frame(row.names =OTU_BL_Stool$patientID,OTU_BL_Stool[,-1])%>%as.matrix()
    OTU_Treat_Saliva_mat<-data.frame(row.names =OTU_Treat_Saliva$patientID,OTU_Treat_Saliva[,-1])%>%as.matrix()
    OTU_Treat_Stool_mat<-data.frame(row.names =OTU_Treat_Stool$patientID,OTU_Treat_Stool[,-1])%>%as.matrix()
    clinical_mat<-subset(clinical,patientID%in%paired_stool_saliva_TB$patientID)
    clinical_mat<-data.frame(row.names =clinical_mat$patientID,clinical_mat[,-1])%>%as.matrix()
    res.hc1 <- eclust(OTU_BL_Saliva_mat, "hclust", k = 2,
                      method = "ward.D2", graph =T) 
    res.hc2 <- eclust(OTU_BL_Stool_mat, "hclust", k = 2,
                      method = "ward.D2", graph =T) 
    res.hc3 <- eclust(OTU_Treat_Saliva_mat, "hclust", k = 2,
                      method = "ward.D2", graph =T) 
    res.hc4 <- eclust(OTU_Treat_Stool_mat, "hclust", k = 2,
                      method = "ward.D2", graph =T) 
    res.hc5 <- eclust(clinical_mat, "hclust", k = 2,
                      method = "ward.D2", graph =T) 

    p1<-fviz_dend(res.hc1, k = 2,main = "BL_Saliva",ggtheme = theme_few(base_size = 6),palette = "lancet",cex = 0.5)

    p2<-fviz_dend(res.hc2, k = 2, main = "BL_Stool",ggtheme = theme_few(base_size = 6),palette = "lancet", cex = 0.5)

    p3<-fviz_dend(res.hc3, k = 2,main = "Treat_Saliva",ggtheme = theme_few(base_size = 6),palette = "lancet", cex = 0.5)

    p4<-fviz_dend(res.hc4, k = 2, main = "BL_Treat",ggtheme = theme_few(base_size = 6),palette ="lancet", cex = 0.5)

    p5<-fviz_dend(res.hc5, k = 2,main = "Clinical", ggtheme = theme_few(base_size = 6),palette = "lancet", cex = 0.5)

    ## [1] "BL_Saliva"    "BL_Stool"     "Treat_Saliva" "Treat_Stool"

2.1 Cluser Dendrograms
----------------------

    plot_grid(p1,p2,p3,p4,p5, labels = c("A","B","C","D","E"), ncol =5 ,nrow = 1)

<img src="AutoDecoder_files/figure-markdown_strict/unnamed-chunk-3-1.png" width="80%" style="display: block; margin: auto;" />

2.2 Phylum and Family composition
---------------------------------

    load("../Data/Data/regaMicrcobiome.RData")
    patients8<-fread("../Data/Data/paired_stool_saliva_TB_8patients.csv",data.table = F)
    phylum_stool<-regaMicrobiome$StoolMicrobiome$TaxonomyComposition$Phylum
    phylum_saliva<-regaMicrobiome$SalivaMicrobiome$TaxonomyComposition$Phylum
    family_stool<-regaMicrobiome$StoolMicrobiome$TaxonomyComposition$Family
    family_saliva<-regaMicrobiome$SalivaMicrobiome$TaxonomyComposition$Family

    data_list<-list(phylum_stool=phylum_stool,phylum_saliva=phylum_saliva,
                    family_stool=family_stool,family_saliva=family_saliva)
    data_list<-lapply(data_list, function(x){
      x<-x[,which(colnames(x)%in%patients8$Samples)]
      x<-data.frame(Samples=colnames(x),t(x))
      x<-merge(patients8,x,by="Samples")[,-1]
      x<-melt(x,id.vars = c("patientID" ,"Response","Group","Site" ),
              variable.name = "MicroName",
              value.name = "Abundance")
    })

    data_list[[1]]$patientID<-factor(data_list[[1]]$patientID,
                                     levels = c("Patient38","Patient41","Patient10",
                                                "Patient13","Patient16","Patient26","Patient32",
                                                "Patient35"))
    data_list[[2]]$patientID<-factor(data_list[[2]]$patientID,
                                     levels = c("Patient38","Patient41","Patient10",
                                                "Patient13","Patient16","Patient26","Patient32",
                                                "Patient35"))
    data_list[[3]]$patientID<-factor(data_list[[3]]$patientID,
                                     levels = c("Patient38","Patient41","Patient10",
                                                "Patient13","Patient16","Patient26","Patient32",
                                                "Patient35"))
    data_list[[4]]$patientID<-factor(data_list[[4]]$patientID,
                                     levels = c("Patient38","Patient41","Patient10",
                                                "Patient13","Patient16","Patient26","Patient32",
                                                "Patient35"))
    p1<-ggplot(data_list$phylum_saliva,aes(Group,Abundance,fill=MicroName))+
      geom_bar(stat = "identity", width=1)+
      facet_grid(Site~patientID,space="free",scales = "free")+
      theme_few(base_size = 5)+
      scale_fill_manual(name="Phylum",values = col16)+
      theme(legend.box.just="top",
            legend.spacing = unit(0.1,"cm"),
            legend.spacing.y = unit(0.1,"cm"),
            legend.spacing.x =unit(0.1,"cm"),
            legend.box.spacing = unit(0.1,"cm"),
            legend.justification=c(.4,.4),
            legend.position="top",legend.key.size=unit(.1,"inches"),axis.text.x = element_text(size=5,angle = 90,vjust = 1,hjust = 1),
            axis.title.x = element_blank())
    p2<-ggplot(data_list$phylum_stool,aes(Group,Abundance,fill=MicroName))+
      geom_bar(stat = "identity", width=1)+
      facet_grid(Site~patientID,space="free",scales = "free")+
      theme_few(base_size = 5)+
      scale_fill_manual(name="Phylum",values = col16)+
      theme(legend.box.just="top",
            legend.spacing = unit(0.1,"cm"),
            legend.spacing.y = unit(0.1,"cm"),
            legend.spacing.x =unit(0.1,"cm"),
            legend.box.spacing = unit(0.1,"cm"),
            legend.justification=c(.4,.4),
            legend.position="top",legend.key.size=unit(.1,"inches"),axis.text.x = element_text(size=5,angle = 90,vjust = 1,hjust = 1),
            axis.title.x = element_blank())

    p3<-ggplot(data_list$family_saliva,aes(Group,Abundance,fill=MicroName))+
      geom_bar(stat = "identity", width=1)+
      facet_grid(Site~patientID,space="free",scales = "free")+
      theme_few(base_size = 5)+
      scale_fill_manual(name="Family",values = col31[c(1:21,24)])+
      theme(legend.box.just="top",
            legend.spacing = unit(0.1,"cm"),
            legend.spacing.y = unit(0.1,"cm"),
            legend.spacing.x =unit(0.1,"cm"),
            legend.box.spacing = unit(0.1,"cm"),
            legend.justification=c(.4,.4),
            legend.position="top",legend.key.size=unit(.1,"inches"), axis.text.x = element_text(size=5,angle = 90,vjust = 1,hjust = 1),
            axis.title.x = element_blank())
    p4<-ggplot(data_list$family_stool,aes(Group,Abundance,fill=MicroName))+
      geom_bar(stat = "identity", width=1)+
      facet_grid(Site~patientID,space="free",scales = "free")+
      theme_few(base_size = 5)+
      scale_fill_manual(name="Family",values =col31[c(1:21,24)])+
      theme(legend.box.just="top",
            legend.spacing = unit(0.1,"cm"),
            legend.spacing.y = unit(0.1,"cm"),
            legend.spacing.x =unit(0.1,"cm"),
            legend.box.spacing = unit(0.1,"cm"),
            legend.justification=c(.4,.4),
            legend.position="top",legend.key.size=unit(.1,"inches"),axis.text.x = element_text(size=5,angle = 90,vjust = 1,hjust = 1),
            axis.title.x = element_blank())

    plot_grid(p1,p3,p2,p4, labels = c("A","B","C","D"), ncol =2 ,nrow = 2)

<img src="AutoDecoder_files/figure-markdown_strict/unnamed-chunk-5-1.png" width="100%" style="display: block; margin: auto;" />

3 Autoencoder for the patients with completion data
===================================================

    clinical<-fread("../Data/Data/clinical_adver_41p.csv",data.table = F)
    cd<-fread("../Data/Data/cd3cd8.csv",data.table = F)
    clinical_cd<-merge(clinical,cd,by="patientID")
    treat_BL<-fread("../Data/Data/16pt_BL_treat_pairs_stool.csv",data.table = F)
    treat_BL<-subset(treat_BL,patientID%in%clinical_cd$patientID)
    treat_BL_BL<-subset(treat_BL,Group=="BL")
    treat_BL_Treat<-subset(treat_BL,Group=="Treat")
    clinical_cd<-subset(clinical_cd,patientID%in%treat_BL$patientID)
    clinical_cd_mat<-data.frame(row.names = clinical_cd$patientID,clinical_cd[,-1])
    otu<-regaMicrobiome$StoolMicrobiome$TaxonomyReads$OTU
    colnames(otu)[1]="OTUid"
    otu_BL<-otu[,which(colnames(otu)%in%c("OTUid",treat_BL_BL$Samples))]
    otu_Treat<-otu[,which(colnames(otu)%in%c("OTUid",treat_BL_Treat$Samples))]
    otu_BL_stat<-data.frame(OTUid=otu_BL$OTUid,num=apply(otu_BL[,-1],1,sum))%>%
      filter(.,num>10)
    otu_Treat_stat<-data.frame(OTUid=otu_Treat$OTUid,num=apply(otu_Treat[,-1],1,sum))%>%
      filter(.,num>10)
    otu_BL<-subset(otu_BL,OTUid%in%otu_BL_stat$OTUid)
    otu_Treat<-subset(otu_Treat,OTUid%in%otu_Treat_stat$OTUid)

    otu_BL_mat<-data.frame(row.names = otu_BL$OTUid,otu_BL[,-1])
    otu_Treat_mat<-data.frame(row.names =otu_Treat$OTUid,otu_Treat[,-1])
    otu_BL_mat<-data.frame(Samples=colnames(otu_BL_mat),t(otu_BL_mat))
    otu_Treat_mat<-data.frame(Samples=colnames(otu_Treat_mat),t(otu_Treat_mat))
    otu_BL_mat<-merge(dplyr::select(treat_BL,c(Samples,patientID)),otu_BL_mat,by="Samples")[,-1]
    otu_Treat_mat<-merge(dplyr::select(treat_BL,c(Samples,patientID)),otu_Treat_mat,by="Samples")[,-1]
    otu_BL_mat<-data.frame(row.names = otu_BL_mat$patientID,otu_BL_mat[,-1])%>%as.matrix()
    otu_Treat_mat<-data.frame(row.names = otu_Treat_mat$patientID,otu_Treat_mat[,-1])%>%as.matrix()

    colnames(otu_BL_mat)<-paste0(colnames(otu_BL_mat),"_BL")
    colnames(otu_Treat_mat)<-paste0(colnames(otu_Treat_mat),"_Treat")
    merged_cli_otu<-do.call(cbind,args = list(clinical_cd_mat,otu_BL_mat,otu_Treat_mat))
    merged_otu<-do.call(cbind,args = list(otu_BL_mat,otu_Treat_mat))
    cd_13<-data.frame(row.names = cd$patientID,cd[,-1])
    cd_13<-cd_13[which(rownames(cd_13)%in%rownames(otu_BL_mat)),]

    dfcol<-treat_BL[-which(duplicated(treat_BL$patientID)),]
    dfcol<-data.frame(row.names = dfcol$patientID,Group=dfcol$Response)
    df<-cbind(dfcol,merged_cli_otu)
    df$Group<-ifelse(df$Group=="R",1,2)
    X<-df[,-c(1:5)]%>%as.matrix()
    X<-data.frame(row.names =rownames(X),apply(X, 2, as.numeric))

    AE <- autoencoder(X, c(100,10,100), random.seed=1234,
                      loss.type = 'huber',drop.last = F,
                      activ.functions = c('tanh','linear','tanh'),
                      batch.size =3, optim.type = 'rmsprop',
                      n.epochs = 1000, val.prop = 0)

    recX <- reconstruct(AE, X)
    sort(recX$anomaly_scores, decreasing = TRUE)[1:5]
    AE_df<-recX$reconstructed
    rownames(AE_df)<-rownames(X)

    ## Artificial Neural Network: 
    ##   Layer - 1842 nodes - input 
    ##   Layer - 100 nodes - tanh 
    ##   Layer - 10 nodes - linear 
    ##   Layer - 100 nodes - tanh 
    ##   Layer - 1842 nodes - linear 
    ## With huber loss and RMSprop optimizer 
    ## Training progress:
    ## [|-------------------------------------------------] 0% - Training loss: 1025.27[|-------------------------------------------------] 1% - Training loss: 1034.31[|-------------------------------------------------] 1% - Training loss: 668.349[|-------------------------------------------------] 1% - Training loss: 644.843[|-------------------------------------------------] 1% - Training loss: 595.502[|-------------------------------------------------] 1% - Training loss: 587.038[+|------------------------------------------------] 2% - Training loss: 619.841[+|------------------------------------------------] 2% - Training loss: 631.189[+|------------------------------------------------] 2% - Training loss: 548.483[+|------------------------------------------------] 2% - Training loss: 717.025[+|------------------------------------------------] 2% - Training loss: 637.72[+|------------------------------------------------] 3% - Training loss: 564.308[+|------------------------------------------------] 3% - Training loss: 749.362[+|------------------------------------------------] 3% - Training loss: 603.041[+|------------------------------------------------] 3% - Training loss: 626.58[+|------------------------------------------------] 3% - Training loss: 552.343[++|-----------------------------------------------] 4% - Training loss: 575.603[++|-----------------------------------------------] 4% - Training loss: 550.654[++|-----------------------------------------------] 4% - Training loss: 647.091[++|-----------------------------------------------] 4% - Training loss: 551.128[++|-----------------------------------------------] 4% - Training loss: 626.721[++|-----------------------------------------------] 5% - Training loss: 717.774[++|-----------------------------------------------] 5% - Training loss: 531.11[++|-----------------------------------------------] 5% - Training loss: 512.017[++|-----------------------------------------------] 5% - Training loss: 505.012[++|-----------------------------------------------] 5% - Training loss: 511.269[+++|----------------------------------------------] 6% - Training loss: 566.548[+++|----------------------------------------------] 6% - Training loss: 508.907[+++|----------------------------------------------] 6% - Training loss: 640.721[+++|----------------------------------------------] 6% - Training loss: 760.594[+++|----------------------------------------------] 6% - Training loss: 469.851[+++|----------------------------------------------] 7% - Training loss: 605.276[+++|----------------------------------------------] 7% - Training loss: 659.755[+++|----------------------------------------------] 7% - Training loss: 505.925[+++|----------------------------------------------] 7% - Training loss: 533.187[+++|----------------------------------------------] 7% - Training loss: 462.351[++++|---------------------------------------------] 8% - Training loss: 589.429[++++|---------------------------------------------] 8% - Training loss: 526.871[++++|---------------------------------------------] 8% - Training loss: 560.092[++++|---------------------------------------------] 8% - Training loss: 535.164[++++|---------------------------------------------] 8% - Training loss: 524.665[++++|---------------------------------------------] 9% - Training loss: 570.58[++++|---------------------------------------------] 9% - Training loss: 503.023[++++|---------------------------------------------] 9% - Training loss: 473.65[++++|---------------------------------------------] 9% - Training loss: 540.612[++++|---------------------------------------------] 9% - Training loss: 535.134[+++++|--------------------------------------------] 10% - Training loss: 529.583[+++++|--------------------------------------------] 10% - Training loss: 540.425[+++++|--------------------------------------------] 10% - Training loss: 498.882[+++++|--------------------------------------------] 10% - Training loss: 505.548[+++++|--------------------------------------------] 10% - Training loss: 514.791[+++++|--------------------------------------------] 11% - Training loss: 571.953[+++++|--------------------------------------------] 11% - Training loss: 473.514[+++++|--------------------------------------------] 11% - Training loss: 492.206[+++++|--------------------------------------------] 11% - Training loss: 478.229[+++++|--------------------------------------------] 11% - Training loss: 522.755[++++++|-------------------------------------------] 12% - Training loss: 511.931[++++++|-------------------------------------------] 12% - Training loss: 520.19[++++++|-------------------------------------------] 12% - Training loss: 555.243[++++++|-------------------------------------------] 12% - Training loss: 521.611[++++++|-------------------------------------------] 12% - Training loss: 467.199[++++++|-------------------------------------------] 13% - Training loss: 494.793[++++++|-------------------------------------------] 13% - Training loss: 530.375[++++++|-------------------------------------------] 13% - Training loss: 452.258[++++++|-------------------------------------------] 13% - Training loss: 493.034[++++++|-------------------------------------------] 13% - Training loss: 436.233[+++++++|------------------------------------------] 14% - Training loss: 481.335[+++++++|------------------------------------------] 14% - Training loss: 412.63[+++++++|------------------------------------------] 14% - Training loss: 437.734[+++++++|------------------------------------------] 14% - Training loss: 513.531[+++++++|------------------------------------------] 14% - Training loss: 467.725[+++++++|------------------------------------------] 15% - Training loss: 483.148[+++++++|------------------------------------------] 15% - Training loss: 450.432[+++++++|------------------------------------------] 15% - Training loss: 454.075[+++++++|------------------------------------------] 15% - Training loss: 425.79[+++++++|------------------------------------------] 15% - Training loss: 391.579[++++++++|-----------------------------------------] 16% - Training loss: 422.526[++++++++|-----------------------------------------] 16% - Training loss: 458.47[++++++++|-----------------------------------------] 16% - Training loss: 470.08[++++++++|-----------------------------------------] 16% - Training loss: 350.489[++++++++|-----------------------------------------] 16% - Training loss: 409.133[++++++++|-----------------------------------------] 17% - Training loss: 443.787[++++++++|-----------------------------------------] 17% - Training loss: 313.49[++++++++|-----------------------------------------] 17% - Training loss: 419.159[++++++++|-----------------------------------------] 17% - Training loss: 345.316[++++++++|-----------------------------------------] 17% - Training loss: 492.9[+++++++++|----------------------------------------] 18% - Training loss: 391.358[+++++++++|----------------------------------------] 18% - Training loss: 409.817[+++++++++|----------------------------------------] 18% - Training loss: 386.469[+++++++++|----------------------------------------] 18% - Training loss: 393.461[+++++++++|----------------------------------------] 18% - Training loss: 313.063[+++++++++|----------------------------------------] 19% - Training loss: 404.768[+++++++++|----------------------------------------] 19% - Training loss: 443.543[+++++++++|----------------------------------------] 19% - Training loss: 344.074[+++++++++|----------------------------------------] 19% - Training loss: 435.828[+++++++++|----------------------------------------] 19% - Training loss: 415.443[++++++++++|---------------------------------------] 20% - Training loss: 334.648[++++++++++|---------------------------------------] 20% - Training loss: 328.711[++++++++++|---------------------------------------] 20% - Training loss: 371.651[++++++++++|---------------------------------------] 20% - Training loss: 402.971[++++++++++|---------------------------------------] 20% - Training loss: 311.372[++++++++++|---------------------------------------] 21% - Training loss: 305.507[++++++++++|---------------------------------------] 21% - Training loss: 326.427[++++++++++|---------------------------------------] 21% - Training loss: 378.418[++++++++++|---------------------------------------] 21% - Training loss: 345.642[++++++++++|---------------------------------------] 21% - Training loss: 357.133[+++++++++++|--------------------------------------] 22% - Training loss: 414.63[+++++++++++|--------------------------------------] 22% - Training loss: 390.511[+++++++++++|--------------------------------------] 22% - Training loss: 348.29[+++++++++++|--------------------------------------] 22% - Training loss: 273.643[+++++++++++|--------------------------------------] 22% - Training loss: 315.804[+++++++++++|--------------------------------------] 23% - Training loss: 390.819[+++++++++++|--------------------------------------] 23% - Training loss: 260.544[+++++++++++|--------------------------------------] 23% - Training loss: 295.347[+++++++++++|--------------------------------------] 23% - Training loss: 258.274[+++++++++++|--------------------------------------] 23% - Training loss: 392.595[++++++++++++|-------------------------------------] 24% - Training loss: 221.186[++++++++++++|-------------------------------------] 24% - Training loss: 381.174[++++++++++++|-------------------------------------] 24% - Training loss: 292.964[++++++++++++|-------------------------------------] 24% - Training loss: 274.416[++++++++++++|-------------------------------------] 24% - Training loss: 372.435[++++++++++++|-------------------------------------] 25% - Training loss: 312.594[++++++++++++|-------------------------------------] 25% - Training loss: 305.727[++++++++++++|-------------------------------------] 25% - Training loss: 349.658[++++++++++++|-------------------------------------] 25% - Training loss: 344.504[++++++++++++|-------------------------------------] 25% - Training loss: 269.886[+++++++++++++|------------------------------------] 26% - Training loss: 338.366[+++++++++++++|------------------------------------] 26% - Training loss: 234.814[+++++++++++++|------------------------------------] 26% - Training loss: 265.163[+++++++++++++|------------------------------------] 26% - Training loss: 310.502[+++++++++++++|------------------------------------] 26% - Training loss: 363.73[+++++++++++++|------------------------------------] 27% - Training loss: 347.353[+++++++++++++|------------------------------------] 27% - Training loss: 256.573[+++++++++++++|------------------------------------] 27% - Training loss: 344.339[+++++++++++++|------------------------------------] 27% - Training loss: 334.98[+++++++++++++|------------------------------------] 27% - Training loss: 230.268[++++++++++++++|-----------------------------------] 28% - Training loss: 287.583[++++++++++++++|-----------------------------------] 28% - Training loss: 371.094[++++++++++++++|-----------------------------------] 28% - Training loss: 276.085[++++++++++++++|-----------------------------------] 28% - Training loss: 238.322[++++++++++++++|-----------------------------------] 28% - Training loss: 340.293[++++++++++++++|-----------------------------------] 29% - Training loss: 270.293[++++++++++++++|-----------------------------------] 29% - Training loss: 191.721[++++++++++++++|-----------------------------------] 29% - Training loss: 262.632[++++++++++++++|-----------------------------------] 29% - Training loss: 337.277[++++++++++++++|-----------------------------------] 29% - Training loss: 203.723[+++++++++++++++|----------------------------------] 30% - Training loss: 267.245[+++++++++++++++|----------------------------------] 30% - Training loss: 255.607[+++++++++++++++|----------------------------------] 30% - Training loss: 288.414[+++++++++++++++|----------------------------------] 30% - Training loss: 227.188[+++++++++++++++|----------------------------------] 30% - Training loss: 154.374[+++++++++++++++|----------------------------------] 31% - Training loss: 108.564[+++++++++++++++|----------------------------------] 31% - Training loss: 323.032[+++++++++++++++|----------------------------------] 31% - Training loss: 64.1306[+++++++++++++++|----------------------------------] 31% - Training loss: 193.766[+++++++++++++++|----------------------------------] 31% - Training loss: 340.756[++++++++++++++++|---------------------------------] 32% - Training loss: 252.552[++++++++++++++++|---------------------------------] 32% - Training loss: 117.21[++++++++++++++++|---------------------------------] 32% - Training loss: 134.794[++++++++++++++++|---------------------------------] 32% - Training loss: 275.739[++++++++++++++++|---------------------------------] 32% - Training loss: 288.892[++++++++++++++++|---------------------------------] 33% - Training loss: 269.453[++++++++++++++++|---------------------------------] 33% - Training loss: 169.448[++++++++++++++++|---------------------------------] 33% - Training loss: 149.895[++++++++++++++++|---------------------------------] 33% - Training loss: 201.661[++++++++++++++++|---------------------------------] 33% - Training loss: 332.972[+++++++++++++++++|--------------------------------] 34% - Training loss: 276.197[+++++++++++++++++|--------------------------------] 34% - Training loss: 215.793[+++++++++++++++++|--------------------------------] 34% - Training loss: 299.403[+++++++++++++++++|--------------------------------] 34% - Training loss: 201.819[+++++++++++++++++|--------------------------------] 34% - Training loss: 200.113[+++++++++++++++++|--------------------------------] 35% - Training loss: 165.412[+++++++++++++++++|--------------------------------] 35% - Training loss: 177.856[+++++++++++++++++|--------------------------------] 35% - Training loss: 357.373[+++++++++++++++++|--------------------------------] 35% - Training loss: 99.4932[+++++++++++++++++|--------------------------------] 35% - Training loss: 257.492[++++++++++++++++++|-------------------------------] 36% - Training loss: 244.675[++++++++++++++++++|-------------------------------] 36% - Training loss: 246.404[++++++++++++++++++|-------------------------------] 36% - Training loss: 326.236[++++++++++++++++++|-------------------------------] 36% - Training loss: 203.303[++++++++++++++++++|-------------------------------] 36% - Training loss: 157.16[++++++++++++++++++|-------------------------------] 37% - Training loss: 186.659[++++++++++++++++++|-------------------------------] 37% - Training loss: 198.356[++++++++++++++++++|-------------------------------] 37% - Training loss: 194.532[++++++++++++++++++|-------------------------------] 37% - Training loss: 154.441[++++++++++++++++++|-------------------------------] 37% - Training loss: 331.74[+++++++++++++++++++|------------------------------] 38% - Training loss: 124.804[+++++++++++++++++++|------------------------------] 38% - Training loss: 312.955[+++++++++++++++++++|------------------------------] 38% - Training loss: 49.3305[+++++++++++++++++++|------------------------------] 38% - Training loss: 157.274[+++++++++++++++++++|------------------------------] 38% - Training loss: 103.36[+++++++++++++++++++|------------------------------] 39% - Training loss: 96.6496[+++++++++++++++++++|------------------------------] 39% - Training loss: 109.749[+++++++++++++++++++|------------------------------] 39% - Training loss: 107.131[+++++++++++++++++++|------------------------------] 39% - Training loss: 176.448[+++++++++++++++++++|------------------------------] 39% - Training loss: 86.0809[++++++++++++++++++++|-----------------------------] 40% - Training loss: 206.662[++++++++++++++++++++|-----------------------------] 40% - Training loss: 164.784[++++++++++++++++++++|-----------------------------] 40% - Training loss: 149.493[++++++++++++++++++++|-----------------------------] 40% - Training loss: 283.796[++++++++++++++++++++|-----------------------------] 40% - Training loss: 136.899[++++++++++++++++++++|-----------------------------] 41% - Training loss: 123.103[++++++++++++++++++++|-----------------------------] 41% - Training loss: 191.125[++++++++++++++++++++|-----------------------------] 41% - Training loss: 187.94[++++++++++++++++++++|-----------------------------] 41% - Training loss: 78.5216[++++++++++++++++++++|-----------------------------] 41% - Training loss: 38.5104[+++++++++++++++++++++|----------------------------] 42% - Training loss: 235.359[+++++++++++++++++++++|----------------------------] 42% - Training loss: 101.877[+++++++++++++++++++++|----------------------------] 42% - Training loss: 187.511[+++++++++++++++++++++|----------------------------] 42% - Training loss: 208.811[+++++++++++++++++++++|----------------------------] 42% - Training loss: 163.793[+++++++++++++++++++++|----------------------------] 43% - Training loss: 291.816[+++++++++++++++++++++|----------------------------] 43% - Training loss: 297.636[+++++++++++++++++++++|----------------------------] 43% - Training loss: 153.392[+++++++++++++++++++++|----------------------------] 43% - Training loss: 127.895[+++++++++++++++++++++|----------------------------] 43% - Training loss: 178.237[++++++++++++++++++++++|---------------------------] 44% - Training loss: 163.362[++++++++++++++++++++++|---------------------------] 44% - Training loss: 101.008[++++++++++++++++++++++|---------------------------] 44% - Training loss: 186.733[++++++++++++++++++++++|---------------------------] 44% - Training loss: 46.0741[++++++++++++++++++++++|---------------------------] 44% - Training loss: 293.772[++++++++++++++++++++++|---------------------------] 45% - Training loss: 39.4864[++++++++++++++++++++++|---------------------------] 45% - Training loss: 41.7005[++++++++++++++++++++++|---------------------------] 45% - Training loss: 65.1735[++++++++++++++++++++++|---------------------------] 45% - Training loss: 73.8696[++++++++++++++++++++++|---------------------------] 45% - Training loss: 134.954[+++++++++++++++++++++++|--------------------------] 46% - Training loss: 176.301[+++++++++++++++++++++++|--------------------------] 46% - Training loss: 73.5575[+++++++++++++++++++++++|--------------------------] 46% - Training loss: 170.215[+++++++++++++++++++++++|--------------------------] 46% - Training loss: 54.0268[+++++++++++++++++++++++|--------------------------] 46% - Training loss: 27.1519[+++++++++++++++++++++++|--------------------------] 47% - Training loss: 178.435[+++++++++++++++++++++++|--------------------------] 47% - Training loss: 173.114[+++++++++++++++++++++++|--------------------------] 47% - Training loss: 29.4907[+++++++++++++++++++++++|--------------------------] 47% - Training loss: 31.3243[+++++++++++++++++++++++|--------------------------] 47% - Training loss: 177.388[++++++++++++++++++++++++|-------------------------] 48% - Training loss: 21.3263[++++++++++++++++++++++++|-------------------------] 48% - Training loss: 186.834[++++++++++++++++++++++++|-------------------------] 48% - Training loss: 127.414[++++++++++++++++++++++++|-------------------------] 48% - Training loss: 148.887[++++++++++++++++++++++++|-------------------------] 48% - Training loss: 45.8133[++++++++++++++++++++++++|-------------------------] 49% - Training loss: 152.163[++++++++++++++++++++++++|-------------------------] 49% - Training loss: 157.739[++++++++++++++++++++++++|-------------------------] 49% - Training loss: 168.162[++++++++++++++++++++++++|-------------------------] 49% - Training loss: 35.5936[++++++++++++++++++++++++|-------------------------] 49% - Training loss: 13.2308[+++++++++++++++++++++++++|------------------------] 50% - Training loss: 26.0673[+++++++++++++++++++++++++|------------------------] 50% - Training loss: 155.52[+++++++++++++++++++++++++|------------------------] 50% - Training loss: 135.77[+++++++++++++++++++++++++|------------------------] 50% - Training loss: 167.951[+++++++++++++++++++++++++|------------------------] 50% - Training loss: 11.6213[+++++++++++++++++++++++++|------------------------] 51% - Training loss: 33.8934[+++++++++++++++++++++++++|------------------------] 51% - Training loss: 152.599[+++++++++++++++++++++++++|------------------------] 51% - Training loss: 164.08[+++++++++++++++++++++++++|------------------------] 51% - Training loss: 14.435[+++++++++++++++++++++++++|------------------------] 51% - Training loss: 151.148[++++++++++++++++++++++++++|-----------------------] 52% - Training loss: 162.477[++++++++++++++++++++++++++|-----------------------] 52% - Training loss: 16.3068[++++++++++++++++++++++++++|-----------------------] 52% - Training loss: 144.556[++++++++++++++++++++++++++|-----------------------] 52% - Training loss: 147.188[++++++++++++++++++++++++++|-----------------------] 52% - Training loss: 148.669[++++++++++++++++++++++++++|-----------------------] 53% - Training loss: 16.1742[++++++++++++++++++++++++++|-----------------------] 53% - Training loss: 19.8634[++++++++++++++++++++++++++|-----------------------] 53% - Training loss: 20.9455[++++++++++++++++++++++++++|-----------------------] 53% - Training loss: 13.8196[++++++++++++++++++++++++++|-----------------------] 53% - Training loss: 271.515[+++++++++++++++++++++++++++|----------------------] 54% - Training loss: 148.105[+++++++++++++++++++++++++++|----------------------] 54% - Training loss: 124.273[+++++++++++++++++++++++++++|----------------------] 54% - Training loss: 144.736[+++++++++++++++++++++++++++|----------------------] 54% - Training loss: 128.955[+++++++++++++++++++++++++++|----------------------] 54% - Training loss: 155.626[+++++++++++++++++++++++++++|----------------------] 55% - Training loss: 145.48[+++++++++++++++++++++++++++|----------------------] 55% - Training loss: 134.124[+++++++++++++++++++++++++++|----------------------] 55% - Training loss: 147.24[+++++++++++++++++++++++++++|----------------------] 55% - Training loss: 137.364[+++++++++++++++++++++++++++|----------------------] 55% - Training loss: 7.05745[++++++++++++++++++++++++++++|---------------------] 56% - Training loss: 17.0802[++++++++++++++++++++++++++++|---------------------] 56% - Training loss: 8.77955[++++++++++++++++++++++++++++|---------------------] 56% - Training loss: 131.891[++++++++++++++++++++++++++++|---------------------] 56% - Training loss: 12.4034[++++++++++++++++++++++++++++|---------------------] 56% - Training loss: 9.70541[++++++++++++++++++++++++++++|---------------------] 57% - Training loss: 129.456[++++++++++++++++++++++++++++|---------------------] 57% - Training loss: 8.26006[++++++++++++++++++++++++++++|---------------------] 57% - Training loss: 130.331[++++++++++++++++++++++++++++|---------------------] 57% - Training loss: 132.884[++++++++++++++++++++++++++++|---------------------] 57% - Training loss: 132.82[+++++++++++++++++++++++++++++|--------------------] 58% - Training loss: 135.188[+++++++++++++++++++++++++++++|--------------------] 58% - Training loss: 140.396[+++++++++++++++++++++++++++++|--------------------] 58% - Training loss: 121.899[+++++++++++++++++++++++++++++|--------------------] 58% - Training loss: 8.29209[+++++++++++++++++++++++++++++|--------------------] 58% - Training loss: 13.1184[+++++++++++++++++++++++++++++|--------------------] 59% - Training loss: 6.62837[+++++++++++++++++++++++++++++|--------------------] 59% - Training loss: 124.007[+++++++++++++++++++++++++++++|--------------------] 59% - Training loss: 123.159[+++++++++++++++++++++++++++++|--------------------] 59% - Training loss: 6.23702[+++++++++++++++++++++++++++++|--------------------] 59% - Training loss: 4.82803[++++++++++++++++++++++++++++++|-------------------] 60% - Training loss: 16.053[++++++++++++++++++++++++++++++|-------------------] 60% - Training loss: 5.80455[++++++++++++++++++++++++++++++|-------------------] 60% - Training loss: 6.71821[++++++++++++++++++++++++++++++|-------------------] 60% - Training loss: 237.987[++++++++++++++++++++++++++++++|-------------------] 60% - Training loss: 121.746[++++++++++++++++++++++++++++++|-------------------] 61% - Training loss: 6.15383[++++++++++++++++++++++++++++++|-------------------] 61% - Training loss: 3.79168[++++++++++++++++++++++++++++++|-------------------] 61% - Training loss: 6.96344[++++++++++++++++++++++++++++++|-------------------] 61% - Training loss: 118.034[++++++++++++++++++++++++++++++|-------------------] 61% - Training loss: 119.667[+++++++++++++++++++++++++++++++|------------------] 62% - Training loss: 16.4797[+++++++++++++++++++++++++++++++|------------------] 62% - Training loss: 9.62117[+++++++++++++++++++++++++++++++|------------------] 62% - Training loss: 119.3[+++++++++++++++++++++++++++++++|------------------] 62% - Training loss: 119.737[+++++++++++++++++++++++++++++++|------------------] 62% - Training loss: 231.132[+++++++++++++++++++++++++++++++|------------------] 63% - Training loss: 4.34861[+++++++++++++++++++++++++++++++|------------------] 63% - Training loss: 10.6091[+++++++++++++++++++++++++++++++|------------------] 63% - Training loss: 2.74788[+++++++++++++++++++++++++++++++|------------------] 63% - Training loss: 118.386[+++++++++++++++++++++++++++++++|------------------] 63% - Training loss: 4.70794[++++++++++++++++++++++++++++++++|-----------------] 64% - Training loss: 118.51[++++++++++++++++++++++++++++++++|-----------------] 64% - Training loss: 223.958[++++++++++++++++++++++++++++++++|-----------------] 64% - Training loss: 2.54312[++++++++++++++++++++++++++++++++|-----------------] 64% - Training loss: 126.488[++++++++++++++++++++++++++++++++|-----------------] 64% - Training loss: 6.80608[++++++++++++++++++++++++++++++++|-----------------] 65% - Training loss: 4.66251[++++++++++++++++++++++++++++++++|-----------------] 65% - Training loss: 7.1503[++++++++++++++++++++++++++++++++|-----------------] 65% - Training loss: 215.062[++++++++++++++++++++++++++++++++|-----------------] 65% - Training loss: 119.673[++++++++++++++++++++++++++++++++|-----------------] 65% - Training loss: 3.72621[+++++++++++++++++++++++++++++++++|----------------] 66% - Training loss: 3.87819[+++++++++++++++++++++++++++++++++|----------------] 66% - Training loss: 5.00474[+++++++++++++++++++++++++++++++++|----------------] 66% - Training loss: 207.829[+++++++++++++++++++++++++++++++++|----------------] 66% - Training loss: 88.1458[+++++++++++++++++++++++++++++++++|----------------] 66% - Training loss: 1.68651[+++++++++++++++++++++++++++++++++|----------------] 67% - Training loss: 119.529[+++++++++++++++++++++++++++++++++|----------------] 67% - Training loss: 5.19894[+++++++++++++++++++++++++++++++++|----------------] 67% - Training loss: 122.04[+++++++++++++++++++++++++++++++++|----------------] 67% - Training loss: 85.5819[+++++++++++++++++++++++++++++++++|----------------] 67% - Training loss: 2.06475[++++++++++++++++++++++++++++++++++|---------------] 68% - Training loss: 120.879[++++++++++++++++++++++++++++++++++|---------------] 68% - Training loss: 4.8402[++++++++++++++++++++++++++++++++++|---------------] 68% - Training loss: 117.901[++++++++++++++++++++++++++++++++++|---------------] 68% - Training loss: 74.3885[++++++++++++++++++++++++++++++++++|---------------] 68% - Training loss: 117.345[++++++++++++++++++++++++++++++++++|---------------] 69% - Training loss: 122.33[++++++++++++++++++++++++++++++++++|---------------] 69% - Training loss: 2.59963[++++++++++++++++++++++++++++++++++|---------------] 69% - Training loss: 3.88795[++++++++++++++++++++++++++++++++++|---------------] 69% - Training loss: 5.07086[++++++++++++++++++++++++++++++++++|---------------] 69% - Training loss: 2.75841[+++++++++++++++++++++++++++++++++++|--------------] 70% - Training loss: 3.50614[+++++++++++++++++++++++++++++++++++|--------------] 70% - Training loss: 118.319[+++++++++++++++++++++++++++++++++++|--------------] 70% - Training loss: 120.619[+++++++++++++++++++++++++++++++++++|--------------] 70% - Training loss: 61.3241[+++++++++++++++++++++++++++++++++++|--------------] 70% - Training loss: 4.67172[+++++++++++++++++++++++++++++++++++|--------------] 71% - Training loss: 5.70963[+++++++++++++++++++++++++++++++++++|--------------] 71% - Training loss: 3.84929[+++++++++++++++++++++++++++++++++++|--------------] 71% - Training loss: 6.88515[+++++++++++++++++++++++++++++++++++|--------------] 71% - Training loss: 5.80001[+++++++++++++++++++++++++++++++++++|--------------] 71% - Training loss: 5.26713[++++++++++++++++++++++++++++++++++++|-------------] 72% - Training loss: 51.731[++++++++++++++++++++++++++++++++++++|-------------] 72% - Training loss: 113.605[++++++++++++++++++++++++++++++++++++|-------------] 72% - Training loss: 2.53777[++++++++++++++++++++++++++++++++++++|-------------] 72% - Training loss: 46.9808[++++++++++++++++++++++++++++++++++++|-------------] 72% - Training loss: 6.9013[++++++++++++++++++++++++++++++++++++|-------------] 73% - Training loss: 4.18959[++++++++++++++++++++++++++++++++++++|-------------] 73% - Training loss: 41.1686[++++++++++++++++++++++++++++++++++++|-------------] 73% - Training loss: 112.166[++++++++++++++++++++++++++++++++++++|-------------] 73% - Training loss: 37.5704[++++++++++++++++++++++++++++++++++++|-------------] 73% - Training loss: 112.525[+++++++++++++++++++++++++++++++++++++|------------] 74% - Training loss: 2.08163[+++++++++++++++++++++++++++++++++++++|------------] 74% - Training loss: 2.36504[+++++++++++++++++++++++++++++++++++++|------------] 74% - Training loss: 1.69734[+++++++++++++++++++++++++++++++++++++|------------] 74% - Training loss: 34.4197[+++++++++++++++++++++++++++++++++++++|------------] 74% - Training loss: 6.18053[+++++++++++++++++++++++++++++++++++++|------------] 75% - Training loss: 6.30497[+++++++++++++++++++++++++++++++++++++|------------] 75% - Training loss: 1.33707[+++++++++++++++++++++++++++++++++++++|------------] 75% - Training loss: 5.68108[+++++++++++++++++++++++++++++++++++++|------------] 75% - Training loss: 108.561[+++++++++++++++++++++++++++++++++++++|------------] 75% - Training loss: 24.2783[++++++++++++++++++++++++++++++++++++++|-----------] 76% - Training loss: 4.44008[++++++++++++++++++++++++++++++++++++++|-----------] 76% - Training loss: 4.27104[++++++++++++++++++++++++++++++++++++++|-----------] 76% - Training loss: 2.62892[++++++++++++++++++++++++++++++++++++++|-----------] 76% - Training loss: 18.9342[++++++++++++++++++++++++++++++++++++++|-----------] 76% - Training loss: 120.852[++++++++++++++++++++++++++++++++++++++|-----------] 77% - Training loss: 2.61904[++++++++++++++++++++++++++++++++++++++|-----------] 77% - Training loss: 4.60423[++++++++++++++++++++++++++++++++++++++|-----------] 77% - Training loss: 1.97257[++++++++++++++++++++++++++++++++++++++|-----------] 77% - Training loss: 1.86942[++++++++++++++++++++++++++++++++++++++|-----------] 77% - Training loss: 2.49726[+++++++++++++++++++++++++++++++++++++++|----------] 78% - Training loss: 3.15326[+++++++++++++++++++++++++++++++++++++++|----------] 78% - Training loss: 4.39018[+++++++++++++++++++++++++++++++++++++++|----------] 78% - Training loss: 2.43831[+++++++++++++++++++++++++++++++++++++++|----------] 78% - Training loss: 2.1044[+++++++++++++++++++++++++++++++++++++++|----------] 78% - Training loss: 3.98181[+++++++++++++++++++++++++++++++++++++++|----------] 79% - Training loss: 14.4678[+++++++++++++++++++++++++++++++++++++++|----------] 79% - Training loss: 14.089[+++++++++++++++++++++++++++++++++++++++|----------] 79% - Training loss: 91.7111[+++++++++++++++++++++++++++++++++++++++|----------] 79% - Training loss: 10.4581[+++++++++++++++++++++++++++++++++++++++|----------] 79% - Training loss: 93.556[++++++++++++++++++++++++++++++++++++++++|---------] 80% - Training loss: 2.39286[++++++++++++++++++++++++++++++++++++++++|---------] 80% - Training loss: 8.07278[++++++++++++++++++++++++++++++++++++++++|---------] 80% - Training loss: 1.44608[++++++++++++++++++++++++++++++++++++++++|---------] 80% - Training loss: 4.77497[++++++++++++++++++++++++++++++++++++++++|---------] 80% - Training loss: 7.37108[++++++++++++++++++++++++++++++++++++++++|---------] 81% - Training loss: 2.35897[++++++++++++++++++++++++++++++++++++++++|---------] 81% - Training loss: 1.22853[++++++++++++++++++++++++++++++++++++++++|---------] 81% - Training loss: 2.34204[++++++++++++++++++++++++++++++++++++++++|---------] 81% - Training loss: 2.82745[++++++++++++++++++++++++++++++++++++++++|---------] 81% - Training loss: 1.77402[+++++++++++++++++++++++++++++++++++++++++|--------] 82% - Training loss: 6.30888[+++++++++++++++++++++++++++++++++++++++++|--------] 82% - Training loss: 5.43112[+++++++++++++++++++++++++++++++++++++++++|--------] 82% - Training loss: 1.13853[+++++++++++++++++++++++++++++++++++++++++|--------] 82% - Training loss: 80.7382[+++++++++++++++++++++++++++++++++++++++++|--------] 82% - Training loss: 1.01559[+++++++++++++++++++++++++++++++++++++++++|--------] 83% - Training loss: 1.06918[+++++++++++++++++++++++++++++++++++++++++|--------] 83% - Training loss: 2.96332[+++++++++++++++++++++++++++++++++++++++++|--------] 83% - Training loss: 7.48202[+++++++++++++++++++++++++++++++++++++++++|--------] 83% - Training loss: 2.83849[+++++++++++++++++++++++++++++++++++++++++|--------] 83% - Training loss: 7.01956[++++++++++++++++++++++++++++++++++++++++++|-------] 84% - Training loss: 0.750244[++++++++++++++++++++++++++++++++++++++++++|-------] 84% - Training loss: 0.778491[++++++++++++++++++++++++++++++++++++++++++|-------] 84% - Training loss: 1.34313[++++++++++++++++++++++++++++++++++++++++++|-------] 84% - Training loss: 68.2549[++++++++++++++++++++++++++++++++++++++++++|-------] 84% - Training loss: 5.7234[++++++++++++++++++++++++++++++++++++++++++|-------] 85% - Training loss: 3.21164[++++++++++++++++++++++++++++++++++++++++++|-------] 85% - Training loss: 6.92968[++++++++++++++++++++++++++++++++++++++++++|-------] 85% - Training loss: 5.90788[++++++++++++++++++++++++++++++++++++++++++|-------] 85% - Training loss: 0.866231[++++++++++++++++++++++++++++++++++++++++++|-------] 85% - Training loss: 60.9461[+++++++++++++++++++++++++++++++++++++++++++|------] 86% - Training loss: 3.27731[+++++++++++++++++++++++++++++++++++++++++++|------] 86% - Training loss: 1.1649[+++++++++++++++++++++++++++++++++++++++++++|------] 86% - Training loss: 0.669345[+++++++++++++++++++++++++++++++++++++++++++|------] 86% - Training loss: 55.6195[+++++++++++++++++++++++++++++++++++++++++++|------] 86% - Training loss: 4.27687[+++++++++++++++++++++++++++++++++++++++++++|------] 87% - Training loss: 2.51772[+++++++++++++++++++++++++++++++++++++++++++|------] 87% - Training loss: 52.0593[+++++++++++++++++++++++++++++++++++++++++++|------] 87% - Training loss: 51.6988[+++++++++++++++++++++++++++++++++++++++++++|------] 87% - Training loss: 3.10242[+++++++++++++++++++++++++++++++++++++++++++|------] 87% - Training loss: 50.0777[++++++++++++++++++++++++++++++++++++++++++++|-----] 88% - Training loss: 1.01781[++++++++++++++++++++++++++++++++++++++++++++|-----] 88% - Training loss: 2.14905[++++++++++++++++++++++++++++++++++++++++++++|-----] 88% - Training loss: 4.07732[++++++++++++++++++++++++++++++++++++++++++++|-----] 88% - Training loss: 3.81711[++++++++++++++++++++++++++++++++++++++++++++|-----] 88% - Training loss: 2.59069[++++++++++++++++++++++++++++++++++++++++++++|-----] 89% - Training loss: 44.1278[++++++++++++++++++++++++++++++++++++++++++++|-----] 89% - Training loss: 1.0759[++++++++++++++++++++++++++++++++++++++++++++|-----] 89% - Training loss: 2.57205[++++++++++++++++++++++++++++++++++++++++++++|-----] 89% - Training loss: 3.47654[++++++++++++++++++++++++++++++++++++++++++++|-----] 89% - Training loss: 1.62835[+++++++++++++++++++++++++++++++++++++++++++++|----] 90% - Training loss: 1.78939[+++++++++++++++++++++++++++++++++++++++++++++|----] 90% - Training loss: 2.77181[+++++++++++++++++++++++++++++++++++++++++++++|----] 90% - Training loss: 2.74509[+++++++++++++++++++++++++++++++++++++++++++++|----] 90% - Training loss: 3.9206[+++++++++++++++++++++++++++++++++++++++++++++|----] 90% - Training loss: 3.13931[+++++++++++++++++++++++++++++++++++++++++++++|----] 91% - Training loss: 1.28613[+++++++++++++++++++++++++++++++++++++++++++++|----] 91% - Training loss: 3.18916[+++++++++++++++++++++++++++++++++++++++++++++|----] 91% - Training loss: 1.46758[+++++++++++++++++++++++++++++++++++++++++++++|----] 91% - Training loss: 1.18631[+++++++++++++++++++++++++++++++++++++++++++++|----] 91% - Training loss: 1.62753[++++++++++++++++++++++++++++++++++++++++++++++|---] 92% - Training loss: 4.91494[++++++++++++++++++++++++++++++++++++++++++++++|---] 92% - Training loss: 2.02316[++++++++++++++++++++++++++++++++++++++++++++++|---] 92% - Training loss: 27.0669[++++++++++++++++++++++++++++++++++++++++++++++|---] 92% - Training loss: 1.54283[++++++++++++++++++++++++++++++++++++++++++++++|---] 92% - Training loss: 1.60811[++++++++++++++++++++++++++++++++++++++++++++++|---] 93% - Training loss: 1.03913[++++++++++++++++++++++++++++++++++++++++++++++|---] 93% - Training loss: 23.3595[++++++++++++++++++++++++++++++++++++++++++++++|---] 93% - Training loss: 23.1683[++++++++++++++++++++++++++++++++++++++++++++++|---] 93% - Training loss: 1.99476[++++++++++++++++++++++++++++++++++++++++++++++|---] 93% - Training loss: 22.9524[+++++++++++++++++++++++++++++++++++++++++++++++|--] 94% - Training loss: 4.15343[+++++++++++++++++++++++++++++++++++++++++++++++|--] 94% - Training loss: 2.46859[+++++++++++++++++++++++++++++++++++++++++++++++|--] 94% - Training loss: 18.2364[+++++++++++++++++++++++++++++++++++++++++++++++|--] 94% - Training loss: 2.9187[+++++++++++++++++++++++++++++++++++++++++++++++|--] 94% - Training loss: 3.34719[+++++++++++++++++++++++++++++++++++++++++++++++|--] 95% - Training loss: 1.91453[+++++++++++++++++++++++++++++++++++++++++++++++|--] 95% - Training loss: 4.98428[+++++++++++++++++++++++++++++++++++++++++++++++|--] 95% - Training loss: 1.01747[+++++++++++++++++++++++++++++++++++++++++++++++|--] 95% - Training loss: 2.32459[+++++++++++++++++++++++++++++++++++++++++++++++|--] 95% - Training loss: 1.9664[++++++++++++++++++++++++++++++++++++++++++++++++|-] 96% - Training loss: 2.28582[++++++++++++++++++++++++++++++++++++++++++++++++|-] 96% - Training loss: 2.05882[++++++++++++++++++++++++++++++++++++++++++++++++|-] 96% - Training loss: 12.4562[++++++++++++++++++++++++++++++++++++++++++++++++|-] 96% - Training loss: 1.63428[++++++++++++++++++++++++++++++++++++++++++++++++|-] 96% - Training loss: 2.50487[++++++++++++++++++++++++++++++++++++++++++++++++|-] 97% - Training loss: 1.87695[++++++++++++++++++++++++++++++++++++++++++++++++|-] 97% - Training loss: 10.4076[++++++++++++++++++++++++++++++++++++++++++++++++|-] 97% - Training loss: 10.0416[++++++++++++++++++++++++++++++++++++++++++++++++|-] 97% - Training loss: 2.00251[++++++++++++++++++++++++++++++++++++++++++++++++|-] 97% - Training loss: 0.886002[+++++++++++++++++++++++++++++++++++++++++++++++++|] 98% - Training loss: 7.73678[+++++++++++++++++++++++++++++++++++++++++++++++++|] 98% - Training loss: 1.2078[+++++++++++++++++++++++++++++++++++++++++++++++++|] 98% - Training loss: 9.18935[+++++++++++++++++++++++++++++++++++++++++++++++++|] 98% - Training loss: 1.69549[+++++++++++++++++++++++++++++++++++++++++++++++++|] 98% - Training loss: 7.07462[+++++++++++++++++++++++++++++++++++++++++++++++++|] 99% - Training loss: 2.27363[+++++++++++++++++++++++++++++++++++++++++++++++++|] 99% - Training loss: 1.83837[+++++++++++++++++++++++++++++++++++++++++++++++++|] 99% - Training loss: 0.713675[+++++++++++++++++++++++++++++++++++++++++++++++++|] 99% - Training loss: 0.998619[+++++++++++++++++++++++++++++++++++++++++++++++++|] 99% - Training loss: 5.24098[++++++++++++++++++++++++++++++++++++++++++++++++++] 100% - Training loss: 1.46765[++++++++++++++++++++++++++++++++++++++++++++++++++] 100% - Training loss: 5.01243[++++++++++++++++++++++++++++++++++++++++++++++++++] 100% - Training loss: 4.83362[++++++++++++++++++++++++++++++++++++++++++++++++++] 100% - Training loss: 2.846[++++++++++++++++++++++++++++++++++++++++++++++++++] 100% - Training loss: 4.59047
    ## [1] 7195.8612 1982.1491 1710.7365  670.8924  660.1566

3.1 K-means clustering
----------------------

### 3.1.1 Choose k value

    fviz_nbclust(
      AE_df, 
      kmeans, 
      k.max = 10,
      method = "wss",
      verbose = FALSE)

<img src="AutoDecoder_files/figure-markdown_strict/unnamed-chunk-7-1.png" width="30%" style="display: block; margin: auto;" />

### 3.1.2 Plot K-means clustering

    km.res <- kmeans(scale(AE_df),6, nstart = 25)
    fviz_cluster(km.res, data = AE_df,repel = T,labelsize = 12,
                # palette = c("#00AFBB", "#E7B800", "black","red","blue","gray"),
                 ggtheme = theme_minimal(),
                 main = "Kmeans Clustering Plot")

<img src="AutoDecoder_files/figure-markdown_strict/unnamed-chunk-8-1.png" width="30%" style="display: block; margin: auto;" />

### 3.1.3 PFS survival based on Autoencoder data

    cluster<-data.frame(km.res$cluster)
    cluster$patientID<-rownames(cluster)
    cluster$km.res.cluster<-paste0("km",cluster$km.res.cluster)
    stat_num<-cluster%>%group_by(km.res.cluster)%>%summarise(Num=n())
    cluster$km.res.cluster<-ifelse(cluster$km.res.cluster==stat_num$km.res.cluster[which(stat_num$Num>5)],"km1","km2")
    new<-merge(cluster,clinical,by="patientID")

    fit1<-survfit(Surv(PFStime,PFS) ~ km.res.cluster,
                  data = new)
    fit1
    p1<-ggsurvplot(fit1, data=new,pval.method = T,add.all = T,
                   tables.theme = theme_classic2(base_size = 8),
                   palette  = c("black", "#00AFBB","#E7B800"),
                   risk.table = T,
                   pval = TRUE,
                   legend.title="K-means",
                   risk.table.col = "strata",
                   surv.median.line = "hv",
                   risk.table.y.text.col = T,
                   risk.table.y.text = F )
    p1

<img src="AutoDecoder_files/figure-markdown_strict/unnamed-chunk-10-1.png" width="30%" style="display: block; margin: auto;" />

    ## Call: survfit(formula = Surv(PFStime, PFS) ~ km.res.cluster, data = new)
    ## 
    ##                    n events median 0.95LCL 0.95UCL
    ## km.res.cluster=km1 7      7   1.97    1.67      NA
    ## km.res.cluster=km2 6      4   5.27    2.30      NA

3.2 Representative features of the k-means clusters
---------------------------------------------------

    otu_info<-fread("../Data/Data/OTUtabale_regaStool.csv")[,c(1:7)]
    colnames(otu_info)[7]="OTUid"
    otu_info$Taxonomy<-paste(otu_info$Phylum,otu_info$Class,otu_info$Order,otu_info$Family,otu_info$Genus,otu_info$Species,sep = "|")
    otu_info<-otu_info[,7:8]
    ano_data<-X
    ano_data$patientID=rownames(ano_data)
    ano_data<-merge(cluster,ano_data,by="patientID")
    data<-ano_data
    data<-data[,-1]
    lev<-unique(data$km.res.cluster)
    f <- factor(data$km.res.cluster, levels=lev) 
    design <- model.matrix(~0+f)
    colnames(design) <- lev
    eset<-dplyr::select(data,-km.res.cluster)
    eset<-data.frame(t(eset))
    #eset<-data.frame(apply(eset, 2, av))
    cont.wt <- makeContrasts("km1-km2",
                             levels=design) 
    fit <- lmFit(eset, design)
    fit2 <- contrasts.fit(fit, cont.wt) 
    fit2 <- eBayes(fit2) 
    tT=topTable(fit2, adjust="BH",sort.by="logFC",n=Inf)
    tT = subset(tT, select=c("adj.P.Val","P.Value","logFC"))
    colnames(tT)=c("FDR","P.Value","logFC")
    range(tT$logFC)
    limma_res<-filter(tT,P.Value<=0.05&abs(logFC)>1)
    limma_res$Factor<-rownames(limma_res)
    limma_res$OTUid<-rownames(limma_res)
    limma_res$OTUid<-gsub("_BL","",limma_res$OTUid)
    limma_res$OTUid<-gsub("_Treat","",limma_res$OTUid)
    limma_res$Group<-rownames(limma_res)
    limma_res$Group<-gsub("^.*_","",limma_res$Group)
    limma_res<-merge(limma_res,otu_info,by="OTUid")
    limma_res$Taxonomy<-gsub("d__Bacteria;k__norank_d__Bacteria;p__","",limma_res$Taxonomy)

    ## [1] -5276.452  5715.119

3.3 Table of representative OTUs for the k-means clusters
---------------------------------------------------------

    knitr::kable(limma_res)

<table>
<colgroup>
<col style="width: 3%" />
<col style="width: 4%" />
<col style="width: 4%" />
<col style="width: 5%" />
<col style="width: 6%" />
<col style="width: 2%" />
<col style="width: 73%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">OTUid</th>
<th style="text-align: right;">FDR</th>
<th style="text-align: right;">P.Value</th>
<th style="text-align: right;">logFC</th>
<th style="text-align: left;">Factor</th>
<th style="text-align: left;">Group</th>
<th style="text-align: left;">Taxonomy</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">OTU1050</td>
<td style="text-align: right;">0.6525563</td>
<td style="text-align: right;">0.0465581</td>
<td style="text-align: right;">-19.333333</td>
<td style="text-align: left;">OTU1050_BL</td>
<td style="text-align: left;">BL</td>
<td style="text-align: left;">Firmicutes|Clostridia|Oscillospirales|Oscillospiraceae|UCG-005|uncultured_organism_g__UCG-005</td>
</tr>
<tr class="even">
<td style="text-align: left;">OTU11</td>
<td style="text-align: right;">0.6525563</td>
<td style="text-align: right;">0.0338668</td>
<td style="text-align: right;">-21.523809</td>
<td style="text-align: left;">OTU11_BL</td>
<td style="text-align: left;">BL</td>
<td style="text-align: left;">Firmicutes|Bacilli|Erysipelotrichales|Erysipelatoclostridiaceae|Erysipelotrichaceae_UCG-003|unclassified_g__Erysipelotrichaceae_UCG-003</td>
</tr>
<tr class="odd">
<td style="text-align: left;">OTU1303</td>
<td style="text-align: right;">0.6525563</td>
<td style="text-align: right;">0.0410576</td>
<td style="text-align: right;">-35.619048</td>
<td style="text-align: left;">OTU1303_BL</td>
<td style="text-align: left;">BL</td>
<td style="text-align: left;">Firmicutes|Clostridia|Oscillospirales|Oscillospiraceae|UCG-005|unclassified_g__UCG-005</td>
</tr>
<tr class="even">
<td style="text-align: left;">OTU1363</td>
<td style="text-align: right;">0.6525563</td>
<td style="text-align: right;">0.0459235</td>
<td style="text-align: right;">-4.000000</td>
<td style="text-align: left;">OTU1363_BL</td>
<td style="text-align: left;">BL</td>
<td style="text-align: left;">Firmicutes|Clostridia|Lachnospirales|Lachnospiraceae|Moryella|human_gut_metagenome_g__Moryella</td>
</tr>
<tr class="odd">
<td style="text-align: left;">OTU1502</td>
<td style="text-align: right;">0.6525563</td>
<td style="text-align: right;">0.0461802</td>
<td style="text-align: right;">56.666667</td>
<td style="text-align: left;">OTU1502_Treat</td>
<td style="text-align: left;">Treat</td>
<td style="text-align: left;">Firmicutes|Clostridia|Lachnospirales|Lachnospiraceae|Lachnospiraceae_NK4A136_group|uncultured_organism_g__Lachnospiraceae_NK4A136_group</td>
</tr>
<tr class="even">
<td style="text-align: left;">OTU1571</td>
<td style="text-align: right;">0.6525563</td>
<td style="text-align: right;">0.0256103</td>
<td style="text-align: right;">-33.642857</td>
<td style="text-align: left;">OTU1571_BL</td>
<td style="text-align: left;">BL</td>
<td style="text-align: left;">Firmicutes|Clostridia|Oscillospirales|Ruminococcaceae|UBA1819|uncultured_organism_g__UBA1819</td>
</tr>
<tr class="odd">
<td style="text-align: left;">OTU4886</td>
<td style="text-align: right;">0.6525563</td>
<td style="text-align: right;">0.0286119</td>
<td style="text-align: right;">38.452381</td>
<td style="text-align: left;">OTU4886_Treat</td>
<td style="text-align: left;">Treat</td>
<td style="text-align: left;">Bacteroidetes|Bacteroidia|Bacteroidales|Bacteroidaceae|Bacteroides|unclassified_g__Bacteroides</td>
</tr>
<tr class="even">
<td style="text-align: left;">OTU4961</td>
<td style="text-align: right;">0.6525563</td>
<td style="text-align: right;">0.0498959</td>
<td style="text-align: right;">-5.214286</td>
<td style="text-align: left;">OTU4961_BL</td>
<td style="text-align: left;">BL</td>
<td style="text-align: left;">Firmicutes|Clostridia|Lachnospirales|Lachnospiraceae|unclassified_f__Lachnospiraceae|unclassified_f__Lachnospiraceae</td>
</tr>
<tr class="odd">
<td style="text-align: left;">OTU4981</td>
<td style="text-align: right;">0.6525563</td>
<td style="text-align: right;">0.0278053</td>
<td style="text-align: right;">-7.547619</td>
<td style="text-align: left;">OTU4981_Treat</td>
<td style="text-align: left;">Treat</td>
<td style="text-align: left;">Firmicutes|Clostridia|Lachnospirales|Lachnospiraceae|unclassified_f__Lachnospiraceae|unclassified_f__Lachnospiraceae</td>
</tr>
<tr class="even">
<td style="text-align: left;">OTU5109</td>
<td style="text-align: right;">0.6525563</td>
<td style="text-align: right;">0.0328646</td>
<td style="text-align: right;">116.119048</td>
<td style="text-align: left;">OTU5109_Treat</td>
<td style="text-align: left;">Treat</td>
<td style="text-align: left;">Bacteroidetes|Bacteroidia|Bacteroidales|Bacteroidaceae|Bacteroides|unclassified_g__Bacteroides</td>
</tr>
<tr class="odd">
<td style="text-align: left;">OTU5205</td>
<td style="text-align: right;">0.6525563</td>
<td style="text-align: right;">0.0278150</td>
<td style="text-align: right;">514.928571</td>
<td style="text-align: left;">OTU5205_Treat</td>
<td style="text-align: left;">Treat</td>
<td style="text-align: left;">Bacteroidetes|Bacteroidia|Bacteroidales|Bacteroidaceae|Bacteroides|uncultured_organism_g__Bacteroides</td>
</tr>
<tr class="even">
<td style="text-align: left;">OTU5272</td>
<td style="text-align: right;">0.6525563</td>
<td style="text-align: right;">0.0384252</td>
<td style="text-align: right;">5.119048</td>
<td style="text-align: left;">OTU5272_BL</td>
<td style="text-align: left;">BL</td>
<td style="text-align: left;">Proteobacteria|Gammaproteobacteria|Enterobacterales|Enterobacteriaceae|Escherichia-Shigella|unclassified_g__Escherichia-Shigella</td>
</tr>
<tr class="odd">
<td style="text-align: left;">OTU5828</td>
<td style="text-align: right;">0.6525563</td>
<td style="text-align: right;">0.0494315</td>
<td style="text-align: right;">-52.166667</td>
<td style="text-align: left;">OTU5828_BL</td>
<td style="text-align: left;">BL</td>
<td style="text-align: left;">Firmicutes|Clostridia|Oscillospirales|Oscillospiraceae|UCG-003|uncultured_organism_g__UCG-003</td>
</tr>
<tr class="even">
<td style="text-align: left;">OTU6004</td>
<td style="text-align: right;">0.6525563</td>
<td style="text-align: right;">0.0408601</td>
<td style="text-align: right;">-27.952381</td>
<td style="text-align: left;">OTU6004_Treat</td>
<td style="text-align: left;">Treat</td>
<td style="text-align: left;">Firmicutes|Clostridia|Oscillospirales|Ruminococcaceae|norank_f__Ruminococcaceae|unclassified_g__norank_f__Ruminococcaceae</td>
</tr>
<tr class="odd">
<td style="text-align: left;">OTU6239</td>
<td style="text-align: right;">0.6525563</td>
<td style="text-align: right;">0.0475586</td>
<td style="text-align: right;">-18.238095</td>
<td style="text-align: left;">OTU6239_Treat</td>
<td style="text-align: left;">Treat</td>
<td style="text-align: left;">Firmicutes|Negativicutes|Veillonellales-Selenomonadales|Veillonellaceae|Allisonella|uncultured_bacterium_g__Allisonella</td>
</tr>
<tr class="even">
<td style="text-align: left;">OTU6472</td>
<td style="text-align: right;">0.6525563</td>
<td style="text-align: right;">0.0377707</td>
<td style="text-align: right;">-40.428571</td>
<td style="text-align: left;">OTU6472_Treat</td>
<td style="text-align: left;">Treat</td>
<td style="text-align: left;">Firmicutes|Clostridia|Peptostreptococcales-Tissierellales|Peptostreptococcaceae|Intestinibacter|uncultured_bacterium_g__Intestinibacter</td>
</tr>
<tr class="odd">
<td style="text-align: left;">OTU6659</td>
<td style="text-align: right;">0.6525563</td>
<td style="text-align: right;">0.0367304</td>
<td style="text-align: right;">-9.000000</td>
<td style="text-align: left;">OTU6659_BL</td>
<td style="text-align: left;">BL</td>
<td style="text-align: left;">Firmicutes|Negativicutes|Acidaminococcales|Acidaminococcaceae|Acidaminococcus|unclassified_g__Acidaminococcus</td>
</tr>
</tbody>
</table>

    final_res<-dplyr::select(ano_data,c(patientID,km.res.cluster,limma_res$Factor))
    final_mat<-data.frame(patientID = final_res$patientID,final_res[,-c(1,2)])
    final_col<-data.frame(patientID =ano_data$patientID,Class=ano_data$km.res.cluster)
    final_mat<-merge(final_col,final_mat,by="patientID")
    limma_pca <- PCA(final_mat[,-c(1,2)], graph = FALSE)

3.4 PCA plot based on the representative OTUs
---------------------------------------------

    fviz_pca_ind(limma_pca,legend.title = "k-means cluster",
                 label = "none", # hide individual labels
                 col.ind = final_mat$Class, # color by groups
                 palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                 addEllipses = TRUE # Concentration ellipses
    )

<img src="AutoDecoder_files/figure-markdown_strict/unnamed-chunk-14-1.png" width="30%" style="display: block; margin: auto;" />

4 Reclassify the patients with incompletion data by the representative OTUs
===========================================================================

    cli<-fread("../Data/Data/final_clinical_40pt.csv")
    meta<-fread("../Data/Data/meta.csv",data.table = F)%>%subset(.,Site=="Stool")
    meta<-merge(dplyr::select(meta,c("patientID","Samples","Cycle")),dplyr::select(cli,c("patientID","PFS","PFStime")),by="patientID")
    colnames(meta)[2]="Samples"
    meta$Group<-ifelse(meta$Cycle=="BL","BL","Treat")
    meta_bl<-subset(meta,Group=="BL")
    meta_treat<-subset(meta,Group=="Treat")
    otu<-regaMicrobiome$StoolMicrobiome$TaxonomyReads$OTU
    colnames(otu)[1]="OTUid"
    limma_res_BL<-subset(limma_res,Group=="BL")
    limma_res_treat<-subset(limma_res,Group=="Treat")
    otu<-subset(otu,OTUid%in%limma_res$OTUid)

    otu_bl<-otu[,which(colnames(otu)%in%c("OTUid",meta_bl$Samples))]
    otu_treat<-otu[,which(colnames(otu)%in%c("OTUid",meta_treat$Samples))]
    otu_bl_limma<-subset(otu_bl,OTUid%in%limma_res_BL$OTUid)
    otu_treat_limma<-subset(otu_treat,OTUid%in%limma_res$OTUid)
    otu_bl_limma<-data.frame(row.names = otu_bl_limma$OTUid,otu_bl_limma[,-1])
    otu_bl_limma<-data.frame(Samples=colnames(otu_bl_limma),t(otu_bl_limma))

    otu_treat_limma<-data.frame(row.names = otu_treat_limma$OTUid,otu_treat_limma[,-1])
    otu_treat_limma<-data.frame(Samples=colnames(otu_treat_limma),t(otu_treat_limma))


    df_bl<-merge(meta_bl,otu_bl_limma,by="Samples")
    df_treat<-merge(meta_treat,otu_treat_limma,by="Samples")

    otu_bl_limma_mat<-data.frame(row.names = df_bl$patientID,df_bl[,-c(1:6)])
    otu_treat_limma_mat<-data.frame(row.names = paste0(df_treat$patientID,"_",df_treat$Cycle),df_treat[,-c(1:6)])

4.1 All the patients with baseline samples
------------------------------------------

### 4.1.1 Choose k value

    fviz_nbclust(
      otu_bl_limma_mat, 
      kmeans, 
      k.max = 10,
      method = "wss",
      verbose = FALSE)

<img src="AutoDecoder_files/figure-markdown_strict/unnamed-chunk-16-1.png" width="30%" style="display: block; margin: auto;" />

### 4.1.2 Plot K-means clustering

    km.res <- kmeans(scale(otu_bl_limma_mat),3, nstart = 25)
    fviz_cluster(km.res, data = otu_bl_limma_mat,repel = T,labelsize = 12,
                 # palette = c("#00AFBB", "#E7B800", "black","red","blue","gray"),
                 ggtheme = theme_minimal(),
                 main = "Kmeans Clustering Plot")

<img src="AutoDecoder_files/figure-markdown_strict/unnamed-chunk-17-1.png" width="30%" style="display: block; margin: auto;" />

4.2 All the patients with treated samples
-----------------------------------------

    cluster_BL<-data.frame(km.res$cluster)
    cluster_BL$patientID<-rownames(cluster_BL)
    cluster_BL$km.res.cluster<-paste0("km",cluster_BL$km.res.cluster)
    stat_num<-cluster_BL%>%group_by(km.res.cluster)%>%summarise(Num=n())
    cluster_BL$km.res.cluster<-ifelse(cluster_BL$km.res.cluster==stat_num$km.res.cluster[which(stat_num$Num>5)],"km1","km2")

### 4.2.1 Choose k value

    fviz_nbclust(
      otu_treat_limma_mat, 
      kmeans, 
      k.max = 10,
      method = "wss",
      verbose = FALSE)

<img src="AutoDecoder_files/figure-markdown_strict/unnamed-chunk-19-1.png" width="30%" style="display: block; margin: auto;" />

### 4.2.2 Plot K-means clustering

    set.seed(123)
    km.res <- kmeans(scale(otu_treat_limma_mat),4, nstart = 25)
    fviz_cluster(km.res, data = otu_treat_limma_mat,repel = T,labelsize = 12,
                  palette = c("#00AFBB", "#E7B800", "black","red","blue","gray"),
                 ggtheme = theme_minimal(),
                 main = "Kmeans Clustering Plot")

![](AutoDecoder_files/figure-markdown_strict/unnamed-chunk-20-1.png)

4.3 PFS survival of the Patients with incompletion data
-------------------------------------------------------

    cluster_treat<-data.frame(km.res$cluster)
    cluster_treat$patientID<-rownames(cluster_treat)
    cluster_treat$km.res.cluster<-paste0("km",cluster_treat$km.res.cluster)
    stat_num<-cluster_treat%>%group_by(km.res.cluster)%>%summarise(Num=n())
    knitr::kable(stat_num)
    cluster_treat$km.res.cluster<-gsub("km4","km2",cluster_treat$km.res.cluster)
    df_bl<-merge(cluster_BL,df_bl,by="patientID")
    df_treat$patientID<-paste0(df_treat$patientID,"_",df_treat$Cycle)
    df_treat<-merge(cluster_treat,df_treat,by="patientID")


    fit1<-survfit(Surv(PFStime,PFS) ~ km.res.cluster,
                  data = df_bl)
    fit1
    fit2<-survfit(Surv(PFStime,PFS) ~ km.res.cluster,
                  data = df_treat)
    fit2

<table>
<thead>
<tr class="header">
<th style="text-align: left;">km.res.cluster</th>
<th style="text-align: right;">Num</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">km1</td>
<td style="text-align: right;">6</td>
</tr>
<tr class="even">
<td style="text-align: left;">km2</td>
<td style="text-align: right;">4</td>
</tr>
<tr class="odd">
<td style="text-align: left;">km3</td>
<td style="text-align: right;">20</td>
</tr>
<tr class="even">
<td style="text-align: left;">km4</td>
<td style="text-align: right;">1</td>
</tr>
</tbody>
</table>

    ## Call: survfit(formula = Surv(PFStime, PFS) ~ km.res.cluster, data = df_bl)
    ## 
    ##                     n events median 0.95LCL 0.95UCL
    ## km.res.cluster=km1 31     27   3.80    2.03     5.6
    ## km.res.cluster=km2  7      7   1.67    1.23      NA
    ## Call: survfit(formula = Surv(PFStime, PFS) ~ km.res.cluster, data = df_treat)
    ## 
    ##                     n events median 0.95LCL 0.95UCL
    ## km.res.cluster=km1  6      6  11.40    2.30      NA
    ## km.res.cluster=km2  5      5   2.20    2.03      NA
    ## km.res.cluster=km3 20     16   6.33    4.27      NA

    p1<-ggsurvplot(fit1, data=df_bl,pval.method = T,combine = T, 
                   palette  = c("black", "#00AFBB","#E7B800","red","green"),
                   tables.theme = theme_bw(base_size = 8),
                   risk.table = T,
                   pval = TRUE,
                   ggtheme = theme_survminer(),
                   legend.title="k-means",
                   risk.table.col = "strata",
                   surv.median.line = "hv",
                   risk.table.y.text.col = T,
                   risk.table.y.text = F )
    p1

<img src="AutoDecoder_files/figure-markdown_strict/unnamed-chunk-22-1.png" width="50%" style="display: block; margin: auto;" />

    p2<-ggsurvplot(fit2, data=df_treat,pval.method = T,combine = T, 
                   palette  = c("black","#00AFBB","#E7B800","red"),
                   tables.theme = theme_bw(base_size = 8),
                   risk.table = T,
                   pval = TRUE,
                   ggtheme = theme_survminer(),
                   legend.title="k-means",
                   risk.table.col = "strata",
                   surv.median.line = "hv",
                   risk.table.y.text.col = T,
                   risk.table.y.text = F )
    p2

<img src="AutoDecoder_files/figure-markdown_strict/unnamed-chunk-23-1.png" width="50%" style="display: block; margin: auto;" />
