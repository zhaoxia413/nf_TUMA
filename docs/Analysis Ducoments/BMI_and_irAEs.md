-   [1 BMI related features](#bmi-related-features)
    -   [1.1 BMI and FBratio of gut
        microbiome](#bmi-and-fbratio-of-gut-microbiome)
        -   [1.1.1 Histgrams for BMI and
            FBartio](#histgrams-for-bmi-and-fbartio)
        -   [1.1.2 PFS survival plot for
            BMI](#pfs-survival-plot-for-bmi)
        -   [1.1.3 Correlation of BMI and gut
            FBratio](#correlation-of-bmi-and-gut-fbratio)
    -   [1.2 Correlation of BMI and stroma CD3
        expression](#correlation-of-bmi-and-stroma-cd3-expression)
    -   [1.3 BMI ralated KEGG pathways](#bmi-ralated-kegg-pathways)
    -   [1.4 Validate the clinical value of gut FBratio with public
        datasets](#validate-the-clinical-value-of-gut-fbratio-with-public-datasets)
        -   [1.4.1 Correlation Firmucutes and
            Bacteroidetes](#correlation-firmucutes-and-bacteroidetes)
        -   [1.4.2 Diagnostic value of Gut FBratio in CRC
            cohorts](#diagnostic-value-of-gut-fbratio-in-crc-cohorts)
-   [2 Treated-related aderse events](#treated-related-aderse-events)
    -   [2.1 PFS survival curves for each
        events](#pfs-survival-curves-for-each-events)
    -   [2.2 Abundance of Diarrhea-related
        *Desulfovibrionaceae*](#abundance-of-diarrhea-related-desulfovibrionaceae)

[`Return`](./)

1 BMI related features
======================

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
    library(ComplexHeatmap)
    library(scico)
    library(colorspace)
    library(RColorBrewer)
    library(lubridate)
    library(tableone)
    library(kableExtra)
    library(ROCR)
    library(caret)
    library(ggrepel)
    library(randomForest)
    library(lavaan)
    library(Hmisc)
    library(corrplot)
    library(colortools)
    library(visibly)
    library(plotly)
    source("../R_function/colors.R")
    source("../R_function/surv_plot.R")
    source("../R_function/geom_flat_violin.R")
    source("../R_function/summarySE.R")
    theme_set(theme_cowplot())
    "%ni%" <- Negate("%in%")
    options(stringsAsFactors = F)

</details>

1.1 BMI and FBratio of gut microbiome
-------------------------------------

### 1.1.1 Histgrams for BMI and FBartio

    df<-fread("../Data/Data/Phylum_cli_111samples.csv",data.table = F)
    df$FBratio<-df$Firmicutes/df$Bacteroidetes
    df$FBratio_g<-ifelse(df$FBratio>=median(df$FBratio),"High","Low")
    data<-subset(df,Site=="Stool"&Response!="NE"&Cycle=="BL")
    par(mfrow=c(1,2))
    hist(data$BMI,main="Frequence of BMI",xlab = "BMI")
    hist(data$FBratio,main="Frequence of FBratio",xlab = "FBratio")

![](BMI_and_irAEs_files/figure-markdown_strict/unnamed-chunk-2-1.png)

    data$BMI_g<-ifelse(data$BMI>25,"High","Low")
    data$FBratio_g<-ifelse(data$FBratio>median(data$FBratio),"High","Low")
    fit<-survfit(Surv(PFStime,PFS) ~ BMI_g,
                       data = data)
    fit

    ## Call: survfit(formula = Surv(PFStime, PFS) ~ BMI_g, data = data)
    ## 
    ##             n events median 0.95LCL 0.95UCL
    ## BMI_g=High 10      8   3.28    2.20      NA
    ## BMI_g=Low  22     21   1.97    1.87     4.2

### 1.1.2 PFS survival plot for BMI

    ggsurvplot(fit, data=data,xlab = "Time(months)",
               censor.size=0.5, size = 0.5,
               tables.theme = theme_few(base_size = 5),
               legend.labs = c(">25", "<=25"),
                    legend.title = "BMI",palette = c("red","black"),
                    risk.table = T,
                    pval = TRUE,pval.size = 3, 
                    pval.coord=c(0.8,0.2),pval.method=F,
                    pval.method.coord=c(0.05,0.3), 
                    ggtheme = theme_minimal() + 
                      theme(line = element_line(size = 0.1),
                            text  = element_text(size = 6)),
                    risk.table.col = "strata",
                    surv.median.line = "hv",
                    risk.table.y.text.col = T,
                    risk.table.y.text = FALSE )

<img src="BMI_and_irAEs_files/figure-markdown_strict/unnamed-chunk-3-1.png" width="40%" style="display: block; margin: auto;" />

### 1.1.3 Correlation of BMI and gut FBratio

    p1<-ggscatter(subset(df,Cycle=="BL"&Response!="NE"&FBratio<10), x = "FBratio", y = "BMI",size=0.5,mean.point = T,
              color = "Site", add.params = list(c(size=0.5,color="Site")),
              add = "reg.line", conf.int = TRUE)+
      stat_cor(label.x = 0.3,aes(color=Site))+
      theme_few(base_size = 8)+
      scale_color_aaas()

    p2<-ggscatter(subset(df,Cycle=="BL"&Response!="NE"), x = "Firmicutes", y = "Bacteroidetes",size=0.5,mean.point = T,
              color = "Site", add.params = list(c(size=0.5,color="Site")),
              add = "reg.line", conf.int = TRUE)+
      stat_cor(label.x = 0.2,aes(color=Site))+
      theme_few(base_size = 8)+
      scale_color_aaas()

    p3<-ggstatsplot::ggbarstats(data = data,x=Response,ggtheme = ggplot2::theme_bw(base_size=8),bias.correct = T,
                            y =FBratio_g,subtitle = F,results.subtitle=T,
                            ggstatsplot.layer = FALSE,
                            legend.position="right",
                            messages = FALSE,
                            package = "ggsci",
                            palette = "default_nejm",
                            main = Response, nboot = 100,
                            legend.title = "Response")

    ## Registered S3 methods overwritten by 'lme4':
    ##   method                          from
    ##   cooks.distance.influence.merMod car 
    ##   influence.merMod                car 
    ##   dfbeta.influence.merMod         car 
    ##   dfbetas.influence.merMod        car

    plot_grid(p1,p2, p3,labels = c("A","B","C"), ncol =3, nrow = 1)

    ## `geom_smooth()` using formula 'y ~ x'

    ## `geom_smooth()` using formula 'y ~ x'

<img src="BMI_and_irAEs_files/figure-markdown_strict/unnamed-chunk-4-1.png" width="90%" style="display: block; margin: auto;" />

1.2 Correlation of BMI and stroma CD3 expression
------------------------------------------------

    bmi<-fread("../Data/Data/BMI_CD3.csv",data.table = F)
    ggscatter(bmi,x = "value", y = "BMI",size=1,mean.point = T,
              add.params = list(c(size=0.5,color="variable")),color="variable",
              add = "reg.line", conf.int = TRUE)+
      stat_cor(aes(color=variable),label.x = 3)+
      facet_wrap(~variable,scales = "free")+
      theme_few(base_size = 8)+
      scale_color_aaas()

<img src="BMI_and_irAEs_files/figure-markdown_strict/unnamed-chunk-5-1.png" width="70%" style="display: block; margin: auto;" />

1.3 BMI ralated KEGG pathways
-----------------------------

    kegg_mat<-read.csv("../Data/Data/BMI_Gut_kegg.csv",header = T,row.names = 1)%>%as.matrix()
    kegg_info<-fread("../Data/Data/Gut_kegg_description.csv",data.table = F)
    CorMatrix <- function(cor,p){
      ut <- upper.tri(cor) 
      data.frame(row = rownames(cor)[row(cor)[ut]],
                 column = rownames(cor)[col(cor)[ut]], 
                 cor =(cor)[ut],
                 p = p[ut])}

    cor_res  <- rcorr(as.matrix(kegg_mat),type ="spearman")
    corMat<-CorMatrix(cor_res$r, cor_res$P)%>%filter(p<=0.05)%>%
      filter(row%in%c("BMI"))
    colnames(corMat)[2]="Description"
    corMat<-merge(kegg_info,corMat,by="Description")
    ggplot(corMat,aes(cor,reorder(Description,cor),color=kegg_group))+
      geom_point(aes(size=abs(cor)))+
      scale_color_manual(values = col11)+
      theme_bw()+
      xlab("corValue (spearman)")+
      ylab("KEGG")

<img src="BMI_and_irAEs_files/figure-markdown_strict/unnamed-chunk-6-1.png" width="80%" style="display: block; margin: auto;" />

1.4 Validate the clinical value of gut FBratio with public datasets
-------------------------------------------------------------------

### 1.4.1 Correlation Firmucutes and Bacteroidetes

    df1<-fread("../Data/Data/publicData/PRJNA541981_phylum.csv",data.table = F)
    df2<-fread("../Data/Data/publicData/PRJEB22863_phylum.csv",data.table = F)
    df3<-fread("../Data/Data/publicData/PRJNA399742_phylum.csv",data.table = F)
    df4<-fread("../Data/Data/publicData/gFBratio_subsetCRC_810samples.csv",data.table = F)[,-c(2:4)]
    colnames(df4)[c(2,3)]=c("Cancer","datasets")
    colnames(df3)[1]="samples"
    colnames(df4)[1]="samples"
    index<-colnames(df4)[-9]
    data<-bind_rows(select(df1,index),select(df2,index),select(df3,index),select(df4,index))
    p1<-ggscatter(data, x = "Firmicutes", 
              y = "Bacteroidetes",size=1,alpha=0.4,color="datasets",cor.method = "spearman",mean.point = T,
              palette = "jco",add.params = list(alpha=0.5,size=0.5), ggtheme = theme_few(base_size = 7),
              add = "reg.line", conf.int = TRUE)+
      stat_cor(aes(color=datasets),label.x = 0.55,size=2)+
      theme(legend.text = element_text(size=7))

    p2<-ggscatter(data, x = "Firmicutes", mean.point = T,cor.method = "spearman",
              ggtheme = theme_few(base_size = 7),
              y = "Bacteroidetes",size=1,alpha=0.4,color="Cancer",
              palette = "jco",add.params = list(alpha=0.5,size=0.5),
              add = "reg.line", conf.int = TRUE)+
      stat_cor(aes(color=Cancer),label.x = 0.55,size=2)+
      theme(legend.text = element_text(size=7))

    plot_grid(p1,p2,labels = c("A","B"), ncol =2, nrow = 1)

<img src="BMI_and_irAEs_files/figure-markdown_strict/unnamed-chunk-7-1.png" width="80%" style="display: block; margin: auto;" />

### 1.4.2 Diagnostic value of Gut FBratio in CRC cohorts

#### 1.4.2.1 Datasets infomation

    data<-fread("../Data/Data/publicData/gFBratio_subsetCRC_810samples.csv",data.table = F)
    data%>%group_by(dataset_name,Disease)%>%dplyr::summarise(Samples_num=n())%>%
      knitr::kable(caption = "Datasets and Disease") 

    ## `summarise()` has grouped output by 'dataset_name'. You can override using the `.groups` argument.

<table>
<caption>
Datasets and Disease
</caption>
<thead>
<tr>
<th style="text-align:left;">
dataset\_name
</th>
<th style="text-align:left;">
Disease
</th>
<th style="text-align:right;">
Samples\_num
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
FengQ\_2015
</td>
<td style="text-align:left;">
Adenoma
</td>
<td style="text-align:right;">
47
</td>
</tr>
<tr>
<td style="text-align:left;">
FengQ\_2015
</td>
<td style="text-align:left;">
CRC
</td>
<td style="text-align:right;">
46
</td>
</tr>
<tr>
<td style="text-align:left;">
FengQ\_2015
</td>
<td style="text-align:left;">
Normal
</td>
<td style="text-align:right;">
61
</td>
</tr>
<tr>
<td style="text-align:left;">
PRJEB6070
</td>
<td style="text-align:left;">
Adenoma
</td>
<td style="text-align:right;">
38
</td>
</tr>
<tr>
<td style="text-align:left;">
PRJEB6070
</td>
<td style="text-align:left;">
CRC
</td>
<td style="text-align:right;">
41
</td>
</tr>
<tr>
<td style="text-align:left;">
PRJEB6070
</td>
<td style="text-align:left;">
Normal
</td>
<td style="text-align:right;">
50
</td>
</tr>
<tr>
<td style="text-align:left;">
PRJNA290926
</td>
<td style="text-align:left;">
Adenoma
</td>
<td style="text-align:right;">
68
</td>
</tr>
<tr>
<td style="text-align:left;">
PRJNA290926
</td>
<td style="text-align:left;">
CRC
</td>
<td style="text-align:right;">
90
</td>
</tr>
<tr>
<td style="text-align:left;">
PRJNA290926
</td>
<td style="text-align:left;">
Normal
</td>
<td style="text-align:right;">
92
</td>
</tr>
<tr>
<td style="text-align:left;">
ThomasAM\_2018a
</td>
<td style="text-align:left;">
Adenoma
</td>
<td style="text-align:right;">
26
</td>
</tr>
<tr>
<td style="text-align:left;">
ThomasAM\_2018a
</td>
<td style="text-align:left;">
CRC
</td>
<td style="text-align:right;">
28
</td>
</tr>
<tr>
<td style="text-align:left;">
ThomasAM\_2018a
</td>
<td style="text-align:left;">
Normal
</td>
<td style="text-align:right;">
24
</td>
</tr>
<tr>
<td style="text-align:left;">
ZellerG\_2014
</td>
<td style="text-align:left;">
Adenoma
</td>
<td style="text-align:right;">
42
</td>
</tr>
<tr>
<td style="text-align:left;">
ZellerG\_2014
</td>
<td style="text-align:left;">
CRC
</td>
<td style="text-align:right;">
91
</td>
</tr>
<tr>
<td style="text-align:left;">
ZellerG\_2014
</td>
<td style="text-align:left;">
Normal
</td>
<td style="text-align:right;">
66
</td>
</tr>
</tbody>
</table>

#### 1.4.2.2 Gut Phylum composition

    data$Disease<-factor(data$Disease,levels = c("Normal","Adenoma","CRC"))
    comp<-list(c("Normal","Adenoma"),
               c("Normal","CRC"),
               c("Adenoma","CRC"))
    ggplot(data,aes(Disease,log10(FBratio),fill=Disease))+
      geom_jitter(size=0.5,alpha=0.5)+
      geom_boxplot(color="black",outlier.size = 0.2,outlier.color = "gray")+
      theme_few()+
      stat_compare_means(comparisons = comp,label = "p.signif")+
      theme(axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,size = 6))+
      scale_fill_uchicago()+
      facet_grid(~dataset_name,scales = "free",space = "free_x")+
      theme(strip.text = element_text(size = 3))

![](BMI_and_irAEs_files/figure-markdown_strict/unnamed-chunk-9-1.png)

    compdata<-melt(data[,-c(1:4,6:7)],
                   id.vars = "Disease",
                   variable.name = "Phylum",
                   value.name = "Abundance")%>%
      group_by(Disease,Phylum)%>%
      dplyr::summarise_each(mean)

    compdata$Abundance<-round(compdata$Abundance,3)
    compdata$Phylum<-factor(compdata$Phylum,levels = c("Firmicutes","Bacteroidetes","Actinobacteria","Proteobacteria","Other"))
    df<-split.data.frame(compdata,f = compdata$Disease,drop = T)
    levels(factor(compdata$Disease))
    fb<-dcast(compdata,Disease~Phylum)

    ## Using 'Abundance' as value column. Use 'value.var' to override

    fb$FBratio=fb$Firmicutes/fb$Bacteroidetes
    p<-list()

    for (i in seq_along(df)) {
      df[[i]]$fraction<-df[[i]]$Abundance/sum(df[[i]]$Abundance)
      df[[i]]$ymax<-cumsum(df[[i]]$fraction)
      df[[i]]$ymin<-c(0,head(df[[i]]$ymax,n=-1))
      df[[i]]$labelPosition<-(df[[i]]$ymax + df[[i]]$ymin)/2
      df[[i]]$label<-paste0(df[[i]]$Abundance)
      p[[i]]<-ggplot(df[[i]],aes(ymax=ymax,ymin=ymin,
                                 xmax=4,xmin=3))+
        geom_rect(aes(fill=Phylum))+
        geom_text_repel(x=3.5,box.padding = 0.5, max.overlaps = Inf,
                        aes(y=labelPosition,label=label),size=3)+
        scale_fill_d3()+
        coord_polar(theta = "y")+
        xlim(2,4)+
        #facet_grid(~Disease)+
        theme_void()+
        ggtitle(paste0(names(df)[i]))+
        theme(legend.position = "none",
            legend.key.size=unit(.1,"inches"),
            legend.text.align=0,
            legend.title=element_text(colour="black",size=8,face = "bold"),
            legend.direction ="vertical",
              legend.spacing = unit(0.1,"cm"),
            legend.spacing.y = unit(0.1,"cm"),
            legend.spacing.x =unit(0.1,"cm"),
            legend.box.spacing = unit(0.1,"cm"),
            legend.justification=c(.4,.4),
              plot.title = element_text(vjust = -65,hjust = +0.5))+
        ggtitle(paste0(names(df)[i]))
      }

    ## [1] "Normal"  "Adenoma" "CRC"

    plot_grid(p[[1]],p[[2]],p[[3]]+theme(legend.position = "right"),
              labels = c("A","B","C"), ncol =2, nrow = 2)

<img src="BMI_and_irAEs_files/figure-markdown_strict/unnamed-chunk-10-1.png" width="80%" style="display: block; margin: auto;" />

#### 1.4.2.3 RF classifiers for each dataset

    data<-fread("../Data/Data/publicData/gFBratio_subsetCRC_810samples.csv",data.table = F)
    dflist<-split.data.frame(data,f=data$dataset_name,drop = F)
    data<-list();data1<-list();data2<-list();data3<-list()
    trainData<-list();testData<-list();trainData<-list();testData<-list()
    re_rf<-list();pred_rf<-list();pred<-list();perf<-list();auc<-list()
    trainData1<-list();testData1<-list();trainData1<-list();testData1<-list()
    re_rf1<-list();pred_rf1<-list();pred1<-list();perf1<-list();auc1<-list()
    trainData2<-list();testData2<-list();trainData2<-list();testData2<-list()
    re_rf2<-list();pred_rf2<-list();pred2<-list();perf2<-list();auc2<-list()
    trainData3<-list();testData3<-list();trainData3<-list();testData3<-list()
    re_rf3<-list();pred_rf3<-list();pred3<-list();perf3<-list();auc3<-list()
    NC_res<-list();NA_res<-list();CA_res<-list();NCA_res<-list()
    predictor<-c("Normal vs CRC","Normal vs Adenoma",
                 "CRC vs Adenoma","Normal vs CRC/Adenoma")
    for (i in seq_along(dflist)) {
      data[[i]]<-subset(dflist[[i]],Disease%in%c("Normal","CRC"))
      data1[[i]]<-subset(dflist[[i]],Disease%in%c("Normal","Adenoma"))
      data2[[i]]<-subset(dflist[[i]],Disease%in%c("CRC","Adenoma"))
      data3[[i]]<-subset(dflist[[i]],Disease%in%c("Normal","CRC","Adenoma"))
      
      data[[i]]$Disease<-ifelse(data[[i]]$Disease=="Normal","Normal","CRC")
      data1[[i]]$Disease<-ifelse(data1[[i]]$Disease=="Normal","Normal","Adenoma")
      data2[[i]]$Disease<-ifelse(data2[[i]]$Disease=="CRC","CRC","Adenoma")
      data3[[i]]$Disease<-ifelse(data3[[i]]$Disease=="Normal","Normal","Other")
      set.seed(1000)
      trainIndex<-sample(nrow(data[[i]]),nrow(data[[i]])*0.8)
      trainData[[i]]<-data[[i]][trainIndex,]
      testData[[i]]<-data[[i]][-trainIndex,]
      trainData[[i]]$Disease = as.factor(trainData[[i]]$Disease)
      testData[[i]]$Disease = as.factor(testData[[i]]$Disease)
      re_rf[[i]] = randomForest(Disease~FBratio+Age,
                                data = trainData[[i]],ntree=500)
      pred_rf[[i]]=predict(re_rf[[i]],newdata=testData[[i]],type="prob") 
      pred[[i]]<-prediction(pred_rf[[i]][,2],testData[[i]]$Disease)
      perf[[i]]<-performance(pred[[i]],"tpr","fpr")
      auc[[i]] <- performance(pred[[i]],'auc')
      auc[[i]] = unlist(slot(auc[[i]],"y.values"))
      auc[[i]]=round(auc[[i]],digits = 3)
      
      set.seed(1234)
      trainIndex1<-sample(nrow(data1[[i]]),nrow(data1[[i]])*0.8)
      trainData1[[i]]<-data1[[i]][trainIndex1,]
      testData1[[i]]<-data1[[i]][-trainIndex1,]
      trainData1[[i]]$Disease = as.factor(trainData1[[i]]$Disease)
      testData1[[i]]$Disease = as.factor(testData1[[i]]$Disease)
      re_rf1[[i]] = randomForest(Disease~FBratio+Age,
                                 data = trainData1[[i]],ntree=500)
      pred_rf1[[i]]=predict(re_rf1[[i]],newdata=testData1[[i]],type="prob") 
      pred1[[i]]<-prediction(pred_rf1[[i]][,2],testData1[[i]]$Disease)
      perf1[[i]]<-performance(pred1[[i]],"tpr","fpr")
      auc1[[i]] <- performance(pred1[[i]],'auc')
      auc1[[i]] = unlist(slot(auc1[[i]],"y.values"))
      auc1[[i]]=round(auc1[[i]],digits = 3)
      
      set.seed(123467)
      trainIndex2<-sample(nrow(data2[[i]]),nrow(data2[[i]])*0.8)
      trainData2[[i]]<-data2[[i]][trainIndex2,]
      testData2[[i]]<-data2[[i]][-trainIndex2,]
      trainData2[[i]]$Disease = as.factor(trainData2[[i]]$Disease)
      testData2[[i]]$Disease = as.factor(testData2[[i]]$Disease)
      re_rf2[[i]] = randomForest(Disease~FBratio+Age,
                                 data = trainData2[[i]],ntree=500)
      pred_rf2[[i]]=predict(re_rf2[[i]],newdata=testData2[[i]],type="prob") 
      pred2[[i]]<-prediction(pred_rf2[[i]][,2],testData2[[i]]$Disease)
      perf2[[i]]<-performance(pred2[[i]],"tpr","fpr")
      auc2[[i]] <- performance(pred2[[i]],'auc')
      auc2[[i]] = unlist(slot(auc2[[i]],"y.values"))
      auc2[[i]]=round(auc2[[i]],digits = 3)
      
      set.seed(1001)
      trainIndex3<-sample(nrow(data3[[i]]),nrow(data3[[i]])*0.8)
      trainData3[[i]]<-data3[[i]][trainIndex3,]
      testData3[[i]]<-data3[[i]][-trainIndex3,]
      trainData3[[i]]$Disease = as.factor(trainData3[[i]]$Disease)
      testData3[[i]]$Disease = as.factor(testData3[[i]]$Disease)
      re_rf3[[i]] = randomForest(Disease~FBratio+Age,
                                 data = trainData3[[i]],ntree=500)
      pred_rf3[[i]]=predict(re_rf3[[i]],newdata=testData3[[i]],type="prob") 
      pred3[[i]]<-prediction(pred_rf3[[i]][,2],testData3[[i]]$Disease)
      perf3[[i]]<-performance(pred3[[i]],"tpr","fpr")
      auc3[[i]] <- performance(pred3[[i]],'auc')
      auc3[[i]] = unlist(slot(auc3[[i]],"y.values"))
      auc3[[i]]=round(auc3[[i]],digits = 3)
      NC_res[[i]]<-list(PerRF=perf[[i]],AUC=auc[[i]])
      NA_res[[i]]<-list(PerRF=perf1[[i]],AUC=auc1[[i]])
      CA_res[[i]]<-list(PerRF=perf2[[i]],AUC=auc2[[i]])
      NCA_res[[i]]<-list(PerRF=perf3[[i]],AUC=auc3[[i]])
      names(NC_res)[i]=names(dflist)[i]
      names(NA_res)[i]=names(dflist)[i]
      names(CA_res)[i]=names(dflist)[i]
      names( NCA_res)[i]=names(dflist)[i]
      print(plot(perf[[i]],col='blue',colorize=F,fontsize=8,
                 xlim=c(0,1), ylim=c(0,1),
                 main=paste(names(dflist)[i]))+
              plot(perf1[[i]],col='red',colorize=F,add=T,
                   xlim=c(0,1), ylim=c(0,1))+
              plot(perf2[[i]],col='black',colorize=F,add=T,
                   xlim=c(0,1), ylim=c(0,1))+
              plot(perf3[[i]],col='orange',colorize=F,add=T,
                   xlim=c(0,1), ylim=c(0,1))+
              abline(0,1))
      legend(0.4,0.2,text.col=c('blue','red','black','orange'),cex=0.7,
             legend = c(paste0(predictor[1],":",auc[[i]]), 
                        paste0(predictor[2],":",auc1[[i]]),
                        paste0(predictor[3],":",auc2[[i]]),
                        paste0(predictor[4],":",auc3[[i]])),
             lty = 1, lwd=2,col = c('blue','red','black','orange'))

    }

<img src="BMI_and_irAEs_files/figure-markdown_strict/unnamed-chunk-11-1.png" width="80%" style="display: block; margin: auto;" /><img src="BMI_and_irAEs_files/figure-markdown_strict/unnamed-chunk-11-2.png" width="80%" style="display: block; margin: auto;" /><img src="BMI_and_irAEs_files/figure-markdown_strict/unnamed-chunk-11-3.png" width="80%" style="display: block; margin: auto;" /><img src="BMI_and_irAEs_files/figure-markdown_strict/unnamed-chunk-11-4.png" width="80%" style="display: block; margin: auto;" /><img src="BMI_and_irAEs_files/figure-markdown_strict/unnamed-chunk-11-5.png" width="80%" style="display: block; margin: auto;" />

    dev.off()

    ## integer(0)
    ## integer(0)
    ## integer(0)
    ## integer(0)
    ## integer(0)
    ## null device 
    ##           1

#### 1.4.2.4 AUC results

    ML_res<-list(Normal_CRC=NC_res,
                 Normal_Adenoma=NA_res,
                 CRC_Adenoma=CA_res,
                 Normal_CRCAdenoma=NCA_res)

    auc_res<-data.frame(dataset=names(dflist),"Normal vs CRC"=unlist(auc),
                        "Normal vs Adenoma" =unlist(auc1),
                        "Adenoma vs CRC" =unlist(auc2),
                        "Normal vs CRC/Adenoma" =unlist(auc3))
    auc_res%>%
      knitr::kable(caption = "auc_res") 
    auc_res<-melt(auc_res,id.vars = "dataset",value.name = "AUC",variable.name = "Group")

    ## Warning in melt(auc_res, id.vars = "dataset", value.name = "AUC", variable.name
    ## = "Group"): The melt generic in data.table has been passed a data.frame and will
    ## attempt to redirect to the relevant reshape2 method; please note that reshape2
    ## is deprecated, and this redirection is now deprecated as well. To continue using
    ## melt methods from reshape2 while both libraries are attached, e.g. melt.list,
    ## you can prepend the namespace like reshape2::melt(auc_res). In the next version,
    ## this warning will become an error.

    ggplot(auc_res,aes(Group,AUC,fill=Group))+
      geom_boxplot(outlier.size = 0.5,outlier.alpha = 0.5)+
      geom_jitter(aes(color=dataset))+
      theme_few(base_size = 8)+
      theme(axis.text.x = element_text(angle = 60,vjust = 1,hjust = 1),
            axis.title.x = element_blank())+
      scale_fill_jama()+
      scale_color_lancet()

<img src="BMI_and_irAEs_files/figure-markdown_strict/unnamed-chunk-12-1.png" width="50%" style="display: block; margin: auto;" />
<table>
<caption>
auc\_res
</caption>
<thead>
<tr>
<th style="text-align:left;">
dataset
</th>
<th style="text-align:right;">
Normal.vs.CRC
</th>
<th style="text-align:right;">
Normal.vs.Adenoma
</th>
<th style="text-align:right;">
Adenoma.vs.CRC
</th>
<th style="text-align:right;">
Normal.vs.CRC.Adenoma
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
FengQ\_2015
</td>
<td style="text-align:right;">
0.667
</td>
<td style="text-align:right;">
0.752
</td>
<td style="text-align:right;">
0.567
</td>
<td style="text-align:right;">
0.598
</td>
</tr>
<tr>
<td style="text-align:left;">
PRJEB6070
</td>
<td style="text-align:right;">
0.702
</td>
<td style="text-align:right;">
0.416
</td>
<td style="text-align:right;">
0.508
</td>
<td style="text-align:right;">
0.481
</td>
</tr>
<tr>
<td style="text-align:left;">
PRJNA290926
</td>
<td style="text-align:right;">
0.561
</td>
<td style="text-align:right;">
0.617
</td>
<td style="text-align:right;">
0.512
</td>
<td style="text-align:right;">
0.614
</td>
</tr>
<tr>
<td style="text-align:left;">
ThomasAM\_2018a
</td>
<td style="text-align:right;">
0.571
</td>
<td style="text-align:right;">
0.571
</td>
<td style="text-align:right;">
0.433
</td>
<td style="text-align:right;">
0.675
</td>
</tr>
<tr>
<td style="text-align:left;">
ZellerG\_2014
</td>
<td style="text-align:right;">
0.700
</td>
<td style="text-align:right;">
0.339
</td>
<td style="text-align:right;">
0.648
</td>
<td style="text-align:right;">
0.466
</td>
</tr>
</tbody>
</table>


    ### Gut FBratio in immunotherapy

    #### Gut FBratio and PFS 


    ```r
    ggplot(subset(df1,PFS_3_month!="NA"),
           aes(PFS_3_month,log10(FBratio),fill=PFS_3_month))+
      geom_boxplot(outlier.size = 0.2,outlier.alpha = 0.5)+
      geom_jitter(color="black",alpha=0.8,size=1)+
      stat_compare_means(bracket.size = 0.1,ref.group = "S_PFS",
                         label = "p.signif",vjust = 0.5)+
      theme_few()+
      facet_grid(~Time,scales = "free",space = "free_x")+
      scale_fill_manual(values = col31[c(10,18)])+
      theme(axis.title.x = element_blank(),
            strip.text.x = element_text(size = 8),
            legend.position = "top",legend.text = element_text(size = 8))

<img src="BMI_and_irAEs_files/figure-markdown_strict/unnamed-chunk-13-1.png" width="50%" style="display: block; margin: auto;" />

    surv_cut_PFS <- surv_cutpoint(
      df1,
      time = "PFSmonth",
      event = "PFS",
      variables = c("FBratio"))
    df1$FB<-ifelse(df1$FBratio>=summary(surv_cut_PFS)[,1],"High","Low")
    fit1<-survfit(Surv(PFSmonth,PFS) ~ FB,
                  data = df1)
    fit1
    ggsurvplot(fit1, df1,xlab = "Time(months)",
               censor.size=0.5, size = 0.5,
               tables.theme = theme_few(base_size = 5),
               legend.labs = c("High","Low"),
                    legend.title = "FBratio",
                    risk.table = T,
                    pval = TRUE,pval.size = 3, 
                    pval.coord=c(0.8,0.2),pval.method=F,
                    pval.method.coord=c(0.05,0.3), 
                    ggtheme = theme_minimal() + 
                      theme(line = element_line(size = 0.1),
                            text  = element_text(size = 6)),
                    risk.table.col = "strata",
                    surv.median.line = "hv",
                    risk.table.y.text.col = T,
                    risk.table.y.text = FALSE )

<img src="BMI_and_irAEs_files/figure-markdown_strict/unnamed-chunk-13-2.png" width="50%" style="display: block; margin: auto;" />

    ## Call: survfit(formula = Surv(PFSmonth, PFS) ~ FB, data = df1)
    ## 
    ##    25 observations deleted due to missingness 
    ##          n events median 0.95LCL 0.95UCL
    ## FB=High 10     10   5.68    3.32      NA
    ## FB=Low  13     13   2.63    2.56      NA

#### 1.4.2.5 Gut FBratio and Response

    ICI_stat<-fread("../Data/Data/publicData/ICI_FBratio_stat.csv",data.table = F)
    ICI_stat$Cancer<-factor(ICI_stat$Cancer,levels = c("Melanoma","RCC","NSCLC"))
    sumrepdat <- summarySE(ICI_stat, measurevar = "FBratio", groupvars=c("Cancer", "Response"))
    ggplot(ICI_stat, aes(x = Cancer, y =  FBratio, fill = Response)) +
      stat_compare_means(comparisons = list(c("Melanoma","RCC"),
                                            c("NSCLC","RCC"),c("Melanoma","NSCLC")),label = "p.signif")+
      geom_flat_violin(aes(fill =  Response),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
      geom_point(aes(x = as.numeric(Cancer)-.15, y = FBratio, colour =  Response),position = position_jitter(width = .05), size = .25, shape = 20)+
      geom_boxplot(aes(x = Cancer, y = FBratio, fill =  Response),outlier.shape = NA, alpha = .2, width = .4, colour = "black")+
      geom_line(data = sumrepdat, aes(x = as.numeric(Cancer)+.1, y = FBratio_mean, group =  Response, colour =  Response), linetype = 3)+
      geom_point(data = sumrepdat, aes(x = as.numeric(Cancer)+.1, y = FBratio_mean, group =  Response, colour =  Response),shape = 18) +
      geom_errorbar(data = sumrepdat, aes(x = as.numeric(Cancer)+.1, y = FBratio_mean, group =  Response, colour =  Response, ymin = FBratio_mean-se, ymax = FBratio_mean+se), width = .05)+
      scale_colour_aaas()+
      scale_fill_aaas()+
     # ggtitle("Figure 12: Repeated Measures - Factorial (Extended)")+
      theme_few()+
      ylab("Log10(FBratio)")+
      theme(axis.title.x = element_blank())

<img src="BMI_and_irAEs_files/figure-markdown_strict/unnamed-chunk-14-1.png" width="50%" style="display: block; margin: auto;" />

    ICI_stat$Response<-factor(ICI_stat$Response,levels = c("R","NR"))
    ggstatsplot::grouped_ggbarstats(data = ICI_stat,x=Response,
                                    ggtheme= theme_few(base_size = 6),
                                    y = FBratio_group,results.subtitle=F,
                                    grouping.var = Cancer,
                                    ggstatsplot.layer = FALSE,
                                    messages = FALSE,
                                    main = Response, nboot = 10,
                                    legend.title = "Response")

<img src="BMI_and_irAEs_files/figure-markdown_strict/unnamed-chunk-14-2.png" width="50%" style="display: block; margin: auto;" />

    ggstatsplot::ggbarstats(data = ICI_stat,x=Response,ggtheme = ggplot2::theme_bw(base_size=6),
                                    y = FBratio_group,
                                    ggstatsplot.layer = FALSE,
                                    messages = FALSE,
                            package = "ggsci",
                            palette = "default_aaas",
                                    main = Response, nboot = 10,
                                    legend.title = "Response")

<img src="BMI_and_irAEs_files/figure-markdown_strict/unnamed-chunk-14-3.png" width="50%" style="display: block; margin: auto;" />

2 Treated-related aderse events
===============================

2.1 PFS survival curves for each events
---------------------------------------


    df <- fread("../Data/Data/paired_BL_treat_16patients.csv", data.table = F)
    df$Hand_food_syndrom <- as.factor(df$Hand_food_syndrom)
    df$Hand_food_syndrom_g <- ifelse(df$Hand_food_syndrom %in% c("0", "1"), "no", "yes")
    df$Rash <- as.factor(df$Rash)
    df$Rash_g <- ifelse(df$Rash == "0", "no", "yes")
    df$Fever <- as.factor(df$Fever)
    df$Fever_g <- ifelse(df$Fever == "0", "no", "yes")
    df$Diarrhea <- as.factor(df$Diarrhea)
    df$Diarrhea_g <- ifelse(df$Diarrhea == "0", "no", "yes")

    df_treat <- subset(df, Group == "Treat")
    # List of ggsurvplots

    splots <- list()
    fit_PFS <- survfit(Surv(PFStime, PFS) ~ Hand_food_syndrom_g, data = df_treat)
    fit_PFS
    fit_OS <- survfit(Surv(OStime, OS) ~ Hand_food_syndrom_g, data = df_treat)
    fit_OS
    splots[[1]] <- surv_plot(fit_PFS, df_treat, colors = c("darkgreen", "darkorange"), 
        title = "HandFoodSyndrom_PFS")
    ## Loading required package: prodlim
    splots[[2]] <- surv_plot(fit_OS, df_treat, colors = c("black", "red"), title = "HandFoodSyndrom_OS")

    fit_PFS <- survfit(Surv(PFStime, PFS) ~ Rash_g, data = df_treat)
    fit_PFS
    fit_OS <- survfit(Surv(OStime, OS) ~ Rash_g, data = df_treat)
    fit_OS
    splots[[3]] <- surv_plot(fit_PFS, df_treat, colors = c("darkgreen", "darkorange"), 
        title = "Rash_PFS")
    splots[[4]] <- surv_plot(fit_OS, df_treat, colors = c("black", "red"), title = "Rash_OS")

    fit_PFS <- survfit(Surv(PFStime, PFS) ~ Fever_g, data = df_treat)
    fit_PFS
    fit_OS <- survfit(Surv(OStime, OS) ~ Fever_g, data = df_treat)
    fit_OS
    splots[[5]] <- surv_plot(fit_PFS, df_treat, colors = c("darkgreen", "darkorange"), 
        title = "Fever_PFS")
    splots[[6]] <- surv_plot(fit_OS, df_treat, colors = c("black", "red"), title = "Fever_OS")
    fit_PFS <- survfit(Surv(PFStime, PFS) ~ Diarrhea_g, data = df_treat)
    fit_PFS
    fit_OS <- survfit(Surv(OStime, OS) ~ Diarrhea_g, data = df_treat)
    fit_OS
    splots[[7]] <- surv_plot(fit_PFS, df_treat, colors = c("darkgreen", "darkorange"), 
        title = "Diarrhea_PFS")
    splots[[8]] <- surv_plot(fit_OS, df_treat, colors = c("black", "red"), title = "Diarrhea_OS")
    require(survminer)
    arrange_ggsurvplots(x = splots, print = TRUE, ncol = 4, nrow = 2)

<img src="BMI_and_irAEs_files/figure-markdown_strict/unnamed-chunk-15-1.png" style="display: block; margin: auto;" />

    ## Call: survfit(formula = Surv(PFStime, PFS) ~ Hand_food_syndrom_g, data = df_treat)
    ## 
    ##                         n events median 0.95LCL 0.95UCL
    ## Hand_food_syndrom_g=no  7      7    2.3    1.97      NA
    ## Hand_food_syndrom_g=yes 9      7    4.2    2.20      NA
    ## Call: survfit(formula = Surv(OStime, OS) ~ Hand_food_syndrom_g, data = df_treat)
    ## 
    ##                         n events median 0.95LCL 0.95UCL
    ## Hand_food_syndrom_g=no  7      2     NA    5.17      NA
    ## Hand_food_syndrom_g=yes 9      3   15.5   10.33      NA
    ## Call: survfit(formula = Surv(PFStime, PFS) ~ Rash_g, data = df_treat)
    ## 
    ##             n events median 0.95LCL 0.95UCL
    ## Rash_g=no  11     10   3.80    2.03      NA
    ## Rash_g=yes  5      4   4.23    2.20      NA
    ## Call: survfit(formula = Surv(OStime, OS) ~ Rash_g, data = df_treat)
    ## 
    ##             n events median 0.95LCL 0.95UCL
    ## Rash_g=no  11      3   10.3    10.3      NA
    ## Rash_g=yes  5      2   15.5    15.5      NA
    ## Call: survfit(formula = Surv(PFStime, PFS) ~ Fever_g, data = df_treat)
    ## 
    ##              n events median 0.95LCL 0.95UCL
    ## Fever_g=no  14     12   3.05    2.03      NA
    ## Fever_g=yes  2      2   5.28    4.23      NA
    ## Call: survfit(formula = Surv(OStime, OS) ~ Fever_g, data = df_treat)
    ## 
    ##              n events median 0.95LCL 0.95UCL
    ## Fever_g=no  14      4   15.5    15.5      NA
    ## Fever_g=yes  2      1   10.3      NA      NA
    ## Call: survfit(formula = Surv(PFStime, PFS) ~ Diarrhea_g, data = df_treat)
    ## 
    ##                 n events median 0.95LCL 0.95UCL
    ## Diarrhea_g=no  13     11   4.23    2.30      NA
    ## Diarrhea_g=yes  3      3   1.97    1.87      NA
    ## Call: survfit(formula = Surv(OStime, OS) ~ Diarrhea_g, data = df_treat)
    ## 
    ##                 n events median 0.95LCL 0.95UCL
    ## Diarrhea_g=no  13      4   15.5    10.3      NA
    ## Diarrhea_g=yes  3      1     NA     3.9      NA

2.2 Abundance of Diarrhea-related *Desulfovibrionaceae*
-------------------------------------------------------

    bar1<-ggplot(df,aes(Group,Desulfovibrionaceae,fill=Response))+
      geom_boxplot()+
      geom_line(aes(group=patientID,color=Response,size=Desulfovibrionaceae),alpha=0.5)+
      geom_point(aes(size=Desulfovibrionaceae),color="darkblue",alpha=0.5)+
      theme_few(base_size = 8)+
      stat_compare_means(label = "p.signif")+
      scale_fill_d3()+
      theme(legend.key = element_blank(),
            axis.title.x = element_blank())
    bar2<-ggplot(df,aes(Group,Desulfovibrionaceae,fill=Diarrhea_g))+
      geom_boxplot()+
      geom_line(aes(group=patientID,color=Diarrhea_g,size=Desulfovibrionaceae),alpha=0.5)+
      geom_point(aes(size=Desulfovibrionaceae),color="darkblue",alpha=0.5)+
      theme_few(base_size = 8)+
      stat_compare_means(label = "p.signif")+
      scale_fill_jama()+
      theme(legend.key = element_blank(),
            axis.title.x = element_blank())

    plot_grid(bar1, bar2, labels = c("A", "B"), ncol = 2, nrow = 1)

<img src="BMI_and_irAEs_files/figure-markdown_strict/unnamed-chunk-17-1.png" width="50%" style="display: block; margin: auto;" />
