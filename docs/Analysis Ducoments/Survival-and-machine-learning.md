---
sort: 5
---
-   [1 Requires](#requires)
-   [2 Risk prediction model for PFS](#risk-prediction-model-for-pfs)
    -   [2.1 univariable CoxPH screening](#univariable-coxph-screening)
    -   [2.2 multivariable CoxPH
        screening](#multivariable-coxph-screening)
    -   [2.3 Final CoxPH modle](#final-coxph-modle)
    -   [2.4 Reclassified patients based on CoxPH
        modle](#reclassified-patients-based-on-coxph-modle)
    -   [2.5 PFS survival curve](#pfs-survival-curve)
    -   [2.6 ROC and Precision-Recall
        curves](#roc-and-precision-recall-curves)
    -   [2.7 Time‐dependent ROC curves](#timedependent-roc-curves)

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
    library(forestmodel)
    library(precrec)
    library(timeROC)
    library(survival)
    library(survivalROC)
    library("scales")
    pal_nejm("default")(8)
    "%ni%" <- Negate("%in%")
    theme_set(theme_few())
    options(stringsAsFactors = F)

    ## [1] "#BC3C29FF" "#0072B5FF" "#E18727FF" "#20854EFF" "#7876B1FF" "#6F99ADFF"
    ## [7] "#FFDC91FF" "#EE4C97FF"

</details>

2 Risk prediction model for PFS
===============================

2.1 univariable CoxPH screening
-------------------------------

    PFSdata<-fread("../Data/Data/coxModle_rawdata.csv",data.table = F)
    expr<-list()
    surv<-list()
    group_data<-list()
    survival_dat<-list()
    survival_dat_merge<-list()
    checkGroup<-list()
    covariates<-list()
    colnames(PFSdata)
    expr<-PFSdata[,-c(1:3)]
    rownames(expr)<-PFSdata$patientID
    surv<-PFSdata[,c(2,3)]
    colnames(surv)<-c("event","time")
    rownames(surv)<-PFSdata$patientID
    group_data <- apply(expr, 2 , function(genus){
      name <- colnames(genus)
      genus <- unlist(genus)
      group <- ifelse(genus >= median(genus), 'high', 'low')
      names(group) <- name
      return(group)
    })
    group_data <- as.data.frame(group_data, stringsAsFactors = F)

    survival_dat <- data.frame(row.names = rownames(surv),status = surv$event,
                               time = surv$time,
                               stringsAsFactors = F)
    survival_dat<-survival_dat[which(rownames(survival_dat)%in%rownames(expr)),]
    survival_dat_merge <- cbind(survival_dat,group_data)

    checkGroup<-apply(survival_dat_merge,2 , function(genus){
      facter_lenth<-length(levels(factor(genus)))
      check<-ifelse(facter_lenth==1,facter_lenth,"Yes")
      names(check)<-colnames(genus)
      return(check)
    })

    checkGroup <- as.data.frame(checkGroup, stringsAsFactors = F)
    covariates <- as.character(rownames(subset(checkGroup,checkGroup!=1)))[-c(1,2)]

    univ_formulas <- sapply(covariates,
                            function(x){
                              ##print(x)
                              as.formula(paste('Surv(time, status)~', x))
                            })
    univ_formulas[1:3]

    univ_models<-list()
    for (j in 1:length(univ_formulas)){
      print(paste0("j=",j))
      univ_models[[j]]<-coxph(univ_formulas[[j]], data = survival_dat_merge)  
    }

    univ_results<-list()

    univ_results<- lapply(univ_models,
                          function(x){
                            x <- summary(x)
                            p.value <- signif(x$wald["pvalue"], digits = 2)
                            beta <- signif(x$coef[1], digits = 2)
                            HR <- signif(x$coef[2], digits = 2)
                            HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                            HR.confint.upper <- signif(x$conf.int[,"upper .95"], 2)
                            HR <- paste0(HR, " (",
                                         HR.confint.lower, "-", HR.confint.upper, ")")
                            res <- c(beta, HR, p.value)
                            names(res) <- c("coef", "HR (95% CI for HR)", "p.value")
                            return(res)
                          })
    res_single<- as.data.frame(t(do.call(cbind, univ_results)))
    rownames(res_single)<-covariates
    res_single$p.value=round(as.numeric(res_single$p.value),4)
    res_single <- res_single[res_single$p.value <= 0.2, ]
    res_single <- res_single[order(res_single$p.value), ]
    knitr::kable(res_single)
    single_pick <-rownames(res_single)

    ##  [1] "patientID"       "PFS"             "PFStime"         "Gender"         
    ##  [5] "Age"             "BMI"             "History"         "Smoking"        
    ##  [9] "Dringking"       "anitEGFR"        "antiVEGF"        "LiverM"         
    ## [13] "LungM"           "LymphM"          "PeritonealM"     "OtherM"         
    ## [17] "MetastasisNum"   "Location"        "Treatment_lines" "Fusobacteriota" 
    ## [21] "Bacteroidetes"   "FBratio"         "Firmicutes"      "Proteobacteria" 
    ## [25] "shannon"         "simpson"         "Alistipes"       "Fusobacterium"  
    ## $Age
    ## Surv(time, status) ~ Age
    ## <environment: 0x7fa1ada526d8>
    ## 
    ## $BMI
    ## Surv(time, status) ~ BMI
    ## <environment: 0x7fa1ada4f938>
    ## 
    ## $antiVEGF
    ## Surv(time, status) ~ antiVEGF
    ## <environment: 0x7fa1ad93dba0>
    ## 
    ## [1] "j=1"
    ## [1] "j=2"
    ## [1] "j=3"
    ## [1] "j=4"
    ## [1] "j=5"
    ## [1] "j=6"
    ## [1] "j=7"
    ## [1] "j=8"
    ## [1] "j=9"
    ## [1] "j=10"
    ## [1] "j=11"
    ## [1] "j=12"
    ## [1] "j=13"
    ## [1] "j=14"
    ## [1] "j=15"
    ## [1] "j=16"

<table>
<thead>
<tr class="header">
<th style="text-align: left;"></th>
<th style="text-align: left;">coef</th>
<th style="text-align: left;">HR (95% CI for HR)</th>
<th style="text-align: right;">p.value</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">MetastasisNum</td>
<td style="text-align: left;">-1</td>
<td style="text-align: left;">0.36 (0.13-0.95)</td>
<td style="text-align: right;">0.039</td>
</tr>
<tr class="even">
<td style="text-align: left;">BMI</td>
<td style="text-align: left;">0.67</td>
<td style="text-align: left;">2 (0.91-4.2)</td>
<td style="text-align: right;">0.085</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Fusobacteriota</td>
<td style="text-align: left;">-0.65</td>
<td style="text-align: left;">0.52 (0.24-1.2)</td>
<td style="text-align: right;">0.110</td>
</tr>
<tr class="even">
<td style="text-align: left;">Fusobacterium</td>
<td style="text-align: left;">-0.65</td>
<td style="text-align: left;">0.52 (0.24-1.2)</td>
<td style="text-align: right;">0.110</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Proteobacteria</td>
<td style="text-align: left;">-0.62</td>
<td style="text-align: left;">0.54 (0.25-1.2)</td>
<td style="text-align: right;">0.120</td>
</tr>
<tr class="even">
<td style="text-align: left;">shannon</td>
<td style="text-align: left;">-0.54</td>
<td style="text-align: left;">0.58 (0.28-1.2)</td>
<td style="text-align: right;">0.160</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Alistipes</td>
<td style="text-align: left;">-0.49</td>
<td style="text-align: left;">0.61 (0.29-1.3)</td>
<td style="text-align: right;">0.200</td>
</tr>
</tbody>
</table>

2.2 multivariable CoxPH screening
---------------------------------

    fmla <- as.formula(paste0("Surv(time, status) ~",paste0(single_pick,collapse = '+')))
    colnames(PFSdata)[c(2,3)]=c("status","time")
    survival_dat_merge$BMI<-factor(survival_dat_merge$BMI,levels = c("low","high"))
    survival_dat_merge$Age<-factor(survival_dat_merge$Age,levels = c("low","high"))
    survival_dat_merge$Proteobacteria<-factor(survival_dat_merge$Proteobacteria,levels = c("low","high"))
    survival_dat_merge$shannon<-factor(survival_dat_merge$shannon,levels = c("low","high"))
    survival_dat_merge$Fusobacteriota<-factor(survival_dat_merge$Fusobacteriota,levels = c("low","high"))
    cox <- coxph(fmla, data = survival_dat_merge)
    cox=step(cox,direction = "both")
    coVar<-gsub("high","",names(cox$coefficients))
    final_fmla <- as.formula(paste0("Surv(time, status) ~",paste0(coVar,collapse = '+')))
    model <- coxph( final_fmla ,data = survival_dat_merge )

    ## Start:  AIC=151.14
    ## Surv(time, status) ~ MetastasisNum + BMI + Fusobacteriota + Fusobacterium + 
    ##     Proteobacteria + shannon + Alistipes
    ## 
    ## 
    ## Step:  AIC=151.14
    ## Surv(time, status) ~ MetastasisNum + BMI + Fusobacteriota + Proteobacteria + 
    ##     shannon + Alistipes
    ## 
    ##                  Df    AIC
    ## - Alistipes       1 149.77
    ## - shannon         1 149.94
    ## - MetastasisNum   1 150.30
    ## - BMI             1 150.56
    ## <none>              151.14
    ## - Proteobacteria  1 151.76
    ## - Fusobacteriota  1 153.27
    ## 
    ## Step:  AIC=149.77
    ## Surv(time, status) ~ MetastasisNum + BMI + Fusobacteriota + Proteobacteria + 
    ##     shannon
    ## 
    ##                  Df    AIC
    ## - MetastasisNum   1 148.75
    ## - BMI             1 149.33
    ## <none>              149.77
    ## - Proteobacteria  1 150.45
    ## + Alistipes       1 151.14
    ## - Fusobacteriota  1 151.47
    ## - shannon         1 153.75
    ## 
    ## Step:  AIC=148.75
    ## Surv(time, status) ~ BMI + Fusobacteriota + Proteobacteria + 
    ##     shannon
    ## 
    ##                  Df    AIC
    ## <none>              148.75
    ## + MetastasisNum   1 149.77
    ## - Fusobacteriota  1 150.20
    ## + Alistipes       1 150.30
    ## - Proteobacteria  1 150.70
    ## - BMI             1 151.31
    ## - shannon         1 155.28

2.3 Final CoxPH modle
---------------------

    ggforest(model,data = survival_dat_merge,
             noDigits = 2,fontsize =0.5)+
      theme_few()

<img src="Survival-and-machine-learning_files/figure-markdown_strict/unnamed-chunk-4-1.png" width="50%" style="display: block; margin: auto;" />

2.4 Reclassified patients based on CoxPH modle
----------------------------------------------

    riskScore=predict(cox,type="risk",newdata=survival_dat_merge)
    risk=as.vector(ifelse(riskScore>1,"high","low"))
    multiCOX_risk_result<-cbind(id=rownames(cbind(survival_dat_merge[,1:2],riskScore,risk)),
                                cbind(survival_dat_merge[,1:2],riskScore,risk))

2.5 PFS survival curve
----------------------

    rt1<-multiCOX_risk_result
    colnames(rt1)[2:3]=c("PFS","PFStime")
    #rt1<-merge(dplyr::select(data,c(PFStime,PFS,SangerID,patientID)),rt1,by="SangerID")

    fit1<-survfit(Surv(PFStime,PFS) ~ risk,
                  data = rt1)
    fit1

    p1<-ggsurvplot(fit1, data=rt1,pval.method = T,
                   tables.theme = theme_bw(base_size = 3),
                   risk.table = F,
                   pval = TRUE,
                   ggtheme = theme_survminer(),
                   palette = c("#FFDC91FF","#20854EFF"),
                   legend.title="coxph_risk",
                   risk.table.col = "strata",
                   surv.median.line = "hv",
                   risk.table.y.text.col = T,
                   risk.table.y.text = T )

    ## Call: survfit(formula = Surv(PFStime, PFS) ~ risk, data = rt1)
    ## 
    ##            n events median 0.95LCL 0.95UCL
    ## risk=high 17     17   1.97    1.87     2.3
    ## risk=low  15     12   4.20    2.27      NA

    p1

<img src="Survival-and-machine-learning_files/figure-markdown_strict/unnamed-chunk-7-1.png" width="30%" style="display: block; margin: auto;" />

2.6 ROC and Precision-Recall curves
-----------------------------------

    rt<-multiCOX_risk_result
    sscurves <- evalmod(scores = rt$time, labels = rt$risk)
    knitr::kable(auc(sscurves))
    precrec_obj <- evalmod(scores = rt$time, labels = rt$risk)
    sspoints <- evalmod(mode = "basic", scores = rt$time, labels = rt$risk)
    autoplot(sspoints)

![](Survival-and-machine-learning_files/figure-markdown_strict/unnamed-chunk-8-1.png)

<table>
<thead>
<tr class="header">
<th style="text-align: left;">modnames</th>
<th style="text-align: right;">dsids</th>
<th style="text-align: left;">curvetypes</th>
<th style="text-align: right;">aucs</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">m1</td>
<td style="text-align: right;">1</td>
<td style="text-align: left;">ROC</td>
<td style="text-align: right;">0.7843137</td>
</tr>
<tr class="even">
<td style="text-align: left;">m1</td>
<td style="text-align: right;">1</td>
<td style="text-align: left;">PRC</td>
<td style="text-align: right;">0.8041761</td>
</tr>
</tbody>
</table>

    autoplot(sscurves,size=6)

<img src="Survival-and-machine-learning_files/figure-markdown_strict/unnamed-chunk-9-1.png" width="30%" style="display: block; margin: auto;" />

2.7 Time‐dependent ROC curves
-----------------------------

    rt<-multiCOX_risk_result
    rt<-rt[,-2]
    rt$risk_binary<-ifelse(rt$risk=="high",1,0)
    ## Define a helper functio nto evaluate at various t
    survivalROC_helper <- function(t) {
      survivalROC(Stime        = rt$time,
                  status       = rt$risk_binary,
                  marker       = rt$riskScore,
                  predict.time = t,
                  method       = "NNE",
                  span = 0.25 * nrow(rt)^(-0.20))
    }
    ## Evaluate every 2.5 month
    survivalROC_data <- data_frame(t = 3 * c(1,2,3,3.66)) %>%
      mutate(survivalROC = map(t, survivalROC_helper),
             ## Extract scalar AUC
             auc = map_dbl(survivalROC, magrittr::extract2, "AUC"),
             ## Put cut off dependent values in a data_frame
             df_survivalROC = map(survivalROC, function(obj) {
               as_data_frame(obj[c("cut.values","TP","FP")])
             })) %>%
      dplyr::select(-survivalROC) %>%
      unnest(cols = c(df_survivalROC)) %>%
      arrange(t, FP, TP)

    survivalROC_data %>%
      ggplot(mapping = aes(x = FP, y = TP)) +
      geom_point() +
      geom_line() +
      geom_label(data = survivalROC_data %>% dplyr::select(t,auc) %>% unique,
                 mapping = aes(label = paste0("AUC :",sprintf("%.3f", auc))), x = 0.5, y = 0.5) +
      facet_wrap( ~ t) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
            legend.key = element_blank(),
            plot.title = element_text(hjust = 0.5),
            strip.background = element_blank())

![](Survival-and-machine-learning_files/figure-markdown_strict/unnamed-chunk-10-1.png)

    SROC= survivalROC(Stime = rt$time, status = rt$risk_binary,
                      marker = rt$riskScore,    
                      predict.time = 11, method= "KM" ) 

    SROC1= survivalROC(Stime = rt$time, status = rt$risk_binary,
                       marker = rt$riskScore,    
                       predict.time =11,  span = 0.25*NROW(rt)^(-0.20)) 

    plot(SROC$FP,SROC$TP, type="l", xlim=c(0,1), ylim=c(0,1),  
         ylab = "TP",main = "11-month PFS ROC", col="#0072B5FF",xlab="FP")
    lines(SROC1$FP, SROC1$TP, type="l",col="#E18727FF",xlim=c(0,1), ylim=c(0,1))
    legend(0.3,0.2,c(paste("AUC of KM =",round(SROC$AUC,3)),
                     paste("AUC of NNE =",round(SROC1$AUC,3))),
           x.intersp=1, y.intersp=1,
           lty= 1 ,lwd= 2,col=c( "#0072B5FF","#E18727FF" ),
           bty = "n",
           seg.len=1,cex=1)
    abline(0,1)

<img src="Survival-and-machine-learning_files/figure-markdown_strict/unnamed-chunk-11-1.png" width="50%" style="display: block; margin: auto;" />
